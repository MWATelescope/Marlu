use flate2::read::GzDecoder;
use ndarray::{Array2, Axis};
use rubbl_casatables::{GlueDataType, Table, TableOpenMode};
use std::{
    fs::create_dir_all,
    path::{Path, PathBuf},
};
use tar::Archive;

use lazy_static::lazy_static;

use crate::io::error::MeasurementSetWriteError;

lazy_static! {
    static ref DEFAULT_TABLES_GZ: &'static [u8] =
        include_bytes!("../../data/default_tables.tar.gz");
    static ref SOURCE_TABLE_GZ: &'static [u8] = include_bytes!("../../data/source_table.tar.gz");
}

const PKG_VERSION: &'static str = env!("CARGO_PKG_VERSION");
const PKG_NAME: &'static str = env!("CARGO_PKG_NAME");

/// A helper struct to write out a uvfits file.
///
/// TODO: reading
// pub struct MeasurementSetWriter<'a> {
pub struct MeasurementSetWriter {
    /// The path to the root of the measurement set (should end in .ms)
    path: PathBuf,
}

// impl<'a> MeasurementSetWriter<'a> {
impl MeasurementSetWriter {
    pub fn new<T: AsRef<Path>>(path: T) -> Self {
        // let
        MeasurementSetWriter {
            path: path.as_ref().to_path_buf(),
        }
    }

    /// Create the default measurement set tables from a compressed archive
    pub fn decompress_default_tables(&self) -> Result<(), std::io::Error> {
        let tar = GzDecoder::new(&DEFAULT_TABLES_GZ[..]);
        let mut archive = Archive::new(tar);
        if !(self.path.exists() && self.path.is_dir()) {
            create_dir_all(&self.path)?;
        }
        archive.unpack(&self.path)?;
        Ok(())
    }

    /// Create the SOURCE table, as described in `casacore::MSSource`
    pub fn decompress_source_table(&self) -> Result<(), std::io::Error> {
        let tar = GzDecoder::new(&SOURCE_TABLE_GZ[..]);
        let mut archive = Archive::new(tar);
        let source_table_path = self.path.join("SOURCE");
        if !(source_table_path.exists() && source_table_path.is_dir()) {
            create_dir_all(&source_table_path)?;
        }
        archive.unpack(&source_table_path)?;
        Ok(())
    }

    /// Add additional columns / tables which are created in `cotter::MSWriter::initialize()`
    pub fn add_cotter_mods(&self, num_channels: usize) {
        let comment = format!(
            "added by {} {}, emulating cotter::MSWriter::initialize()",
            PKG_VERSION, PKG_NAME
        );
        let mut main_table = Table::open(self.path.clone(), TableOpenMode::ReadWrite).unwrap();
        // TODO: why isn't it let data_shape = [4, num_channels as _];
        let data_shape = [num_channels as _, 4];
        main_table
            .add_array_column(
                GlueDataType::TpComplex,
                "DATA",
                Some(comment.as_str()),
                Some(&data_shape),
                false,
                false,
            )
            .unwrap();
        main_table
            .add_array_column(
                GlueDataType::TpFloat,
                "WEIGHT_SPECTRUM",
                Some(comment.as_str()),
                Some(&data_shape),
                false,
                false,
            )
            .unwrap();

        let source_table_path = self.path.join("SOURCE");
        let mut source_table = Table::open(&source_table_path, TableOpenMode::ReadWrite).unwrap();
        source_table
            .add_array_column(
                GlueDataType::TpDouble,
                "REST_FREQUENCY",
                Some(comment.as_str()),
                None,
                false,
                false,
            )
            .unwrap();

        // TODO: whatever the hell this is
        // TableMeasRefDesc measRef(MFrequency::DEFAULT);
        // TableMeasValueDesc measVal(sourceTableDesc, MSSource::columnName(MSSourceEnums::REST_FREQUENCY));
        // TableMeasDesc<MFrequency> restFreqColMeas(measVal, measRef);
        // // write makes the Measure column persistent.
        // restFreqColMeas.write(sourceTableDesc);

        main_table
            .put_table_keyword("SOURCE", source_table)
            .unwrap();
    }

    /// Write a row into the spectral window table. Return the row index.
    ///
    /// - `name` - Spectral Window name (`NAME` column)
    /// - `ref_freq` - Reference frequency (`REF_FREQUENCY` column)
    /// - `chan_info` - A two-dimensional array of shape (n, 4), containing the
    ///     following for each channel:
    ///     - `CHAN_FREQ` - the center frequencies
    ///     - `CHAN_WIDTH` - channel widths,
    ///     - `EFFECTIVE_BW` - effective noise bandwidths
    ///     - `RESOLUTION` - resolutions.
    /// - `total_bw` - Total bandwidth (`TOTAL_BANDWIDTH` column)
    /// - `flag` - Row flag (`FLAG_ROW` column)
    pub fn write_spectral_window_row(
        &self,
        name: &str,
        ref_freq: f64,
        chan_info: Array2<f64>,
        total_bw: f64,
        flag: bool,
    ) -> Result<u64, MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        let spw_table_path = self.path.join("SPECTRAL_WINDOW");
        let mut spw_table = Table::open(&spw_table_path, TableOpenMode::ReadWrite).unwrap();
        let spw_idx = spw_table.n_rows();

        match chan_info.shape() {
            [num_chans, 4] => {
                spw_table.add_rows(1).unwrap();
                spw_table
                    .put_cell("NUM_CHAN", spw_idx, &(*num_chans as i32))
                    .unwrap();
            }
            _ => {
                return Err(MeasurementSetWriteError::BadSpwChannelInfoShape {
                    expected: "[n, 4]".into(),
                    received: format!("{:?}", chan_info.shape()).into(),
                })
            }
        }

        spw_table
            .put_cell("NAME", spw_idx, &name.to_string())
            .unwrap();
        spw_table
            .put_cell("REF_FREQUENCY", spw_idx, &ref_freq)
            .unwrap();

        let col_names = ["CHAN_FREQ", "CHAN_WIDTH", "EFFECTIVE_BW", "RESOLUTION"];
        for (value, &col_name) in chan_info.lanes(Axis(0)).into_iter().zip(col_names.iter()) {
            spw_table
                .put_cell(col_name, spw_idx, &value.to_owned())
                .unwrap();
        }

        spw_table.put_cell("MEAS_FREQ_REF", spw_idx, &5).unwrap(); // 5 means "TOPO"
        spw_table
            .put_cell("TOTAL_BANDWIDTH", spw_idx, &total_bw)
            .unwrap();
        spw_table.put_cell("FLAG_ROW", spw_idx, &flag).unwrap();

        Ok(spw_idx)
    }

    /// Write a row into the data description table. Return the row index.
    ///
    /// - `spectral_window_id` - Pointer to spectralwindow table
    /// - `polarization_id` - Pointer to polarization table
    /// - `flag_row` - Flag this row
    pub fn write_data_description_row(
        &self,
        spectral_window_id: i32,
        polarization_id: i32,
        flag_row: bool,
    ) -> Result<u64, MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        let ddesc_table_path = self.path.join("DATA_DESCRIPTION");
        let mut ddesc_table = Table::open(&ddesc_table_path, TableOpenMode::ReadWrite).unwrap();
        let ddesc_idx = ddesc_table.n_rows();
        ddesc_table.add_rows(1).unwrap();

        ddesc_table
            .put_cell("SPECTRAL_WINDOW_ID", ddesc_idx, &spectral_window_id)
            .unwrap();
        ddesc_table
            .put_cell("POLARIZATION_ID", ddesc_idx, &polarization_id)
            .unwrap();
        ddesc_table
            .put_cell("FLAG_ROW", ddesc_idx, &flag_row)
            .unwrap();
        Ok(ddesc_idx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    use approx::abs_diff_eq;
    use itertools::izip;
    use tempfile::tempdir;

    use rubbl_core::Complex;

    lazy_static! {
        static ref PATH_1254670392: PathBuf =
            "tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms".into();
    }

    macro_rules! assert_table_column_names_match {
        ( $left:expr, $right:expr ) => {
            match (&$left.column_names(), &$right.column_names()) {
                (Ok(left_columns), Ok(right_columns)) => {
                    for (col_idx, (left_col, right_col)) in
                        izip!(left_columns, right_columns).enumerate()
                    {
                        assert_eq!(
                            left_col, right_col,
                            "column names at index {} do not match. {} != {}",
                            col_idx, left_col, right_col
                        );
                    }
                }
                (Err(left_err), _) => {
                    panic!("could not read left table column names. {:?}", left_err);
                }
                (_, Err(right_err)) => {
                    panic!("could not read right table column names. {:?}", right_err);
                }
            }
        };
    }

    macro_rules! assert_table_nrows_match {
        ( $left:expr, $right:expr ) => {
            assert_eq!(
                $left.n_rows(),
                $right.n_rows(),
                "row counts do not match. {} != {}",
                $left.n_rows(),
                $right.n_rows(),
            )
        };
    }

    macro_rules! assert_table_column_descriptions_match {
        ( $left:expr, $right:expr, $col_name:expr ) => {
            match (
                &$left.get_col_desc($col_name),
                &$right.get_col_desc($col_name),
            ) {
                (Ok(left_desc), Ok(right_desc)) => {
                    assert_eq!(
                        left_desc, right_desc,
                        "column descriptions at column {} do not match.",
                        $col_name
                    );
                }
                (Err(left_err), _) => {
                    panic!(
                        "could not read left column description for column {}. {:?}",
                        $col_name, left_err
                    );
                }
                (_, Err(right_err)) => {
                    panic!(
                        "could not read right column description for column {}. {:?}",
                        $col_name, right_err
                    );
                }
            }
        };
    }

    macro_rules! assert_table_column_values_match {
        ( $left:expr, $right:expr, $col_name:expr, $col_desc:expr, $type_:ty ) => {
            if $col_desc.is_scalar() {
                match (
                    &$left.get_col_as_vec::<$type_>($col_name),
                    &$right.get_col_as_vec::<$type_>($col_name),
                ) {
                    (Ok(left_col), Ok(right_col)) => {
                        assert_eq!(left_col, right_col, "column {} does not match.", $col_name);
                    }
                    (Err(left_err), _) => {
                        panic!(
                            "could not read left column description for column {}. {:?}",
                            $col_name, left_err
                        );
                    }
                    (_, Err(right_err)) => {
                        panic!(
                            "could not read right column description for column {}. {:?}",
                            $col_name, right_err
                        );
                    }
                }
            } else {
                for row_idx in 0..($left.n_rows()) {
                    match (
                        &$left.get_cell_as_vec::<$type_>($col_name, row_idx),
                        &$right.get_cell_as_vec::<$type_>($col_name, row_idx),
                    ) {
                        (Ok(left_cell), Ok(right_cell)) => {
                            for (vec_idx, (left_value, right_value)) in
                                izip!(left_cell, right_cell).enumerate()
                            {
                                assert_eq!(
                                    left_value, right_value,
                                    "cells don't match at index {} in column {}.",
                                    vec_idx, $col_name
                                );
                            }
                        }
                        (Err(left_err), _) => {
                            panic!(
                                "could not read left cell value for column {}, row {}. {:?}",
                                $col_name, row_idx, left_err
                            );
                        }
                        (_, Err(right_err)) => {
                            panic!(
                                "could not read right cell value for column {}, row {}. {:?}",
                                $col_name, row_idx, right_err
                            );
                        }
                    };
                }
            }
        };
    }

    macro_rules! assert_table_column_values_match_approx {
        ( $left:expr, $right:expr, $col_name:expr, $col_desc:expr, $type_:ty ) => {
            if $col_desc.is_scalar() {
                match (
                    &$left.get_col_as_vec::<$type_>($col_name),
                    &$right.get_col_as_vec::<$type_>($col_name),
                ) {
                    (Ok(left_col), Ok(right_col)) => {
                        for (row_idx, (left_cell, right_cell)) in
                            izip!(left_col, right_col).enumerate()
                        {
                            assert!(
                                abs_diff_eq!(left_cell, right_cell, epsilon = 1e-7),
                                "cells don't match in column {}, row {}.",
                                $col_name,
                                row_idx,
                            );
                        }
                    }
                    (Err(left_err), _) => {
                        panic!(
                            "could not read left column description for column {}. {:?}",
                            $col_name, left_err
                        );
                    }
                    (_, Err(right_err)) => {
                        panic!(
                            "could not read right column description for column {}. {:?}",
                            $col_name, right_err
                        );
                    }
                }
            } else {
                for row_idx in 0..($left.n_rows()) {
                    match (
                        &$left.get_cell_as_vec::<$type_>($col_name, row_idx),
                        &$right.get_cell_as_vec::<$type_>($col_name, row_idx),
                    ) {
                        (Ok(left_cell), Ok(right_cell)) => {
                            for (vec_idx, (left_value, right_value)) in
                                izip!(left_cell, right_cell).enumerate()
                            {
                                assert!(
                                    abs_diff_eq!(left_value, right_value, epsilon = 1e-7),
                                    "cells don't match at index {} in column {}, row {}.",
                                    vec_idx,
                                    $col_name,
                                    row_idx,
                                );
                            }
                        }
                        (Err(left_err), _) => {
                            panic!(
                                "could not read left cell value for column {}, row {}. {:?}",
                                $col_name, row_idx, left_err
                            );
                        }
                        (_, Err(right_err)) => {
                            panic!(
                                "could not read right cell value for column {}, row {}. {:?}",
                                $col_name, row_idx, right_err
                            );
                        }
                    };
                }
            }
        };
    }

    macro_rules! assert_tables_match {
        ( $left:expr, $right:expr ) => {
            assert_table_column_names_match!($left, $right);
            assert_table_nrows_match!($left, $right);
            for col_name in $left.column_names().unwrap().iter() {
                assert_table_column_descriptions_match!($left, $right, col_name);
                let col_desc = $left.get_col_desc(col_name).unwrap();
                match col_desc.data_type() {
                    GlueDataType::TpBool => {
                        assert_table_column_values_match!($left, $right, col_name, col_desc, bool)
                    }
                    GlueDataType::TpChar => {
                        assert_table_column_values_match!($left, $right, col_name, col_desc, i8)
                    }
                    GlueDataType::TpUChar => {
                        assert_table_column_values_match!($left, $right, col_name, col_desc, u8)
                    }
                    GlueDataType::TpShort => {
                        assert_table_column_values_match!($left, $right, col_name, col_desc, i16)
                    }
                    GlueDataType::TpUShort => {
                        assert_table_column_values_match!($left, $right, col_name, col_desc, u16)
                    }
                    GlueDataType::TpInt => {
                        assert_table_column_values_match!($left, $right, col_name, col_desc, i32)
                    }
                    GlueDataType::TpUInt => {
                        assert_table_column_values_match!($left, $right, col_name, col_desc, u32)
                    }
                    GlueDataType::TpInt64 => {
                        assert_table_column_values_match!($left, $right, col_name, col_desc, i64)
                    }
                    GlueDataType::TpString => {
                        assert_table_column_values_match!($left, $right, col_name, col_desc, String)
                    }
                    GlueDataType::TpFloat => assert_table_column_values_match_approx!(
                        $left, $right, col_name, col_desc, f32
                    ),
                    GlueDataType::TpDouble => assert_table_column_values_match_approx!(
                        $left, $right, col_name, col_desc, f64
                    ),
                    GlueDataType::TpComplex => assert_table_column_values_match_approx!(
                        $left,
                        $right,
                        col_name,
                        col_desc,
                        Complex<f32>
                    ),
                    GlueDataType::TpDComplex => assert_table_column_values_match_approx!(
                        $left,
                        $right,
                        col_name,
                        col_desc,
                        Complex<f64>
                    ),
                    x => println!("unhandled data type in column {}: {:?}", col_name, x),
                }
            }
        };
    }

    #[test]
    fn test_decompress_default_tables() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        drop(ms_writer);

        assert!(table_path.exists());

        let mut main_table = Table::open(&table_path.clone(), TableOpenMode::Read).unwrap();
        let main_table_keywords = main_table.table_keyword_names().unwrap().clone();
        drop(main_table);

        for (table_name, col_names) in [
            (
                "",
                vec![
                    "TIME",
                    "TIME_CENTROID",
                    "ANTENNA1",
                    "ANTENNA2",
                    "DATA_DESC_ID",
                    "UVW",
                    "INTERVAL",
                    "EXPOSURE",
                    "PROCESSOR_ID",
                    "SCAN_NUMBER",
                    "STATE_ID",
                    "SIGMA",
                    "WEIGHT",
                    "FLAG",
                ],
            ),
            (
                "ANTENNA",
                vec![
                    "OFFSET",
                    "POSITION",
                    "TYPE",
                    "DISH_DIAMETER",
                    "FLAG_ROW",
                    "MOUNT",
                    "NAME",
                    "STATION",
                ],
            ),
            (
                "DATA_DESCRIPTION",
                vec!["FLAG_ROW", "POLARIZATION_ID", "SPECTRAL_WINDOW_ID"],
            ),
            (
                "FEED",
                vec![
                    "POSITION",
                    "BEAM_OFFSET",
                    "POLARIZATION_TYPE",
                    "POL_RESPONSE",
                    "RECEPTOR_ANGLE",
                    "ANTENNA_ID",
                    "BEAM_ID",
                    "FEED_ID",
                    "INTERVAL",
                    "NUM_RECEPTORS",
                    "SPECTRAL_WINDOW_ID",
                    "TIME",
                ],
            ),
            (
                "FIELD",
                vec![
                    "DELAY_DIR",
                    "PHASE_DIR",
                    "REFERENCE_DIR",
                    "CODE",
                    "FLAG_ROW",
                    "NAME",
                    "NUM_POLY",
                    "SOURCE_ID",
                    "TIME",
                ],
            ),
            (
                "FLAG_CMD",
                vec![
                    "APPLIED", "COMMAND", "INTERVAL", "LEVEL", "REASON", "SEVERITY", "TIME", "TYPE",
                ],
            ),
            (
                "HISTORY",
                vec![
                    "APP_PARAMS",
                    "CLI_COMMAND",
                    "APPLICATION",
                    "MESSAGE",
                    "OBJECT_ID",
                    "OBSERVATION_ID",
                    "ORIGIN",
                    "PRIORITY",
                    "TIME",
                ],
            ),
            (
                "OBSERVATION",
                vec![
                    "TIME_RANGE",
                    "LOG",
                    "SCHEDULE",
                    "FLAG_ROW",
                    "OBSERVER",
                    "PROJECT",
                    "RELEASE_DATE",
                    "SCHEDULE_TYPE",
                    "TELESCOPE_NAME",
                ],
            ),
            (
                "POINTING",
                vec![
                    "DIRECTION",
                    "ANTENNA_ID",
                    "INTERVAL",
                    "NAME",
                    "NUM_POLY",
                    "TARGET",
                    "TIME",
                    "TIME_ORIGIN",
                    "TRACKING",
                ],
            ),
            (
                "POLARIZATION",
                vec!["CORR_TYPE", "CORR_PRODUCT", "FLAG_ROW", "NUM_CORR"],
            ),
            (
                "PROCESSOR",
                vec!["FLAG_ROW", "MODE_ID", "TYPE", "TYPE_ID", "SUB_TYPE"],
            ),
            (
                "SPECTRAL_WINDOW",
                vec![
                    "MEAS_FREQ_REF",
                    "CHAN_FREQ",
                    "REF_FREQUENCY",
                    "CHAN_WIDTH",
                    "EFFECTIVE_BW",
                    "RESOLUTION",
                    "FLAG_ROW",
                    "FREQ_GROUP",
                    "FREQ_GROUP_NAME",
                    "IF_CONV_CHAIN",
                    "NAME",
                    "NET_SIDEBAND",
                    "NUM_CHAN",
                    "TOTAL_BANDWIDTH",
                ],
            ),
            (
                "STATE",
                vec![
                    "CAL", "FLAG_ROW", "LOAD", "OBS_MODE", "REF", "SIG", "SUB_SCAN",
                ],
            ),
        ] {
            let mut table = Table::open(&table_path.join(table_name), TableOpenMode::Read).unwrap();
            let mut exp_table =
                Table::open(PATH_1254670392.join(table_name), TableOpenMode::Read).unwrap();
            for col_name in col_names {
                assert_table_column_descriptions_match!(table, exp_table, col_name);
            }
            if table_name != "" {
                assert!(main_table_keywords.contains(&table_name.into()));
            }
        }
    }

    #[test]
    fn test_add_source_table() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_source_table().unwrap();
        drop(ms_writer);

        let mut table = Table::open(&table_path.join("SOURCE"), TableOpenMode::Read).unwrap();
        let mut exp_table =
            Table::open(PATH_1254670392.join("SOURCE"), TableOpenMode::Read).unwrap();
        for col_name in [
            "SOURCE_ID",
            "TIME",
            "INTERVAL",
            "SPECTRAL_WINDOW_ID",
            "NUM_LINES",
            "NAME",
            "CALIBRATION_GROUP",
            "CODE",
            "DIRECTION",
            "PROPER_MOTION",
        ] {
            assert_table_column_descriptions_match!(table, exp_table, col_name);
        }
    }

    #[test]
    fn test_add_cotter_mods() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        // let table_path: PathBuf = "/tmp/marlu.ms".into();
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);
        drop(ms_writer);

        for (table_name, col_names) in [
            ("", vec!["DATA", "WEIGHT_SPECTRUM"]),
            ("SOURCE", vec!["REST_FREQUENCY"]),
        ] {
            let mut table = Table::open(&table_path.join(table_name), TableOpenMode::Read).unwrap();
            let mut exp_table =
                Table::open(PATH_1254670392.join(table_name), TableOpenMode::Read).unwrap();
            for col_name in col_names {
                assert_table_column_descriptions_match!(table, exp_table, col_name);
            }
        }

        let mut main_table = Table::open(&table_path.clone(), TableOpenMode::Read).unwrap();
        let main_table_keywords = main_table.table_keyword_names().unwrap();
        assert!(main_table_keywords.contains(&"SOURCE".into()));
    }

    #[test]
    fn test_write_spectral_window_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        // let table_path: PathBuf = "/tmp/marlu.ms".into();
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let chan_info = Array2::from_shape_fn((768, 4), |(c, i)| {
            (if i == 0 {
                167055000. + (c as f64) * 40000.
            } else {
                40000.
            }) as f64
        });
        let result = ms_writer
            .write_spectral_window_row("MWA_BAND_182.4", 182395000., chan_info, 30720000., false)
            .unwrap();
        drop(ms_writer);

        assert_eq!(result, 0);

        let spw_table_path = table_path.join("SPECTRAL_WINDOW");
        let mut spw_table = Table::open(&spw_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("SPECTRAL_WINDOW"), TableOpenMode::Read).unwrap();

        assert_tables_match!(spw_table, expected_table);
    }

    #[test]
    fn handle_bad_spw_chan_info() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        // let table_path: PathBuf = "/tmp/marlu.ms".into();
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let chan_info = Array2::from_shape_fn((768, 3), |(_, _)| 40000.);
        let result = ms_writer.write_spectral_window_row(
            "MWA_BAND_182.4",
            182395000.,
            chan_info,
            30720000.,
            false,
        );

        assert!(matches!(
            result,
            Err(MeasurementSetWriteError::BadSpwChannelInfoShape { .. })
        ))
    }

    
    #[test]
    fn test_write_data_description_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        // let table_path: PathBuf = "/tmp/marlu.ms".into();
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let result = ms_writer
            .write_data_description_row(0, 0, false)
            .unwrap();
        drop(ms_writer);

        assert_eq!(result, 0);

        let ddesc_table_path = table_path.join("DATA_DESCRIPTION");
        let mut ddesc_table = Table::open(&ddesc_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("DATA_DESCRIPTION"), TableOpenMode::Read).unwrap();

        assert_tables_match!(ddesc_table, expected_table);
    }

}
