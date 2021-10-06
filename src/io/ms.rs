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

    /// Write a row into the SPECTRAL_WINDOW table, unless another table is provided.
    /// Return the row index.
    ///
    /// - `table` - optional [`rubbl_casatables::Table`] object.
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
        table: Option<&mut Table>,
        name: &str,
        ref_freq: f64,
        chan_info: Array2<f64>,
        total_bw: f64,
        flag: bool,
    ) -> Result<u64, MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        // This is to set the default table's lifetime
        let mut table_deref: Table;
        let table = match table {
            Some(table) => table,
            None => {
                let table_path = &self.path.join("SPECTRAL_WINDOW");
                table_deref = Table::open(table_path, TableOpenMode::ReadWrite).unwrap();
                &mut table_deref
            }
        };
        let spw_idx = table.n_rows();

        match chan_info.shape() {
            [num_chans, 4] => {
                table.add_rows(1).unwrap();
                table
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

        table.put_cell("NAME", spw_idx, &name.to_string()).unwrap();
        table.put_cell("REF_FREQUENCY", spw_idx, &ref_freq).unwrap();

        let col_names = ["CHAN_FREQ", "CHAN_WIDTH", "EFFECTIVE_BW", "RESOLUTION"];
        for (value, &col_name) in chan_info.lanes(Axis(0)).into_iter().zip(col_names.iter()) {
            table
                .put_cell(col_name, spw_idx, &value.to_owned())
                .unwrap();
        }

        table.put_cell("MEAS_FREQ_REF", spw_idx, &5).unwrap(); // 5 means "TOPO"
        table
            .put_cell("TOTAL_BANDWIDTH", spw_idx, &total_bw)
            .unwrap();
        table.put_cell("FLAG_ROW", spw_idx, &flag).unwrap();

        Ok(spw_idx)
    }

    /// Write a row into the DATA_DESCRIPTION table, unless another table is provided.
    /// Return the row index.
    ///
    /// - `table` - optional [`rubbl_casatables::Table`] object.
    /// - `spectral_window_id` - Pointer to spectralwindow table
    /// - `polarization_id` - Pointer to polarization table
    /// - `flag_row` - Flag this row
    pub fn write_data_description_row(
        &self,
        table: Option<&mut Table>,
        spectral_window_id: i32,
        polarization_id: i32,
        flag_row: bool,
    ) -> Result<u64, MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        // This is to set the default table's lifetime
        let mut table_deref: Table;
        let table = match table {
            Some(table) => table,
            None => {
                let table_path = &self.path.join("DATA_DESCRIPTION");
                table_deref = Table::open(table_path, TableOpenMode::ReadWrite).unwrap();
                &mut table_deref
            }
        };
        let ddesc_idx = table.n_rows();
        table.add_rows(1).unwrap();

        table
            .put_cell("SPECTRAL_WINDOW_ID", ddesc_idx, &spectral_window_id)
            .unwrap();
        table
            .put_cell("POLARIZATION_ID", ddesc_idx, &polarization_id)
            .unwrap();
        table.put_cell("FLAG_ROW", ddesc_idx, &flag_row).unwrap();
        Ok(ddesc_idx)
    }

    /// Write a row into the `ANTENNA` table, unless another table is provided.
    /// Return the row index.
    ///
    /// - `table` - optional [`rubbl_casatables::Table`] object.
    /// - `name` - Antenna name, e.g. VLA22, CA03
    /// - `station` - Station (antenna pad) name
    /// - `ant_type` - Antenna type (e.g. SPACE-BASED)
    /// - `mount` - Mount type e.g. alt-az, equatorial, etc.
    /// - `position` - Antenna X,Y,Z phase reference position
    /// - `dish_diameter` - Physical diameter of dish
    pub fn write_antenna_row(
        &self,
        table: Option<&mut Table>,
        name: &str,
        station: &str,
        ant_type: &str,
        mount: &str,
        position: &Vec<f64>,
        dish_diameter: f64,
    ) -> Result<u64, MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        // This is to set the default table's lifetime
        let mut table_deref: Table;
        let table = match table {
            Some(table) => table,
            None => {
                let table_path = &self.path.join("ANTENNA");
                table_deref = Table::open(table_path, TableOpenMode::ReadWrite).unwrap();
                &mut table_deref
            }
        };
        let ant_idx = table.n_rows();
        table.add_rows(1).unwrap();

        table.put_cell("NAME", ant_idx, &name.to_string()).unwrap();
        table.put_cell("STATION", ant_idx, &station.to_string()).unwrap();
        table.put_cell("TYPE", ant_idx, &ant_type.to_string()).unwrap();
        table.put_cell("MOUNT", ant_idx, &mount.to_string()).unwrap();
        table.put_cell("POSITION", ant_idx, position).unwrap();
        table
            .put_cell("DISH_DIAMETER", ant_idx, &dish_diameter)
            .unwrap();
        Ok(ant_idx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    use approx::abs_diff_eq;
    use itertools::{Itertools, izip};
    use ndarray::array;
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

    macro_rules! assert_table_columns_match {
        ( $left:expr, $right:expr, $col_name:expr ) => {
            assert_table_column_descriptions_match!($left, $right, $col_name);
            let col_desc = $left.get_col_desc($col_name).unwrap();
            match col_desc.data_type() {
                GlueDataType::TpBool => {
                    assert_table_column_values_match!($left, $right, $col_name, col_desc, bool)
                }
                GlueDataType::TpChar => {
                    assert_table_column_values_match!($left, $right, $col_name, col_desc, i8)
                }
                GlueDataType::TpUChar => {
                    assert_table_column_values_match!($left, $right, $col_name, col_desc, u8)
                }
                GlueDataType::TpShort => {
                    assert_table_column_values_match!($left, $right, $col_name, col_desc, i16)
                }
                GlueDataType::TpUShort => {
                    assert_table_column_values_match!($left, $right, $col_name, col_desc, u16)
                }
                GlueDataType::TpInt => {
                    assert_table_column_values_match!($left, $right, $col_name, col_desc, i32)
                }
                GlueDataType::TpUInt => {
                    assert_table_column_values_match!($left, $right, $col_name, col_desc, u32)
                }
                GlueDataType::TpInt64 => {
                    assert_table_column_values_match!($left, $right, $col_name, col_desc, i64)
                }
                GlueDataType::TpString => {
                    assert_table_column_values_match!($left, $right, $col_name, col_desc, String)
                }
                GlueDataType::TpFloat => assert_table_column_values_match_approx!(
                    $left, $right, $col_name, col_desc, f32
                ),
                GlueDataType::TpDouble => assert_table_column_values_match_approx!(
                    $left, $right, $col_name, col_desc, f64
                ),
                GlueDataType::TpComplex => assert_table_column_values_match_approx!(
                    $left,
                    $right,
                    $col_name,
                    col_desc,
                    Complex<f32>
                ),
                GlueDataType::TpDComplex => assert_table_column_values_match_approx!(
                    $left,
                    $right,
                    $col_name,
                    col_desc,
                    Complex<f64>
                ),
                x => println!("unhandled data type in column {}: {:?}", $col_name, x),
            }
        };
    }
    macro_rules! assert_tables_match {
        ( $left:expr, $right:expr ) => {
            assert_table_column_names_match!($left, $right);
            assert_table_nrows_match!($left, $right);
            for col_name in $left.column_names().unwrap().iter() {
                assert_table_columns_match!($left, $right, col_name);
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
            .write_spectral_window_row(
                None,
                "MWA_BAND_182.4",
                182395000.,
                chan_info,
                30720000.,
                false,
            )
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
            None,
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
            .write_data_description_row(None, 0, 0, false)
            .unwrap();
        drop(ms_writer);

        assert_eq!(result, 0);

        let ddesc_table_path = table_path.join("DATA_DESCRIPTION");
        let mut ddesc_table = Table::open(&ddesc_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table = Table::open(
            PATH_1254670392.join("DATA_DESCRIPTION"),
            TableOpenMode::Read,
        )
        .unwrap();

        assert_tables_match!(ddesc_table, expected_table);
    }

    /// Test data:
    /// ```python
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms/ANTENNA')
    /// tb.getcol("POSITION").transpose()
    ///
    /// ```
    #[test]
    fn test_write_antenna_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        // let table_path: PathBuf = "/tmp/marlu.ms".into();
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let positions = array![
            [-2559524.23682043, 5095846.67363471, -2848988.72758185],
            [-2559573.85868766, 5095824.22162118, -2848984.94323095],
            [-2559579.15163694, 5095819.49141935, -2848988.67657305],
            [-2559586.38326204, 5095812.02407298, -2848995.44896997],
            [-2559577.02462702, 5095811.75478512, -2849004.38523846],
            [-2559562.36386962, 5095816.84349801, -2849008.24740001],
            [-2559570.02742687, 5095812.09631166, -2849009.86884274],
            [-2559516.43210438, 5095827.89444716, -2849029.50443245],
            [-2559542.67487828, 5095845.84987386, -2848973.71803064],
            [-2559535.05568249, 5095846.11608886, -2848980.07958805],
            [-2559551.43200163, 5095842.62700645, -2848971.7358047],
            [-2559578.31186847, 5095827.57088939, -2848974.95061817],
            [-2559592.06244703, 5095821.40759147, -2848973.62713226],
            [-2559585.97629236, 5095821.85852986, -2848978.36154236],
            [-2559590.77956999, 5095816.96664687, -2848982.76126562],
            [-2559505.50905793, 5095867.80668442, -2848967.03341176],
            [-2559623.93203249, 5095764.07097682, -2849046.04852502],
            [-2559627.99317149, 5095772.56137801, -2849027.49081096],
            [-2559615.70261996, 5095785.26937374, -2849016.08079147],
            [-2559603.6968219, 5095792.64269869, -2849014.05600249],
            [-2559597.19950629, 5095802.70810172, -2849002.23111709],
            [-2559594.27461354, 5095808.5234834, -2848994.57029268],
            [-2559604.05430348, 5095805.67436744, -2848990.71350117],
            [-2559612.58057688, 5095801.74907642, -2848989.90035145],
            [-2559598.4859026, 5095812.26227806, -2848984.1037881],
            [-2559605.72288271, 5095808.9984095, -2848983.37295141],
            [-2559597.46151388, 5095816.67254755, -2848977.21883162],
            [-2559618.91111694, 5095802.9097269, -2848982.11385852],
            [-2559612.72483892, 5095807.56133422, -2848979.42404531],
            [-2559637.38376624, 5095788.70978988, -2848990.59316456],
            [-2559612.64392954, 5095811.72451522, -2848972.30299137],
            [-2559646.35453727, 5095789.28297456, -2848981.4854132],
            [-2559554.29367401, 5095807.36283376, -2849032.28419724],
            [-2559562.94692796, 5095809.54957696, -2849020.689028],
            [-2559577.69423937, 5095801.73696795, -2849021.43051958],
            [-2559581.52629478, 5095804.23287999, -2849013.5771202],
            [-2559587.26083052, 5095805.79991338, -2849005.64284638],
            [-2559588.68737068, 5095799.49900897, -2849015.52767181],
            [-2559594.82640672, 5095792.89532829, -2849021.71623476],
            [-2559603.95335699, 5095784.60739105, -2849028.00146619],
            [-2559575.79885824, 5095973.11697647, -2848713.11381342],
            [-2559560.56644989, 5095870.56766467, -2848912.88094033],
            [-2559601.50147232, 5095822.23349762, -2848963.60138247],
            [-2559594.14753799, 5095827.5547599, -2848960.7194767],
            [-2559572.29855907, 5095840.22173704, -2848957.38868042],
            [-2559573.28535158, 5095833.7894811, -2848968.17300371],
            [-2559324.9552274, 5095997.25143554, -2848893.50274043],
            [-2559446.80446483, 5095945.96528192, -2848877.67618636],
            [-2559935.74760999, 5095605.62083096, -2849042.9212563],
            [-2559797.6967689, 5095671.01646134, -2849052.35479234],
            [-2559661.24112631, 5095791.32737906, -2848964.52965989],
            [-2559620.89371688, 5095809.11573327, -2848969.42892551],
            [-2559614.04824423, 5095816.84610126, -2848961.89590579],
            [-2559645.1441408, 5095828.76798455, -2848912.56845977],
            [-2559757.07568572, 5095864.83693965, -2848746.84315802],
            [-2559784.75166311, 5095842.51003371, -2848762.10587847],
            [-2559738.65514646, 5095747.41468543, -2848972.03113726],
            [-2559751.10939783, 5095741.07673235, -2848972.02113379],
            [-2559763.50964093, 5095734.57555975, -2848971.90425245],
            [-2559775.96095263, 5095728.18720491, -2848971.83061971],
            [-2559729.94465145, 5095745.67762816, -2848982.86739382],
            [-2559742.35806589, 5095739.26730021, -2848982.79992339],
            [-2559754.83975112, 5095732.86811016, -2848982.76298578],
            [-2559767.29995078, 5095726.506339, -2848982.67867908],
            [-2559779.74305434, 5095720.0615301, -2848982.62369524],
            [-2559721.25096366, 5095743.94063141, -2848993.71082022],
            [-2559733.67672752, 5095737.55709839, -2848993.65148687],
            [-2559746.15303707, 5095731.15165965, -2848993.58589663],
            [-2559758.59139801, 5095724.75978859, -2848993.5301347],
            [-2559771.04566174, 5095718.38846077, -2848993.46816989],
            [-2559783.48858163, 5095711.938816, -2848993.37739059],
            [-2559712.52999265, 5095742.21390819, -2849004.53813469],
            [-2559724.98764116, 5095735.81588544, -2849004.45374705],
            [-2559737.43228728, 5095729.42537685, -2849004.38995624],
            [-2559762.37758752, 5095716.67568517, -2849004.29019445],
            [-2559774.81691468, 5095710.28575466, -2849004.25233016],
            [-2559787.2616333, 5095703.85974742, -2849004.16250454],
            [-2559716.28667322, 5095734.09556809, -2849015.29631442],
            [-2559728.73117662, 5095727.72260534, -2849015.24596354],
            [-2559741.21602658, 5095721.33641253, -2849015.1831127],
            [-2559753.65009199, 5095714.97163288, -2849015.14974617],
            [-2559766.10324871, 5095708.54238463, -2849015.05814847],
            [-2559778.55264822, 5095702.15690283, -2849015.02566832],
            [-2559720.06369651, 5095725.99544314, -2849026.0930153],
            [-2559732.50110758, 5095719.63509742, -2849026.0265794],
            [-2559744.97867445, 5095713.254445, -2849026.00486744],
            [-2559757.42753429, 5095706.81441589, -2849025.89177366],
            [-2559769.88357391, 5095700.46664613, -2849025.85313128],
            [-2559723.8167694, 5095717.89888394, -2849036.88254641],
            [-2559736.27622115, 5095711.50890434, -2849036.8268114],
            [-2559748.72503752, 5095705.11333816, -2849036.74600889],
            [-2559761.17582237, 5095698.74174742, -2849036.68761551],
            [-2559627.31768437, 5095740.43444833, -2849084.4272822],
            [-2559639.76871994, 5095734.0522183, -2849084.36529039],
            [-2559652.22984472, 5095727.6343786, -2849084.31307286],
            [-2559664.68469533, 5095721.25974039, -2849084.25646519],
            [-2559618.60818901, 5095738.75285625, -2849095.29762911],
            [-2559631.07686103, 5095732.3233046, -2849095.2150006],
            [-2559643.54839723, 5095725.94624065, -2849095.17002014],
            [-2559655.97458262, 5095719.53904331, -2849095.08559228],
            [-2559668.45237752, 5095713.16330236, -2849095.05495177],
            [-2559609.94471255, 5095736.99578644, -2849106.10527376],
            [-2559622.37654902, 5095730.5931573, -2849106.0521702],
            [-2559634.83233422, 5095724.22706593, -2849105.97947743],
            [-2559647.26434808, 5095717.82033456, -2849105.93350307],
            [-2559659.75202994, 5095711.44868177, -2849105.88677754],
            [-2559672.22402978, 5095705.07254134, -2849105.79532843],
            [-2559601.23330926, 5095735.30370979, -2849116.95146649],
            [-2559613.68834554, 5095728.90493985, -2849116.89301227],
            [-2559626.15538672, 5095722.4966482, -2849116.84082857],
            [-2559651.07173741, 5095709.74947699, -2849116.74014692],
            [-2559663.53217349, 5095703.36813688, -2849116.64597861],
            [-2559675.98323678, 5095696.95476912, -2849116.61700223],
            [-2559604.99340349, 5095727.20107635, -2849127.74545497],
            [-2559617.45356918, 5095720.78578159, -2849127.66729721],
            [-2559629.89983307, 5095714.42747003, -2849127.61162303],
            [-2559642.37761784, 5095708.01383455, -2849127.53707051],
            [-2559654.83998482, 5095701.63634271, -2849127.50103223],
            [-2559667.28946303, 5095695.23541424, -2849127.42916509],
            [-2559608.75663761, 5095719.08241484, -2849138.52063277],
            [-2559621.20690764, 5095712.70311669, -2849138.45418358],
            [-2559633.69923758, 5095706.34517109, -2849138.40123508],
            [-2559646.14470992, 5095699.96523292, -2849138.34013613],
            [-2559658.62296149, 5095693.5792595, -2849138.3049907],
            [-2559612.52547674, 5095710.9905088, -2849149.28782889],
            [-2559625.00345224, 5095704.63071907, -2849149.2509654],
            [-2559637.46936929, 5095698.25806493, -2849149.18095809],
            [-2559649.92256756, 5095691.87122651, -2849149.12522308]
        ];

        let names = vec![
            "Tile011", "Tile012", "Tile013", "Tile014", "Tile015", "Tile016", "Tile017", "Tile018",
            "Tile021", "Tile022", "Tile023", "Tile024", "Tile025", "Tile026", "Tile027", "Tile028",
            "Tile031", "Tile032", "Tile033", "Tile034", "Tile035", "Tile036", "Tile037", "Tile038",
            "Tile041", "Tile042", "Tile043", "Tile044", "Tile045", "Tile046", "Tile047", "Tile048",
            "Tile061", "Tile062", "Tile063", "Tile064", "Tile065", "Tile066", "Tile067", "Tile068",
            "Tile081", "Tile082", "Tile083", "Tile084", "Tile085", "Tile086", "Tile087", "Tile088",
            "Tile091", "Tile092", "Tile093", "Tile094", "Tile095", "Tile096", "Tile097", "Tile098",
            "HexE1", "HexE2", "HexE3", "HexE4", "HexE5", "HexE6", "HexE7", "HexE8", "HexE9",
            "HexE10", "HexE11", "HexE12", "HexE13", "HexE14", "HexE15", "HexE16", "HexE17",
            "HexE18", "HexE19", "HexE20", "HexE21", "HexE22", "HexE23", "HexE24", "HexE25",
            "HexE26", "HexE27", "HexE28", "HexE29", "HexE30", "HexE31", "HexE32", "HexE33",
            "HexE34", "HexE35", "HexE36", "HexS1", "HexS2", "HexS3", "HexS4", "HexS5", "HexS6",
            "HexS7", "HexS8", "HexS9", "HexS10", "HexS11", "HexS12", "HexS13", "HexS14", "HexS15",
            "HexS16", "HexS17", "HexS18", "HexS19", "HexS20", "HexS21", "HexS22", "HexS23",
            "HexS24", "HexS25", "HexS26", "HexS27", "HexS28", "HexS29", "HexS30", "HexS31",
            "HexS32", "HexS33", "HexS34", "HexS35", "HexS36",
        ];

        let ddesc_table_path = ms_writer.path.join("ANTENNA");
        let mut ddesc_table = Table::open(ddesc_table_path, TableOpenMode::ReadWrite).unwrap();

        for (idx, (name, position)) in izip!(names, positions.outer_iter()).enumerate() {
            let position = position.iter().cloned().collect();

            let result = ms_writer
            .write_antenna_row(Some(&mut ddesc_table), name, "MWA", "GROUND-BASED", "ALT-AZ", &position, 4.0)
            .unwrap();
            assert_eq!(result, idx as _);
        }

        drop(ms_writer);

        let mut expected_table =
            Table::open(PATH_1254670392.join("ANTENNA"), TableOpenMode::Read).unwrap();

        assert_tables_match!(ddesc_table, expected_table);
    }
}
