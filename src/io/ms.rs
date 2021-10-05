use flate2::read::GzDecoder;
use ndarray::{Array2, Axis};
use rubbl_casatables::{GlueDataType, Table, TableOpenMode};
use std::{
    fs::create_dir_all,
    path::{Path, PathBuf},
};
use tar::Archive;

use lazy_static::lazy_static;


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

}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    use itertools::izip;
    use tempfile::tempdir;
    // use float_cmp::{approx_eq, F32Margin, F64Margin};
    use approx::abs_diff_eq;

    use num_complex::Complex;

    type c32 = Complex<f32>;
    type c64 = Complex<f64>;

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

}
