use crate::{
    average_chunk_f64, c32,
    io::error::IOError,
    ndarray::{array, Array2, Array3, ArrayView, ArrayView3, Axis},
    num_complex::Complex,
    LatLngHeight, RADec, VisContext,
};

use std::{
    fs::create_dir_all,
    path::{Path, PathBuf},
};

use flate2::read::GzDecoder;
use lazy_static::lazy_static;
use rubbl_casatables::{
    GlueDataType, Table, TableCreateMode, TableDesc, TableDescCreateMode, TableOpenMode,
    TableRecord,
};
use tar::Archive;

use super::VisWritable;

use indicatif::{ProgressDrawTarget, ProgressStyle};

cfg_if::cfg_if! {
    if #[cfg(feature = "mwalib")] {
        use std::{
            ops::Range,
            time::SystemTime,
            f64::consts::FRAC_PI_2,
        };

        use log::trace;
        use itertools::izip;
        use mwalib::CorrelatorContext;

        use crate::{precession::precess_time,
            Jones, ENH, UVW, hifitime::Epoch
        };
    }
}

lazy_static! {
    static ref DEFAULT_TABLES_GZ: &'static [u8] =
        include_bytes!("../../data/default_tables.tar.gz");
    static ref SOURCE_TABLE_GZ: &'static [u8] = include_bytes!("../../data/source_table.tar.gz");
}

const PKG_VERSION: &str = env!("CARGO_PKG_VERSION");
const PKG_NAME: &str = env!("CARGO_PKG_NAME");

/// A helper struct to write out a CASA Measurement Set.
pub struct MeasurementSetWriter {
    /// The path to the root of the measurement set (typically ends in .ms)
    path: PathBuf,

    /// The RA/Dec where this observation is phased to
    phase_centre: RADec,

    /// Array Position [Latitude (radians), Longitude (radians), Height (m)]
    array_pos: LatLngHeight,

    /// The next row to write to in the main table
    pub main_row_idx: usize,
}

impl MeasurementSetWriter {
    pub fn new<T: AsRef<Path>>(
        path: T,
        phase_centre: RADec,
        array_pos: Option<LatLngHeight>,
    ) -> Self {
        let array_pos = match array_pos {
            Some(pos) => pos,
            None => {
                // The results here are slightly different to those given by cotter.
                // This is at least partly due to different constants (the altitude is
                // definitely slightly different), but possibly also because ERFA is
                // more accurate than cotter's "homebrewed" Geodetic2XYZ.
                LatLngHeight::new_mwa()
            }
        };

        MeasurementSetWriter {
            path: path.as_ref().to_path_buf(),
            phase_centre,
            array_pos,
            main_row_idx: 0,
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

    /// Add additional columns / tables / keywords from `cotter::MSWriter::initialize()`
    pub fn add_cotter_mods(&self, num_channels: usize) {
        let comment = format!(
            "added by {} {}, emulating cotter::MSWriter::initialize()",
            PKG_VERSION, PKG_NAME
        );
        let mut main_table = Table::open(&self.path, TableOpenMode::ReadWrite).unwrap();
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

        source_table
            .put_column_keyword("REST_FREQUENCY", "QuantumUnits", &vec!["s".to_string()])
            .unwrap();

        let mut meas_info = TableRecord::new().unwrap();
        meas_info
            .put_field("type", &"frequency".to_string())
            .unwrap();
        meas_info.put_field("Ref", &"LSRK".to_string()).unwrap();

        source_table
            .put_column_keyword("REST_FREQUENCY", "MEASINFO", &meas_info)
            .unwrap();

        main_table
            .put_table_keyword("SOURCE", source_table)
            .unwrap();
    }

    /// Add additional columns / tables / keywords from `cotter::MWAMS::addMWAAntennaFields()`
    pub fn add_mwa_ant_mods(&self) {
        let comment = format!(
            "added by {} {}, emulating cotter::MWAMS::addMWAAntennaFields()",
            PKG_VERSION, PKG_NAME
        );

        let ant_table_path = self.path.join("ANTENNA");
        let mut ant_table = Table::open(&ant_table_path, TableOpenMode::ReadWrite).unwrap();
        ant_table
            .add_array_column(
                GlueDataType::TpInt,
                "MWA_INPUT",
                Some(comment.as_str()),
                None,
                false,
                false,
            )
            .unwrap();
        ant_table
            .add_scalar_column(
                GlueDataType::TpInt,
                "MWA_TILE_NR",
                Some(comment.as_str()),
                false,
                false,
            )
            .unwrap();
        ant_table
            .add_scalar_column(
                GlueDataType::TpInt,
                "MWA_RECEIVER",
                Some(comment.as_str()),
                false,
                false,
            )
            .unwrap();
        ant_table
            .add_array_column(
                GlueDataType::TpInt,
                "MWA_SLOT",
                Some(comment.as_str()),
                None,
                false,
                false,
            )
            .unwrap();
        ant_table
            .add_array_column(
                GlueDataType::TpDouble,
                "MWA_CABLE_LENGTH",
                Some(comment.as_str()),
                None,
                false,
                false,
            )
            .unwrap();
    }

    /// Add additional columns / tables / keywords from `cotter::MWAMS::addMWAFieldFields()`
    pub fn add_mwa_field_mods(&self) {
        let comment = format!(
            "added by {} {}, emulating cotter::MWAMS::addMWAFieldFields()",
            PKG_VERSION, PKG_NAME
        );

        let field_table_path = self.path.join("FIELD");
        let mut field_table = Table::open(&field_table_path, TableOpenMode::ReadWrite).unwrap();
        field_table
            .add_scalar_column(
                GlueDataType::TpBool,
                "MWA_HAS_CALIBRATOR",
                Some(comment.as_str()),
                false,
                false,
            )
            .unwrap();
    }

    /// Add additional columns / tables / keywords from `cotter::MWAMS::addMWAObservationFields()`
    pub fn add_mwa_obs_mods(&self) {
        let comment = format!(
            "added by {} {}, emulating cotter::MWAMS::addMWAObservationFields()",
            PKG_VERSION, PKG_NAME
        );

        let obs_table_path = self.path.join("OBSERVATION");
        let mut obs_table = Table::open(&obs_table_path, TableOpenMode::ReadWrite).unwrap();
        obs_table
            .add_scalar_column(
                GlueDataType::TpDouble,
                "MWA_GPS_TIME",
                Some(comment.as_str()),
                false,
                false,
            )
            .unwrap();
        obs_table
            .add_scalar_column(
                GlueDataType::TpString,
                "MWA_FILENAME",
                Some(comment.as_str()),
                false,
                false,
            )
            .unwrap();
        obs_table
            .add_scalar_column(
                GlueDataType::TpString,
                "MWA_OBSERVATION_MODE",
                Some(comment.as_str()),
                false,
                false,
            )
            .unwrap();
        obs_table
            .add_scalar_column(
                GlueDataType::TpInt,
                "MWA_FLAG_WINDOW_SIZE",
                Some(comment.as_str()),
                false,
                false,
            )
            .unwrap();
        obs_table
            .add_scalar_column(
                GlueDataType::TpDouble,
                "MWA_DATE_REQUESTED",
                Some(comment.as_str()),
                false,
                false,
            )
            .unwrap();

        obs_table
            .put_column_keyword("MWA_DATE_REQUESTED", "QuantumUnits", &vec!["s".to_string()])
            .unwrap();

        let mut meas_info = TableRecord::new().unwrap();
        meas_info.put_field("type", &"epoch".to_string()).unwrap();
        meas_info.put_field("Ref", &"UTC".to_string()).unwrap();

        obs_table
            .put_column_keyword("MWA_DATE_REQUESTED", "MEASINFO", &meas_info)
            .unwrap();
    }

    /// Add additional columns / tables / keywords from `cotter::MWAMS::addMWASpectralWindowFields()`
    pub fn add_mwa_spw_mods(&self) {
        let comment = format!(
            "added by {} {}, emulating cotter::MWAMS::addMWASpectralWindowFields()",
            PKG_VERSION, PKG_NAME
        );

        let spw_table_path = self.path.join("SPECTRAL_WINDOW");
        let mut spw_table = Table::open(&spw_table_path, TableOpenMode::ReadWrite).unwrap();
        spw_table
            .add_scalar_column(
                GlueDataType::TpInt,
                "MWA_CENTRE_SUBBAND_NR",
                Some(comment.as_str()),
                false,
                false,
            )
            .unwrap();
    }

    /// Add additional columns / tables / keywords from `cotter::MWAMS::addMWATilePointingFields()`
    pub fn add_mwa_pointing_mods(&self) {
        let comment = format!(
            "added by {} {}, emulating cotter::MWAMS::addMWATilePointingFields()",
            PKG_VERSION, PKG_NAME
        );

        let mut pointing_table_desc =
            TableDesc::new("MWA_TILE_POINTING", TableDescCreateMode::TDM_SCRATCH).unwrap();

        // let mut pointing_table = Table::open(&pointing_table_path, TableOpenMode::ReadWrite).unwrap();
        pointing_table_desc
            .add_array_column(
                GlueDataType::TpDouble,
                "INTERVAL",
                Some(comment.as_str()),
                None,
                false,
                false,
            )
            .unwrap();
        pointing_table_desc
            .add_array_column(
                GlueDataType::TpInt,
                "DELAYS",
                Some(comment.as_str()),
                None,
                false,
                false,
            )
            .unwrap();
        pointing_table_desc
            .add_array_column(
                GlueDataType::TpDouble,
                "DIRECTION",
                Some(comment.as_str()),
                None,
                false,
                false,
            )
            .unwrap();

        pointing_table_desc
            .put_column_keyword("INTERVAL", "QuantumUnits", &vec!["s".to_string()])
            .unwrap();

        let mut meas_info = TableRecord::new().unwrap();
        meas_info.put_field("type", &"epoch".to_string()).unwrap();
        meas_info.put_field("Ref", &"UTC".to_string()).unwrap();

        pointing_table_desc
            .put_column_keyword("INTERVAL", "MEASINFO", &meas_info)
            .unwrap();

        let pointing_table_path = self.path.join("MWA_TILE_POINTING");
        let pointing_table = Table::new(
            pointing_table_path,
            pointing_table_desc,
            0,
            TableCreateMode::New,
        )
        .unwrap();

        let mut main_table = Table::open(&self.path, TableOpenMode::ReadWrite).unwrap();
        main_table
            .put_table_keyword("MWA_TILE_POINTING", pointing_table)
            .unwrap();
    }

    /// Add additional columns / tables / keywords from `cotter::MWAMS::addMWASubbandFields()`
    pub fn add_mwa_subband_mods(&self) {
        let comment = format!(
            "added by {} {}, emulating cotter::MWAMS::addMWASubbandFields()",
            PKG_VERSION, PKG_NAME
        );

        let mut subband_table_desc =
            TableDesc::new("MWA_SUBBAND", TableDescCreateMode::TDM_SCRATCH).unwrap();

        subband_table_desc
            .add_scalar_column(
                GlueDataType::TpInt,
                "NUMBER",
                Some(comment.as_str()),
                false,
                false,
            )
            .unwrap();
        subband_table_desc
            .add_scalar_column(
                GlueDataType::TpDouble,
                "GAIN",
                Some(comment.as_str()),
                false,
                false,
            )
            .unwrap();
        subband_table_desc
            .add_scalar_column(
                GlueDataType::TpBool,
                "FLAG_ROW",
                Some(comment.as_str()),
                false,
                false,
            )
            .unwrap();

        let subband_table_path = self.path.join("MWA_SUBBAND");
        let subband_table = Table::new(
            subband_table_path,
            subband_table_desc,
            0,
            TableCreateMode::New,
        )
        .unwrap();

        let mut main_table = Table::open(&self.path, TableOpenMode::ReadWrite).unwrap();
        main_table
            .put_table_keyword("MWA_SUBBAND", subband_table)
            .unwrap();
    }

    /// Add additional columns / tables / keywords from `cotter::MWAMS::InitializeMWAFields()`
    pub fn add_mwa_mods(&self) {
        self.add_mwa_ant_mods();
        self.add_mwa_field_mods();
        self.add_mwa_obs_mods();
        self.add_mwa_spw_mods();
        self.add_mwa_pointing_mods();
        self.add_mwa_subband_mods();
    }

    /// Write a row into the SPECTRAL_WINDOW table. Remember to also write to
    /// the DATA_DESCRIPTION table.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `idx` - row index to write to (ensure enough rows have been added)
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
    #[allow(clippy::too_many_arguments)]
    pub fn write_spectral_window_row(
        &self,
        table: &mut Table,
        idx: u64,
        name: &str,
        ref_freq: f64,
        chan_info: Array2<f64>,
        total_bw: f64,
        flag: bool,
    ) -> Result<(), IOError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        match chan_info.shape() {
            [num_chans, 4] => {
                table
                    .put_cell("NUM_CHAN", idx, &(*num_chans as i32))
                    .unwrap();
            }
            sh => {
                return Err(IOError::BadArrayShape {
                    argument: "chan_info".into(),
                    function: "write_spectral_window_row".into(),
                    expected: "[n, 4]".into(),
                    received: format!("{:?}", sh),
                })
            }
        }

        table.put_cell("NAME", idx, &name.to_string()).unwrap();
        table.put_cell("REF_FREQUENCY", idx, &ref_freq).unwrap();

        let col_names = ["CHAN_FREQ", "CHAN_WIDTH", "EFFECTIVE_BW", "RESOLUTION"];
        for (value, &col_name) in chan_info.lanes(Axis(0)).into_iter().zip(col_names.iter()) {
            table.put_cell(col_name, idx, &value.to_vec()).unwrap();
        }

        table.put_cell("MEAS_FREQ_REF", idx, &5).unwrap(); // 5 means "TOPO"
        table.put_cell("TOTAL_BANDWIDTH", idx, &total_bw).unwrap();
        table.put_cell("FLAG_ROW", idx, &flag).unwrap();

        Ok(())
    }

    /// Write a row into the SPECTRAL_WINDOW table with extra mwa columns enabled.
    /// Remember to also write to the DATA_DESCRIPTION table.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `idx` - row index to write to (ensure enough rows have been added)
    /// - `name` - Spectral Window name (`NAME` column)
    /// - `ref_freq` - Reference frequency (`REF_FREQUENCY` column)
    /// - `chan_info` - A two-dimensional array of shape (n, 4), containing the
    ///     following for each channel:
    ///     - `CHAN_FREQ` - the center frequencies
    ///     - `CHAN_WIDTH` - channel widths,
    ///     - `EFFECTIVE_BW` - effective noise bandwidths
    ///     - `RESOLUTION` - resolutions.
    /// - `total_bw` - Total bandwidth (`TOTAL_BANDWIDTH` column)
    /// - `centre_subband_nr` - This is the "sky" channel number of the center coarse channel in
    ///     the spectral window.
    /// - `flag` - Row flag (`FLAG_ROW` column)
    #[allow(clippy::too_many_arguments)]
    pub fn write_spectral_window_row_mwa(
        &self,
        table: &mut Table,
        idx: u64,
        name: &str,
        ref_freq: f64,
        chan_info: Array2<f64>,
        total_bw: f64,
        centre_subband_nr: i32,
        flag: bool,
    ) -> Result<(), IOError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        self.write_spectral_window_row(table, idx, name, ref_freq, chan_info, total_bw, flag)
            .unwrap();

        table
            .put_cell("MWA_CENTRE_SUBBAND_NR", idx, &centre_subband_nr)
            .unwrap();

        Ok(())
    }

    /// Write a row into the DATA_DESCRIPTION table.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `idx` - row index to write to (ensure enough rows have been added)
    /// - `spectral_window_id` - Pointer to spectralwindow table
    /// - `polarization_id` - Pointer to polarization table
    /// - `flag_row` - Flag this row
    pub fn write_data_description_row(
        &self,
        table: &mut Table,
        idx: u64,
        spectral_window_id: i32,
        polarization_id: i32,
        flag_row: bool,
    ) -> Result<(), IOError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        table
            .put_cell("SPECTRAL_WINDOW_ID", idx, &spectral_window_id)
            .unwrap();
        table
            .put_cell("POLARIZATION_ID", idx, &polarization_id)
            .unwrap();
        table.put_cell("FLAG_ROW", idx, &flag_row).unwrap();
        Ok(())
    }

    /// Write a row into the `ANTENNA` table.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `idx` - row index to write to (ensure enough rows have been added)
    /// - `name` - Antenna name, e.g. VLA22, CA03
    /// - `station` - Station (antenna pad) name
    /// - `ant_type` - Antenna type (e.g. SPACE-BASED)
    /// - `mount` - Mount type e.g. alt-az, equatorial, etc.
    /// - `position` - Antenna X,Y,Z phase reference position
    /// - `dish_diameter` - Physical diameter of dish
    /// - `flag_row` - Row flag
    #[allow(clippy::ptr_arg)]
    #[allow(clippy::too_many_arguments)]
    pub fn write_antenna_row(
        &self,
        table: &mut Table,
        idx: u64,
        name: &str,
        station: &str,
        ant_type: &str,
        mount: &str,
        // TODO: should this be an XyzGeodetic/XyzGeocentric?
        position: &Vec<f64>,
        dish_diameter: f64,
        flag_row: bool,
    ) -> Result<(), IOError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        table.put_cell("NAME", idx, &name.to_string()).unwrap();
        table
            .put_cell("STATION", idx, &station.to_string())
            .unwrap();
        table.put_cell("TYPE", idx, &ant_type.to_string()).unwrap();
        table.put_cell("MOUNT", idx, &mount.to_string()).unwrap();
        table.put_cell("POSITION", idx, position).unwrap();
        table
            .put_cell("DISH_DIAMETER", idx, &dish_diameter)
            .unwrap();
        table.put_cell("FLAG_ROW", idx, &flag_row).unwrap();
        Ok(())
    }

    /// Write a row into the `ANTENNA` table with extra mwa columns enabled.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `idx` - row index to write to (ensure enough rows have been added)
    /// - `name` - Antenna name, e.g. VLA22, CA03
    /// - `station` - Station (antenna pad) name
    /// - `ant_type` - Antenna type (e.g. SPACE-BASED)
    /// - `mount` - Mount type e.g. alt-az, equatorial, etc.
    /// - `position` - Antenna X,Y,Z phase reference position
    /// - `dish_diameter` - Physical diameter of dish
    /// - `input` - A vector containing the (legacy) correlator input number for each polarization
    /// - `tile_nr` - The MWA tile ID number
    /// - `receiver` - Receiver number
    /// - `slot` - A vector containing the physical receiver slot number for each polarization
    /// - `cable_length` - A vector containing the electrical length for each polarization
    /// - `flag_row` - Row flag
    #[allow(clippy::ptr_arg)]
    #[allow(clippy::too_many_arguments)]
    pub fn write_antenna_row_mwa(
        &self,
        table: &mut Table,
        idx: u64,
        name: &str,
        station: &str,
        ant_type: &str,
        mount: &str,
        // TODO: should this be an XyzGeodetic/XyzGeocentric?
        position: &Vec<f64>,
        dish_diameter: f64,
        input: &Vec<i32>,
        tile_nr: i32,
        receiver: i32,
        slot: &Vec<i32>,
        cable_length: &Vec<f64>,
        flag_row: bool,
    ) -> Result<(), IOError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        self.write_antenna_row(
            table,
            idx,
            name,
            station,
            ant_type,
            mount,
            position,
            dish_diameter,
            flag_row,
        )
        .unwrap();

        table.put_cell("MWA_INPUT", idx, input).unwrap();
        table.put_cell("MWA_TILE_NR", idx, &tile_nr).unwrap();
        table.put_cell("MWA_RECEIVER", idx, &receiver).unwrap();
        table.put_cell("MWA_SLOT", idx, slot).unwrap();
        table
            .put_cell("MWA_CABLE_LENGTH", idx, cable_length)
            .unwrap();

        Ok(())
    }

    /// Write a row into the `POLARIZATION` table.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `idx` - row index to write to (ensure enough rows have been added)
    /// - `corr_type` - The polarization type for each correlation product, as a Stokes enum.
    /// - `corr_product` - Indices describing receptors of feed going into correlation.
    ///     Shape should be [n, 2] where n is the length of `corr_type`
    /// - `flag_row` - Row flag
    #[allow(clippy::ptr_arg)]
    pub fn write_polarization_row(
        &self,
        table: &mut Table,
        idx: u64,
        corr_type: &Vec<i32>,
        corr_product: &Array2<i32>,
        flag_row: bool,
    ) -> Result<(), IOError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        let num_corr_type = corr_type.len();

        match corr_product.shape() {
            [num_corr, 2] if *num_corr == num_corr_type => {
                table
                    .put_cell("NUM_CORR", idx, &(*num_corr as i32))
                    .unwrap();
            }
            sh => {
                return Err(IOError::BadArrayShape {
                    argument: "corr_product".into(),
                    function: "write_polarization_row".into(),
                    expected: format!("[n, 2] (where n = corr_type.len() = {})", num_corr_type),
                    received: format!("{:?}", sh),
                })
            }
        }

        table.put_cell("CORR_TYPE", idx, corr_type).unwrap();
        table.put_cell("CORR_PRODUCT", idx, corr_product).unwrap();
        table.put_cell("FLAG_ROW", idx, &flag_row).unwrap();
        Ok(())
    }

    /// Write a row into the `SOURCE` table.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `idx` - row index to write to (ensure enough rows have been added)
    /// - `source_id` - Source id
    /// - `time` - Midpoint of time for which this set of parameters is accurate
    /// - `interval` - Interval of time for which this set of parameters is accurate
    /// - `spw_idx` - ID for this spectral window setup
    /// - `name` - Name of this source
    /// - `calibration_group` - Number of grouping for calibration purpose.
    /// - `code` - Special characteristics of source, e.g. Bandpass calibrator
    /// - `direction` - Direction (RA, DEC) [Rad, J2000].
    /// - `proper_motion` - [rad/s].
    #[allow(clippy::too_many_arguments)]
    pub fn write_source_row(
        &self,
        table: &mut Table,
        idx: u64,
        source_id: i32,
        time: f64,
        interval: f64,
        spw_idx: i32,
        num_lines: i32,
        name: &str,
        calibration_group: i32,
        code: &str,
        // TODO: should this be an `RADec`?
        direction: Vec<f64>,
        proper_motion: Vec<f64>,
    ) -> Result<(), IOError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        match (direction.len(), proper_motion.len()) {
            (2, 2) => {
                table.put_cell("DIRECTION", idx, &direction).unwrap();
                table
                    .put_cell("PROPER_MOTION", idx, &proper_motion)
                    .unwrap();
            }
            sh => {
                return Err(IOError::BadArrayShape {
                    argument: "direction|proper_motion".into(),
                    function: "write_source_row".into(),
                    expected: "(2, 2)".into(),
                    received: format!("{:?}", sh),
                })
            }
        }

        table.put_cell("SOURCE_ID", idx, &source_id).unwrap();
        table.put_cell("TIME", idx, &time).unwrap();
        table.put_cell("INTERVAL", idx, &interval).unwrap();
        table.put_cell("SPECTRAL_WINDOW_ID", idx, &spw_idx).unwrap();
        table.put_cell("NUM_LINES", idx, &num_lines).unwrap();
        table.put_cell("NAME", idx, &name.to_string()).unwrap();
        table
            .put_cell("CALIBRATION_GROUP", idx, &calibration_group)
            .unwrap();
        table.put_cell("CODE", idx, &code.to_string()).unwrap();
        Ok(())
    }

    /// Write a row into the `FIELD` table.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `idx` - row index to write to (ensure enough rows have been added)
    /// - `name` - Name of this field
    /// - `code` - Special characteristics of field, e.g. Bandpass calibrator
    /// - `time` - Time origin for direction and rate
    /// - `dir_info` - An array of polynomial coefficients to calculate a direction
    ///     (RA, DEC) relative to `time`. The shape is  [3, p, 2], where p is the maximum
    ///     order of all polynomials, and there are three direction polynomials:
    ///     - `DELAY_DIR` - Direction of delay center (e.g. RA, DEC) in time
    ///     - `PHASE_DIR` - Direction of phase center (e.g. RA, DEC) in time
    ///     - `REFERENCE_DIR` - Direction of reference center (e.g. RA, DEC) in time
    /// - `source_id` - Source id
    /// - `flag_row` - Row Flag
    ///
    /// For the MWA, Phase, Reference and Delay directions are basically the same
    #[allow(clippy::too_many_arguments)]
    pub fn write_field_row(
        &self,
        table: &mut Table,
        idx: u64,
        name: &str,
        code: &str,
        time: f64,
        dir_info: &Array3<f64>,
        source_id: i32,
        flag_row: bool,
    ) -> Result<(), IOError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        match dir_info.shape() {
            [3, p, 2] if *p > 0 => {
                table.put_cell("NUM_POLY", idx, &((*p - 1) as i32)).unwrap();
            }
            sh => {
                return Err(IOError::BadArrayShape {
                    argument: "dir_info".into(),
                    function: "write_field_row".into(),
                    expected: "[3, p, 2] (where p is highest polynomial order)".into(),
                    received: format!("{:?}", sh),
                })
            }
        }

        table.put_cell("NAME", idx, &name.to_string()).unwrap();
        table.put_cell("CODE", idx, &code.to_string()).unwrap();
        table.put_cell("TIME", idx, &time).unwrap();

        let col_names = ["DELAY_DIR", "PHASE_DIR", "REFERENCE_DIR"];
        for (value, &col_name) in dir_info.outer_iter().zip(col_names.iter()) {
            table.put_cell(col_name, idx, &value.to_owned()).unwrap();
        }

        table.put_cell("SOURCE_ID", idx, &source_id).unwrap();
        table.put_cell("FLAG_ROW", idx, &flag_row).unwrap();
        Ok(())
    }

    /// Write a row into the `FIELD` table with extra mwa columns enabled.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `idx` - row index to write to (ensure enough rows have been added)
    /// - `name` - Name of this field
    /// - `code` - Special characteristics of field, e.g. Bandpass calibrator
    /// - `time` - Time origin for direction and rate
    /// - `dir_info` - An array of polynomial coefficients to calculate a direction
    ///     (RA, DEC) relative to `time`. The shape is  [3, p, 2], where p is the maximum
    ///     order of all polynomials, and there are three direction polynomials:
    ///     - `DELAY_DIR` - Direction of delay center (e.g. RA, DEC) in time
    ///     - `PHASE_DIR` - Direction of phase center (e.g. RA, DEC) in time
    ///     - `REFERENCE_DIR` - Direction of reference center (e.g. RA, DEC) in time
    /// - `source_id` - Source id
    /// - `has_calibrator` - whether this observation is intended for calibration (from metafits:CALIBRAT)
    /// - `flag_row` - Row Flag
    #[allow(clippy::too_many_arguments)]
    pub fn write_field_row_mwa(
        &self,
        table: &mut Table,
        idx: u64,
        name: &str,
        code: &str,
        time: f64,
        dir_info: &Array3<f64>,
        source_id: i32,
        has_calibrator: bool,
        flag_row: bool,
    ) -> Result<(), IOError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        self.write_field_row(table, idx, name, code, time, dir_info, source_id, flag_row)
            .unwrap();

        table
            .put_cell("MWA_HAS_CALIBRATOR", idx, &has_calibrator)
            .unwrap();

        Ok(())
    }

    /// Write a row into the `OBSERVATION` table.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `idx` - row index to write to (ensure enough rows have been added)
    /// - `telescope_name` - Telescope Name (e.g. WSRT, VLBA)
    /// - `time_range` - Start and end of observation
    /// - `observer` - Antenna type (e.g. SPACE-BASED)
    /// - `schedule_type` - Observing schedule type
    /// - `project` - Project identification string
    /// - `release_date` - Release date when data becomes public
    /// - `flag_row` - Row flag
    #[allow(clippy::too_many_arguments)]
    pub fn write_observation_row(
        &self,
        table: &mut Table,
        idx: u64,
        telescope_name: &str,
        time_range: (f64, f64),
        observer: &str,
        schedule_type: &str,
        project: &str,
        release_date: f64,
        flag_row: bool,
    ) -> Result<(), IOError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        table
            .put_cell("TELESCOPE_NAME", idx, &telescope_name.to_string())
            .unwrap();
        let time_range = vec![time_range.0, time_range.1];
        table.put_cell("TIME_RANGE", idx, &time_range).unwrap();
        table
            .put_cell("OBSERVER", idx, &observer.to_string())
            .unwrap();
        table
            .put_cell("SCHEDULE_TYPE", idx, &schedule_type.to_string())
            .unwrap();
        table
            .put_cell("PROJECT", idx, &project.to_string())
            .unwrap();
        table.put_cell("RELEASE_DATE", idx, &release_date).unwrap();
        table.put_cell("FLAG_ROW", idx, &flag_row).unwrap();
        Ok(())
    }

    /// Write a row into the `OBSERVATION` table with extra mwa columns enabled.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `idx` - row index to write to (ensure enough rows have been added)
    /// - `telescope_name` - Telescope Name (e.g. WSRT, VLBA)
    /// - `time_range` - Start and end of observation
    /// - `observer` - Antenna type (e.g. SPACE-BASED)
    /// - `schedule_type` - Observing schedule type
    /// - `project` - Project identification string
    /// - `release_date` - Release date when data becomes public
    /// - `gps_time` - Scheduled start (gps time) of observation (obsid, from metafits:GPSTIME)
    /// - `filename` - Name of observation (from metafits:FILENAME)
    /// - `observation_mode` - Observation mode (from metafits:MODE)
    /// - `flag_window_size` - Number of scans in this partition
    /// - `date_requested` - from metafits:DATE-OBS
    /// - `flag_row` - Row flag
    #[allow(clippy::too_many_arguments)]
    pub fn write_observation_row_mwa(
        &self,
        table: &mut Table,
        idx: u64,
        telescope_name: &str,
        time_range: (f64, f64),
        observer: &str,
        schedule_type: &str,
        project: &str,
        release_date: f64,
        gps_time: f64,
        filename: &str,
        observation_mode: &str,
        flag_window_size: i32,
        date_requested: f64,
        flag_row: bool,
    ) -> Result<(), IOError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        self.write_observation_row(
            table,
            idx,
            telescope_name,
            time_range,
            observer,
            schedule_type,
            project,
            release_date,
            flag_row,
        )
        .unwrap();

        table.put_cell("MWA_GPS_TIME", idx, &gps_time).unwrap();
        table
            .put_cell("MWA_FILENAME", idx, &filename.to_string())
            .unwrap();
        table
            .put_cell("MWA_OBSERVATION_MODE", idx, &observation_mode.to_string())
            .unwrap();
        table
            .put_cell("MWA_FLAG_WINDOW_SIZE", idx, &flag_window_size)
            .unwrap();
        table
            .put_cell("MWA_DATE_REQUESTED", idx, &date_requested)
            .unwrap();
        Ok(())
    }
    /// Write a row into the `HISTORY_ITERM` table.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `idx` - row index to write to (ensure enough rows have been added)
    /// - `time` - Time of message
    /// - `cmd_line` - CLI command sequence
    /// - `message` - Log message
    /// - `application` - Application name
    /// - `params` - Application parameters
    #[allow(clippy::too_many_arguments)]
    pub fn write_history_row(
        &self,
        table: &mut Table,
        idx: u64,
        time: f64,
        cmd_line: &str,
        message: &str,
        application: &str,
        params: &str,
    ) -> Result<(), IOError> {
        let cmd_line: Vec<String> = vec![cmd_line.to_string()];
        let params: Vec<String> = vec![params.to_string()];

        table.put_cell("TIME", idx, &time).unwrap();
        table.put_cell("OBSERVATION_ID", idx, &0).unwrap();
        table
            .put_cell("MESSAGE", idx, &message.to_string())
            .unwrap();
        table
            .put_cell("APPLICATION", idx, &application.to_string())
            .unwrap();
        table
            .put_cell("PRIORITY", idx, &"NORMAL".to_string())
            .unwrap();
        table
            .put_cell("ORIGIN", idx, &"standalone".to_string())
            .unwrap();
        table.put_cell("APP_PARAMS", idx, &params).unwrap();
        table.put_cell("CLI_COMMAND", idx, &cmd_line).unwrap();

        Ok(())
    }

    /// Write a row into the `FEED` table.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `idx` - row index to write to (ensure enough rows have been added)
    /// - `ant_idx` - ID of antenna in this array
    /// - `feed_idx` - Feed ID
    /// - `spw_idx` - Spectral window ID for this setup
    /// - `time` - Midpoint of time for which this set of parameters is accurate
    /// - `interval` - Interval for which this set of parameters is accurate
    /// - `num_receptors` - Number of receptors on this feed (probably 1 or 2)
    /// - `beam_idx` - Beam model ID
    /// - `beam_offset` - Beam position offset (on sky but in antennareference frame)
    /// - `pol_type` - Type of polarization to which a given RECEPTOR responds
    /// - `pol_response` - D-matrix i.e. leakage between two receptors
    /// - `position` - Position of feed relative to feed reference position
    /// - `receptor_angle` - The reference angle for polarization
    #[allow(clippy::ptr_arg)]
    #[allow(clippy::too_many_arguments)]
    pub fn write_feed_row(
        &self,
        table: &mut Table,
        idx: u64,
        ant_idx: i32,
        feed_idx: i32,
        spw_idx: i32,
        time: f64,
        interval: f64,
        num_receptors: i32,
        beam_idx: i32,
        beam_offset: &Array2<f64>,
        pol_type: &Vec<String>,
        pol_response: &Array2<c32>,
        // TODO: should this be an XyzGeodetic/XyzGeocentric?
        position: &Vec<f64>,
        receptor_angle: &Vec<f64>,
    ) -> Result<(), IOError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        if beam_offset.shape() != [num_receptors as usize, 2] {
            return Err(IOError::BadArrayShape {
                argument: "beam_offset".into(),
                function: "write_feed_row".into(),
                expected: "[n, 2]".into(),
                received: format!("{:?}", beam_offset.shape()),
            });
        }
        if pol_type.len() != num_receptors as usize {
            return Err(IOError::BadArrayShape {
                argument: "pol_type".into(),
                function: "write_feed_row".into(),
                expected: "n".into(),
                received: format!("{:?}", pol_type.len()),
            });
        }
        if pol_response.shape() != [num_receptors as usize, num_receptors as usize] {
            return Err(IOError::BadArrayShape {
                argument: "pol_response".into(),
                function: "write_feed_row".into(),
                expected: "[n, n]".into(),
                received: format!("{:?}", pol_response.shape()),
            });
        }
        if position.len() != 3 {
            return Err(IOError::BadArrayShape {
                argument: "position".into(),
                function: "write_feed_row".into(),
                expected: "3".into(),
                received: format!("{:?}", position.len()),
            });
        }
        if receptor_angle.len() != num_receptors as usize {
            return Err(IOError::BadArrayShape {
                argument: "receptor_angle".into(),
                function: "write_feed_row".into(),
                expected: "n".into(),
                received: format!("{:?}", receptor_angle.len()),
            });
        }

        table.put_cell("ANTENNA_ID", idx, &ant_idx).unwrap();
        table.put_cell("FEED_ID", idx, &feed_idx).unwrap();
        table.put_cell("SPECTRAL_WINDOW_ID", idx, &spw_idx).unwrap();
        table.put_cell("TIME", idx, &time).unwrap();
        table.put_cell("INTERVAL", idx, &interval).unwrap();
        table
            .put_cell("NUM_RECEPTORS", idx, &num_receptors)
            .unwrap();
        table.put_cell("BEAM_ID", idx, &beam_idx).unwrap();
        table.put_cell("BEAM_OFFSET", idx, beam_offset).unwrap();
        table.put_cell("POLARIZATION_TYPE", idx, pol_type).unwrap();
        table.put_cell("POL_RESPONSE", idx, pol_response).unwrap();
        table.put_cell("POSITION", idx, position).unwrap();
        table
            .put_cell("RECEPTOR_ANGLE", idx, receptor_angle)
            .unwrap();
        Ok(())
    }

    /// Write a row into the `MWA_TILE_POINTING` table.
    ///
    /// - `start` - start MJD of observation
    /// - `end` - end MJD of observation
    /// - `delays` - beamformer delays, from metafits:DELAYS
    /// - `direction_{ra|dec}` - pointing direction [Ra/Dec]
    #[allow(clippy::ptr_arg)]
    #[allow(clippy::too_many_arguments)]
    pub fn write_mwa_tile_pointing_row(
        &self,
        table: &mut Table,
        idx: u64,
        start: f64,
        end: f64,
        delays: &Vec<i32>,
        direction_ra: f64,
        direction_dec: f64,
    ) -> Result<(), IOError> {
        table.put_cell("INTERVAL", idx, &vec![start, end]).unwrap();
        table.put_cell("DELAYS", idx, delays).unwrap();
        table
            .put_cell("DIRECTION", idx, &vec![direction_ra, direction_dec])
            .unwrap();

        Ok(())
    }

    /// Write a row into the `MWA_SUBBAND` table.
    ///
    /// - `number` - Subband (coarse channel) index
    /// - `gain` - (deprecated) - from metafits:CHANGAIN, use 0.
    /// - `flag_row` - flag this subband
    pub fn write_mwa_subband_row(
        &self,
        table: &mut Table,
        idx: u64,
        number: i32,
        gain: f64,
        flag_row: bool,
    ) -> Result<(), IOError> {
        table.put_cell("NUMBER", idx, &number).unwrap();
        table.put_cell("GAIN", idx, &gain).unwrap();
        table.put_cell("FLAG_ROW", idx, &flag_row).unwrap();

        Ok(())
    }

    fn get_centre_freq(freqs: &[f64]) -> f64 {
        let len = freqs.len();
        if len % 2 == 0 {
            (freqs[len / 2] + freqs[len / 2 - 1]) * 0.5
        } else {
            freqs[len / 2]
        }
    }

    /// Create an MWA measurement set, with all tables (except the main visibility table)
    /// prefilled with metadata from a [`mwalib::CorrelatorContext`]
    ///
    /// `mwalib_timestep_range` the range of timestep indices (according to mwalib)
    /// of the current chunk being written to the measurement set.
    ///
    /// `mwalib_coarse_chan_range` the range of coarse channel indices (according to mwalib)
    /// of the current chunk being written to the measurement set.
    ///
    /// `baseline_idxs` - the range of indices into `CorrelatorContext.metafits_context.baselines`
    ///     corresponding to the third dimension of the jones array.
    ///
    /// `avg_time` - the temporal averaging factor which determines the number of timesteps that
    ///     will be written
    ///
    /// `avg_freq` - the frequency averaging factor which determines the number of frequencies that
    ///     will be written
    #[cfg(feature = "mwalib")]
    pub fn initialize_from_mwalib(
        &self,
        corr_ctx: &CorrelatorContext,
        timestep_range: &Range<usize>,
        coarse_chan_range: &Range<usize>,
        baseline_idxs: &[usize],
        avg_time: usize,
        avg_freq: usize,
    ) -> Result<(), IOError> {
        // TODO: initialise from marlu
        // pub fn initialize(
        //     &self,
        //     vis_ctx: &VisContext,
        //     timestep_range: &Range<usize>,
        //     coarse_chan_range: &Range<usize>,
        //     baseline_idxs: &[usize],
        //     avg_time: usize,
        //     avg_freq: usize,
        // ) -> Result<(), IOError> {
        trace!("initialize_from_mwalib");

        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let num_sel_coarse_chans = coarse_chan_range.len();
        let num_sel_chans = fine_chans_per_coarse * num_sel_coarse_chans;
        let num_avg_chans = (num_sel_chans as f64 / avg_freq as f64).ceil() as usize;
        trace!(
            "fine_chans_per_coarse={}, num_sel_coarse_chans={}, num_sel_chans={}, num_avg_chans={}",
            fine_chans_per_coarse,
            num_sel_coarse_chans,
            num_sel_chans,
            num_avg_chans
        );

        let ants = &corr_ctx.metafits_context.antennas;
        let num_ant_pols = corr_ctx.metafits_context.num_ant_pols;

        // ////// //
        // Cutoff //
        // ////// //

        self.decompress_default_tables().unwrap();
        self.decompress_source_table().unwrap();
        self.add_cotter_mods(num_avg_chans);
        self.add_mwa_mods();

        // //// //
        // Main //
        // //// //

        let num_sel_timesteps = timestep_range.len();
        let num_avg_timesteps = (num_sel_timesteps as f64 / avg_time as f64).ceil() as usize;

        trace!(
            "num_sel_timesteps={}, num_avg_timesteps={}",
            num_sel_timesteps,
            num_avg_timesteps
        );

        let num_baselines = baseline_idxs.len();
        let total_num_rows = num_avg_timesteps * num_baselines;
        let mut main_table = Table::open(&self.path, TableOpenMode::ReadWrite).unwrap();

        main_table.add_rows(total_num_rows).unwrap();

        // /////////////// //
        // Spectral Window //
        // /////////////// //

        let mut spw_table =
            Table::open(&self.path.join("SPECTRAL_WINDOW"), TableOpenMode::ReadWrite).unwrap();

        let mwalib_centre_coarse_chan_idx = coarse_chan_range.start + (num_sel_coarse_chans / 2);
        let centre_coarse_chan =
            &corr_ctx.metafits_context.metafits_coarse_chans[mwalib_centre_coarse_chan_idx];
        let mwalib_start_fine_chan_idx = coarse_chan_range.start * fine_chans_per_coarse;
        let fine_chan_width_hz =
            avg_freq as u32 * corr_ctx.metafits_context.corr_fine_chan_width_hz;

        let avg_fine_chan_freqs_hz: Vec<f64> =
            corr_ctx.metafits_context.metafits_fine_chan_freqs_hz[mwalib_start_fine_chan_idx..]
                .chunks(avg_freq)
                .map(Self::get_centre_freq)
                .collect();

        let chan_info = Array2::from_shape_fn((num_avg_chans, 4), |(c, i)| {
            if i == 0 {
                avg_fine_chan_freqs_hz[c]
            } else {
                fine_chan_width_hz as f64
            }
        });

        let center_freq_hz = Self::get_centre_freq(&avg_fine_chan_freqs_hz);

        spw_table.add_rows(1).unwrap();
        self.write_spectral_window_row_mwa(
            &mut spw_table,
            0,
            format!("MWA_BAND_{:.1}", center_freq_hz / 1_000_000.).as_str(),
            center_freq_hz,
            chan_info,
            fine_chan_width_hz as f64 * num_avg_chans as f64,
            centre_coarse_chan.rec_chan_number as _,
            false,
        )
        .unwrap();

        // //////////////// //
        // Data Description //
        // //////////////// //

        let mut ddesc_table = Table::open(
            &self.path.join("DATA_DESCRIPTION"),
            TableOpenMode::ReadWrite,
        )
        .unwrap();

        ddesc_table.add_rows(1).unwrap();
        self.write_data_description_row(&mut ddesc_table, 0, 0, 0, false)
            .unwrap();

        // //////// //
        // Antennae //
        // //////// //

        let mut ant_table =
            Table::open(&self.path.join("ANTENNA"), TableOpenMode::ReadWrite).unwrap();

        let num_ants = ants.len();
        ant_table.add_rows(num_ants).unwrap();
        for (idx, antenna) in ants.iter().enumerate() {
            let position_enh = ENH {
                e: antenna.east_m,
                n: antenna.north_m,
                h: antenna.height_m,
            };
            let position = position_enh
                .to_xyz(self.array_pos.latitude_rad)
                .to_geocentric(self.array_pos)
                .unwrap();
            self.write_antenna_row_mwa(
                &mut ant_table,
                idx as _,
                &antenna.tile_name,
                "MWA",
                "GROUND-BASED",
                "ALT-AZ",
                &vec![position.x, position.y, position.z],
                4.0,
                &vec![
                    antenna.rfinput_x.input as i32,
                    antenna.rfinput_y.input as i32,
                ],
                antenna.tile_id as i32,
                antenna.rfinput_x.rec_number as i32,
                &vec![
                    antenna.rfinput_x.rec_slot_number as i32,
                    antenna.rfinput_y.rec_slot_number as i32,
                ],
                &vec![
                    antenna.rfinput_x.electrical_length_m,
                    antenna.rfinput_y.electrical_length_m,
                ],
                false,
            )
            .unwrap();
        }

        // //////////// //
        // Polarization //
        // //////////// //
        //
        // MWA always has the following polarizations (feeds):
        // - XX (0, 0)
        // - XY (0, 1)
        // - YX (1, 0)
        // - YY (1, 1)

        let mut pol_table =
            Table::open(&self.path.join("POLARIZATION"), TableOpenMode::ReadWrite).unwrap();

        let corr_product = array![[0, 0], [0, 1], [1, 0], [1, 1]];
        let corr_type = vec![9, 10, 11, 12];
        pol_table.add_rows(1).unwrap();

        self.write_polarization_row(&mut pol_table, 0, &corr_type, &corr_product, false)
            .unwrap();

        // ///// //
        // Field //
        // ///// //

        let mut field_table =
            Table::open(&self.path.join("FIELD"), TableOpenMode::ReadWrite).unwrap();

        let ra_phase_rad = corr_ctx
            .metafits_context
            .ra_phase_center_degrees
            .unwrap()
            .to_radians();
        let dec_phase_rad = corr_ctx
            .metafits_context
            .dec_phase_center_degrees
            .unwrap()
            .to_radians();

        // TODO: get phase centre from self.phase_centre
        // TODO: is dir_info right?
        //  - `DELAY_DIR` - Direction of delay center (e.g. RA, DEC) in time
        //  - `PHASE_DIR` - Direction of phase center (e.g. RA, DEC) in time
        //  - `REFERENCE_DIR` - Direction of reference center (e.g. RA, DEC) in time

        let dir_info = array![
            [[ra_phase_rad, dec_phase_rad]],
            [[ra_phase_rad, dec_phase_rad]],
            [[ra_phase_rad, dec_phase_rad]],
        ];

        let obs_name = &corr_ctx.metafits_context.obs_name;
        let field_name = obs_name
            .rsplit_once('_')
            .unwrap_or((obs_name.as_str(), ""))
            .0;

        let sched_start_time_mjd_utc_s = Epoch::from_gpst_seconds(
            corr_ctx.metafits_context.sched_start_gps_time_ms as f64 / 1e3,
        )
        .as_mjd_utc_seconds();

        field_table.add_rows(1).unwrap();
        self.write_field_row_mwa(
            &mut field_table,
            0,
            field_name,
            "",
            sched_start_time_mjd_utc_s,
            &dir_info,
            -1,
            corr_ctx.metafits_context.calibrator,
            false,
        )
        .unwrap();

        // ////// //
        // Source //
        // ////// //

        let mut source_table =
            Table::open(&self.path.join("SOURCE"), TableOpenMode::ReadWrite).unwrap();

        let duration_ms = corr_ctx.metafits_context.sched_duration_ms;
        let int_time_ms = corr_ctx.metafits_context.corr_int_time_ms;

        // interval is from start of first scan to end of last scan.
        let source_interval = (duration_ms + int_time_ms) as f64 / 1000.0;
        let source_time = sched_start_time_mjd_utc_s + source_interval / 2.;

        source_table.add_rows(1).unwrap();
        self.write_source_row(
            &mut source_table,
            0,
            0,
            source_time,
            source_interval,
            0,
            0,
            field_name,
            0,
            "",
            vec![ra_phase_rad, dec_phase_rad],
            vec![0., 0.],
        )
        .unwrap();

        // /////////// //
        // Observation //
        // /////////// //

        let mut obs_table =
            Table::open(&self.path.join("OBSERVATION"), TableOpenMode::ReadWrite).unwrap();
        obs_table.add_rows(1).unwrap();

        let start_time_centroid_mjd_utc_s = Epoch::from_gpst_seconds(
            (corr_ctx.metafits_context.sched_start_gps_time_ms + int_time_ms / 2) as f64 / 1e3,
        )
        .as_mjd_utc_seconds();
        let end_time_centroid_mjd_utc_s = Epoch::from_gpst_seconds(
            (corr_ctx.metafits_context.sched_end_gps_time_ms + int_time_ms / 2) as f64 / 1e3,
        )
        .as_mjd_utc_seconds();

        self.write_observation_row_mwa(
            &mut obs_table,
            0,
            "MWA",
            (start_time_centroid_mjd_utc_s, end_time_centroid_mjd_utc_s),
            &corr_ctx.metafits_context.creator,
            "MWA",
            &corr_ctx.metafits_context.project_id,
            0.,
            corr_ctx.metafits_context.obs_id as _,
            &corr_ctx.metafits_context.obs_name,
            &corr_ctx.metafits_context.mode.to_string(),
            (timestep_range.len() + 1) as _,
            sched_start_time_mjd_utc_s,
            false,
        )
        .unwrap();

        // /////// //
        // History //
        // /////// //

        let mut hist_table =
            Table::open(&self.path.join("HISTORY"), TableOpenMode::ReadWrite).unwrap();

        hist_table.add_rows(1).unwrap();

        let time = SystemTime::now()
            .duration_since(SystemTime::UNIX_EPOCH)
            .unwrap()
            .as_millis() as f64
            / 1000.;
        // TODO:
        let cmd_line = "TODO";
        let message = "TODO";
        let params = "TODO";
        self.write_history_row(
            &mut hist_table,
            0,
            time,
            cmd_line,
            message,
            &format!("Birli {}", PKG_VERSION),
            params,
        )
        .unwrap();

        // //// //
        // Feed //
        // //// //

        let mut feed_table =
            Table::open(&self.path.join("FEED"), TableOpenMode::ReadWrite).unwrap();

        // all of this assumes num_pols = 2
        assert_eq!(num_ant_pols, 2);
        feed_table.add_rows(num_ants).unwrap();

        for idx in 0..num_ants {
            self.write_feed_row(
                &mut feed_table,
                idx as _,
                idx as _,
                0,
                -1,
                source_time,
                source_interval,
                num_ant_pols as _,
                -1,
                &array![[0., 0.], [0., 0.]],
                &vec!["X".into(), "Y".into()],
                &array![
                    [c32::new(1., 0.), c32::new(0., 0.)],
                    [c32::new(0., 0.), c32::new(1., 0.)]
                ],
                &vec![0., 0., 0.],
                &vec![0., FRAC_PI_2],
            )
            .unwrap();
        }

        // ///////////////// //
        // MWA Tile Pointing //
        // ///////////////// //

        let mut point_table = Table::open(
            &self.path.join("MWA_TILE_POINTING"),
            TableOpenMode::ReadWrite,
        )
        .unwrap();

        point_table.add_rows(1).unwrap();

        let delays = corr_ctx
            .metafits_context
            .delays
            .iter()
            .map(|&x| x as _)
            .collect();

        self.write_mwa_tile_pointing_row(
            &mut point_table,
            0,
            start_time_centroid_mjd_utc_s,
            end_time_centroid_mjd_utc_s,
            &delays,
            ra_phase_rad,
            dec_phase_rad,
        )
        .unwrap();

        // /////////// //
        // MWA Subband //
        // /////////// //

        let mut subband_table =
            Table::open(&self.path.join("MWA_SUBBAND"), TableOpenMode::ReadWrite).unwrap();

        subband_table.add_rows(num_sel_coarse_chans).unwrap();

        for i in 0..num_sel_coarse_chans {
            self.write_mwa_subband_row(&mut subband_table, i as _, i as _, 0 as _, false)
                .unwrap();
        }

        Ok(())
    }

    /// Write a row into the main table.
    ///
    /// The main table holds measurements from a Telescope
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `idx` - row index to write to (ensure enough rows have been added)
    /// - `time` - Modified Julian Day, at midpoint of scan
    /// - `time_centroid` - Modified Julian Day, at centroid of scan
    /// - `antenna1` - ID of first antenna in interferometer
    /// - `antenna2` - ID of second antenna in interferometer
    /// - `data_desc_id` - The data description table index
    /// - `uvw` - Vector with uvw coordinates (in meters)
    /// - `interval` - The sampling interval
    /// - `processor_id` - Id for backend processor, index in PROCESSOR table
    /// - `scan_number` - Sequential scan number from on-line system
    /// - `state_id` - ID for this observing state
    /// - `sigma` - Estimated rms noise for channel with unity bandpass response
    /// - `data` - an `[n, p]` shaped ndarray of complex visibilities, where `n`
    ///     is the number of channels, and p is the number of polarizations
    /// - `flags` - an `[n, p]` shaped ndarray of boolean flags.
    /// - `weights` - a `[p]` shaped ndarray of weights for each polarization
    ///
    /// # Gorey details
    ///
    /// According to <https://casa.nrao.edu/Memos/229.html#SECTION00061000000000000000>,
    /// midpoint and centroid time are different things? It doesn't explain how, and
    /// Cotter doesn't treat them differently either.
    #[allow(clippy::ptr_arg)]
    #[allow(clippy::too_many_arguments)]
    pub fn write_main_row(
        &self,
        table: &mut Table,
        idx: u64,
        time: f64,
        time_centroid: f64,
        antenna1: i32,
        antenna2: i32,
        data_desc_id: i32,
        // TODO: take UVW
        uvw: &Vec<f64>,
        interval: f64,
        // TODO: is this not just interval?
        // exposure: f64,
        processor_id: i32,
        scan_number: i32,
        state_id: i32,
        sigma: &Vec<f32>,
        data: &Array2<c32>,
        flags: &Array2<bool>,
        weights: &Array2<f32>,
        flag_row: bool,
    ) -> Result<(), IOError> {
        let num_pols = 4;

        if uvw.len() != 3 {
            return Err(IOError::BadArrayShape {
                argument: "uvw".into(),
                function: "write_main_row".into(),
                expected: "3".into(),
                received: format!("{:?}", uvw.len()),
            });
        }

        if sigma.len() != num_pols {
            return Err(IOError::BadArrayShape {
                argument: "sigma".into(),
                function: "write_main_row".into(),
                expected: format!("{}", num_pols),
                received: format!("{:?}", sigma.len()),
            });
        }

        match (data.shape(), flags.shape(), weights.shape()) {
            ([d0, d1], [f0, f1], [w0, w1])
                if d0 == f0
                    && f0 == w0
                    && d1 == &num_pols
                    && f1 == &num_pols
                    && w1 == &num_pols => {}
            (dsh, fsh, wsh) => {
                return Err(IOError::BadArrayShape {
                    argument: "data|flags|weights".into(),
                    function: "write_main_row".into(),
                    expected: format!(
                        "[n, p]|[n, p]|[n, p] where n=num_chans, p=num_pols({})",
                        num_pols
                    ),
                    received: format!("{:?}|{:?}|{:?}", dsh, fsh, wsh),
                })
            }
        }

        let weight_pol = weights
            .axis_iter(Axis(1))
            .map(|weights_pol_view| weights_pol_view.sum())
            .collect::<Vec<f32>>();

        table.put_cell("TIME", idx, &time).unwrap();
        table
            .put_cell("TIME_CENTROID", idx, &time_centroid)
            .unwrap();
        table.put_cell("ANTENNA1", idx, &antenna1).unwrap();
        table.put_cell("ANTENNA2", idx, &antenna2).unwrap();
        table.put_cell("DATA_DESC_ID", idx, &data_desc_id).unwrap();
        table.put_cell("UVW", idx, uvw).unwrap();
        table.put_cell("INTERVAL", idx, &interval).unwrap();
        // TODO: really?
        table.put_cell("EXPOSURE", idx, &interval).unwrap();
        table.put_cell("PROCESSOR_ID", idx, &processor_id).unwrap();
        table.put_cell("SCAN_NUMBER", idx, &scan_number).unwrap();
        table.put_cell("STATE_ID", idx, &state_id).unwrap();
        table.put_cell("SIGMA", idx, sigma).unwrap();
        table.put_cell("DATA", idx, data).unwrap();
        table.put_cell("WEIGHT_SPECTRUM", idx, weights).unwrap();
        table.put_cell("WEIGHT", idx, &weight_pol).unwrap();
        table.put_cell("FLAG", idx, flags).unwrap();
        table.put_cell("FLAG_ROW", idx, &flag_row).unwrap();

        Ok(())
    }
}

impl VisWritable for MeasurementSetWriter {
    fn write_vis_marlu(
        &mut self,
        vis: ArrayView3<Jones<f32>>,
        weights: ArrayView3<f32>,
        flags: ArrayView3<bool>,
        vis_ctx: &VisContext,
        draw_progress: bool,
    ) -> Result<(), IOError> {
        let num_avg_timesteps = vis_ctx.num_avg_timesteps();
        let num_avg_chans = vis_ctx.num_avg_chans();
        let num_vis_pols = vis_ctx.num_vis_pols;
        let num_avg_rows = num_avg_timesteps * vis_ctx.sel_baselines.len();

        // Progress bars
        let draw_target = if draw_progress {
            ProgressDrawTarget::stderr()
        } else {
            ProgressDrawTarget::hidden()
        };
        let write_progress =
            indicatif::ProgressBar::with_draw_target(num_avg_rows as u64, draw_target);
        write_progress.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
                )
                .progress_chars("=> "),
        );
        write_progress.set_message("write ms vis");

        // Open the table for writing
        let mut main_table = Table::open(&self.path, TableOpenMode::ReadWrite).unwrap();
        let num_main_rows = main_table.n_rows();
        trace!(
            "num_main_rows={}, self.main_row_idx={}, num_avg_rows (selected)={}",
            num_main_rows,
            self.main_row_idx,
            num_avg_rows
        );
        assert!(num_main_rows - self.main_row_idx as u64 >= num_avg_rows as u64);

        let mut uvw_tmp = Vec::with_capacity(3);
        let sigma_tmp = vec![1., 1., 1., 1.];
        let mut data_tmp = Array2::zeros((num_avg_chans, num_vis_pols));
        let mut weights_tmp = Array2::zeros((num_avg_chans, num_vis_pols));
        let mut flags_tmp = Array2::from_elem((num_avg_chans, num_vis_pols), false);

        for (avg_centroid_timestamp, vis_chunk, weight_chunk, flag_chunk) in izip!(
            vis_ctx.timeseries(true, true),
            vis.axis_chunks_iter(Axis(0), vis_ctx.avg_time),
            weights.axis_chunks_iter(Axis(0), vis_ctx.avg_time),
            flags.axis_chunks_iter(Axis(0), vis_ctx.avg_time),
        ) {
            let scan_centroid_mjd_utc_s = avg_centroid_timestamp.as_mjd_utc_seconds();

            let prec_info = precess_time(
                self.phase_centre,
                avg_centroid_timestamp,
                self.array_pos.longitude_rad,
                self.array_pos.latitude_rad,
            );

            let tiles_xyz_precessed = prec_info.precess_xyz_parallel(&vis_ctx.tiles_xyz_geod);

            for ((ant1_idx, ant2_idx), vis_chunk, weight_chunk, flag_chunk) in izip!(
                vis_ctx.sel_baselines.clone().into_iter(),
                vis_chunk.axis_iter(Axis(2)),
                weight_chunk.axis_iter(Axis(2)),
                flag_chunk.axis_iter(Axis(2)),
            ) {
                let baseline_xyz_precessed =
                    tiles_xyz_precessed[ant1_idx] - tiles_xyz_precessed[ant2_idx];
                let uvw = UVW::from_xyz(baseline_xyz_precessed, prec_info.hadec_j2000);

                // copy values into temporary arrays to avoid heap allocs.
                uvw_tmp.clear();
                uvw_tmp.extend_from_slice(&[uvw.u, uvw.v, uvw.w]);

                data_tmp.fill(Complex::default());
                weights_tmp.fill(0.);
                flags_tmp.fill(false);

                // iterate through the channel dimension of the arrays in chunks of size `avg_freq`,
                // averaging the chunks into the tmp arrays.
                for (
                    vis_chunk,
                    weight_chunk,
                    flag_chunk,
                    mut data_tmp_view,
                    mut weights_tmp_view,
                    mut flags_tmp_view,
                ) in izip!(
                    vis_chunk.axis_chunks_iter(Axis(1), vis_ctx.avg_freq),
                    weight_chunk.axis_chunks_iter(Axis(1), vis_ctx.avg_freq),
                    flag_chunk.axis_chunks_iter(Axis(1), vis_ctx.avg_freq),
                    data_tmp.outer_iter_mut(),
                    weights_tmp.outer_iter_mut(),
                    flags_tmp.outer_iter_mut()
                ) {
                    let mut avg_weight: f32 = weight_chunk[[0, 0]];
                    let mut avg_flag: bool = flag_chunk[[0, 0]];
                    if vis_ctx.trivial_averaging() {
                        data_tmp_view.assign(&ArrayView::from(vis_chunk[[0, 0]].as_slice()));
                    } else {
                        average_chunk_f64!(
                            vis_chunk,
                            weight_chunk,
                            flag_chunk,
                            data_tmp_view,
                            avg_weight,
                            avg_flag
                        );
                    }
                    weights_tmp_view.fill(avg_weight);
                    flags_tmp_view.fill(avg_flag);
                }

                let flag_row = flags_tmp.iter().all(|&x| x);
                self.write_main_row(
                    &mut main_table,
                    self.main_row_idx as _,
                    scan_centroid_mjd_utc_s,
                    scan_centroid_mjd_utc_s,
                    ant1_idx as _,
                    ant2_idx as _,
                    0,
                    &uvw_tmp,
                    vis_ctx.avg_int_time().in_seconds(),
                    -1,
                    1,
                    -1,
                    &sigma_tmp,
                    &data_tmp,
                    &flags_tmp,
                    &weights_tmp,
                    flag_row,
                )?;

                self.main_row_idx += 1;

                write_progress.inc(1);
            }
        }
        write_progress.finish();
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::{
        collections::{BTreeMap, HashSet},
        f64::consts::FRAC_PI_2,
        path::PathBuf,
    };

    use super::*;

    use approx::abs_diff_eq;
    use itertools::izip;
    use lexical::parse;
    use regex::Regex;
    use serial_test::serial;
    use tempfile::tempdir;

    use crate::{
        c64,
        ndarray::{s, Array, Array4},
    };

    cfg_if::cfg_if! {
        if #[cfg(feature = "mwalib")] {
            use crate::{
                c32,
                constants::{
                    COTTER_MWA_HEIGHT_METRES, COTTER_MWA_LATITUDE_RADIANS, COTTER_MWA_LONGITUDE_RADIANS,
                },
                ndarray::array,
            };
        }
    }

    lazy_static! {
        static ref PATH_1254670392: PathBuf =
            "tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms".into();
        static ref PATH_1254670392_AVG_4S_80KHZ: PathBuf =
            "tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_80khz.trunc.ms".into();
    }

    macro_rules! assert_table_column_names_match {
        ( $left:expr, $right:expr ) => {
            match (&$left.column_names(), &$right.column_names()) {
                (Ok(left_columns), Ok(right_columns)) => {
                    assert_eq!(left_columns, right_columns, "column names do not match");
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
                "row counts do not match. {} != {} for tables {:?} and {:?}",
                $left.n_rows(),
                $right.n_rows(),
                &$left,
                &$right
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
            assert_table_column_values_match_approx!($left, $right, $col_name, $col_desc, $type_, 1e-7)
        };
        ( $left:expr, $right:expr, $col_name:expr, $col_desc:expr, $type_:ty, $epsilon:expr ) => {
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
                                abs_diff_eq!(left_cell, right_cell, epsilon = $epsilon),
                                "cells don't match in column {}, row {}. {:?} != {:?}",
                                $col_name,
                                row_idx,
                                left_cell,
                                right_cell
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
                                    "cells don't match at index {} in column {}, row {}. {:?} != {:?}",
                                    vec_idx,
                                    $col_name,
                                    row_idx,
                                    left_cell,
                                    right_cell
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
            assert_table_columns_match!($left, $right, $col_name, 1e-7)
        };
        ( $left:expr, $right:expr, $col_name:expr, $epsilon:expr ) => {
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
                    $left,
                    $right,
                    $col_name,
                    col_desc,
                    f32,
                    $epsilon as f32
                ),
                GlueDataType::TpDouble => assert_table_column_values_match_approx!(
                    $left,
                    $right,
                    $col_name,
                    col_desc,
                    f64,
                    $epsilon as f64
                ),
                GlueDataType::TpComplex => assert_table_column_values_match_approx!(
                    $left,
                    $right,
                    $col_name,
                    col_desc,
                    c32,
                    $epsilon as f32
                ),
                GlueDataType::TpDComplex => assert_table_column_values_match_approx!(
                    $left,
                    $right,
                    $col_name,
                    col_desc,
                    c64,
                    $epsilon as f64
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
    #[serial]
    fn test_decompress_default_tables() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        drop(ms_writer);

        assert!(table_path.exists());

        let mut main_table = Table::open(&table_path, TableOpenMode::Read).unwrap();
        let main_table_keywords = main_table.table_keyword_names().unwrap();
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
            if !table_name.is_empty() {
                assert!(main_table_keywords.contains(&table_name.into()));
            }
        }
    }

    #[test]
    #[serial]
    fn test_add_source_table() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
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
    #[serial]
    fn test_add_cotter_mods() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
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

        let mut main_table = Table::open(&table_path, TableOpenMode::Read).unwrap();
        let main_table_keywords = main_table.table_keyword_names().unwrap();
        assert!(main_table_keywords.contains(&"SOURCE".into()));
    }

    #[test]
    #[serial]
    fn test_add_mwa_mods() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.add_mwa_mods();
        drop(ms_writer);

        for (table_name, col_names) in [
            (
                "ANTENNA",
                vec![
                    "MWA_INPUT",
                    "MWA_TILE_NR",
                    "MWA_RECEIVER",
                    "MWA_SLOT",
                    "MWA_CABLE_LENGTH",
                ],
            ),
            ("FIELD", vec!["MWA_HAS_CALIBRATOR"]),
            (
                "OBSERVATION",
                vec![
                    "MWA_GPS_TIME",
                    "MWA_FILENAME",
                    "MWA_OBSERVATION_MODE",
                    "MWA_FLAG_WINDOW_SIZE",
                    "MWA_DATE_REQUESTED",
                ],
            ),
            ("SPECTRAL_WINDOW", vec!["MWA_CENTRE_SUBBAND_NR"]),
            ("MWA_TILE_POINTING", vec!["INTERVAL", "DELAYS", "DIRECTION"]),
            ("MWA_SUBBAND", vec!["NUMBER", "GAIN", "FLAG_ROW"]),
        ] {
            let mut table = Table::open(&table_path.join(table_name), TableOpenMode::Read).unwrap();
            let mut exp_table =
                Table::open(PATH_1254670392.join(table_name), TableOpenMode::Read).unwrap();
            for col_name in col_names {
                assert_table_column_descriptions_match!(table, exp_table, col_name);
            }
        }

        let mut main_table = Table::open(&table_path, TableOpenMode::Read).unwrap();
        let main_table_keywords = main_table.table_keyword_names().unwrap();
        assert!(main_table_keywords.contains(&"MWA_TILE_POINTING".into()));
        assert!(main_table_keywords.contains(&"MWA_SUBBAND".into()));
    }

    #[test]
    #[serial]
    fn test_write_spectral_window_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
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
        let spw_table_path = table_path.join("SPECTRAL_WINDOW");
        let mut spw_table = Table::open(&spw_table_path, TableOpenMode::ReadWrite).unwrap();

        spw_table.add_rows(1).unwrap();

        ms_writer
            .write_spectral_window_row(
                &mut spw_table,
                0,
                "MWA_BAND_182.4",
                182395000.,
                chan_info,
                30720000.,
                false,
            )
            .unwrap();
        drop(ms_writer);

        let mut spw_table = Table::open(&spw_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("SPECTRAL_WINDOW"), TableOpenMode::Read).unwrap();

        assert_table_nrows_match!(spw_table, expected_table);
        for col_name in [
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
        ] {
            assert_table_columns_match!(spw_table, expected_table, col_name);
        }
    }

    #[test]
    #[serial]
    fn test_write_spectral_window_row_mwa() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);
        ms_writer.add_mwa_mods();

        let chan_info = Array2::from_shape_fn((768, 4), |(c, i)| {
            (if i == 0 {
                167055000. + (c as f64) * 40000.
            } else {
                40000.
            }) as f64
        });
        let spw_table_path = table_path.join("SPECTRAL_WINDOW");
        let mut spw_table = Table::open(&spw_table_path, TableOpenMode::ReadWrite).unwrap();

        spw_table.add_rows(1).unwrap();

        ms_writer
            .write_spectral_window_row_mwa(
                &mut spw_table,
                0,
                "MWA_BAND_182.4",
                182395000.,
                chan_info,
                30720000.,
                143,
                false,
            )
            .unwrap();
        drop(ms_writer);

        let mut spw_table = Table::open(&spw_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("SPECTRAL_WINDOW"), TableOpenMode::Read).unwrap();

        assert_tables_match!(spw_table, expected_table);
    }

    #[test]
    #[serial]
    fn handle_bad_spw_chan_info() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let chan_info = Array2::from_shape_fn((768, 3), |(_, _)| 40000.);

        let spw_table_path = table_path.join("SPECTRAL_WINDOW");
        let mut spw_table = Table::open(&spw_table_path, TableOpenMode::ReadWrite).unwrap();

        spw_table.add_rows(1).unwrap();

        let result = ms_writer.write_spectral_window_row(
            &mut spw_table,
            0,
            "MWA_BAND_182.4",
            182395000.,
            chan_info,
            30720000.,
            false,
        );

        assert!(matches!(result, Err(IOError::BadArrayShape { .. })))
    }

    #[test]
    #[serial]
    fn test_write_data_description_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let ddesc_table_path = table_path.join("DATA_DESCRIPTION");
        let mut ddesc_table = Table::open(&ddesc_table_path, TableOpenMode::ReadWrite).unwrap();

        ddesc_table.add_rows(1).unwrap();

        ms_writer
            .write_data_description_row(&mut ddesc_table, 0, 0, 0, false)
            .unwrap();
        drop(ms_writer);

        let mut ddesc_table = Table::open(&ddesc_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table = Table::open(
            PATH_1254670392.join("DATA_DESCRIPTION"),
            TableOpenMode::Read,
        )
        .unwrap();

        assert_tables_match!(ddesc_table, expected_table);
    }

    const ANT_POSITIONS: &[[f64; 3]] = &[
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
        [-2559649.92256756, 5095691.87122651, -2849149.12522308],
    ];

    const ANT_NAMES: &[&str] = &[
        "Tile011", "Tile012", "Tile013", "Tile014", "Tile015", "Tile016", "Tile017", "Tile018",
        "Tile021", "Tile022", "Tile023", "Tile024", "Tile025", "Tile026", "Tile027", "Tile028",
        "Tile031", "Tile032", "Tile033", "Tile034", "Tile035", "Tile036", "Tile037", "Tile038",
        "Tile041", "Tile042", "Tile043", "Tile044", "Tile045", "Tile046", "Tile047", "Tile048",
        "Tile061", "Tile062", "Tile063", "Tile064", "Tile065", "Tile066", "Tile067", "Tile068",
        "Tile081", "Tile082", "Tile083", "Tile084", "Tile085", "Tile086", "Tile087", "Tile088",
        "Tile091", "Tile092", "Tile093", "Tile094", "Tile095", "Tile096", "Tile097", "Tile098",
        "HexE1", "HexE2", "HexE3", "HexE4", "HexE5", "HexE6", "HexE7", "HexE8", "HexE9", "HexE10",
        "HexE11", "HexE12", "HexE13", "HexE14", "HexE15", "HexE16", "HexE17", "HexE18", "HexE19",
        "HexE20", "HexE21", "HexE22", "HexE23", "HexE24", "HexE25", "HexE26", "HexE27", "HexE28",
        "HexE29", "HexE30", "HexE31", "HexE32", "HexE33", "HexE34", "HexE35", "HexE36", "HexS1",
        "HexS2", "HexS3", "HexS4", "HexS5", "HexS6", "HexS7", "HexS8", "HexS9", "HexS10", "HexS11",
        "HexS12", "HexS13", "HexS14", "HexS15", "HexS16", "HexS17", "HexS18", "HexS19", "HexS20",
        "HexS21", "HexS22", "HexS23", "HexS24", "HexS25", "HexS26", "HexS27", "HexS28", "HexS29",
        "HexS30", "HexS31", "HexS32", "HexS33", "HexS34", "HexS35", "HexS36",
    ];

    const ANT_INPUTS: &[[i32; 2]] = &[
        [87, 86],
        [85, 84],
        [83, 82],
        [81, 80],
        [95, 94],
        [93, 92],
        [91, 90],
        [89, 88],
        [167, 166],
        [165, 164],
        [163, 162],
        [161, 160],
        [175, 174],
        [173, 172],
        [139, 138],
        [169, 168],
        [71, 70],
        [69, 68],
        [67, 66],
        [65, 64],
        [79, 78],
        [77, 76],
        [75, 74],
        [73, 72],
        [135, 134],
        [133, 132],
        [131, 130],
        [129, 128],
        [143, 142],
        [141, 140],
        [171, 170],
        [137, 136],
        [151, 150],
        [149, 148],
        [147, 146],
        [145, 144],
        [159, 158],
        [157, 156],
        [155, 154],
        [153, 152],
        [119, 118],
        [117, 116],
        [115, 114],
        [113, 112],
        [127, 126],
        [125, 124],
        [123, 122],
        [121, 120],
        [103, 102],
        [101, 100],
        [99, 98],
        [97, 96],
        [111, 110],
        [109, 108],
        [107, 106],
        [105, 104],
        [183, 182],
        [181, 180],
        [179, 178],
        [177, 176],
        [191, 190],
        [189, 188],
        [187, 186],
        [185, 184],
        [23, 22],
        [21, 20],
        [19, 18],
        [17, 16],
        [31, 30],
        [29, 28],
        [27, 26],
        [25, 24],
        [199, 198],
        [197, 196],
        [195, 194],
        [193, 192],
        [207, 206],
        [205, 204],
        [203, 202],
        [201, 200],
        [215, 214],
        [213, 212],
        [211, 210],
        [209, 208],
        [223, 222],
        [221, 220],
        [219, 218],
        [217, 216],
        [39, 38],
        [37, 36],
        [35, 34],
        [33, 32],
        [47, 46],
        [45, 44],
        [43, 42],
        [41, 40],
        [55, 54],
        [53, 52],
        [51, 50],
        [49, 48],
        [63, 62],
        [61, 60],
        [59, 58],
        [57, 56],
        [231, 230],
        [229, 228],
        [227, 226],
        [225, 224],
        [239, 238],
        [237, 236],
        [235, 234],
        [233, 232],
        [247, 246],
        [245, 244],
        [243, 242],
        [241, 240],
        [255, 254],
        [253, 252],
        [251, 250],
        [249, 248],
        [7, 6],
        [5, 4],
        [3, 2],
        [1, 0],
        [15, 14],
        [13, 12],
        [11, 10],
        [9, 8],
    ];

    const ANT_TILE_NRS: &[i32] = &[
        11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 31, 32, 33, 34, 35, 36, 37,
        38, 41, 42, 43, 44, 45, 46, 47, 48, 61, 62, 63, 64, 65, 66, 67, 68, 81, 82, 83, 84, 85, 86,
        87, 88, 91, 92, 93, 94, 95, 96, 97, 98, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008,
        1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023,
        1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038,
        1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053,
        1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068,
        1069, 1070, 1071, 1072,
    ];

    const ANT_CABLE_LENGTHS: &[[f64; 2]] = &[
        [-656.14, -656.14],
        [-655.66, -655.66],
        [-582.91, -582.91],
        [-583., -583.],
        [-583.02, -583.02],
        [-655.79, -655.79],
        [-583.32, -583.32],
        [-655.39, -655.39],
        [-485.73, -485.73],
        [-487.73, -487.73],
        [-485.73, -485.73],
        [-485.73, -485.73],
        [-488.23, -488.23],
        [-486.53, -486.53],
        [-582.7, -582.7],
        [-485.73, -485.73],
        [-499.13, -499.13],
        [-498.93, -498.93],
        [-595.63, -595.63],
        [-596.38, -596.38],
        [-596.33, -596.33],
        [-595.63, -595.63],
        [-598.63, -598.63],
        [-595.63, -595.63],
        [-582.7, -582.7],
        [-582.7, -582.7],
        [-582.7, -582.7],
        [-582.7, -582.7],
        [-582.7, -582.7],
        [-485.7, -485.7],
        [-387.53, -387.53],
        [-483.9, -483.9],
        [-580.62, -580.62],
        [-580.62, -580.62],
        [-580.62, -580.62],
        [-579.92, -579.92],
        [-484.62, -484.62],
        [-580.62, -580.62],
        [-580.62, -580.62],
        [-580.62, -580.62],
        [-482.93, -482.93],
        [-578.99, -578.99],
        [-483.04, -483.04],
        [-481.13, -481.13],
        [-482.99, -482.99],
        [-482.98, -482.98],
        [-482.33, -482.33],
        [-651.32, -651.32],
        [-398.6, -398.6],
        [-497.18, -497.18],
        [-595.63, -595.63],
        [-497.08, -497.08],
        [-497.36, -497.36],
        [-497.59, -497.59],
        [-497.7, -497.7],
        [-497.4, -497.4],
        [-655.85, -655.85],
        [-655.85, -655.85],
        [-654.85, -654.85],
        [-582.85, -582.85],
        [-582.85, -582.85],
        [-582.85, -582.85],
        [-582.85, -582.85],
        [-582.85, -582.85],
        [-216.43, -216.43],
        [-289.43, -289.43],
        [-216.43, -216.43],
        [-216.43, -216.43],
        [-216.43, -216.43],
        [-216.43, -216.43],
        [-216.43, -216.43],
        [-289.43, -289.43],
        [-291.22, -291.22],
        [-218.22, -218.22],
        [-218.22, -218.22],
        [-221.72, -221.72],
        [-121.22, -121.22],
        [-291.22, -291.22],
        [-218.22, -218.22],
        [-218.22, -218.22],
        [-214.06, -214.06],
        [-215.46, -215.46],
        [-214.06, -214.06],
        [-290.46, -290.46],
        [-214.06, -214.06],
        [-213.96, -213.96],
        [-219.06, -219.06],
        [-218.06, -218.06],
        [-220.85, -220.85],
        [-219.85, -219.85],
        [-219.85, -219.85],
        [-223.85, -223.85],
        [-292.85, -292.85],
        [-293.85, -293.85],
        [-292.35, -292.35],
        [-219.85, -219.85],
        [-211.98, -211.98],
        [-211.98, -211.98],
        [-211.98, -211.98],
        [-211.98, -211.98],
        [-211.23, -211.23],
        [-284.18, -284.18],
        [-211.98, -211.98],
        [-211.98, -211.98],
        [-216.87, -216.87],
        [-210.87, -210.87],
        [-114.72, -114.72],
        [-285.87, -285.87],
        [-285.87, -285.87],
        [-212.12, -212.12],
        [-211.47, -211.47],
        [-212.87, -212.87],
        [-117.84, -117.84],
        [-288.84, -288.84],
        [-215.84, -215.84],
        [-215.84, -215.84],
        [-215.84, -215.84],
        [-215.09, -215.09],
        [-115.84, -115.84],
        [-290.14, -290.14],
        [-217.53, -217.53],
        [-211.53, -211.53],
        [-210.53, -210.53],
        [-211.53, -211.53],
        [-211.53, -211.53],
        [-211.53, -211.53],
        [-211.53, -211.53],
        [-211.53, -211.53],
    ];

    /// Test data:
    /// ```python
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms/ANTENNA')
    /// tb.getcol("POSITION").transpose()
    /// tb.getcol("NAME")
    /// ```
    #[test]
    #[serial]
    fn test_write_antenna_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let ant_table_path = ms_writer.path.join("ANTENNA");
        let mut ant_table = Table::open(ant_table_path, TableOpenMode::ReadWrite).unwrap();

        ant_table.add_rows(ANT_NAMES.len()).unwrap();

        for (idx, (name, position)) in izip!(ANT_NAMES, ANT_POSITIONS).enumerate() {
            let position = position.to_vec();

            ms_writer
                .write_antenna_row(
                    &mut ant_table,
                    idx as _,
                    name,
                    "MWA",
                    "GROUND-BASED",
                    "ALT-AZ",
                    &position,
                    4.0,
                    false,
                )
                .unwrap();
        }

        drop(ms_writer);

        let mut expected_table =
            Table::open(PATH_1254670392.join("ANTENNA"), TableOpenMode::Read).unwrap();

        assert_table_nrows_match!(ant_table, expected_table);
        for col_name in [
            "OFFSET",
            "POSITION",
            "TYPE",
            "DISH_DIAMETER",
            "FLAG_ROW",
            "MOUNT",
            "NAME",
            "STATION",
        ] {
            assert_table_columns_match!(ant_table, expected_table, col_name);
        }
    }

    #[test]
    #[serial]
    fn test_write_antenna_row_mwa() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);
        ms_writer.add_mwa_mods();

        let ant_table_path = ms_writer.path.join("ANTENNA");
        let mut ant_table = Table::open(&ant_table_path, TableOpenMode::ReadWrite).unwrap();

        ant_table.add_rows(ANT_NAMES.len()).unwrap();

        for (idx, (name, position, input, tile_nr, cable_length)) in izip!(
            ANT_NAMES,
            ANT_POSITIONS,
            ANT_INPUTS,
            ANT_TILE_NRS,
            ANT_CABLE_LENGTHS
        )
        .enumerate()
        {
            let position = position.to_vec();

            ms_writer
                .write_antenna_row_mwa(
                    &mut ant_table,
                    idx as _,
                    name,
                    "MWA",
                    "GROUND-BASED",
                    "ALT-AZ",
                    &position,
                    4.0,
                    &input.to_vec(),
                    *tile_nr,
                    0,
                    &vec![0, 0],
                    &cable_length.to_vec(),
                    false,
                )
                .unwrap();
        }

        drop(ms_writer);

        let mut ant_table = Table::open(&ant_table_path, TableOpenMode::Read).unwrap();
        let mut expected_table =
            Table::open(PATH_1254670392.join("ANTENNA"), TableOpenMode::Read).unwrap();

        assert_tables_match!(ant_table, expected_table);
    }

    #[test]
    #[serial]
    fn test_write_polarization_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let corr_product = array![[0, 0], [0, 1], [1, 0], [1, 1]];
        let corr_type = vec![9, 10, 11, 12];

        let pol_table_path = table_path.join("POLARIZATION");
        let mut pol_table = Table::open(&pol_table_path, TableOpenMode::ReadWrite).unwrap();

        pol_table.add_rows(1).unwrap();

        ms_writer
            .write_polarization_row(&mut pol_table, 0, &corr_type, &corr_product, false)
            .unwrap();
        drop(ms_writer);

        let mut pol_table = Table::open(&pol_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("POLARIZATION"), TableOpenMode::Read).unwrap();

        assert_tables_match!(pol_table, expected_table);
    }

    #[test]
    #[serial]
    fn handle_bad_pol_small_corr_type() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let corr_product = array![[0, 0], [0, 1], [1, 0], [1, 1]];
        let corr_type = vec![9, 10, 11];

        let pol_table_path = table_path.join("POLARIZATION");
        let mut pol_table = Table::open(&pol_table_path, TableOpenMode::ReadWrite).unwrap();

        pol_table.add_rows(1).unwrap();

        let result =
            ms_writer.write_polarization_row(&mut pol_table, 0, &corr_type, &corr_product, false);

        assert!(matches!(result, Err(IOError::BadArrayShape { .. })))
    }

    #[test]
    #[serial]
    fn handle_bad_pol_big_corr_product() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let corr_product = array![[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]];
        let corr_type = vec![9, 10, 11, 12];

        let pol_table_path = table_path.join("POLARIZATION");
        let mut pol_table = Table::open(&pol_table_path, TableOpenMode::ReadWrite).unwrap();

        pol_table.add_rows(1).unwrap();

        let result =
            ms_writer.write_polarization_row(&mut pol_table, 0, &corr_type, &corr_product, false);

        assert!(matches!(result, Err(IOError::BadArrayShape { .. })))
    }

    #[test]
    #[serial]
    fn test_write_source_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let source_table_path = table_path.join("SOURCE");
        let mut source_table = Table::open(&source_table_path, TableOpenMode::ReadWrite).unwrap();

        source_table.add_rows(1).unwrap();

        ms_writer
            .write_source_row(
                &mut source_table,
                0,
                0,
                5077351979.000001,
                0.00009259259240934625,
                0,
                0,
                "high_2019B_2458765_EOR0_RADec0.0,-27.0",
                0,
                "",
                vec![phase_centre.ra, phase_centre.dec],
                vec![0., 0.],
            )
            .unwrap();
        drop(ms_writer);

        let mut source_table = Table::open(&source_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("SOURCE"), TableOpenMode::Read).unwrap();

        assert_table_nrows_match!(source_table, expected_table);
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
            assert_table_columns_match!(source_table, expected_table, col_name);
        }
    }

    #[test]
    #[serial]
    fn test_write_field_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let field_table_path = table_path.join("FIELD");
        let mut field_table = Table::open(&field_table_path, TableOpenMode::ReadWrite).unwrap();

        field_table.add_rows(1).unwrap();

        let dir_info = array![
            [[phase_centre.ra, phase_centre.dec]],
            [[phase_centre.ra, phase_centre.dec]],
            [[phase_centre.ra, phase_centre.dec]]
        ];
        ms_writer
            .write_field_row(
                &mut field_table,
                0,
                "high_2019B_2458765_EOR0_RADec0.0,-27.0",
                "",
                5077351974.,
                &dir_info,
                -1,
                false,
            )
            .unwrap();
        drop(ms_writer);

        let mut field_table = Table::open(&field_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("FIELD"), TableOpenMode::Read).unwrap();

        assert_table_nrows_match!(field_table, expected_table);
        for col_name in [
            "DELAY_DIR",
            "PHASE_DIR",
            "REFERENCE_DIR",
            "CODE",
            "FLAG_ROW",
            "NAME",
            "NUM_POLY",
            "SOURCE_ID",
            "TIME",
        ] {
            assert_table_columns_match!(field_table, expected_table, col_name);
        }
    }

    #[test]
    #[serial]
    fn test_write_field_row_mwa() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);
        ms_writer.add_mwa_mods();

        let field_table_path = table_path.join("FIELD");
        let mut field_table = Table::open(&field_table_path, TableOpenMode::ReadWrite).unwrap();

        field_table.add_rows(1).unwrap();

        let dir_info = array![
            [[phase_centre.ra, phase_centre.dec]],
            [[phase_centre.ra, phase_centre.dec]],
            [[phase_centre.ra, phase_centre.dec]]
        ];
        ms_writer
            .write_field_row_mwa(
                &mut field_table,
                0,
                "high_2019B_2458765_EOR0_RADec0.0,-27.0",
                "",
                5077351974.,
                &dir_info,
                -1,
                false,
                false,
            )
            .unwrap();
        drop(ms_writer);

        let mut field_table = Table::open(&field_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("FIELD"), TableOpenMode::Read).unwrap();

        assert_tables_match!(field_table, expected_table);
    }

    #[test]
    #[serial]
    fn handle_bad_field_shape() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let field_table_path = table_path.join("FIELD");
        let mut field_table = Table::open(&field_table_path, TableOpenMode::ReadWrite).unwrap();

        let dir_info = array![
            [[phase_centre.ra, phase_centre.dec]],
            [[phase_centre.ra, phase_centre.dec]],
            [[phase_centre.ra, phase_centre.dec]],
            [[0., 1.]]
        ];
        let result = ms_writer.write_field_row(
            &mut field_table,
            0,
            "high_2019B_2458765_EOR0_RADec0.0,-27.0",
            "",
            5077351974.,
            &dir_info,
            -1,
            false,
        );

        assert!(matches!(result, Err(IOError::BadArrayShape { .. })))
    }

    #[test]
    #[serial]
    fn test_write_observation_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let obs_table_path = table_path.join("OBSERVATION");
        let mut obs_table = Table::open(&obs_table_path, TableOpenMode::ReadWrite).unwrap();

        obs_table.add_rows(1).unwrap();

        ms_writer
            .write_observation_row(
                &mut obs_table,
                0,
                "MWA",
                (5077351975.0, 5077351983.0),
                "andrew",
                "MWA",
                "G0009",
                0.,
                false,
            )
            .unwrap();
        drop(ms_writer);

        let mut obs_table = Table::open(&obs_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("OBSERVATION"), TableOpenMode::Read).unwrap();

        assert_table_nrows_match!(obs_table, expected_table);
        for col_name in [
            "TIME_RANGE",
            // These are never written
            // "LOG",
            // "SCHEDULE",
            "FLAG_ROW",
            "OBSERVER",
            "PROJECT",
            "RELEASE_DATE",
            "SCHEDULE_TYPE",
            "TELESCOPE_NAME",
        ] {
            assert_table_columns_match!(obs_table, expected_table, col_name);
        }
    }

    #[test]
    #[serial]
    fn test_write_observation_row_mwa() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);
        ms_writer.add_mwa_mods();

        let obs_table_path = table_path.join("OBSERVATION");
        let mut obs_table = Table::open(&obs_table_path, TableOpenMode::ReadWrite).unwrap();

        obs_table.add_rows(1).unwrap();

        ms_writer
            .write_observation_row_mwa(
                &mut obs_table,
                0,
                "MWA",
                (5077351975.0, 5077351983.0),
                "andrew",
                "MWA",
                "G0009",
                0.,
                1254670392.,
                "high_2019B_2458765_EOR0_RADec0.0,-27.0_143",
                "HW_LFILES",
                4,
                5077351974.,
                false,
            )
            .unwrap();
        drop(ms_writer);

        let mut obs_table = Table::open(&obs_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("OBSERVATION"), TableOpenMode::Read).unwrap();

        assert_table_nrows_match!(obs_table, expected_table);
        for col_name in [
            "MWA_GPS_TIME",
            "MWA_FILENAME",
            "MWA_OBSERVATION_MODE",
            "MWA_FLAG_WINDOW_SIZE",
            "MWA_DATE_REQUESTED",
        ] {
            assert_table_columns_match!(obs_table, expected_table, col_name);
        }
    }

    #[test]
    #[serial]
    fn test_write_history_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);
        ms_writer.add_mwa_mods();

        let hist_table_path = table_path.join("HISTORY");
        let mut hist_table = Table::open(&hist_table_path, TableOpenMode::ReadWrite).unwrap();

        hist_table.add_rows(1).unwrap();

        ms_writer
            .write_history_row(
                &mut hist_table,
                0,
                5149221788.625,
                "cotter \"-m\" \"tests/data/1254670392_avg/1254670392.metafits\" \"-o\" \"tests/data/1254670392_avg/1254670392.cotter.none.ms\" \"-allowmissing\" \"-nostats\" \"-nogeom\" \"-noantennapruning\" \"-nosbgains\" \"-noflagautos\" \"-noflagdcchannels\" \"-nocablelength\" \"-edgewidth\" \"0\" \"-initflag\" \"0\" \"-endflag\" \"0\" \"-sbpassband\" \"tests/data/subband-passband-32ch-unitary.txt\" \"-nostats\" \"-flag-strategy\" \"/usr/share/aoflagger/strategies/mwa-default.lua\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox01_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox02_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox03_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox04_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox05_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox06_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox07_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox08_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox09_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox10_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox11_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox12_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox13_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox14_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox15_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox16_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox17_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox18_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox19_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox20_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox21_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox22_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox23_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox24_00.fits\"",
                "Preprocessed & AOFlagged",
                "Cotter MWA preprocessor",
                "timeavg=1,freqavg=1,windowSize=4",
            )
            .unwrap();
        drop(ms_writer);

        let mut hist_table = Table::open(&hist_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("HISTORY"), TableOpenMode::Read).unwrap();

        assert_table_nrows_match!(hist_table, expected_table);
        for col_name in [
            "TIME",
            "OBSERVATION_ID",
            "MESSAGE",
            "APPLICATION",
            "PRIORITY",
            "ORIGIN",
            "APP_PARAMS",
            "CLI_COMMAND",
        ] {
            assert_table_columns_match!(hist_table, expected_table, col_name);
        }
    }

    #[test]
    #[serial]
    fn test_write_feed_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);
        ms_writer.add_mwa_mods();

        let feed_table_path = table_path.join("FEED");
        let mut feed_table = Table::open(&feed_table_path, TableOpenMode::ReadWrite).unwrap();

        feed_table.add_rows(128).unwrap();

        for idx in 0..128 {
            ms_writer
                .write_feed_row(
                    &mut feed_table,
                    idx as _,
                    idx as _,
                    0,
                    -1,
                    5077351975.,
                    0.,
                    2,
                    -1,
                    &array![[0., 0.], [0., 0.]],
                    &vec!["X".into(), "Y".into()],
                    &array![
                        [c32::new(1., 0.), c32::new(0., 0.)],
                        [c32::new(0., 0.), c32::new(1., 0.)]
                    ],
                    &vec![0., 0., 0.],
                    &vec![0., FRAC_PI_2],
                )
                .unwrap();
        }
        drop(ms_writer);

        let mut feed_table = Table::open(&feed_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("FEED"), TableOpenMode::Read).unwrap();

        assert_table_nrows_match!(feed_table, expected_table);
        for col_name in [
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
        ] {
            assert_table_columns_match!(feed_table, expected_table, col_name);
        }
    }

    #[test]
    #[serial]
    fn test_write_mwa_tile_pointing_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);
        ms_writer.add_mwa_mods();

        let point_table_path = table_path.join("MWA_TILE_POINTING");
        let mut point_table = Table::open(&point_table_path, TableOpenMode::ReadWrite).unwrap();

        point_table.add_rows(1).unwrap();

        ms_writer
            .write_mwa_tile_pointing_row(
                &mut point_table,
                0,
                5077351975.,
                5077351984.,
                &vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0],
                6.283106909188959,
                -0.4644033662289352,
            )
            .unwrap();
        drop(ms_writer);

        let mut point_table = Table::open(&point_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table = Table::open(
            PATH_1254670392.join("MWA_TILE_POINTING"),
            TableOpenMode::Read,
        )
        .unwrap();

        assert_tables_match!(point_table, expected_table);
    }

    #[test]
    #[serial]
    fn test_write_mwa_subband_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);
        ms_writer.add_mwa_mods();

        let subband_table_path = table_path.join("MWA_SUBBAND");
        let mut subband_table = Table::open(&subband_table_path, TableOpenMode::ReadWrite).unwrap();

        subband_table.add_rows(24).unwrap();

        for idx in 0..24 {
            ms_writer
                .write_mwa_subband_row(&mut subband_table, idx, idx as _, 0., false)
                .unwrap();
        }

        drop(ms_writer);

        let mut subband_table = Table::open(&subband_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("MWA_SUBBAND"), TableOpenMode::Read).unwrap();

        assert_tables_match!(subband_table, expected_table);
    }

    #[cfg(feature = "mwalib")]
    #[test]
    #[serial]
    fn test_initialize_from_mwalib_all() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let array_pos = Some(LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        });

        let corr_ctx = CorrelatorContext::new(
            &String::from("tests/data/1254670392_avg/1254670392.metafits"),
            &((1..=24)
                .map(|i| {
                    format!(
                        "tests/data/1254670392_avg/1254670392_20191009153257_gpubox{:02}_00.fits",
                        i
                    )
                })
                .collect::<Vec<_>>()),
        )
        .unwrap();

        let phase_centre = RADec::from_mwalib_phase_or_pointing(&corr_ctx.metafits_context);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, array_pos);

        let mwalib_timestep_range = 0..3;
        let mwalib_coarse_chan_range = *corr_ctx.provided_coarse_chan_indices.first().unwrap()
            ..(*corr_ctx.provided_coarse_chan_indices.last().unwrap() + 1);
        // let mwalib_baseline_idxs: Vec<usize> = (0..corr_ctx.metafits_context.num_baselines).collect();
        let mwalib_baseline_idxs: Vec<usize> = vec![0];

        let (avg_time, avg_freq) = (1, 1);

        ms_writer
            .initialize_from_mwalib(
                &corr_ctx,
                &mwalib_timestep_range,
                &mwalib_coarse_chan_range,
                &mwalib_baseline_idxs,
                avg_time,
                avg_freq,
            )
            .unwrap();

        for (table_name, col_names) in [
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
                    "MWA_INPUT",
                    "MWA_TILE_NR",
                    "MWA_CABLE_LENGTH",
                    // These are wrong in Cotter
                    // "MWA_RECEIVER",
                    // "MWA_SLOT",
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
                    // interval is hardcoded to zero in cotter, it should be obs time
                    // "INTERVAL",
                    "NUM_RECEPTORS",
                    "SPECTRAL_WINDOW_ID",
                    // time is also wrong in Cotter, it should be midpoint, not start time.
                    // "TIME",
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
                    "MWA_HAS_CALIBRATOR",
                ],
            ),
            // WONTDO: this is not written in Cotter
            // (
            //     "FLAG_CMD",
            //     vec![
            //         "APPLIED", "COMMAND", "INTERVAL", "LEVEL", "REASON", "SEVERITY", "TIME", "TYPE",
            //     ],
            // ),
            (
                "HISTORY",
                vec![
                    "OBJECT_ID",
                    "OBSERVATION_ID",
                    "ORIGIN",
                    "PRIORITY",
                    // TODO:
                    // "APP_PARAMS",
                    // "CLI_COMMAND",
                    // "MESSAGE",
                    // WONTDO:
                    // "TIME", // this is wrong in Cotter.
                    // "APPLICATION", // Different application so these will never match
                ],
            ),
            (
                "OBSERVATION",
                vec![
                    "TIME_RANGE",
                    "FLAG_ROW",
                    "OBSERVER",
                    "PROJECT",
                    "RELEASE_DATE",
                    "SCHEDULE_TYPE",
                    "TELESCOPE_NAME",
                    "MWA_GPS_TIME",
                    "MWA_FILENAME",
                    "MWA_OBSERVATION_MODE",
                    "MWA_FLAG_WINDOW_SIZE",
                    "MWA_DATE_REQUESTED",
                    // WONTDO: these are never written in Cotter
                    // "LOG",
                    // "SCHEDULE",
                ],
            ),
            // WONTDO: this is not written in Cotter
            // (
            //     "POINTING",
            //     vec![
            //         "DIRECTION",
            //         "ANTENNA_ID",
            //         "INTERVAL",
            //         "NAME",
            //         "NUM_POLY",
            //         "TARGET",
            //         "TIME",
            //         "TIME_ORIGIN",
            //         "TRACKING",
            //     ],
            // ),
            (
                "POLARIZATION",
                vec!["CORR_TYPE", "CORR_PRODUCT", "FLAG_ROW", "NUM_CORR"],
            ),
            // WONTDO: this is not written in Cotter
            // (
            //     "PROCESSOR",
            //     vec!["FLAG_ROW", "MODE_ID", "TYPE", "TYPE_ID", "SUB_TYPE"],
            // ),
            (
                "SOURCE",
                vec![
                    "SOURCE_ID",
                    "TIME",
                    // Interval is wrong in cotter, it's in days instead of seconds.
                    // "INTERVAL",
                    "SPECTRAL_WINDOW_ID",
                    "NUM_LINES",
                    "NAME",
                    "CALIBRATION_GROUP",
                    "CODE",
                    "DIRECTION",
                    "PROPER_MOTION",
                ],
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
                    "MWA_CENTRE_SUBBAND_NR",
                ],
            ),
            // WONTDO: this is not written in Cotter
            // (
            //     "STATE",
            //     vec![
            //         "CAL", "FLAG_ROW", "LOAD", "OBS_MODE", "REF", "SIG", "SUB_SCAN",
            //     ],
            // ),
        ] {
            let mut table = Table::open(&table_path.join(table_name), TableOpenMode::Read).unwrap();
            let mut exp_table =
                Table::open(PATH_1254670392.join(table_name), TableOpenMode::Read).unwrap();
            if col_names.is_empty() {
                assert_tables_match!(table, exp_table);
            } else {
                assert_table_nrows_match!(table, exp_table);
                for col_name in col_names {
                    let epsilon = 1e-6;
                    assert_table_columns_match!(table, exp_table, col_name, epsilon);
                }
            }
        }

        let main_table = Table::open(&table_path, TableOpenMode::Read).unwrap();

        assert_eq!(
            main_table.n_rows() as usize,
            mwalib_timestep_range.len() * mwalib_baseline_idxs.len()
        );

        let mut ant_table = Table::open(&table_path.join("ANTENNA"), TableOpenMode::Read).unwrap();
        let receivers: Vec<i32> = ant_table.get_col_as_vec("MWA_RECEIVER").unwrap();
        assert_eq!(
            receivers,
            [
                1, 1, 1, 1, 1, 1, 1, 1, 14, 14, 14, 14, 14, 14, 12, 14, 2, 2, 2, 2, 2, 2, 2, 2, 12,
                12, 12, 12, 12, 12, 14, 12, 11, 11, 11, 11, 11, 11, 11, 11, 8, 8, 8, 8, 8, 8, 8, 8,
                9, 9, 9, 9, 9, 9, 9, 9, 13, 13, 13, 13, 13, 13, 13, 13, 15, 15, 15, 15, 15, 15, 15,
                15, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 16, 16, 16, 16, 16, 16, 16, 16,
                7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 10, 10, 10,
                10, 10, 10, 10, 10
            ]
        );
        let num_ants = ant_table.n_rows();
        let exp_slots = Array2::from_shape_fn((128, 2), |(i, _)| (i % 8 + 1) as i32);
        assert_eq!(num_ants, exp_slots.shape()[0] as u64);
        for (idx, exp_slot) in exp_slots.outer_iter().enumerate() {
            let slot: Vec<i32> = ant_table.get_cell_as_vec("MWA_SLOT", idx as _).unwrap();
            assert_eq!(slot, exp_slot.to_vec());
        }
    }

    type VisTestData = (
        Array3<Jones<f32>>,
        Array4<f32>,
        Array4<bool>,
        Array3<f64>,
        Vec<f64>,
        Vec<(usize, usize)>,
    );

    fn get_test_data(
        csv_path: PathBuf,
        num_timesteps: usize,
        num_freqs: usize,
        num_baselines: usize,
    ) -> VisTestData {
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .flexible(true)
            .trim(csv::Trim::All)
            .from_path(csv_path)
            .unwrap();
        let headers = reader.headers().unwrap();
        let keys = ["time", "ant1", "ant2", "u", "v", "w", "pol", "type", "0"];
        let indices = parse_csv_headers(headers, &keys);
        let freq_start_header = indices.get("0").unwrap().to_owned();

        let mut jones = Array3::<Jones<f32>>::zeros((num_timesteps, num_freqs, num_baselines));
        let mut weights = Array4::<f32>::zeros((num_timesteps, num_freqs, num_baselines, 4));
        let mut flags =
            Array4::<bool>::from_elem((num_timesteps, num_freqs, num_baselines, 4), false);
        let mut uvws = Array3::<f64>::zeros((num_timesteps, num_baselines, 3));

        let mut times = Vec::<f64>::new();
        let mut baselines = Vec::<(usize, usize)>::new();
        let mut pols: Vec<String> = Vec::new();

        let mut timestep_idx;
        let mut baseline_idx;
        let mut pol_idx;

        for row in reader.records() {
            let record = row.unwrap();
            let time = record[indices["time"]].parse::<f64>().unwrap();
            timestep_idx = if let Some(idx) = times.iter().position(|&x| x == time) {
                idx
            } else {
                times.push(time);
                times.len() - 1
            };

            let ant1 = record[indices["ant1"]].parse::<usize>().unwrap();
            let ant2 = record[indices["ant2"]].parse::<usize>().unwrap();
            let baseline = (ant1, ant2);
            baseline_idx = if let Some(idx) = baselines.iter().position(|&x| x == baseline) {
                idx
            } else {
                baselines.push(baseline);
                baselines.len() - 1
            };

            let pol: String = record[indices["pol"]].into();
            pol_idx = if let Some(idx) = pols.iter().position(|x| x == &pol) {
                idx
            } else {
                pols.push(pol);
                pols.len() - 1
            };

            for (uvw_idx, &uvw) in ["u", "v", "w"].iter().enumerate() {
                uvws[[timestep_idx, baseline_idx, uvw_idx]] = record[indices[uvw]].parse().unwrap();
            }

            let row_type: &str = record[indices["type"]].into();

            for (freq_idx, val) in record.iter().skip(freq_start_header).enumerate() {
                match row_type {
                    "vis" => {
                        jones[(timestep_idx, freq_idx, baseline_idx)][pol_idx] = parse_complex(val);
                    }
                    "weight" => {
                        weights[(timestep_idx, freq_idx, baseline_idx, pol_idx)] =
                            val.parse().unwrap();
                    }
                    "flag" => {
                        flags[(timestep_idx, freq_idx, baseline_idx, pol_idx)] =
                            val.to_lowercase().parse().unwrap();
                    }
                    _ => {
                        panic!("Unknown row type: {}", row_type);
                    }
                }
            }
        }

        (jones, weights, flags, uvws, times, baselines)
    }

    lazy_static! {
        static ref COMPLEX_REGEX: Regex = Regex::new(format!(
                r"^(?P<only_real>{0})$|^(?P<only_imag>{0})j$|^\((?P<complex_real>{0})\+?(?P<complex_imag>{0})j\)$",
                r"-?[\d\.]+(e-?\d+)?"
            ).as_str()
        ).unwrap();
    }

    fn parse_complex(cell: &str) -> Complex<f32> {
        let captures = COMPLEX_REGEX.captures(cell).unwrap();
        let (real, imag) = match (
            captures.name("complex_real"),
            captures.name("complex_imag"),
            captures.name("only_real"),
            captures.name("only_imag"),
        ) {
            (Some(real), Some(imag), _, _) => (
                parse::<f32, _>(real.as_str()).unwrap(),
                parse::<f32, _>(imag.as_str()).unwrap(),
            ),
            (None, None, Some(real), None) => (parse::<f32, _>(real.as_str()).unwrap(), 0.0),
            (None, None, None, Some(imag)) => (0.0, parse::<f32, _>(imag.as_str()).unwrap()),
            _ => panic!("can't parse complex {}", cell),
        };
        Complex::new(real, imag)
    }

    fn parse_csv_headers(headers: &csv::StringRecord, keys: &[&str]) -> BTreeMap<String, usize> {
        let mut remaining_keys: HashSet<_> = keys.iter().map(|x| String::from(*x)).collect();
        let mut indices = BTreeMap::<String, usize>::new();

        for (idx, cell) in headers.iter().enumerate() {
            let mut remove: Option<String> = None;
            for key in remaining_keys.iter() {
                if cell == key {
                    indices.insert(String::from(cell), idx);
                    remove = Some(key.clone());
                    break;
                }
            }
            if let Some(key) = remove {
                remaining_keys.remove(&key);
            }
        }

        if !remaining_keys.is_empty() {
            panic!("not all keys found: {:?}", remaining_keys);
        }

        indices
    }

    #[test]
    #[serial]
    fn test_write_main_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.47123889803846897);
        let ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let mut main_table = Table::open(&table_path, TableOpenMode::ReadWrite).unwrap();

        let num_timesteps = 2;
        let num_freqs = 768;
        let num_baselines = 1;

        main_table.add_rows(num_timesteps * num_baselines).unwrap();

        let (jones, weights, flags, uvws, times, baselines) = get_test_data(
            "tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms.csv".into(),
            num_timesteps,
            num_freqs,
            num_baselines,
        );

        let mut row_flags = Array::from_elem((768, 4), false);
        let mut row_weights = Array::zeros((768, 4));

        let mut row_idx = 0;
        for (timestep_idx, &time) in times.iter().enumerate() {
            for (baseline_idx, &(ant1, ant2)) in baselines.iter().enumerate() {
                let uvw: Vec<f64> = uvws
                    .slice(s![timestep_idx, baseline_idx, ..])
                    .iter()
                    .copied()
                    .collect();
                let data_array =
                    Array2::from_shape_fn((768, 4), |(c, p)| jones[[timestep_idx, c, 0]][p]);

                row_flags.assign(&flags.slice(s![timestep_idx, .., baseline_idx, ..]));
                row_weights.assign(&weights.slice(s![timestep_idx, .., baseline_idx, ..]));

                ms_writer
                    .write_main_row(
                        &mut main_table,
                        row_idx,
                        time,
                        time,
                        ant1 as _,
                        ant2 as _,
                        0,
                        &uvw,
                        2.,
                        -1,
                        1,
                        -1,
                        &vec![1., 1., 1., 1.],
                        &data_array,
                        &row_flags,
                        &row_weights,
                        false,
                    )
                    .unwrap();

                row_idx += 1;
            }
        }

        drop(ms_writer);

        let mut main_table = Table::open(&table_path, TableOpenMode::Read).unwrap();
        let mut expected_table =
            Table::open(PATH_1254670392.join(""), TableOpenMode::Read).unwrap();

        assert_table_nrows_match!(main_table, expected_table);
        for col_name in [
            "ANTENNA1",
            "ANTENNA2",
            "ARRAY_ID",
            "DATA_DESC_ID",
            "DATA",
            "EXPOSURE",
            "FEED1",
            "FEED2",
            "FIELD_ID",
            // TODO
            // "FLAG_CATEGORY",
            "FLAG_ROW",
            "FLAG",
            "INTERVAL",
            "OBSERVATION_ID",
            "PROCESSOR_ID",
            "SCAN_NUMBER",
            "SIGMA",
            "STATE_ID",
            "TIME_CENTROID",
            "TIME",
            "UVW",
            "WEIGHT_SPECTRUM",
            "WEIGHT",
        ] {
            assert_table_columns_match!(main_table, expected_table, col_name);
        }
    }

    const REPRODUCIBLE_TABLE_COLNAMES: &[(&str, &[&str])] = &[
        (
            "",
            &[
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
            &[
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
            &["FLAG_ROW", "POLARIZATION_ID", "SPECTRAL_WINDOW_ID"],
        ),
        (
            "FEED",
            &[
                "POSITION",
                "BEAM_OFFSET",
                "POLARIZATION_TYPE",
                "POL_RESPONSE",
                "RECEPTOR_ANGLE",
                "ANTENNA_ID",
                "BEAM_ID",
                "FEED_ID",
                "NUM_RECEPTORS",
                "SPECTRAL_WINDOW_ID",
                // These are wrong in Cotter, at least according to https://casa.nrao.edu/Memos/229.html#SECTION00066000000000000000
                // "INTERVAL",
                // "TIME",
            ],
        ),
        (
            "FIELD",
            &[
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
        // (
        //     "FLAG_CMD",
        //     &[
        //         "APPLIED", "COMMAND", "INTERVAL", "LEVEL", "REASON", "SEVERITY", "TIME", "TYPE",
        //     ],
        // ),
        // (
        //     "HISTORY",
        //     &[
        //         "APP_PARAMS",
        //         "CLI_COMMAND",
        //         "APPLICATION",
        //         "MESSAGE",
        //         "OBJECT_ID",
        //         "OBSERVATION_ID",
        //         "ORIGIN",
        //         "PRIORITY",
        //         "TIME",
        //     ],
        // ),
        (
            "OBSERVATION",
            &[
                "TIME_RANGE",
                // "LOG",
                // "SCHEDULE",
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
            &[
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
            &["CORR_TYPE", "CORR_PRODUCT", "FLAG_ROW", "NUM_CORR"],
        ),
        (
            "PROCESSOR",
            &["FLAG_ROW", "MODE_ID", "TYPE", "TYPE_ID", "SUB_TYPE"],
        ),
        (
            "SPECTRAL_WINDOW",
            &[
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
            &[
                "CAL", "FLAG_ROW", "LOAD", "OBS_MODE", "REF", "SIG", "SUB_SCAN",
            ],
        ),
    ];

    #[cfg(feature = "mwalib")]
    #[test]
    #[serial]
    fn test_write_vis_from_mwalib() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let array_pos = Some(LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        });

        let corr_ctx = CorrelatorContext::new(
            &String::from("tests/data/1254670392_avg/1254670392.metafits"),
            &((1..=24)
                .map(|i| {
                    format!(
                        "tests/data/1254670392_avg/1254670392_20191009153257_gpubox{:02}_00.fits",
                        i
                    )
                })
                .collect::<Vec<_>>()),
        )
        .unwrap();

        let phase_centre = RADec::from_mwalib_phase_or_pointing(&corr_ctx.metafits_context);
        let mut ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, array_pos);

        let mwalib_timestep_range = 0..2_usize;
        let mwalib_coarse_chan_range = *corr_ctx.provided_coarse_chan_indices.first().unwrap()
            ..(*corr_ctx.provided_coarse_chan_indices.last().unwrap() + 1);
        let mwalib_baseline_idxs = vec![1_usize];

        let (avg_time, avg_freq) = (1, 1);

        ms_writer
            .initialize_from_mwalib(
                &corr_ctx,
                &mwalib_timestep_range,
                &mwalib_coarse_chan_range,
                &mwalib_baseline_idxs,
                avg_time,
                avg_freq,
            )
            .unwrap();

        let (jones_array, weight_array, flag_array, _, _, _) = get_test_data(
            "tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms.csv".into(),
            2,
            768,
            1,
        );

        ms_writer
            .write_vis_mwalib(
                jones_array.view(),
                weight_array.view(),
                flag_array.view(),
                &corr_ctx,
                &mwalib_timestep_range,
                &mwalib_coarse_chan_range,
                &mwalib_baseline_idxs,
                avg_time,
                avg_freq,
                false,
            )
            .unwrap();

        for (table_name, col_names) in REPRODUCIBLE_TABLE_COLNAMES {
            let mut table = Table::open(&table_path.join(table_name), TableOpenMode::Read).unwrap();
            let mut exp_table =
                Table::open(PATH_1254670392.join(table_name), TableOpenMode::Read).unwrap();
            assert_table_nrows_match!(table, exp_table);
            for col_name in col_names.iter() {
                if ["TIME_CENTROID", "TIME"].contains(col_name) {
                    assert_table_columns_match!(table, exp_table, col_name, 5e-6);
                } else {
                    assert_table_columns_match!(table, exp_table, col_name);
                }
            }
        }
    }

    #[cfg(feature = "mwalib")]
    #[test]
    #[serial]
    fn test_write_vis_from_mwalib_averaging() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let array_pos = Some(LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        });

        let corr_ctx = CorrelatorContext::new(
            &String::from("tests/data/1254670392_avg/1254670392.metafits"),
            &((1..=24)
                .map(|i| {
                    format!(
                        "tests/data/1254670392_avg/1254670392_20191009153257_gpubox{:02}_00.fits",
                        i
                    )
                })
                .collect::<Vec<_>>()),
        )
        .unwrap();

        let phase_centre = RADec::from_mwalib_phase_or_pointing(&corr_ctx.metafits_context);
        let mut ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, array_pos);

        let mwalib_timestep_range = 0..2_usize;
        let mwalib_coarse_chan_range = *corr_ctx.provided_coarse_chan_indices.first().unwrap()
            ..(*corr_ctx.provided_coarse_chan_indices.last().unwrap() + 1);
        let mwalib_baseline_idxs = vec![1_usize];

        let (avg_time, avg_freq) = (2, 2);

        ms_writer
            .initialize_from_mwalib(
                &corr_ctx,
                &mwalib_timestep_range,
                &mwalib_coarse_chan_range,
                &mwalib_baseline_idxs,
                avg_time,
                avg_freq,
            )
            .unwrap();

        let (jones_array, weight_array, flag_array, _, _, _) = get_test_data(
            "tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms.csv".into(),
            2,
            768,
            1,
        );

        ms_writer
            .write_vis_mwalib(
                jones_array.view(),
                weight_array.view(),
                flag_array.view(),
                &corr_ctx,
                &mwalib_timestep_range,
                &mwalib_coarse_chan_range,
                &mwalib_baseline_idxs,
                avg_time,
                avg_freq,
                false,
            )
            .unwrap();

        for (table_name, col_names) in REPRODUCIBLE_TABLE_COLNAMES {
            let mut table = Table::open(&table_path.join(table_name), TableOpenMode::Read).unwrap();
            let mut exp_table = Table::open(
                PATH_1254670392_AVG_4S_80KHZ.join(table_name),
                TableOpenMode::Read,
            )
            .unwrap();
            assert_table_nrows_match!(table, exp_table);
            for col_name in col_names.iter() {
                if ["TIME_CENTROID", "TIME"].contains(col_name) {
                    assert_table_columns_match!(table, exp_table, col_name, 5e-6);
                } else {
                    assert_table_columns_match!(table, exp_table, col_name);
                }
            }
        }
    }

    /// as above, but with two consecutive calls to write_vis_mwalib
    #[cfg(feature = "mwalib")]
    #[test]
    #[serial]
    fn test_write_vis_from_mwalib_chunks() {
        use itertools::Itertools;

        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let array_pos = Some(LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        });

        let corr_ctx = CorrelatorContext::new(
            &String::from("tests/data/1254670392_avg/1254670392.metafits"),
            &((1..=24)
                .map(|i| {
                    format!(
                        "tests/data/1254670392_avg/1254670392_20191009153257_gpubox{:02}_00.fits",
                        i
                    )
                })
                .collect::<Vec<_>>()),
        )
        .unwrap();

        let phase_centre = RADec::from_mwalib_phase_or_pointing(&corr_ctx.metafits_context);
        let mut ms_writer = MeasurementSetWriter::new(&table_path, phase_centre, array_pos);

        let mwalib_timestep_range = 0..2_usize;
        let mwalib_coarse_chan_range = *corr_ctx.provided_coarse_chan_indices.first().unwrap()
            ..(*corr_ctx.provided_coarse_chan_indices.last().unwrap() + 1);
        let mwalib_baseline_idxs = vec![1_usize];

        let (avg_time, avg_freq) = (1, 1);

        ms_writer
            .initialize_from_mwalib(
                &corr_ctx,
                &mwalib_timestep_range,
                &mwalib_coarse_chan_range,
                &mwalib_baseline_idxs,
                avg_time,
                avg_freq,
            )
            .unwrap();

        let (jones_array, weight_array, flag_array, _, _, _) = get_test_data(
            "tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms.csv".into(),
            2,
            768,
            1,
        );

        let num_chunk_timesteps = 1;

        for (mut timestep_chunk, jones_array_chunk, weight_array_chunk, flag_array_chunk) in izip!(
            &mwalib_timestep_range.chunks(num_chunk_timesteps),
            jones_array.axis_chunks_iter(Axis(0), num_chunk_timesteps),
            weight_array.axis_chunks_iter(Axis(0), num_chunk_timesteps),
            flag_array.axis_chunks_iter(Axis(0), num_chunk_timesteps)
        ) {
            let first_timestep = timestep_chunk.next().unwrap();
            let last_timestep = timestep_chunk.last().unwrap_or(first_timestep);
            let timestep_range: Range<usize> = first_timestep..last_timestep + 1;
            ms_writer
                .write_vis_mwalib(
                    jones_array_chunk.view(),
                    weight_array_chunk.view(),
                    flag_array_chunk.view(),
                    &corr_ctx,
                    &timestep_range,
                    &mwalib_coarse_chan_range,
                    &mwalib_baseline_idxs,
                    avg_time,
                    avg_freq,
                    false,
                )
                .unwrap();
        }

        for (table_name, col_names) in REPRODUCIBLE_TABLE_COLNAMES {
            let mut table = Table::open(&table_path.join(table_name), TableOpenMode::Read).unwrap();
            let mut exp_table =
                Table::open(PATH_1254670392.join(table_name), TableOpenMode::Read).unwrap();
            assert_table_nrows_match!(table, exp_table);
            for col_name in col_names.iter() {
                if ["TIME_CENTROID", "TIME"].contains(col_name) {
                    assert_table_columns_match!(table, exp_table, col_name, 5e-6);
                } else {
                    assert_table_columns_match!(table, exp_table, col_name);
                }
            }
        }
    }
}
