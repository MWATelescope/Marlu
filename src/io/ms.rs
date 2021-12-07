use crate::{
    c32,
    io::error::MeasurementSetWriteError,
    ndarray::{array, Array2, Array3, ArrayView3, ArrayView4, Axis},
    precession::precess_time,
    time::{gps_millis_to_epoch, gps_to_epoch},
    Jones, LatLngHeight, RADec, RubblArray, XyzGeodetic, ENH, UVW,
};
use flate2::read::GzDecoder;
use itertools::izip;
#[cfg(feature = "mwalib")]
use mwalib::CorrelatorContext;
use rubbl_casatables::{
    GlueDataType, Table, TableCreateMode, TableDesc, TableDescCreateMode, TableOpenMode,
    TableRecord,
};
use std::{
    f64::consts::PI,
    fs::create_dir_all,
    ops::Range,
    path::{Path, PathBuf},
    time::SystemTime,
};
use tar::Archive;

use lazy_static::lazy_static;
use log::trace;

use super::{error::IOError, VisWritable};

lazy_static! {
    static ref DEFAULT_TABLES_GZ: &'static [u8] =
        include_bytes!("../../data/default_tables.tar.gz");
    static ref SOURCE_TABLE_GZ: &'static [u8] = include_bytes!("../../data/source_table.tar.gz");
}

const PKG_VERSION: &'static str = env!("CARGO_PKG_VERSION");
const PKG_NAME: &'static str = env!("CARGO_PKG_NAME");

/// This is a very stupid hack that lets us use rubbl's ndarray. (┛◉Д◉)┛彡┻━┻
/// TODO: reduce double-handling
macro_rules! rubblify_array {
    ($array:expr) => {
        RubblArray::from_shape_vec($array.dim(), $array.to_owned().into_raw_vec()).unwrap()
    };
}

/// A helper struct to write out a uvfits file.
///
// pub struct MeasurementSetWriter<'a> {
pub struct MeasurementSetWriter {
    /// The path to the root of the measurement set (typically ends in .ms)
    path: PathBuf,

    /// The RA/Dec where this observation is phased to
    phase_centre: RADec,

    /// Array Position [Latitude (radians), Longitude (radians), Height (m)]
    array_pos: LatLngHeight,
}

// impl<'a> MeasurementSetWriter<'a> {
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

        let mut main_table = Table::open(self.path.clone(), TableOpenMode::ReadWrite).unwrap();
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

        let mut main_table = Table::open(self.path.clone(), TableOpenMode::ReadWrite).unwrap();
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
    pub fn write_spectral_window_row(
        &self,
        table: &mut Table,
        idx: u64,
        name: &str,
        ref_freq: f64,
        chan_info: Array2<f64>,
        total_bw: f64,
        flag: bool,
    ) -> Result<(), MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        match chan_info.shape() {
            [num_chans, 4] => {
                table
                    .put_cell("NUM_CHAN", idx, &(*num_chans as i32))
                    .unwrap();
            }
            sh => {
                return Err(MeasurementSetWriteError::BadArrayShape {
                    argument: "chan_info".into(),
                    function: "write_spectral_window_row".into(),
                    expected: "[n, 4]".into(),
                    received: format!("{:?}", sh).into(),
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
    ) -> Result<(), MeasurementSetWriteError> {
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
    ) -> Result<(), MeasurementSetWriteError> {
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
    pub fn write_antenna_row(
        &self,
        table: &mut Table,
        idx: u64,
        name: &str,
        station: &str,
        ant_type: &str,
        mount: &str,
        position: &Vec<f64>,
        dish_diameter: f64,
        flag_row: bool,
    ) -> Result<(), MeasurementSetWriteError> {
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
    pub fn write_antenna_row_mwa(
        &self,
        table: &mut Table,
        idx: u64,
        name: &str,
        station: &str,
        ant_type: &str,
        mount: &str,
        position: &Vec<f64>,
        dish_diameter: f64,
        input: &Vec<i32>,
        tile_nr: i32,
        receiver: i32,
        slot: &Vec<i32>,
        cable_length: &Vec<f64>,
        flag_row: bool,
    ) -> Result<(), MeasurementSetWriteError> {
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
    pub fn write_polarization_row(
        &self,
        table: &mut Table,
        idx: u64,
        corr_type: &Vec<i32>,
        corr_product: &Array2<i32>,
        flag_row: bool,
    ) -> Result<(), MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        let num_corr_type = corr_type.len();

        match corr_product.shape() {
            [num_corr, 2] if *num_corr == num_corr_type => {
                table
                    .put_cell("NUM_CORR", idx, &(*num_corr as i32))
                    .unwrap();
            }
            sh => {
                return Err(MeasurementSetWriteError::BadArrayShape {
                    argument: "corr_product".into(),
                    function: "write_polarization_row".into(),
                    expected: format!("[n, 2] (where n = corr_type.len() = {})", num_corr_type)
                        .into(),
                    received: format!("{:?}", sh).into(),
                })
            }
        }

        table.put_cell("CORR_TYPE", idx, corr_type).unwrap();

        let corr_product = rubblify_array!(corr_product);
        table.put_cell("CORR_PRODUCT", idx, &corr_product).unwrap();
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
        direction: Vec<f64>,
        proper_motion: Vec<f64>,
    ) -> Result<(), MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        match (direction.len(), proper_motion.len()) {
            (2, 2) => {
                table.put_cell("DIRECTION", idx, &direction).unwrap();
                table
                    .put_cell("PROPER_MOTION", idx, &proper_motion)
                    .unwrap();
            }
            sh => {
                return Err(MeasurementSetWriteError::BadArrayShape {
                    argument: "direction|proper_motion".into(),
                    function: "write_source_row".into(),
                    expected: format!("(2, 2)").into(),
                    received: format!("{:?}", sh).into(),
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
    ) -> Result<(), MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        let dir_info = rubblify_array!(dir_info);

        match dir_info.shape() {
            [3, p, 2] if *p > 0 => {
                table.put_cell("NUM_POLY", idx, &((*p - 1) as i32)).unwrap();
            }
            sh => {
                return Err(MeasurementSetWriteError::BadArrayShape {
                    argument: "dir_info".into(),
                    function: "write_field_row".into(),
                    expected: format!("[3, p, 2] (where p is highest polynomial order)").into(),
                    received: format!("{:?}", sh).into(),
                })
            }
        }

        table.put_cell("NAME", idx, &name.to_string()).unwrap();
        table.put_cell("CODE", idx, &code.to_string()).unwrap();
        table.put_cell("TIME", idx, &time).unwrap();

        let col_names = ["DELAY_DIR", "PHASE_DIR", "REFERENCE_DIR"];
        for (value, &col_name) in dir_info.outer_iter().zip(col_names.iter()) {
            // println!("{:?}", value.shape());
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
    ) -> Result<(), MeasurementSetWriteError> {
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
    ) -> Result<(), MeasurementSetWriteError> {
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
    ) -> Result<(), MeasurementSetWriteError> {
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
    pub fn write_history_row(
        &self,
        table: &mut Table,
        idx: u64,
        time: f64,
        cmd_line: &str,
        message: &str,
        application: &str,
        params: &str,
    ) -> Result<(), MeasurementSetWriteError> {
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
        position: &Vec<f64>,
        receptor_angle: &Vec<f64>,
    ) -> Result<(), MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        if beam_offset.shape() != &[num_receptors as _, 2] {
            return Err(MeasurementSetWriteError::BadArrayShape {
                argument: "beam_offset".into(),
                function: "write_feed_row".into(),
                expected: "[n, 2]".into(),
                received: format!("{:?}", beam_offset.shape()).into(),
            });
        }
        if pol_type.len() != num_receptors as _ {
            return Err(MeasurementSetWriteError::BadArrayShape {
                argument: "pol_type".into(),
                function: "write_feed_row".into(),
                expected: "n".into(),
                received: format!("{:?}", pol_type.len()).into(),
            });
        }
        if pol_response.shape() != &[num_receptors as _, num_receptors as _] {
            return Err(MeasurementSetWriteError::BadArrayShape {
                argument: "pol_response".into(),
                function: "write_feed_row".into(),
                expected: "[n, n]".into(),
                received: format!("{:?}", pol_response.shape()).into(),
            });
        }
        if position.len() != 3 as _ {
            return Err(MeasurementSetWriteError::BadArrayShape {
                argument: "position".into(),
                function: "write_feed_row".into(),
                expected: "3".into(),
                received: format!("{:?}", position.len()).into(),
            });
        }
        if receptor_angle.len() != num_receptors as _ {
            return Err(MeasurementSetWriteError::BadArrayShape {
                argument: "receptor_angle".into(),
                function: "write_feed_row".into(),
                expected: "n".into(),
                received: format!("{:?}", receptor_angle.len()).into(),
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
    pub fn write_mwa_tile_pointing_row(
        &self,
        table: &mut Table,
        idx: u64,
        start: f64,
        end: f64,
        delays: &Vec<i32>,
        direction_ra: f64,
        direction_dec: f64,
    ) -> Result<(), MeasurementSetWriteError> {
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
    ) -> Result<(), MeasurementSetWriteError> {
        table.put_cell("NUMBER", idx, &number).unwrap();
        table.put_cell("GAIN", idx, &gain).unwrap();
        table.put_cell("FLAG_ROW", idx, &flag_row).unwrap();

        Ok(())
    }

    /// Create an MWA measurement set, with all tables (except the main visibility table)
    /// prefilled with metadata from a [`mwalib::CorrelatorContext`]
    ///
    /// `mwalib_timestep_range` the range of timestep indices (according to mwalib)
    /// of the current chunk being written to the measurement set.
    ///
    /// `mwalib_coarse_chan_range` the range of coarse channel indices (according to mwalib)
    /// of the current chunk being written to the measurement set.
    #[cfg(feature = "mwalib")]
    pub fn initialize_from_mwalib(
        &self,
        context: &CorrelatorContext,
        mwalib_timestep_range: &Range<usize>,
        mwalib_coarse_chan_range: &Range<usize>,
    ) -> Result<(), MeasurementSetWriteError> {
        let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
        let num_img_coarse_chans = mwalib_coarse_chan_range.len();
        let num_img_chans = fine_chans_per_coarse * num_img_coarse_chans;

        self.decompress_default_tables().unwrap();
        self.decompress_source_table().unwrap();
        self.add_cotter_mods(num_img_chans);
        self.add_mwa_mods();

        // /////////////// //
        // Spectral Window //
        // /////////////// //

        let mut spw_table =
            Table::open(&self.path.join("SPECTRAL_WINDOW"), TableOpenMode::ReadWrite).unwrap();

        let mwalib_centre_coarse_chan_idx =
            mwalib_coarse_chan_range.start + (num_img_coarse_chans / 2);
        let centre_coarse_chan =
            context.metafits_context.metafits_coarse_chans[mwalib_centre_coarse_chan_idx].clone();
        let fine_chan_width_hz = context.metafits_context.corr_fine_chan_width_hz;
        let mwalib_start_fine_chan_idx = mwalib_coarse_chan_range.start * fine_chans_per_coarse;

        let chan_info = Array2::from_shape_fn((num_img_chans, 4), |(c, i)| {
            if i == 0 {
                context.metafits_context.metafits_fine_chan_freqs_hz[c + mwalib_start_fine_chan_idx]
            } else {
                fine_chan_width_hz as f64
            }
        });

        let img_center_freq_hz = if num_img_chans % 2 == 0 {
            (context.metafits_context.metafits_fine_chan_freqs_hz
                [mwalib_start_fine_chan_idx + (num_img_chans / 2)]
                + context.metafits_context.metafits_fine_chan_freqs_hz
                    [mwalib_start_fine_chan_idx + (num_img_chans / 2) - 1])
                * 0.5
        } else {
            context.metafits_context.metafits_fine_chan_freqs_hz
                [mwalib_start_fine_chan_idx + (num_img_chans / 2)]
        };

        spw_table.add_rows(1).unwrap();
        self.write_spectral_window_row_mwa(
            &mut spw_table,
            0,
            format!("MWA_BAND_{:.1}", img_center_freq_hz / 1_000_000.).as_str(),
            img_center_freq_hz,
            chan_info,
            fine_chan_width_hz as f64 * num_img_chans as f64,
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

        let antennae = context.metafits_context.antennas.clone();
        ant_table.add_rows(antennae.len()).unwrap();
        for (idx, antenna) in antennae.into_iter().enumerate() {
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
                &antenna.tile_name.as_str(),
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

        let ra_phase_rad = context.metafits_context.ra_phase_center_degrees.unwrap() * (PI / 180.);
        let dec_phase_rad =
            context.metafits_context.dec_phase_center_degrees.unwrap() * (PI / 180.);

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

        let obs_name = context.metafits_context.obs_name.clone();
        let field_name = obs_name
            .rsplit_once("_")
            .unwrap_or((obs_name.as_str(), ""))
            .0;

        let sched_start_time_mjd_utc_s =
            gps_millis_to_epoch(context.metafits_context.sched_start_gps_time_ms)
                .as_mjd_utc_seconds();

        field_table.add_rows(1).unwrap();
        self.write_field_row_mwa(
            &mut field_table,
            0,
            &field_name,
            "",
            sched_start_time_mjd_utc_s,
            &dir_info,
            -1,
            context.metafits_context.calibrator,
            false,
        )
        .unwrap();

        // ////// //
        // Source //
        // ////// //

        let mut source_table =
            Table::open(&self.path.join("SOURCE"), TableOpenMode::ReadWrite).unwrap();

        let duration_ms = context.metafits_context.sched_duration_ms;
        let int_time_ms = context.metafits_context.corr_int_time_ms;

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
            &field_name,
            0,
            &"",
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

        let start_time_centroid_mjd_utc_s =
            gps_millis_to_epoch(context.metafits_context.sched_start_gps_time_ms + int_time_ms / 2)
                .as_mjd_utc_seconds();
        let end_time_centroid_mjd_utc_s =
            gps_millis_to_epoch(context.metafits_context.sched_end_gps_time_ms + int_time_ms / 2)
                .as_mjd_utc_seconds();

        self.write_observation_row_mwa(
            &mut obs_table,
            0,
            "MWA",
            (start_time_centroid_mjd_utc_s, end_time_centroid_mjd_utc_s),
            &context.metafits_context.creator,
            "MWA",
            &context.metafits_context.project_id,
            0.,
            context.metafits_context.obs_id as _,
            &context.metafits_context.obs_name,
            &context.metafits_context.mode.to_string(),
            (mwalib_timestep_range.len() + 1) as _,
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
        assert_eq!(context.metafits_context.num_ant_pols, 2);

        let num_ants = context.metafits_context.num_ants;

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
                context.metafits_context.num_ant_pols as _,
                -1,
                &array![[0., 0.], [0., 0.]],
                &vec!["X".into(), "Y".into()],
                &array![
                    [c32::new(1., 0.), c32::new(0., 0.)],
                    [c32::new(0., 0.), c32::new(1., 0.)]
                ],
                &vec![0., 0., 0.],
                &vec![0., PI / 2.],
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

        let delays = context
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

        subband_table.add_rows(num_img_coarse_chans).unwrap();

        for i in 0..num_img_coarse_chans {
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
    /// - `time` - Modified Julian Day, at start of scan
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
    pub fn write_main_row(
        &self,
        table: &mut Table,
        idx: u64,
        time: f64,
        time_centroid: f64,
        antenna1: i32,
        antenna2: i32,
        data_desc_id: i32,
        uvw: &Vec<f64>,
        interval: f64,
        // TODO: is this not just interval?
        // exposure: f64,
        processor_id: i32,
        scan_number: i32,
        state_id: i32,
        sigma: &Vec<f32>,
        data: Array2<c32>,
        flags: Array2<bool>,
        weights: Array2<f32>,
        flag_row: bool,
    ) -> Result<(), MeasurementSetWriteError> {
        let num_pols = 4;

        if uvw.len() != 3 {
            return Err(MeasurementSetWriteError::BadArrayShape {
                argument: "uvw".into(),
                function: "write_main_row".into(),
                expected: format!("3").into(),
                received: format!("{:?}", uvw.len()).into(),
            });
        }

        if sigma.len() != num_pols {
            return Err(MeasurementSetWriteError::BadArrayShape {
                argument: "sigma".into(),
                function: "write_main_row".into(),
                expected: format!("{}", num_pols).into(),
                received: format!("{:?}", sigma.len()).into(),
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
                return Err(MeasurementSetWriteError::BadArrayShape {
                    argument: "data|flags|weights".into(),
                    function: "write_main_row".into(),
                    expected: format!(
                        "[n, p]|[n, p]|[n, p] where n=num_chans, p=num_pols({})",
                        num_pols
                    )
                    .into(),
                    received: format!("{:?}|{:?}|{:?}", dsh, fsh, wsh).into(),
                })
            }
        }

        let weight_pol = weights
            .axis_iter(Axis(1))
            .map(|weights_pol_view| weights_pol_view.sum())
            .collect::<Vec<f32>>();

        // TODO: get rid of this disgusting and wasteful hack.
        let data = rubblify_array!(data);
        let flags = rubblify_array!(flags);
        let weights = rubblify_array!(weights);

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
        table.put_cell("DATA", idx, &data).unwrap();
        table.put_cell("WEIGHT_SPECTRUM", idx, &weights).unwrap();
        table.put_cell("WEIGHT", idx, &weight_pol).unwrap();
        table.put_cell("FLAG", idx, &flags).unwrap();
        table.put_cell("FLAG_ROW", idx, &flag_row).unwrap();

        Ok(())
    }
}

impl VisWritable for MeasurementSetWriter {
    #[cfg(feature = "mwalib")]
    fn write_vis_mwalib(
        &mut self,
        jones_array: ArrayView3<Jones<f32>>,
        weight_array: ArrayView4<f32>,
        flag_array: ArrayView4<bool>,
        context: &CorrelatorContext,
        timestep_range: &Range<usize>,
        coarse_chan_range: &Range<usize>,
        baseline_idxs: &[usize],
    ) -> Result<(), IOError> {
        trace!(
            "timestep range {:?}, coarse chan range {:?}, baseline idxs {:?}",
            &timestep_range,
            &coarse_chan_range,
            &baseline_idxs
        );

        let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;

        let jones_dims = jones_array.dim();
        let weight_dims = weight_array.dim();
        assert_eq!(jones_dims.0, weight_dims.0);
        assert_eq!(jones_dims.1, weight_dims.1);
        assert_eq!(jones_dims.2, weight_dims.2);

        let num_img_timesteps = timestep_range.len();
        assert_eq!(num_img_timesteps, jones_dims.0);

        let num_img_coarse_chans = coarse_chan_range.len();
        let num_img_chans = fine_chans_per_coarse * num_img_coarse_chans;
        assert_eq!(num_img_coarse_chans * fine_chans_per_coarse, jones_dims.1);

        let num_baselines = baseline_idxs.len();
        assert_eq!(num_baselines, jones_dims.2);

        let total_num_rows = num_img_timesteps * num_baselines;

        // let phase_centre = match phase_centre {
        //     Some(radec) => radec,
        //     None => RADec::from_mwalib_phase_or_pointing(&context.metafits_context),
        // };

        // Weights are normalized so that default res of 10 kHz, 1s has weight of "1" per sample
        let integration_time_s = context.metafits_context.corr_int_time_ms as f64 / 1000.0;

        let tiles_xyz_geod = XyzGeodetic::get_tiles_mwa(&context.metafits_context);
        let img_timesteps = &context.timesteps[timestep_range.clone()];
        let img_baselines = baseline_idxs
            .iter()
            .map(|&idx| {
                context
                    .metafits_context
                    .baselines
                    .get(idx)
                    .unwrap()
                    .to_owned()
            })
            .collect::<Vec<_>>();

        let mut main_table = Table::open(&self.path.clone(), TableOpenMode::ReadWrite).unwrap();

        main_table.add_rows(total_num_rows).unwrap();

        let mut main_idx = 0;

        for (timestep, jones_timestep_view, weight_timestep_view, flag_timestep_view) in izip!(
            img_timesteps.iter(),
            jones_array.outer_iter(),
            weight_array.outer_iter(),
            flag_array.outer_iter()
        ) {
            let gps_time_s = timestep.gps_time_ms as f64 / 1000.0;
            let centroid_epoch = gps_to_epoch(gps_time_s + integration_time_s / 2.0);

            let scan_start_mjd_utc_s =
                gps_millis_to_epoch(timestep.gps_time_ms).as_mjd_utc_seconds();
            let scan_centroid_mjd_utc_s = gps_millis_to_epoch(
                timestep.gps_time_ms + context.metafits_context.corr_int_time_ms / 2,
            )
            .as_mjd_utc_seconds();

            let prec_info = precess_time(
                self.phase_centre,
                centroid_epoch,
                self.array_pos.longitude_rad,
                self.array_pos.latitude_rad,
            );

            let tiles_xyz_precessed = prec_info.precess_xyz_parallel(&tiles_xyz_geod);

            for (baseline, jones_baseline_view, weight_baseline_view, flag_baseline_view) in izip!(
                img_baselines.iter(),
                jones_timestep_view.axis_iter(Axis(1)),
                weight_timestep_view.axis_iter(Axis(1)),
                flag_timestep_view.axis_iter(Axis(1))
            ) {
                let ant1_idx = baseline.ant1_index;
                let ant2_idx = baseline.ant2_index;

                let baseline_xyz_precessed =
                    tiles_xyz_precessed[ant1_idx] - tiles_xyz_precessed[ant2_idx];
                let uvw = UVW::from_xyz(baseline_xyz_precessed, prec_info.hadec_j2000);

                let data: Array2<c32> =
                    Array2::from_shape_fn((num_img_chans, 4), |(c, p)| jones_baseline_view[c][p]);

                let flag_row = flag_baseline_view.iter().all(|&x| x);

                self.write_main_row(
                    &mut main_table,
                    main_idx,
                    scan_start_mjd_utc_s,
                    scan_centroid_mjd_utc_s,
                    ant1_idx as _,
                    ant2_idx as _,
                    0,
                    &vec![uvw.u, uvw.v, uvw.w],
                    integration_time_s,
                    -1,
                    1,
                    -1,
                    // TODO
                    &vec![1., 1., 1., 1.],
                    data,
                    flag_baseline_view.to_owned(),
                    weight_baseline_view.to_owned(),
                    flag_row,
                )?;

                main_idx += 1;
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    use crate::{
        approx::abs_diff_eq,
        c32, c64,
        constants::{
            COTTER_MWA_HEIGHT_METRES, COTTER_MWA_LATITUDE_RADIANS, COTTER_MWA_LONGITUDE_RADIANS,
        },
        ndarray::{arr2, array, Array},
    };
    use itertools::izip;
    use tempfile::tempdir;

    lazy_static! {
        static ref PATH_1254670392: PathBuf =
            "tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms".into();
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

    macro_rules! c32x4 {
        ($re0:expr, $im0:expr, $re1:expr, $im1:expr, $re2:expr, $im2:expr, $re3:expr, $im3:expr) => {
            [
                c32::new($re0, $im0),
                c32::new($re1, $im1),
                c32::new($re2, $im2),
                c32::new($re3, $im3),
            ]
        };
    }

    #[test]
    fn test_decompress_default_tables() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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

        // TODO: check TableMeasDesc

        let mut main_table = Table::open(&table_path.clone(), TableOpenMode::Read).unwrap();
        let main_table_keywords = main_table.table_keyword_names().unwrap();
        assert!(main_table_keywords.contains(&"SOURCE".into()));
    }

    #[test]
    fn test_add_mwa_mods() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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

        // TODO: check TableMeasDesc

        let mut main_table = Table::open(&table_path.clone(), TableOpenMode::Read).unwrap();
        let main_table_keywords = main_table.table_keyword_names().unwrap();
        assert!(main_table_keywords.contains(&"MWA_TILE_POINTING".into()));
        assert!(main_table_keywords.contains(&"MWA_SUBBAND".into()));
    }

    #[test]
    fn test_write_spectral_window_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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
    fn test_write_spectral_window_row_mwa() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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
    fn handle_bad_spw_chan_info() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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

        assert!(matches!(
            result,
            Err(MeasurementSetWriteError::BadArrayShape { .. })
        ))
    }

    #[test]
    fn test_write_data_description_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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

    const ANT_POSITIONS: &'static [[f64; 3]] = &[
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

    const ANT_NAMES: &'static [&'static str] = &[
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

    const ANT_INPUTS: &'static [[i32; 2]] = &[
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

    const ANT_TILE_NRS: &'static [i32] = &[
        11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 31, 32, 33, 34, 35, 36, 37,
        38, 41, 42, 43, 44, 45, 46, 47, 48, 61, 62, 63, 64, 65, 66, 67, 68, 81, 82, 83, 84, 85, 86,
        87, 88, 91, 92, 93, 94, 95, 96, 97, 98, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008,
        1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023,
        1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038,
        1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053,
        1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068,
        1069, 1070, 1071, 1072,
    ];

    const ANT_CABLE_LENGTHS: &'static [[f64; 2]] = &[
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
    fn test_write_antenna_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let ant_table_path = ms_writer.path.join("ANTENNA");
        let mut ant_table = Table::open(ant_table_path, TableOpenMode::ReadWrite).unwrap();

        ant_table.add_rows(ANT_NAMES.len()).unwrap();

        for (idx, (name, position)) in izip!(ANT_NAMES, ANT_POSITIONS).enumerate() {
            let position = position.iter().cloned().collect();

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
    fn test_write_antenna_row_mwa() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);
        ms_writer.add_mwa_mods();

        let ant_table_path = ms_writer.path.join("ANTENNA");
        let mut ant_table = Table::open(ant_table_path.clone(), TableOpenMode::ReadWrite).unwrap();

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
            let position = position.iter().cloned().collect();

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
    fn test_write_polarization_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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
    fn handle_bad_pol_small_corr_type() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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

        assert!(matches!(
            result,
            Err(MeasurementSetWriteError::BadArrayShape { .. })
        ))
    }

    #[test]
    fn handle_bad_pol_big_corr_product() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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

        assert!(matches!(
            result,
            Err(MeasurementSetWriteError::BadArrayShape { .. })
        ))
    }

    #[test]
    fn test_write_source_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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
    fn test_write_field_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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
    fn test_write_field_row_mwa() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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
    fn handle_bad_field_shape() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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

        assert!(matches!(
            result,
            Err(MeasurementSetWriteError::BadArrayShape { .. })
        ))
    }

    #[test]
    fn test_write_observation_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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
    fn test_write_observation_row_mwa() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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
    fn test_write_history_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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
                5140115737.974,
                "cotter \"-m\" \"tests/data/1254670392_avg/1254670392.metafits\" \"-o\" \"tests/data/1254670392_avg/1254670392.cotter.none.ms\" \"-allowmissing\" \"-norfi\" \"-nostats\" \"-nogeom\" \"-noantennapruning\" \"-nosbgains\" \"-noflagautos\" \"-noflagdcchannels\" \"-nocablelength\" \"-edgewidth\" \"0\" \"-initflag\" \"0\" \"-endflag\" \"0\" \"-sbpassband\" \"tests/data/subband-passband-32ch-unitary.txt\" \"-nostats\" \"-flag-strategy\" \"/usr/local/share/aoflagger/strategies/mwa-default.lua\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox01_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox02_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox03_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox04_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox05_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox06_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox07_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox08_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox09_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox10_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox11_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox12_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox13_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox14_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox15_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox16_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox17_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox18_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox19_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox20_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox21_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox22_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox23_00.fits\" \"tests/data/1254670392_avg/1254670392_20191009153257_gpubox24_00.fits\"",
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
            // TODO:
            "APP_PARAMS",
            "CLI_COMMAND",
        ] {
            assert_table_columns_match!(hist_table, expected_table, col_name);
        }
    }

    #[test]
    fn test_write_feed_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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
                    &vec![0., PI / 2.],
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
    fn test_write_mwa_tile_pointing_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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
                6.28310690918895887,
                -0.464403366228935188,
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
    fn test_write_mwa_subband_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
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
    fn test_initialize_from_mwalib_all() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let array_pos = Some(LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        });

        let context = CorrelatorContext::new(
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

        let phase_centre = RADec::from_mwalib_phase_or_pointing(&context.metafits_context);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, array_pos);

        // let mwalib_timestep_range = (*context.common_timestep_indices.first().unwrap())
        //     ..(*context.provided_timestep_indices.last().unwrap() + 1);
        let mwalib_timestep_range = 0..(context.metafits_context.num_metafits_timesteps - 1);
        let mwalib_coarse_chan_range = *context.provided_coarse_chan_indices.first().unwrap()
            ..(*context.provided_coarse_chan_indices.last().unwrap() + 1);

        ms_writer
            .initialize_from_mwalib(&context, &mwalib_timestep_range, &mwalib_coarse_chan_range)
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
            // (
            //     "FLAG_CMD",
            //     vec![
            //         "APPLIED", "COMMAND", "INTERVAL", "LEVEL", "REASON", "SEVERITY", "TIME", "TYPE",
            //     ],
            // ),
            (
                "HISTORY",
                vec![
                    // TODO:
                    // "APP_PARAMS",
                    // TODO:
                    // "CLI_COMMAND",
                    // Different application so these will never match
                    // "APPLICATION",
                    // TODO:
                    // "MESSAGE",
                    "OBJECT_ID",
                    "OBSERVATION_ID",
                    "ORIGIN",
                    "PRIORITY",
                    // "TIME",
                ],
            ),
            (
                "OBSERVATION",
                vec![
                    "TIME_RANGE",
                    // "LOG",
                    // "SCHEDULE",
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
                ],
            ),
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
        let exp_slots = array![
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6],
            [7, 7],
            [8, 8],
        ];
        assert_eq!(num_ants, exp_slots.shape()[0] as _);
        for (idx, exp_slot) in exp_slots.outer_iter().enumerate() {
            let slot: Vec<i32> = ant_table.get_cell_as_vec("MWA_SLOT", idx as _).unwrap();
            assert_eq!(slot, exp_slot.to_vec());
        }
    }

    /// This was dumped with
    /// ```python
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms')
    /// with open('dump.txt', 'w') as dump:
    ///     print(
    ///         "\n".join(
    ///             f"c32x4!({row[0].real}, {row[0].imag}, {row[1].real}, {row[1].imag}, {row[2].real}, {row[2].imag}, {row[3].real}, {row[3].imag}),"
    ///             for row in tb.getcell('DATA', 0).transpose().tolist()
    ///         ),
    ///         file=dump
    ///     )
    /// ```
    const VIS_DATA_1254670392: &'static [[c32; 4]] = &[
        c32x4!(24.25, 1.0, 85.5, 81.75, 35.25, -2.0, 154.5, 9.625),
        c32x4!(58.25, -67.0, 3.875, -12.375, -36.0, 75.75, 17.375, 75.625),
        c32x4!(75.375, -14.75, -90.125, -38.75, 11.0, 48.75, -13.625, -168.125),
        c32x4!(122.125, 75.25, -21.875, -95.375, 41.875, 37.375, 25.125, -160.75),
        c32x4!(-136.625, 4.375, -109.5, 0.0, -243.875, -152.875, 35.25, -26.125),
        c32x4!(99.75, -109.375, 166.25, 21.25, -145.625, 3.375, 64.25, 36.25),
        c32x4!(89.5, -32.0, -47.5, -44.0, -16.0, -233.875, 38.25, 184.5),
        c32x4!(158.875, -20.625, -35.625, -155.5, 126.125, 28.375, -114.125, -124.125),
        c32x4!(117.125, -15.625, 191.0, 26.875, -41.375, 268.125, 374.75, -143.375),
        c32x4!(-30.125, 71.875, -60.625, -99.625, -126.25, 190.625, 118.0, -257.5),
        c32x4!(-46.625, -17.25, -245.5, -57.5, 83.375, -48.125, 166.875, -31.125),
        c32x4!(94.375, -243.875, -96.375, -60.625, -222.0, 52.0, 58.25, -0.25),
        c32x4!(23.625, -42.5, -24.375, -7.375, 69.875, 112.875, 282.625, 57.875),
        c32x4!(242.75, 196.0, -26.5, 78.375, -81.875, -223.375, 53.75, -198.625),
        c32x4!(-6.625, -91.875, -49.5, -153.0, -44.75, -67.25, 64.125, 30.125),
        c32x4!(95.25, 37.875, 87.75, -162.75, 114.75, 94.0, 120.25, -143.0),
        c32x4!(463.875, -101.375, 397.125, -88.875, 498.125, 159.125, 400.0, -48.75),
        c32x4!(-59.0, -90.5, -138.125, -29.625, 133.625, 60.625, -29.0, -22.625),
        c32x4!(266.375, -164.375, 49.75, -125.375, 180.625, -91.375, 137.0, -118.0),
        c32x4!(113.25, -49.0, 137.0, -4.875, 135.75, -68.875, -1.125, -131.875),
        c32x4!(217.5, 36.5, 106.75, -93.375, 45.625, -126.25, -15.25, 62.125),
        c32x4!(86.25, -171.75, 125.5, 59.375, 189.0, 181.375, 206.75, -24.125),
        c32x4!(178.0, -34.75, 151.625, -6.75, 10.75, -112.625, 20.875, 7.875),
        c32x4!(300.75, -187.25, -212.875, -96.125, -209.25, -76.375, 18.75, 35.0),
        c32x4!(-34.875, -232.5, -22.375, -263.5, 195.75, -144.75, -7.875, 98.125),
        c32x4!(25.25, -110.875, -77.625, 53.125, -1.375, 35.75, 19.75, 108.375),
        c32x4!(165.875, -172.375, -103.625, -17.5, 52.75, -181.875, 9.5, -196.875),
        c32x4!(-173.375, 51.0, 102.375, 57.125, -52.875, 81.125, -117.75, 107.25),
        c32x4!(111.0, -107.125, -121.75, 59.0, -17.625, -79.375, 172.0, 255.5),
        c32x4!(45.375, 16.875, 23.875, -120.375, -8.375, -62.0, -81.5, 107.5),
        c32x4!(51.375, -52.375, 129.125, -18.125, -56.25, 36.875, 69.875, -32.75),
        c32x4!(53.25, -16.375, 24.875, -50.125, -41.25, 56.75, 11.0, -123.125),
        c32x4!(55.75, -10.5, 19.25, 52.0, 38.5, -93.25, 163.0, 3.25),
        c32x4!(113.5, 15.5, -39.0, 25.0, 34.625, 17.625, 89.75, -79.5),
        c32x4!(118.375, -20.375, 19.125, -4.25, 57.75, -68.25, -75.75, -21.25),
        c32x4!(22.5, 92.125, 38.5, -174.0, 66.375, -47.125, -154.625, -66.625),
        c32x4!(-42.75, -238.25, -141.375, 21.125, 9.25, 64.25, -103.875, -24.125),
        c32x4!(257.375, -75.5, -172.875, 161.75, -124.875, 76.875, 30.625, 14.0),
        c32x4!(104.125, -190.625, 92.0, 122.75, -127.625, 147.125, 35.25, -154.25),
        c32x4!(1.625, -175.625, 26.375, -2.125, -88.125, -31.25, -134.625, 233.625),
        c32x4!(-32.0, 53.375, -64.625, 77.875, -82.5, -74.625, 38.5, -16.875),
        c32x4!(173.875, 99.625, 262.25, -242.25, 191.125, -190.25, 159.5, -23.625),
        c32x4!(56.125, -65.75, -123.75, 33.5, -55.625, -79.875, -50.5, -191.75),
        c32x4!(-225.75, 11.25, -38.625, 3.25, -7.25, 8.375, 173.0, 6.875),
        c32x4!(251.5, -27.625, 56.125, 192.0, 183.625, -203.75, -111.0, -45.375),
        c32x4!(-190.5, -61.875, 16.625, 90.125, 147.0, 156.75, 215.75, -85.375),
        c32x4!(41.25, -205.5, 133.125, -163.25, -8.625, 55.625, 11.25, -69.875),
        c32x4!(245.875, 29.125, 57.625, -5.125, 69.0, -130.125, -7.875, -48.5),
        c32x4!(511.875, -44.0, 584.625, -254.5, 274.125, 89.0, 429.0, 49.25),
        c32x4!(-149.375, -156.625, 129.75, 79.625, 104.0, 152.625, 2.375, -74.0),
        c32x4!(229.25, -108.375, 149.25, -250.75, 171.75, -113.875, 6.125, -34.5),
        c32x4!(295.375, 90.125, 166.625, -55.875, -126.25, -27.125, 17.125, -278.625),
        c32x4!(135.5, 29.5, -73.125, -100.125, -140.0, 204.25, 108.0, -54.125),
        c32x4!(157.625, -94.5, 4.875, -56.25, 214.5, -19.0, 167.0, 140.375),
        c32x4!(22.25, 69.75, 68.625, 131.875, -25.25, -33.375, 178.25, -126.5),
        c32x4!(-13.875, 20.125, 27.0, 15.25, -67.25, -355.625, 87.0, -31.875),
        c32x4!(17.75, -79.875, -53.75, -15.625, -17.375, -93.25, 121.25, -147.875),
        c32x4!(-134.375, -194.25, 101.75, 23.75, 0.0, -141.125, 6.25, -0.625),
        c32x4!(-67.25, 61.375, -312.875, -68.375, 74.125, -57.0, 120.875, -5.625),
        c32x4!(100.625, -38.875, -31.25, -56.75, -33.125, -0.25, -113.75, -1.5),
        c32x4!(35.875, -36.125, -142.0, -90.5, 17.375, -7.25, -69.75, 45.0),
        c32x4!(56.5, -66.5, -71.5, -48.0, -116.75, 57.625, 142.875, 44.75),
        c32x4!(71.75, -207.0, -46.5, -1.5, 113.0, -9.5, 122.375, -5.375),
        c32x4!(103.5, 70.0, 46.75, 69.25, -44.5, -162.5, -105.875, -42.375),
        c32x4!(134.75, -86.25, -35.375, -57.25, -54.25, 0.0, 140.25, 103.75),
        c32x4!(82.375, -97.5, -100.0, -1.625, 42.25, 22.875, 112.625, -147.125),
        c32x4!(9.5, 71.375, 40.0, -63.5, -122.25, 75.0, -115.625, 128.625),
        c32x4!(39.875, 16.125, 42.625, 96.75, -38.75, 37.75, 84.375, 10.625),
        c32x4!(80.5, -15.75, 71.875, 6.75, 161.75, -45.125, -46.5, -194.875),
        c32x4!(48.125, 140.25, 136.625, 193.875, 116.625, -141.75, 260.375, 61.75),
        c32x4!(1.625, 3.625, 23.5, -72.625, 33.5, -66.5, 26.75, -39.5),
        c32x4!(32.625, -53.5, 272.5, 96.375, 27.5, -34.5, 0.625, -115.5),
        c32x4!(29.375, 157.875, 43.75, -282.5, 46.875, -21.5, -97.375, -2.75),
        c32x4!(-21.375, 65.5, 107.0, -18.5, -59.75, 184.625, -67.125, -23.625),
        c32x4!(-112.625, -18.375, -63.75, -94.0, -27.0, -63.375, -61.25, -186.75),
        c32x4!(95.875, -89.375, 102.75, -101.875, 155.25, 33.75, 103.125, -170.0),
        c32x4!(94.25, -220.375, -98.25, -33.375, -6.25, -127.125, 43.75, 27.625),
        c32x4!(230.375, -109.0, -109.125, 31.0, -391.0, -105.375, 109.25, -73.25),
        c32x4!(89.25, 16.25, 63.25, 169.625, 54.875, 198.0, 221.0, -67.625),
        c32x4!(289.75, -81.875, -60.125, 121.625, -102.875, -7.375, -106.625, 98.25),
        c32x4!(319.5, -79.625, 417.5, -147.5, 503.0, -206.125, 214.375, 186.875),
        c32x4!(33.5, -83.75, 45.75, 178.625, -99.75, -9.5, 131.875, -68.875),
        c32x4!(-85.875, 37.375, -105.125, -85.375, -6.875, 188.625, 227.375, 51.625),
        c32x4!(-178.875, -302.125, 88.375, -99.0, -25.75, 89.25, 78.875, 88.625),
        c32x4!(24.125, -65.375, 282.25, 131.125, 86.375, -9.625, 262.875, -234.625),
        c32x4!(-54.5, -178.0, 107.125, -162.5, 158.0, 104.375, -71.25, 40.25),
        c32x4!(-3.5, 16.25, -88.125, 17.0, -80.375, -42.0, 269.0, 12.0),
        c32x4!(9.5, 102.0, -102.5, 200.25, -330.875, -166.375, 189.375, -55.5),
        c32x4!(275.875, 41.5, 169.75, -64.25, -168.5, 135.25, 273.0, 340.0),
        c32x4!(-39.875, -37.75, 25.25, -204.5, -163.5, -0.5, 11.125, -171.5),
        c32x4!(-108.625, 11.75, -55.375, -296.375, -93.125, -87.125, 131.125, -51.375),
        c32x4!(58.25, -96.75, -218.75, -172.0, 34.5, -78.875, -47.375, 157.5),
        c32x4!(108.75, -39.75, -80.5, 115.25, 203.75, 132.5, -105.375, 106.875),
        c32x4!(115.375, 69.5, 171.25, 161.125, 60.5, 31.875, 114.375, 34.25),
        c32x4!(2.5, -167.75, -64.625, -7.375, -67.0, 73.5, 148.875, 40.125),
        c32x4!(29.25, 34.5, 1.5, -59.0, 88.875, -88.875, 36.875, 62.125),
        c32x4!(11.5, 88.25, -11.125, -36.625, 7.75, 28.5, -14.625, 63.75),
        c32x4!(70.625, -205.625, -23.25, -28.0, -10.5, -17.625, -30.875, 70.875),
        c32x4!(-70.875, 139.25, 14.0, 76.375, 10.875, 47.625, -97.125, 67.25),
        c32x4!(22.75, -10.5, -107.875, 135.625, -26.125, -62.0, 57.875, 307.0),
        c32x4!(-7.375, 114.75, 105.0, -136.125, -35.375, -131.625, 45.25, -203.25),
        c32x4!(-173.625, -177.25, 99.0, 121.0, 29.75, -64.625, 68.625, 25.0),
        c32x4!(-81.625, -89.875, -242.375, 11.625, -126.125, -26.25, 1.5, -102.625),
        c32x4!(62.125, 52.0, 218.0, -304.0, -241.625, -319.75, -135.375, -128.125),
        c32x4!(-18.125, 18.25, 4.125, -101.625, -73.0, -43.25, -11.625, -150.75),
        c32x4!(233.0, -97.0, 3.125, -64.25, 156.25, 56.375, 61.375, -64.125),
        c32x4!(22.375, -134.875, -18.875, 151.875, -71.25, 80.125, 137.75, -17.125),
        c32x4!(54.375, -198.875, -189.5, -79.125, 66.625, -57.625, 83.625, 20.25),
        c32x4!(49.75, -137.75, 118.5, -44.125, 11.75, -75.25, -22.875, -130.875),
        c32x4!(35.375, -203.375, -4.25, -40.0, 30.375, 233.875, 312.25, 226.375),
        c32x4!(133.75, -94.375, -93.125, -118.25, -21.625, -80.625, 155.5, 115.375),
        c32x4!(-166.875, 4.375, -51.125, 17.0, 88.625, 97.5, 137.625, -145.25),
        c32x4!(479.375, -79.5, 228.25, -208.375, 431.125, 86.5, 294.375, 83.875),
        c32x4!(22.875, -156.625, -74.75, -86.125, -80.875, 199.0, 178.75, -89.5),
        c32x4!(60.875, -157.5, 132.5, 108.75, 0.25, 68.875, 216.125, 66.75),
        c32x4!(46.5, -145.875, -19.625, -31.0, 6.875, -189.0, -53.25, 13.625),
        c32x4!(47.25, -152.0, -144.25, -90.625, 182.125, -68.125, -168.875, -120.25),
        c32x4!(104.0, 69.375, 159.5, -3.625, -99.875, -96.625, 91.5, -63.875),
        c32x4!(3.75, -107.75, 138.75, 49.625, 162.25, -232.125, -104.375, -147.25),
        c32x4!(232.25, 50.5, 196.875, 139.875, 14.375, 15.875, -24.0, -236.75),
        c32x4!(112.5, -151.0, 77.75, 51.75, -109.625, 79.625, -178.5, 87.875),
        c32x4!(-117.5, 95.0, -61.75, -70.125, -1.5, -169.625, 14.0, -10.5),
        c32x4!(-331.375, -48.375, 174.125, -239.75, 25.125, 49.875, -44.75, -211.625),
        c32x4!(50.25, -111.0, -112.0, -89.375, 159.5, -168.25, 40.125, -43.75),
        c32x4!(205.75, 198.75, -110.25, 104.875, 40.0, -23.125, -67.125, -4.625),
        c32x4!(-189.75, -62.0, -54.125, -139.75, 78.625, 46.375, 123.5, 3.875),
        c32x4!(53.125, -158.875, -95.5, 91.375, -76.625, -116.0, 100.375, -67.0),
        c32x4!(-29.625, 104.25, 39.5, -32.0, 19.25, 32.5, 26.375, 31.625),
        c32x4!(-51.375, 42.625, 87.75, -10.625, -65.25, 120.5, 5.625, 16.125),
        c32x4!(36.25, -95.5, 5.625, -4.625, -18.5, -20.625, 70.625, 73.375),
        c32x4!(-25.125, -156.75, 15.0, 1.5, -81.125, 2.375, 145.625, -161.875),
        c32x4!(217.875, 58.125, 55.75, -84.25, 157.5, -17.5, 225.25, 69.625),
        c32x4!(-120.5, -103.5, 78.375, 10.875, 48.125, -57.375, 9.0, 0.125),
        c32x4!(121.875, -88.75, -8.25, -203.75, -132.375, -71.75, -3.625, 111.875),
        c32x4!(36.375, -35.5, -24.25, 131.75, -67.625, -87.375, 131.875, -136.375),
        c32x4!(-8.125, -245.75, 161.625, 300.875, -96.375, 64.125, 99.5, 149.375),
        c32x4!(119.0, -195.75, 13.25, -59.5, -76.375, -6.875, -198.75, -163.5),
        c32x4!(70.0, -205.125, 155.75, 198.0, 76.875, -52.875, -163.875, 126.625),
        c32x4!(140.125, -64.375, 72.875, -10.75, -50.5, 52.375, 113.625, 79.125),
        c32x4!(81.25, -70.0, -119.25, -233.375, -82.5, -118.375, 5.5, -50.5),
        c32x4!(57.875, -99.125, -146.0, 34.75, 51.25, -71.25, 65.625, -81.75),
        c32x4!(223.375, -77.5, 99.375, -82.125, 116.5, 38.125, 18.5, -93.875),
        c32x4!(-351.5, 121.25, 88.75, -5.75, 353.25, 82.25, 254.375, 162.0),
        c32x4!(-161.0, -26.75, 176.0, 68.0, 123.75, 62.5, 30.75, -83.0),
        c32x4!(332.25, 20.875, 438.875, -110.0, 241.25, -49.75, 296.875, 86.25),
        c32x4!(201.75, 47.375, 37.25, 243.875, -64.75, 64.25, 125.75, 93.0),
        c32x4!(-200.75, -290.375, -6.375, 4.75, 13.75, -115.25, 60.25, -238.0),
        c32x4!(81.25, 148.625, 20.125, -49.25, 27.5, 1.125, -147.625, -223.375),
        c32x4!(-12.5, -72.25, -179.25, -65.5, 68.125, -7.625, 22.0, -115.75),
        c32x4!(36.0, -113.25, 63.0, 105.375, -81.125, -95.25, 71.0, -145.75),
        c32x4!(18.0, -117.375, 264.0, -169.875, -69.25, 114.375, 16.0, 104.75),
        c32x4!(-71.625, 3.25, 68.25, -69.5, -1.5, 10.875, 137.5, -153.375),
        c32x4!(194.5, 23.25, -60.625, -257.625, -80.25, 10.0, -109.75, 45.75),
        c32x4!(256.625, 17.625, 23.25, -249.25, -263.5, 143.25, 63.0, 0.5),
        c32x4!(100.25, 130.375, -9.375, 260.75, -44.0, 217.875, 225.125, -65.5),
        c32x4!(-161.625, -59.5, 4.375, 71.875, 207.5, 17.125, 47.5, 39.125),
        c32x4!(77.375, -98.5, -28.75, 200.625, 27.625, -36.0, -2.75, -84.75),
        c32x4!(5.625, -161.5, 51.25, 68.5, -50.0, -150.0, 3.5, 11.75),
        c32x4!(-99.125, 89.5, -60.5, 28.75, 6.375, -116.375, 114.125, -65.375),
        c32x4!(83.5, 0.0, -62.875, -42.875, 33.0, -56.875, -62.375, -21.0),
        c32x4!(13.875, -5.0, 19.375, -55.125, 81.75, 53.375, -8.5, 88.375),
        c32x4!(141.125, -43.0, -63.625, -36.125, -51.25, 18.75, -9.875, -115.875),
        c32x4!(3.0, -106.0, 82.875, -17.75, 80.75, -151.5, 29.5, -10.875),
        c32x4!(-107.75, -140.75, 53.0, 247.625, 24.375, 23.125, -32.75, 205.625),
        c32x4!(-123.375, -64.25, 275.75, 74.0, -61.625, 52.125, 48.125, -151.375),
        c32x4!(7.5, -103.375, 153.0, -13.0, -207.0, 53.875, -181.625, 83.125),
        c32x4!(4.25, 14.75, 67.5, -33.875, 12.75, 18.375, -92.25, 103.125),
        c32x4!(42.875, -12.75, 120.0, 167.75, -63.625, -7.125, -81.125, -222.5),
        c32x4!(-119.0, -22.875, 210.125, 14.5, 44.25, -56.875, -211.75, -64.25),
        c32x4!(-73.25, 198.75, -7.75, -11.375, -158.125, 115.125, -137.25, -146.25),
        c32x4!(41.625, 83.5, -278.75, 23.125, 33.875, 108.5, 53.0, 41.625),
        c32x4!(-72.75, -322.875, 118.625, 0.0, 172.625, -92.5, 65.5, -68.25),
        c32x4!(-51.375, 80.5, 30.0, -92.0, -79.0, -0.125, 200.5, -155.5),
        c32x4!(19.25, 106.375, 233.5, 45.375, 46.125, -52.25, -2.0, -22.625),
        c32x4!(96.5, 27.25, 199.5, -11.875, -75.25, -171.5, 9.375, -156.75),
        c32x4!(131.0, 2.5, -41.25, -204.875, 330.5, 101.25, -47.0, -53.625),
        c32x4!(85.75, 78.5, 157.75, 13.375, 350.75, 9.5, 210.625, -107.875),
        c32x4!(51.5, -141.5, -25.0, 80.25, 24.5, 46.5, -98.625, -113.375),
        c32x4!(-15.75, 175.0, -105.875, 14.25, 42.625, 119.625, 194.25, 31.625),
        c32x4!(72.125, 60.25, -19.0, -56.125, -101.875, 184.75, 342.625, -88.125),
        c32x4!(120.25, -136.5, 36.125, -117.5, 36.625, 152.625, 50.75, 169.125),
        c32x4!(112.125, 27.375, -116.5, -172.125, -138.375, 9.375, -248.625, -169.125),
        c32x4!(284.625, -149.75, 184.75, 259.375, -197.0, -66.875, 99.0, 16.75),
        c32x4!(264.375, 45.0, 78.5, -77.625, -59.625, -73.75, -25.75, 162.25),
        c32x4!(47.75, 72.125, 179.375, 77.375, 32.0, -119.875, -43.625, 20.125),
        c32x4!(6.875, 34.75, -70.125, 29.25, -119.0, 9.625, 38.625, 7.0),
        c32x4!(114.625, 91.5, 189.0, 191.5, -30.875, 6.75, 91.5, -65.0),
        c32x4!(-57.5, 245.125, 154.125, 71.125, 29.625, -152.125, 132.5, -48.75),
        c32x4!(3.375, -40.125, -87.0, -69.25, -206.875, -83.25, 32.125, -9.125),
        c32x4!(70.625, -59.5, -91.75, 17.125, 49.25, -76.5, 49.75, -135.25),
        c32x4!(-62.25, -8.375, 58.5, -63.0, 1.25, -58.5, 57.5, -82.75),
        c32x4!(87.75, 33.5, 75.5, -4.875, -4.75, -20.875, 59.75, 4.625),
        c32x4!(-50.375, -54.125, -42.75, 12.75, 2.875, -70.25, 108.5, 181.625),
        c32x4!(44.125, -24.625, -49.375, -44.625, 146.625, 68.0, -105.875, -44.625),
        c32x4!(-4.375, -16.5, -158.5, 43.875, -14.625, 46.75, 78.875, -119.125),
        c32x4!(137.5, -108.75, -123.625, 99.75, 93.25, 42.75, -168.375, -155.375),
        c32x4!(-27.0, -93.5, -47.25, -127.0, 30.75, 160.125, -93.125, 6.0),
        c32x4!(28.75, 105.75, -11.625, 115.75, 108.625, -94.125, -109.125, -171.625),
        c32x4!(40.25, 52.0, -14.375, 39.25, 3.875, -54.5, -152.5, -137.125),
        c32x4!(-95.25, 51.625, 15.375, 73.25, 20.25, 7.875, -57.25, -40.125),
        c32x4!(-104.125, 14.0, 31.0, -104.625, -75.875, -149.625, -156.625, -175.875),
        c32x4!(82.125, -31.875, 93.0, -15.625, -41.125, 69.375, 132.375, -83.125),
        c32x4!(-44.375, -113.25, 187.5, 192.375, 65.75, -22.5, 139.5, 11.875),
        c32x4!(228.5, 117.125, -19.375, 24.0, 126.5, -29.0, 106.375, -81.375),
        c32x4!(-103.375, 163.875, 178.25, 1.0, 42.375, -119.875, -27.5, -90.625),
        c32x4!(-143.0, -12.25, 413.625, 101.0, 11.25, -75.875, 224.625, -251.25),
        c32x4!(-75.125, -21.5, 131.75, -150.875, 187.375, 279.375, 18.375, 21.0),
        c32x4!(-51.875, -171.875, 146.875, -25.25, -184.125, 72.5, 11.75, -10.5),
        c32x4!(132.5, -213.75, 572.625, -222.625, 595.375, -124.125, 223.875, -16.0),
        c32x4!(32.125, -293.25, 139.75, -22.125, -176.125, -54.75, -283.25, -32.625),
        c32x4!(-71.5, 240.0, 74.625, -81.0, 212.75, -0.75, 99.375, -115.125),
        c32x4!(280.875, -30.375, 99.375, -112.125, 132.5, -36.625, -65.875, 125.5),
        c32x4!(62.875, -7.0, 47.875, 245.125, 123.5, 103.125, 129.875, -126.5),
        c32x4!(1.875, 114.875, -63.75, 10.375, 142.375, 164.5, 6.125, 29.125),
        c32x4!(123.125, -96.75, -55.5, 189.0, 43.25, 107.375, 101.5, -142.0),
        c32x4!(38.375, -126.25, -122.625, -2.375, 7.125, 48.875, -99.5, 116.625),
        c32x4!(-131.0, -268.875, -149.25, -4.375, -23.125, -6.125, -90.625, 8.25),
        c32x4!(-241.875, -111.25, 88.5, -59.5, 24.75, 92.75, 179.5, -291.25),
        c32x4!(-31.5, -80.625, 2.5, 116.625, -25.5, 107.875, 83.625, 6.25),
        c32x4!(19.875, 8.5, 131.0, 94.625, -72.125, 154.125, 27.375, 182.5),
        c32x4!(11.25, -51.375, 95.375, -189.0, -127.25, -25.125, 83.125, -72.25),
        c32x4!(35.625, -36.0, -26.125, -17.5, 41.625, 89.625, 8.125, 12.125),
        c32x4!(-152.5, -70.5, -76.625, 92.5, -51.375, -32.0, -22.75, -36.5),
        c32x4!(-76.625, 9.625, 31.0, 56.25, 22.5, -16.125, 92.0, -114.875),
        c32x4!(-42.5, -59.25, 97.875, -38.75, -48.75, 23.375, 13.125, -63.25),
        c32x4!(109.625, 27.875, 61.375, -74.875, 92.875, -63.75, 11.25, 35.125),
        c32x4!(-107.625, -104.0, -120.5, -54.375, -28.375, 83.0, -98.625, 8.25),
        c32x4!(133.0, 93.25, -0.375, 123.0, 117.5, -62.25, -23.25, -46.25),
        c32x4!(-143.875, 117.25, -113.0, -154.75, -14.625, -75.625, -98.875, -5.375),
        c32x4!(95.375, -209.25, -69.875, 35.75, 116.25, -174.125, 104.25, 125.0),
        c32x4!(-132.5, 29.75, -99.5, 53.25, 180.625, 26.125, -103.75, -196.375),
        c32x4!(25.375, 31.875, 262.375, -38.375, 2.125, -136.125, 70.375, -185.5),
        c32x4!(106.125, 92.125, -204.0, -11.0, 90.75, 8.125, -70.5, 6.0),
        c32x4!(-90.5, -85.125, 229.25, 62.875, 43.875, 62.375, 10.5, -105.0),
        c32x4!(-112.625, 87.0, 189.125, -31.625, -66.25, -39.625, -116.75, 55.0),
        c32x4!(0.0, 94.25, -47.625, 79.125, 25.125, 18.0, 153.375, 133.75),
        c32x4!(-85.0, 17.875, 24.875, 216.125, 146.375, -44.375, 13.875, 66.625),
        c32x4!(82.875, 76.0, -12.5, 221.0, 97.875, 119.125, -202.125, -84.375),
        c32x4!(100.25, -92.25, -98.0, -20.375, 132.0, 129.75, -5.875, -16.125),
        c32x4!(-38.625, 23.875, -29.625, 45.25, 34.625, 135.375, -101.375, 81.875),
        c32x4!(218.25, -79.625, 424.75, 70.25, 334.5, -78.875, 318.75, 38.625),
        c32x4!(-138.0, -66.125, 52.375, 92.625, -7.375, -137.5, 169.875, 115.875),
        c32x4!(-44.75, -22.125, -18.25, -107.75, 71.5, 68.875, 203.5, -18.0),
        c32x4!(-74.625, 118.75, 41.25, -19.125, 188.5, -163.75, 87.875, -0.625),
        c32x4!(-152.875, 74.375, -36.125, -204.5, 84.5, 55.375, 47.75, 83.375),
        c32x4!(-269.875, -211.625, -18.0, 40.5, -4.75, -96.5, -57.5, -22.875),
        c32x4!(-159.75, -117.25, -113.25, -29.25, -161.25, -91.25, 215.375, -55.375),
        c32x4!(200.375, -30.25, 35.375, 18.875, -86.25, -76.5, 132.0, 96.0),
        c32x4!(17.25, 97.125, 131.125, -37.75, -93.875, 143.0, 339.0, -119.25),
        c32x4!(-64.0, 54.25, 16.125, -140.875, 47.5, 27.25, 38.375, 82.625),
        c32x4!(-103.625, -173.375, 73.75, 10.875, 24.625, -63.125, -55.0, 141.5),
        c32x4!(-66.0, -231.125, 58.5, 151.875, 24.25, 133.125, -44.625, -327.75),
        c32x4!(-10.0, 114.75, 97.875, -110.5, 95.375, 10.125, 229.125, -148.875),
        c32x4!(-117.0, -156.375, -96.0, -72.375, 15.0, -83.375, -22.125, -57.5),
        c32x4!(10.375, -72.625, -21.125, -34.875, 36.5, -35.375, 88.625, -78.75),
        c32x4!(-45.25, 61.625, -65.75, 32.125, 113.125, -38.625, -18.625, 12.625),
        c32x4!(22.0, 11.125, 59.0, 22.625, -67.125, -3.375, -23.375, 4.125),
        c32x4!(-55.375, 26.5, 82.625, 70.375, 23.0, -88.25, 103.125, 45.75),
        c32x4!(58.125, 148.75, 14.625, 139.5, 1.125, -45.75, -96.125, -133.0),
        c32x4!(-57.25, 6.5, -35.25, -118.0, 11.625, 67.75, 34.0, 94.875),
        c32x4!(48.5, -192.125, -4.625, -117.375, 208.875, -47.75, 10.0, -30.375),
        c32x4!(-30.5, 39.875, 67.25, -22.625, 22.0, -78.75, -11.875, -166.25),
        c32x4!(151.125, 84.125, 32.0, 130.125, 132.125, -133.375, 45.75, -106.875),
        c32x4!(56.875, 44.25, 194.75, -54.125, 50.0, -138.5, 50.5, -114.625),
        c32x4!(-116.0, -49.125, -109.25, 53.125, -20.0, 105.125, 1.875, 241.0),
        c32x4!(112.625, -34.625, -60.375, -37.75, -47.5, -0.75, -131.5, 205.25),
        c32x4!(86.875, -229.875, -166.25, 181.625, 34.0, -124.375, 49.5, 242.0),
        c32x4!(-56.875, 51.25, 18.5, 91.875, -189.125, 66.5, -73.0, 266.125),
        c32x4!(-22.5, -161.25, -28.75, 33.75, -263.125, -132.25, -82.125, 41.625),
        c32x4!(24.875, -122.625, 195.875, 214.375, 258.5, -8.375, 45.0, -31.75),
        c32x4!(-13.625, -164.75, 119.875, -175.125, 0.125, -67.0, -173.25, 246.625),
        c32x4!(-202.375, 93.25, -58.5, 78.875, -42.5, 45.5, -14.75, -138.625),
        c32x4!(392.25, 97.875, 280.25, 20.375, 553.25, -5.625, 336.875, -73.625),
        c32x4!(-88.125, -61.25, 32.5, -158.375, -5.5, 72.125, 149.75, 188.5),
        c32x4!(12.5, -77.875, -26.125, 124.375, 14.5, 67.125, -67.875, -38.375),
        c32x4!(-24.25, -43.375, 110.875, 140.125, 87.0, -43.125, -26.625, -63.0),
        c32x4!(17.75, -14.5, 53.375, 156.375, 38.125, 18.25, 29.875, 0.125),
        c32x4!(-68.875, -176.25, -192.5, -226.125, -40.75, 381.125, 162.25, 178.5),
        c32x4!(71.375, -1.875, -76.625, 20.25, 61.75, -120.625, -166.75, 35.75),
        c32x4!(-63.625, 74.25, 14.375, 220.0, 232.625, 104.0, -55.25, -36.0),
        c32x4!(-173.125, 27.875, -5.0, -174.625, -116.0, 92.75, -83.75, 175.125),
        c32x4!(199.75, -50.625, -78.0, -229.5, 98.875, -60.875, -180.625, 14.375),
        c32x4!(-169.25, 61.625, 85.875, 38.25, 26.375, 1.125, 199.625, -49.75),
        c32x4!(28.875, 137.375, -158.25, 130.75, -51.125, -102.5, -92.25, -93.625),
        c32x4!(-102.375, 37.625, -164.25, 13.875, 53.125, -58.5, -181.125, 69.0),
        c32x4!(-81.625, -102.5, 45.125, 25.0, -61.875, -72.0, -5.75, -83.125),
        c32x4!(21.5, -50.375, -75.25, -13.125, -42.875, 29.0, 8.625, -52.75),
        c32x4!(-64.75, 24.125, 59.25, -54.75, -1.375, 72.625, -137.75, 112.75),
        c32x4!(-13.0, 56.875, -60.875, 94.875, -30.0, 45.25, 24.25, -38.375),
        c32x4!(-54.0, 26.625, 66.0, 46.5, -63.25, 10.0, 177.75, 45.0),
        c32x4!(-52.625, -97.125, -1.125, -93.375, 62.25, -32.5, -71.375, -76.25),
        c32x4!(-90.875, -21.25, -180.625, 24.5, 202.75, 10.375, 164.375, 94.375),
        c32x4!(78.875, -297.125, -59.25, 27.5, 14.875, 21.875, 95.75, -6.875),
        c32x4!(-45.75, -25.125, -0.125, 3.5, 36.625, 177.75, -48.125, -46.25),
        c32x4!(-154.125, -21.75, -77.5, -175.25, -17.25, 95.375, 46.625, 19.125),
        c32x4!(24.625, -191.25, 70.5, 117.75, -59.375, 50.25, 37.875, 96.375),
        c32x4!(66.0, -90.125, 56.5, -16.125, 229.25, 43.25, 193.0, 57.375),
        c32x4!(-174.75, -25.875, 81.875, 20.625, 351.25, 44.5, 14.375, -144.125),
        c32x4!(33.5, -114.0, -22.25, -65.625, -106.0, -84.0, 122.875, -133.25),
        c32x4!(-256.875, -96.0, -21.5, -27.25, -25.0, -63.125, 174.125, -55.75),
        c32x4!(-51.0, -16.75, 112.625, 72.875, 57.375, 145.25, -30.25, -44.5),
        c32x4!(8.875, 44.625, 280.875, 195.75, 22.25, 69.75, 68.875, 23.5),
        c32x4!(-42.5, 49.125, 46.375, 55.375, -28.5, -67.0, -16.625, -77.125),
        c32x4!(63.0, -5.5, -0.625, -122.75, -69.125, 9.0, -20.625, 106.375),
        c32x4!(275.25, 148.875, 300.625, -64.75, 258.375, -12.25, 233.875, 68.75),
        c32x4!(-84.25, 209.125, -34.0, -36.5, 84.625, -79.75, 17.75, -171.625),
        c32x4!(-70.0, 159.5, 95.25, -47.25, -81.25, -71.625, -140.5, -102.875),
        c32x4!(-3.375, 7.875, 121.125, -9.625, -150.625, -202.25, 178.5, -166.125),
        c32x4!(-28.375, 3.0, 96.375, 25.375, -77.875, 53.5, 140.0, -151.75),
        c32x4!(36.375, 202.875, 124.875, 61.375, 49.25, -88.375, -69.25, 26.125),
        c32x4!(-117.125, -70.75, 221.25, -103.125, -18.625, -130.125, -70.5, -52.125),
        c32x4!(-122.25, 33.5, 53.125, 134.5, 191.25, -5.625, -251.125, -20.875),
        c32x4!(-405.75, 257.375, -194.0, -60.0, 75.75, 89.125, -42.125, -66.625),
        c32x4!(-91.0, -3.375, -186.75, 91.875, 45.375, 143.875, -210.75, -125.0),
        c32x4!(89.75, -133.625, -0.125, -150.875, 169.25, 58.875, 53.375, -163.0),
        c32x4!(-105.5, -111.5, 67.75, 53.125, 35.125, 7.125, -115.125, 43.125),
        c32x4!(-15.25, 94.625, 1.25, -9.75, -178.375, -3.125, -16.75, -7.375),
        c32x4!(11.25, -14.25, -41.625, 62.125, 50.625, 117.75, -42.75, -76.0),
        c32x4!(47.125, 41.5, 89.5, 41.875, 10.25, -79.0, 45.5, -24.875),
        c32x4!(-12.625, -27.5, 81.875, -84.375, 8.5, 77.0, -2.0, -28.875),
        c32x4!(31.125, -11.375, -86.0, 68.375, 15.625, 17.0, 73.375, -84.75),
        c32x4!(-82.25, -105.875, 115.0, 7.0, 43.625, 88.0, -87.75, 16.75),
        c32x4!(-77.75, 42.625, 143.25, -30.625, -53.125, 85.75, 29.625, -116.0),
        c32x4!(1.625, 17.875, -67.0, -49.625, 51.875, -21.125, 116.875, -7.5),
        c32x4!(1.5, -37.25, 115.25, 195.125, 101.125, 155.0, -129.75, 42.75),
        c32x4!(106.5, -10.5, -34.25, -43.25, -51.875, -228.25, 239.125, -121.0),
        c32x4!(-4.375, 103.0, 44.75, 219.125, -183.5, -85.25, 141.875, -71.375),
        c32x4!(-55.25, -26.25, 97.375, -5.0, -29.0, -65.75, 18.5, -157.125),
        c32x4!(-49.625, 145.0, -121.375, -58.125, -55.375, 63.375, 64.0, -24.625),
        c32x4!(-103.75, 6.375, -78.75, 3.625, -101.5, -113.0, 122.875, 204.625),
        c32x4!(122.125, -25.5, 2.125, 109.375, 53.0, 51.25, -64.0, 44.125),
        c32x4!(-9.75, 103.0, -74.0, 34.875, -37.25, -84.5, 184.125, -124.625),
        c32x4!(69.0, 128.75, 47.25, -44.0, -48.5, -81.625, -54.25, -323.375),
        c32x4!(-41.0, 8.625, -73.875, 92.875, 1.75, -112.25, -46.75, -128.125),
        c32x4!(-19.0, -90.75, 116.75, -166.875, 88.875, -19.125, -392.25, -103.75),
        c32x4!(-69.75, -45.625, -83.25, 18.625, 25.0, -19.625, -168.0, 27.125),
        c32x4!(195.125, -18.125, 330.375, 100.625, 538.0, 282.125, 459.75, 115.375),
        c32x4!(-42.625, 73.0, 30.75, 93.25, -10.0, -140.375, 21.25, -10.0),
        c32x4!(-174.375, -32.0, -1.125, 282.0, -2.5, -83.75, -19.0, 2.25),
        c32x4!(-361.25, 86.375, 76.875, -43.0, 110.5, -112.125, -338.625, 142.0),
        c32x4!(-55.5, 137.375, -62.0, 39.875, 127.625, 29.5, 29.25, 65.25),
        c32x4!(-115.5, -22.0, -12.5, 145.0, -23.875, -125.5, -19.75, 18.625),
        c32x4!(127.875, -208.0, 134.125, -231.875, -300.125, 80.5, 23.75, 22.75),
        c32x4!(6.75, 28.25, -191.375, 36.125, -80.875, 25.875, 30.5, 67.375),
        c32x4!(-167.375, -87.75, -91.875, -22.625, 96.625, 124.625, -2.125, 63.375),
        c32x4!(38.375, -11.375, -48.375, 34.75, -108.375, 79.875, -65.0, 167.5),
        c32x4!(26.0, 16.125, -32.75, 103.875, 11.25, -165.25, 79.375, -14.25),
        c32x4!(-113.75, -210.125, 0.625, 53.0, 56.375, 162.875, -30.875, 146.125),
        c32x4!(-30.25, -15.875, -73.375, -43.125, 321.625, -9.75, -23.375, 27.625),
        c32x4!(23.875, -12.375, -22.25, 92.25, 23.375, 52.0, -31.5, -11.125),
        c32x4!(134.0, -19.625, 38.125, 127.625, -47.125, -26.625, -60.625, -32.875),
        c32x4!(-48.125, -94.75, -121.75, -31.0, -123.75, 129.5, -81.75, 33.375),
        c32x4!(-33.375, 32.0, -195.625, -56.125, 36.125, 76.75, 11.0, -53.75),
        c32x4!(-159.625, -37.75, 79.25, 2.625, 59.0, -51.375, -128.5, 11.25),
        c32x4!(70.25, -32.875, -7.375, 64.375, 14.625, 129.5, 89.5, 38.25),
        c32x4!(6.25, 15.125, -114.625, 75.75, -76.75, 70.375, -27.0, -147.0),
        c32x4!(-80.75, 25.75, -18.125, 139.0, -27.875, 19.25, -38.625, -92.625),
        c32x4!(-85.375, -146.25, -46.125, -56.0, 78.375, 65.875, -241.5, 32.0),
        c32x4!(-2.625, -59.875, 8.375, 113.25, 38.125, -133.0, 42.375, -6.625),
        c32x4!(-205.0, 19.0, -71.375, 0.5, -2.5, 62.5, -8.625, -207.375),
        c32x4!(-57.0, -343.5, -274.625, 26.0, -36.0, -92.625, 52.875, -133.75),
        c32x4!(-154.125, 22.375, 22.625, -23.375, -100.0, -128.375, 120.625, 70.375),
        c32x4!(-63.875, 27.75, -184.625, 231.375, -79.75, 18.125, 148.5, 105.0),
        c32x4!(-67.625, 413.0, 27.375, -43.875, -142.0, 9.875, 176.25, 35.75),
        c32x4!(-110.25, 61.875, 34.5, -176.375, 20.125, 2.375, -46.0, -139.375),
        c32x4!(-180.5, -177.75, -20.0, 99.75, -177.625, 16.25, 297.125, 82.125),
        c32x4!(84.25, 52.5, -234.125, 33.125, 8.375, 69.625, -35.625, 12.875),
        c32x4!(11.5, -241.125, 81.75, 80.375, -111.375, -70.0, 4.875, -4.5),
        c32x4!(306.25, -75.375, 486.375, 139.875, 130.375, 239.625, 411.125, 85.125),
        c32x4!(-148.125, 176.375, -143.625, -2.375, -72.375, -18.5, -24.25, -7.875),
        c32x4!(47.625, 154.25, 216.375, -83.75, 40.5, 198.875, -20.875, 91.375),
        c32x4!(96.875, 20.125, -66.375, -139.5, -66.875, -43.75, -80.625, 5.0),
        c32x4!(-106.5, -98.5, 121.875, 64.875, -23.125, -17.75, 58.125, 55.75),
        c32x4!(-53.625, 24.625, 126.375, -36.25, 42.75, -75.375, -29.375, -112.875),
        c32x4!(-14.125, 142.125, -48.875, 185.5, -169.75, -68.875, -54.25, 65.125),
        c32x4!(-59.125, -74.5, 107.75, 46.375, -18.5, 56.25, -187.75, -175.75),
        c32x4!(-67.625, -200.875, 73.375, 22.0, 47.625, -135.125, 252.0, -163.625),
        c32x4!(14.0, 40.875, 166.875, 25.375, 105.0, 78.875, -4.75, 23.625),
        c32x4!(14.0, -284.75, -12.25, 73.75, 53.375, 206.0, -18.375, -31.0),
        c32x4!(76.875, -77.125, 25.0, -92.375, -129.0, 66.875, 13.125, -61.5),
        c32x4!(118.125, -128.125, 21.625, 69.75, 26.75, 34.5, 39.0, -123.75),
        c32x4!(-70.125, -96.0, 5.125, 146.625, -47.25, 43.375, 47.25, 86.875),
        c32x4!(84.25, -48.0, -46.125, -47.875, -31.25, -30.25, 2.625, 5.625),
        c32x4!(-57.5, -60.75, 62.75, -68.75, 39.375, 32.625, -32.375, 8.875),
        c32x4!(65.0, -63.75, 10.25, -50.125, -56.5, 93.0, 101.875, -123.5),
        c32x4!(-97.75, -23.75, -73.25, -17.375, 50.875, 2.0, 30.25, 3.75),
        c32x4!(-46.0, -78.625, 5.5, -54.25, 44.125, 25.0, 92.125, 6.125),
        c32x4!(-156.875, -28.0, -20.25, 66.625, -38.125, -42.0, -260.375, 58.0),
        c32x4!(-29.625, -124.875, -18.75, 96.5, -133.75, 72.375, -91.0, 162.625),
        c32x4!(11.125, -130.625, 19.875, 47.375, 92.125, -13.0, -146.125, -244.5),
        c32x4!(128.5, 11.5, -55.25, 179.25, -325.375, 8.875, -1.125, 20.625),
        c32x4!(88.75, -166.875, -103.0, 269.375, 101.5, 2.125, 38.5, 118.25),
        c32x4!(-102.75, -169.0, -216.75, 50.875, -122.25, -69.5, 42.25, 238.0),
        c32x4!(-28.0, 99.0, 56.5, 109.25, 81.5, -76.0, -0.375, -15.625),
        c32x4!(80.125, -145.75, -20.5, -72.375, -153.0, 155.0, -7.75, 3.125),
        c32x4!(23.0, 31.25, 190.25, 25.125, -2.75, -7.125, -251.625, 48.5),
        c32x4!(-49.5, -121.625, 22.125, -144.625, -81.125, -240.75, -107.625, -87.0),
        c32x4!(121.875, -109.625, -31.625, 92.5, -107.75, -59.25, 67.125, -55.5),
        c32x4!(-64.25, -88.625, -108.375, 116.75, -87.125, -55.125, -148.375, 275.375),
        c32x4!(-18.5, -48.875, 56.375, 123.125, 81.25, -129.25, -180.25, -80.25),
        c32x4!(192.375, 29.75, 154.125, -71.5, 217.25, 150.875, 433.0, -93.875),
        c32x4!(-97.375, 27.75, 103.375, -12.25, -5.75, 45.0, 58.75, 61.0),
        c32x4!(11.25, -148.5, -50.0, -30.875, 11.625, -100.25, 100.375, 70.375),
        c32x4!(-67.625, -17.0, -109.625, -74.125, 71.375, -56.875, 87.75, 90.125),
        c32x4!(-25.375, -173.375, 39.125, -137.75, -161.0, 174.0, -24.25, 17.625),
        c32x4!(56.125, -224.875, -107.875, -29.125, -43.625, -18.0, -138.375, -83.125),
        c32x4!(-115.125, -94.625, -16.375, 134.25, 110.125, 85.75, 103.125, -163.625),
        c32x4!(150.5, -64.125, 18.125, -91.125, -60.375, -57.625, 136.5, 158.5),
        c32x4!(-96.625, -28.5, -226.625, -191.25, 35.0, -245.75, -6.625, 129.125),
        c32x4!(9.25, -53.75, 177.375, 102.0, -51.0, 87.875, -12.125, 68.25),
        c32x4!(53.875, 21.0, -117.0, -9.5, -40.0, -240.0, 123.625, -32.25),
        c32x4!(-43.625, 13.0, 11.875, 27.125, -141.125, 24.375, -207.625, -41.75),
        c32x4!(-102.875, -14.0, -140.875, -66.25, -115.5, -44.5, 43.25, 59.5),
        c32x4!(28.375, -70.125, -58.0, 11.5, 13.0, 28.75, -24.625, -4.375),
        c32x4!(-4.75, -157.875, -143.25, 15.625, -29.75, -20.875, 70.75, 78.625),
        c32x4!(66.875, -86.125, -41.5, -85.75, -20.0, 43.5, -6.875, -70.5),
        c32x4!(72.625, -90.75, -102.125, -32.25, -2.875, 23.75, 17.125, -145.625),
        c32x4!(-43.5, 2.875, 12.0, -110.5, 115.875, -22.125, -178.875, -54.875),
        c32x4!(-27.75, -60.0, 85.75, -113.25, -30.125, -10.375, -91.75, -95.5),
        c32x4!(-113.625, -37.375, -114.5, -53.0, 74.375, 85.75, -40.125, -107.375),
        c32x4!(-55.75, -101.625, 19.875, 20.0, -42.25, -71.0, -159.25, 20.25),
        c32x4!(-109.75, -281.5, 108.625, -35.375, -180.875, 2.75, -40.625, 183.625),
        c32x4!(-240.0, 71.5, -48.0, 147.375, -241.25, 43.875, -243.875, -74.75),
        c32x4!(-51.0, -197.75, -274.25, -90.625, 10.5, -57.0, -98.25, -39.125),
        c32x4!(-22.125, -189.125, 27.625, 147.75, -65.0, -80.0, -320.25, -58.125),
        c32x4!(-70.125, 61.0, -6.375, 139.375, 18.5, -30.0, -168.25, 77.25),
        c32x4!(-13.5, 44.875, 0.0, -14.75, -90.875, 23.625, -69.0, 92.75),
        c32x4!(-49.5, 5.75, 80.125, 58.375, 57.5, -97.125, -95.5, -27.875),
        c32x4!(-231.5, -137.75, 86.125, 61.875, -110.25, -115.5, 51.25, -64.75),
        c32x4!(-124.75, -24.375, -60.125, 119.625, -47.75, 48.875, -15.0, 72.125),
        c32x4!(51.0, -27.125, 96.75, -167.0, -101.875, 123.5, -150.0, 94.0),
        c32x4!(-132.5, -198.375, 130.125, 178.125, 94.875, 35.0, 17.375, 0.5),
        c32x4!(365.0, -130.25, 249.25, -10.25, 270.25, -67.875, 494.0, -135.625),
        c32x4!(3.375, -145.75, -177.25, 19.375, 31.875, 134.5, -89.0, 121.625),
        c32x4!(-184.25, 24.375, 171.875, -110.375, 34.0, 93.625, -130.625, -165.25),
        c32x4!(85.25, -55.0, -87.0, 165.125, 34.5, -108.625, -160.625, 136.625),
        c32x4!(-148.75, 25.75, -42.25, -117.0, 99.125, -45.375, -124.625, -96.75),
        c32x4!(-3.25, 138.625, -89.25, -98.0, -26.125, -105.0, -37.875, -137.875),
        c32x4!(62.125, -156.5, 12.125, 21.125, 57.125, 40.625, -80.375, 42.625),
        c32x4!(-250.25, -22.25, -81.125, 58.625, -39.625, -104.25, -78.0, -73.375),
        c32x4!(-71.375, 58.375, 17.75, -124.25, -236.25, 35.375, -85.875, -132.625),
        c32x4!(81.125, -82.375, -81.75, -60.625, -105.5, 11.5, -125.75, 28.625),
        c32x4!(-65.25, -156.125, 70.0, -63.25, -44.75, -67.375, 26.0, -70.5),
        c32x4!(118.375, -165.5, 88.125, 173.0, -22.625, 91.875, 43.0, -158.5),
        c32x4!(102.625, 72.0, 85.25, -12.875, 12.375, 50.75, 15.25, -96.125),
        c32x4!(99.875, -66.125, -48.0, 122.5, -38.75, -59.125, -9.375, 107.375),
        c32x4!(-65.375, -77.875, -118.0, 59.0, 114.0, -93.875, -86.5, -62.125),
        c32x4!(31.625, -88.0, 14.5, -75.125, 29.375, -47.375, -50.25, -60.75),
        c32x4!(15.0, 33.125, -93.875, 68.0, -8.25, 15.625, 6.25, -22.0),
        c32x4!(-4.25, -42.875, -75.875, -85.625, -5.125, -136.0, 23.875, -137.125),
        c32x4!(30.625, -4.625, -62.625, -21.875, 41.0, 49.25, -176.0, -1.625),
        c32x4!(19.0, -45.5, 114.25, 269.5, -112.125, -71.75, -102.0, -59.25),
        c32x4!(-17.875, -96.25, -104.375, -0.5, -118.25, -225.0, -97.125, -126.875),
        c32x4!(-2.875, -133.875, 73.5, 158.875, -46.25, -75.75, 4.625, -339.25),
        c32x4!(-20.125, -22.875, 50.375, -18.125, 73.875, -73.75, -124.75, 96.625),
        c32x4!(-81.25, 11.125, 166.125, -86.0, 255.25, 90.25, 19.125, -255.375),
        c32x4!(-78.0, -95.375, 311.625, 141.375, -80.125, -61.125, 187.875, -175.5),
        c32x4!(88.875, -64.125, 182.625, -40.875, 106.125, -10.0, -101.875, -71.5),
        c32x4!(-0.375, 26.875, -99.5, -64.5, -106.125, 218.25, -172.625, 90.875),
        c32x4!(-45.875, 86.125, -119.5, 60.875, 28.25, 11.0, -201.0, 40.375),
        c32x4!(-6.375, 43.0, 16.375, 10.0, -86.0, 219.25, 72.0, -156.75),
        c32x4!(57.0, -143.0, 28.625, 150.5, -96.625, -51.75, -12.375, -0.25),
        c32x4!(-187.75, -241.125, -1.875, -15.625, -2.625, 61.875, -183.125, -26.875),
        c32x4!(120.5, 91.0, -153.5, 52.125, -148.875, -33.125, 32.125, -81.0),
        c32x4!(363.625, 15.125, 388.25, -142.0, 559.375, 96.75, 398.875, 106.5),
        c32x4!(110.625, -62.125, -30.125, -166.0, -4.125, -125.875, -88.625, -68.0),
        c32x4!(-165.5, 65.375, 0.125, -134.125, -37.5, 110.5, -7.0, -158.375),
        c32x4!(50.875, 76.75, -2.875, 138.125, 100.5, 45.625, -80.375, 108.75),
        c32x4!(-126.0, 28.0, -357.5, -57.125, -101.625, 31.375, 0.25, -79.875),
        c32x4!(-108.625, -76.25, 94.125, -0.5, -28.375, 1.125, -15.625, -99.625),
        c32x4!(-37.125, 7.0, 126.875, 41.625, 177.125, 105.875, -27.0, -378.0),
        c32x4!(-31.125, 36.25, -8.625, -98.5, -3.875, 17.0, 26.125, -224.625),
        c32x4!(-30.875, -109.5, 100.625, 11.625, -72.25, 172.75, -132.875, -87.75),
        c32x4!(-13.0, -54.25, -114.875, 17.625, -12.0, 90.5, 87.875, 90.875),
        c32x4!(-127.5, -57.875, -111.75, 28.0, -33.875, -132.875, -84.25, -82.0),
        c32x4!(33.875, 38.75, -75.625, 14.375, 98.125, -58.875, 6.875, -19.125),
        c32x4!(2.0, -7.625, 65.875, -52.25, -59.25, 74.25, -19.625, 210.0),
        c32x4!(-43.5, 16.0, 140.625, -44.875, 69.375, -17.125, 37.625, -159.125),
        c32x4!(56.375, -62.0, 14.875, 8.0, 49.25, 38.375, -64.125, 36.75),
        c32x4!(116.125, -13.125, -43.875, 57.625, -57.25, -44.125, 163.0, -82.25),
        c32x4!(-83.125, -48.0, 8.625, 71.25, -66.125, -23.875, -22.5, -9.0),
        c32x4!(-119.125, 34.375, 31.75, 121.25, -96.625, -118.75, -57.125, -96.875),
        c32x4!(-12.75, 14.0, -40.75, -20.25, 32.5, -62.25, 101.5, 36.75),
        c32x4!(-105.0, 20.0, 99.625, 27.125, 185.125, -75.125, 2.5, -172.125),
        c32x4!(-181.375, -36.875, -17.875, 10.125, -63.875, -155.5, -84.0, -72.875),
        c32x4!(35.375, -13.125, -232.875, -131.375, 68.5, -7.875, -62.375, 82.875),
        c32x4!(-118.875, -4.75, -67.375, 96.875, -106.5, 27.625, -5.375, -109.0),
        c32x4!(56.5, -452.5, -16.125, 381.25, -149.75, 23.125, -95.5, -76.625),
        c32x4!(-120.125, 0.625, -160.625, -138.5, 93.75, -39.125, -306.875, 52.625),
        c32x4!(-183.0, -92.125, 14.0, -15.375, 0.0, 30.0, -21.625, -141.75),
        c32x4!(64.5, -132.5, 125.5, 24.125, -102.375, -5.0, -197.25, -148.5),
        c32x4!(-19.25, -26.75, 8.25, 51.625, -192.375, -54.75, 87.375, 71.125),
        c32x4!(37.625, -6.25, 0.0, 113.625, 146.0, -71.375, 130.25, -170.25),
        c32x4!(15.75, -31.875, 30.5, 180.25, 156.125, -75.0, -47.875, -214.0),
        c32x4!(-85.0, -43.25, -79.75, 39.25, -38.625, -120.625, -15.125, -173.375),
        c32x4!(166.25, 129.0, 13.125, 72.125, 68.25, -53.0, -9.75, -298.625),
        c32x4!(466.125, 18.5, 282.625, -108.25, 342.375, 132.875, 178.375, -129.625),
        c32x4!(63.375, -14.625, -155.75, 55.25, 5.375, 13.125, 210.75, -130.375),
        c32x4!(-143.25, 77.875, 32.375, 142.25, 91.375, 32.75, 87.875, 15.5),
        c32x4!(-35.375, -201.75, -53.875, 38.125, -39.75, 24.375, -128.25, -235.625),
        c32x4!(6.5, -103.0, -32.625, 70.25, 204.75, -28.625, 7.625, 77.5),
        c32x4!(-52.5, 89.75, -93.875, -34.875, -60.0, -35.375, 140.875, 58.0),
        c32x4!(119.0, -29.875, 26.875, -24.875, -47.875, 87.75, -104.875, 16.625),
        c32x4!(-84.875, -108.75, -61.125, 66.25, -137.25, 56.375, -75.875, 51.0),
        c32x4!(17.0, 62.5, -69.75, 109.0, -152.375, -11.75, 19.125, -10.125),
        c32x4!(-121.375, 11.375, 32.0, 122.375, 147.625, -148.125, -119.0, -45.0),
        c32x4!(-64.375, -213.75, -41.875, 56.75, -9.125, 157.75, -145.0, -5.625),
        c32x4!(160.5, -51.5, 17.0, -122.0, 20.25, -80.875, -79.875, 176.25),
        c32x4!(-27.5, -162.75, -71.25, -74.25, -11.125, -74.875, 43.375, 13.25),
        c32x4!(25.5, -40.25, -101.25, 22.625, 19.375, 69.375, -12.375, 99.375),
        c32x4!(-76.625, -15.5, -40.0, 41.625, -123.5, 8.75, 11.125, -41.25),
        c32x4!(-39.0, 22.625, 27.0, -2.0, 38.125, -29.5, 39.125, -33.625),
        c32x4!(-96.0, -7.25, 65.0, 54.75, -73.0, -24.125, -42.875, 20.875),
        c32x4!(-63.0, 8.625, 8.5, 148.75, -49.875, -69.125, -51.5, 2.25),
        c32x4!(-14.75, -81.0, -46.625, 39.75, -20.125, 78.5, 37.25, 1.625),
        c32x4!(-18.375, -1.0, -1.25, -208.5, 34.75, -87.5, 20.75, 9.125),
        c32x4!(-93.75, -65.125, 78.125, -41.375, 5.75, 44.125, -35.625, -85.5),
        c32x4!(47.0, -44.0, -65.875, -0.875, -119.625, -37.75, -121.5, 16.25),
        c32x4!(9.625, -152.625, 84.0, 40.375, -65.125, 92.75, -21.375, 27.375),
        c32x4!(116.625, 121.5, 7.25, 141.625, -186.0, -173.75, 57.75, -65.125),
        c32x4!(11.25, -3.375, -34.5, 34.0, -66.875, -15.375, 21.125, 60.75),
        c32x4!(41.125, 82.5, -103.0, -93.125, 37.875, 72.25, 65.25, -42.375),
        c32x4!(217.375, -10.0, -39.375, 87.25, -108.625, 53.875, -158.75, -168.75),
        c32x4!(-84.875, -16.5, -137.375, 61.125, -26.0, 67.75, -68.5, -111.625),
        c32x4!(23.125, -40.0, 82.5, 157.375, 91.125, 54.375, -136.125, -58.25),
        c32x4!(-156.875, 120.0, 77.5, -65.75, -133.75, -136.0, 167.125, 122.625),
        c32x4!(-47.0, -198.5, -127.125, -4.25, 108.25, 73.5, 15.75, 115.5),
        c32x4!(-32.0, 43.5, 14.25, -152.5, -88.625, -57.375, 23.75, -151.0),
        c32x4!(255.375, -6.875, 316.25, 20.625, 431.25, 11.25, 306.0, -5.75),
        c32x4!(-2.125, -28.875, 41.625, 10.875, -86.625, -50.875, -18.625, -18.625),
        c32x4!(-129.125, -132.375, 46.5, 112.5, -101.375, 149.75, 335.125, -31.0),
        c32x4!(-158.625, 57.0, -36.0, -17.875, 131.0, -58.5, -89.375, -23.5),
        c32x4!(46.0, -124.0, 119.125, 88.875, -12.625, 23.25, -49.875, 31.875),
        c32x4!(-169.75, 38.75, 93.625, -106.25, -195.125, -96.375, -181.125, 46.875),
        c32x4!(-1.75, -269.125, -40.0, 34.25, -21.25, 154.5, -149.0, -113.0),
        c32x4!(30.75, -53.875, 48.5, -32.625, 89.875, -12.375, 0.625, -57.375),
        c32x4!(8.0, 98.0, 116.375, 48.125, -41.75, -111.625, -130.125, 60.75),
        c32x4!(-238.5, -155.5, 59.0, 114.625, 81.125, -36.0, 67.875, -111.875),
        c32x4!(-183.0, -37.375, -120.25, -127.0, -15.625, 30.25, 69.625, 149.25),
        c32x4!(68.375, 146.625, 55.875, -68.125, 7.0, 98.125, -6.5, -117.875),
        c32x4!(-51.125, 61.25, 112.75, 135.375, -124.375, -51.0, -20.125, 50.0),
        c32x4!(67.375, 41.5, -44.625, -71.75, -63.25, 60.625, 60.25, 18.0),
        c32x4!(-53.5, -37.625, -52.125, 88.875, -39.25, -29.375, -83.25, -58.0),
        c32x4!(-21.75, -13.375, 0.625, -86.875, -6.0, 46.25, 75.125, 2.125),
        c32x4!(-55.125, 149.375, 28.625, 0.875, -102.625, -16.75, -92.625, -58.125),
        c32x4!(-112.875, 126.0, 35.0, 42.0, -41.5, 45.0, 55.75, 1.125),
        c32x4!(-111.5, 6.125, 99.0, 173.0, 66.375, 101.25, 114.125, -166.5),
        c32x4!(23.75, 117.375, -122.75, 22.75, -79.25, -139.5, -95.25, 30.375),
        c32x4!(82.75, -2.875, 100.875, -21.75, -39.5, -98.75, -93.125, -0.5),
        c32x4!(52.625, 32.75, -42.625, 87.25, 30.375, 116.875, -5.5, 46.5),
        c32x4!(-76.75, 103.625, 27.0, -47.875, 38.0, -22.625, 37.875, -81.0),
        c32x4!(-156.75, 84.625, 108.0, -37.75, 56.125, 44.125, -194.75, 12.625),
        c32x4!(84.25, -55.75, -186.875, 47.625, -125.875, 42.75, 44.5, -186.0),
        c32x4!(69.375, -32.75, 77.625, -10.375, -28.375, 72.375, 80.75, 47.5),
        c32x4!(-155.375, -122.0, -131.625, -92.375, 85.25, -106.125, -29.875, -177.625),
        c32x4!(-81.375, 109.875, 253.625, -107.75, 24.625, -111.875, -129.0, 60.5),
        c32x4!(-64.25, 86.125, 93.625, -140.375, -134.875, 83.0, 92.375, -90.125),
        c32x4!(-10.875, -21.0, 58.0, -37.5, -111.625, 30.625, -57.5, -108.875),
        c32x4!(-69.875, -91.25, 46.375, 91.5, 77.625, -10.375, -42.375, 66.25),
        c32x4!(12.0, 21.5, 119.5, 53.5, -40.0, -45.25, 66.5, 27.75),
        c32x4!(308.625, -86.875, 330.875, 127.25, 245.25, 49.5, 209.5, -18.875),
        c32x4!(133.75, 18.875, -169.75, -64.25, -213.0, -245.0, -166.5, -224.375),
        c32x4!(63.75, 57.375, -64.0, -11.875, -102.25, 84.75, -3.875, -57.5),
        c32x4!(-57.75, 37.375, -72.5, -113.75, 3.5, -27.375, -97.25, -196.125),
        c32x4!(-60.625, 181.875, 80.875, 10.375, 200.875, 135.375, 78.25, 19.625),
        c32x4!(-154.0, 84.125, -10.75, -26.625, -58.5, 36.5, -43.375, -131.0),
        c32x4!(-71.625, 195.0, 70.75, -37.625, 186.375, 7.5, -79.5, -104.0),
        c32x4!(-122.625, 53.75, 201.5, 58.625, 32.5, -52.75, -154.625, -69.625),
        c32x4!(35.875, 41.75, 122.0, 8.875, -93.0, 112.25, -24.0, -78.25),
        c32x4!(-85.125, 79.125, -166.375, -79.25, 105.25, 17.0, 152.75, -20.0),
        c32x4!(23.625, -47.5, -53.625, -135.375, -16.375, -46.875, 70.375, -44.75),
        c32x4!(39.875, 116.375, -11.75, 27.125, 29.875, 42.5, 45.625, -61.75),
        c32x4!(67.0, 66.625, 53.75, 40.25, -42.375, 52.5, -54.0, -54.75),
        c32x4!(67.75, -28.375, 110.75, -76.0, 12.0, 81.0, -51.625, -195.125),
        c32x4!(-63.625, -61.875, -50.375, -58.0, -33.375, -45.375, 6.5, -40.25),
        c32x4!(-51.125, -67.0, -47.375, 45.625, 49.875, 70.625, 15.25, -54.125),
        c32x4!(-4.875, 62.625, 7.0, 74.25, -43.0, -6.875, -140.625, 12.25),
        c32x4!(-51.375, 81.125, 43.75, -36.875, -4.875, -60.75, -121.875, -33.875),
        c32x4!(-54.25, 19.75, 2.875, 44.375, 110.625, -98.125, 68.875, -96.0),
        c32x4!(22.25, -96.25, 5.125, -50.125, 67.0, -22.625, 8.875, -51.875),
        c32x4!(-8.5, 56.875, 134.125, -112.0, 60.125, -222.875, -164.25, -153.625),
        c32x4!(-56.625, -236.875, 67.75, -18.25, -155.125, -119.0, -98.75, -109.125),
        c32x4!(-23.75, -80.125, -24.0, 64.125, -174.375, -56.875, -90.75, -196.125),
        c32x4!(-186.5, -130.125, 123.875, -11.0, -53.75, 63.375, -113.625, 101.5),
        c32x4!(45.75, -15.75, -13.25, 12.5, 80.625, 38.625, -1.0, -18.5),
        c32x4!(65.0, -93.375, -80.75, 148.25, -109.625, -21.5, -11.875, -116.625),
        c32x4!(-189.0, 17.75, 20.875, -84.0, -2.625, -28.375, 68.625, -119.5),
        c32x4!(-61.375, -45.0, -32.25, 91.5, -61.125, 35.25, 7.5, -152.125),
        c32x4!(-205.375, 12.125, -56.125, 5.375, 176.25, 45.25, -136.625, -21.375),
        c32x4!(-15.0, -26.75, 195.625, -108.0, -75.75, 37.875, 109.5, -64.875),
        c32x4!(25.25, -214.75, 56.75, 278.5, 70.25, -55.625, -55.375, 49.625),
        c32x4!(37.25, -241.625, -112.0, -39.875, 81.75, -181.375, -228.875, -138.5),
        c32x4!(336.875, -33.5, 305.625, 59.875, 300.375, 59.625, 206.25, -119.375),
        c32x4!(40.375, -109.875, 212.375, 79.5, 71.375, -56.125, -173.75, -158.0),
        c32x4!(-207.875, -104.125, 1.625, -44.25, -12.75, -165.5, -198.875, 3.0),
        c32x4!(-219.625, 19.625, -69.375, 39.875, 59.125, 9.75, -80.875, 78.625),
        c32x4!(-195.625, 11.875, -110.5, -18.75, 47.0, 21.75, 97.625, -37.625),
        c32x4!(-129.375, -42.25, -20.125, -38.5, 34.625, -40.125, 14.0, -94.625),
        c32x4!(-103.125, 77.25, -211.75, -86.625, 54.25, 59.25, -145.25, -202.25),
        c32x4!(-7.25, -25.5, 10.5, -65.25, 34.125, 12.25, 91.5, -91.125),
        c32x4!(-104.25, -69.0, 1.125, -48.125, -91.625, 21.125, -119.0, 41.125),
        c32x4!(-22.875, 47.5, -38.0, -64.75, -69.125, 116.375, 20.5, -32.0),
        c32x4!(-0.5, -35.875, -124.125, 66.625, -83.5, -54.625, -75.375, -118.75),
        c32x4!(-147.125, -35.0, 74.75, -33.625, -29.625, -92.625, 10.25, -310.375),
        c32x4!(93.125, 62.0, 10.625, -24.625, -122.75, 81.0, -104.75, -199.5),
        c32x4!(-61.25, -34.375, -32.0, 50.75, 100.0, 25.375, -156.375, -311.125),
        c32x4!(-55.125, -56.375, -75.0, -13.375, 97.875, -47.0, -74.625, -138.375),
        c32x4!(57.625, 12.375, 34.0, -111.125, 0.5, -40.75, 16.375, -84.5),
        c32x4!(-7.25, 28.75, 59.5, 64.0, 53.75, 5.0, -69.0, -82.5),
        c32x4!(-123.0, -44.625, 69.75, 34.5, 34.75, 51.5, -22.5, -56.75),
        c32x4!(-26.0, -92.75, -95.0, -71.25, -60.0, 6.0, -84.625, -50.0),
        c32x4!(-110.25, 24.0, 98.625, 70.0, -56.75, -40.375, 63.125, -116.0),
        c32x4!(-112.125, 75.875, -11.75, -69.25, -169.875, 116.0, -26.0, -135.125),
        c32x4!(8.0, -42.625, 163.75, 32.0, 221.0, 5.25, -31.125, -237.25),
        c32x4!(-94.5, -13.375, 14.75, 44.375, -109.625, -73.75, 111.75, -53.125),
        c32x4!(-162.0, 125.5, -87.125, -65.625, -127.125, -141.5, 39.125, -87.375),
        c32x4!(84.75, -131.875, 109.0, 84.5, 15.875, 42.0, -142.25, -148.375),
        c32x4!(-205.375, 35.25, 45.875, -302.5, -87.25, 120.875, 20.875, -237.625),
        c32x4!(-83.625, 143.625, -39.875, -134.25, 17.75, -91.25, -23.25, -50.375),
        c32x4!(-74.625, 159.0, -38.75, -8.625, 120.25, 82.375, -38.875, -102.625),
        c32x4!(-61.625, -136.5, -57.75, 98.875, 70.25, 57.875, 14.125, -121.625),
        c32x4!(-104.0, -39.125, -96.25, -102.625, -34.25, -50.625, 68.0, 30.75),
        c32x4!(-93.375, -22.125, -50.875, -178.75, -12.125, -8.375, -15.25, -330.75),
        c32x4!(-54.375, -43.5, 40.125, -83.0, -49.75, 120.375, -209.75, -334.5),
        c32x4!(336.125, 35.0, 357.875, 68.875, 304.5, 188.125, 370.75, -183.0),
        c32x4!(-236.25, -21.5, 47.375, -24.375, -105.125, -194.375, 5.375, -67.25),
        c32x4!(-110.75, 75.0, 83.375, -71.625, 25.375, -89.625, -3.75, -112.5),
        c32x4!(20.375, -67.0, 11.0, 9.375, -20.125, 65.375, -91.625, -51.125),
        c32x4!(-125.125, -133.125, 10.625, -11.125, -91.25, 58.75, 78.75, -82.125),
        c32x4!(18.875, -178.0, -54.625, -84.375, 102.25, -178.125, 74.625, -194.875),
        c32x4!(-218.75, -13.375, -134.0, 109.625, 14.625, 0.375, -4.875, -29.25),
        c32x4!(157.75, 84.375, 162.75, -72.5, -157.75, 96.75, 41.125, -194.0),
        c32x4!(77.125, 22.625, -77.375, -81.875, 60.375, 128.5, 187.625, -164.75),
        c32x4!(-83.75, 176.125, 149.25, 54.625, -27.5, 21.625, -93.875, -31.25),
        c32x4!(-89.75, 47.75, 58.875, -15.875, 49.25, 5.625, -130.25, -192.5),
        c32x4!(-69.125, -184.25, -56.375, -29.0, 93.625, -5.75, 57.625, -17.125),
        c32x4!(1.375, -122.25, 158.5, 45.0, -34.875, -79.5, 61.625, 10.5),
        c32x4!(-23.25, -37.625, 100.75, 7.375, 49.75, -69.875, 98.375, -115.75),
        c32x4!(-34.125, -79.0, 30.125, 12.25, 77.875, -14.375, -107.875, -14.5),
        c32x4!(-49.0, -20.0, -78.5, 3.5, 19.25, 20.5, 88.125, 10.875),
        c32x4!(-31.625, 0.125, 37.625, 55.125, -74.625, 33.0, -55.875, -56.625),
        c32x4!(-28.375, -95.25, -41.125, -27.625, -37.75, 94.625, -83.375, 20.25),
        c32x4!(-69.25, -4.25, -19.75, -4.625, -27.75, 116.25, -78.625, -112.75),
        c32x4!(-183.375, -95.75, 96.875, -95.375, -56.75, -80.625, -94.75, -99.125),
        c32x4!(-29.5, -31.125, -204.375, -10.0, 40.625, 92.75, -33.5, -158.0),
        c32x4!(-78.125, -101.875, 1.375, -123.125, 65.375, -29.875, -117.25, -2.625),
        c32x4!(137.375, -27.25, -3.875, 25.875, -84.875, -100.5, 259.75, -24.25),
        c32x4!(-11.375, -110.25, 9.25, 93.0, 68.625, -18.5, -87.0, -236.375),
        c32x4!(-89.0, 199.75, -41.5, 114.625, -26.375, -36.0, -12.75, -132.0),
        c32x4!(-193.125, -121.125, 19.25, -1.375, 48.375, 184.875, 105.0, 14.5),
        c32x4!(-66.375, -49.875, 93.625, -159.75, -16.875, 121.125, -11.875, -58.5),
        c32x4!(33.625, -107.125, -67.875, -128.125, -54.875, -214.75, 10.25, -93.125),
        c32x4!(-140.75, -10.0, -115.75, 69.375, 7.625, -132.375, -3.75, -20.75),
        c32x4!(-10.625, -157.125, 207.75, -56.25, 104.75, 176.375, -67.75, -104.125),
        c32x4!(66.0, -23.75, -183.125, 256.875, -66.875, 34.75, 11.875, -287.125),
        c32x4!(5.25, 187.5, 147.125, 9.375, -12.0, -69.0, -19.5, -355.875),
        c32x4!(351.0, -89.25, 493.875, 20.25, 327.625, 130.75, 392.125, -262.75),
        c32x4!(21.875, 21.25, -100.75, -219.375, -15.375, -125.125, -82.625, -259.625),
        c32x4!(-149.5, -20.625, -99.625, 46.875, 7.75, -12.875, 53.125, 99.125),
        c32x4!(91.0, -44.5, 131.5, -104.5, -28.375, -133.375, -38.625, -262.0),
        c32x4!(6.75, -123.25, -58.25, -198.75, 22.5, -124.125, 108.5, -286.375),
        c32x4!(18.5, 34.125, 36.25, -48.125, 36.875, -66.625, 139.375, -21.375),
        c32x4!(-124.0, -62.125, -35.625, 76.0, 94.875, 15.625, -64.75, 101.375),
        c32x4!(58.625, -249.875, -35.125, -2.0, 19.875, -73.0, 62.625, -55.625),
        c32x4!(-75.5, -112.75, 17.375, 56.875, -111.875, -51.75, -87.875, 15.875),
        c32x4!(-116.75, 2.125, -7.875, 0.875, -3.25, 63.125, -15.375, -25.375),
        c32x4!(-108.25, -90.375, 96.875, 63.0, -148.0, -24.625, 23.625, -203.0),
        c32x4!(53.875, -155.5, 16.5, 23.0, 12.375, 119.0, 94.625, -69.625),
        c32x4!(-39.0, -17.5, 56.625, 76.125, 72.0, 133.875, 34.0, -124.375),
        c32x4!(-70.625, -39.125, -82.5, -54.625, -119.75, 3.875, 97.75, 50.5),
        c32x4!(40.75, -66.5, -64.25, -3.875, -35.875, -6.5, -114.125, -51.25),
        c32x4!(-16.75, 27.875, -9.875, 39.0, 50.875, 12.875, 95.25, 7.375),
        c32x4!(-28.625, -32.625, -41.625, 22.625, 28.5, 77.625, -39.5, 22.0),
        c32x4!(-112.25, -14.25, -46.25, -35.875, -19.0, -13.75, 107.375, 60.375),
        c32x4!(-100.75, 99.125, -76.375, -133.5, 66.5, -66.875, -72.5, 26.625),
        c32x4!(-51.375, -45.125, 101.75, -71.625, -56.625, -4.125, -42.0, 21.0),
        c32x4!(-72.5, 0.875, 91.25, 83.625, 22.875, -5.25, 61.25, -185.25),
        c32x4!(-64.25, -11.0, -75.75, -72.5, -191.125, -0.5, 94.75, -46.875),
        c32x4!(-45.25, -54.75, -46.5, -73.0, -51.375, -38.75, 47.0, 121.5),
        c32x4!(-287.75, -7.625, -66.25, -76.875, -30.875, -73.0, -8.875, -2.125),
        c32x4!(-95.25, 17.375, 72.375, 222.375, -181.75, 70.75, 95.25, -177.0),
        c32x4!(-0.5, -34.875, -58.625, -62.875, 36.125, 53.625, 81.375, -188.75),
        c32x4!(-137.125, -58.25, -78.75, -4.625, 95.0, 42.375, -20.25, -95.375),
        c32x4!(-59.25, 11.875, 23.0, 28.5, -17.25, 37.5, -42.875, -123.125),
        c32x4!(-31.625, 9.875, -71.875, 211.75, -144.5, -30.625, 44.25, -220.25),
        c32x4!(-164.625, 35.125, 231.75, -202.25, -68.25, -140.125, 85.0, -172.375),
        c32x4!(-198.0, 172.125, 28.0, 343.25, -108.125, -64.625, -40.0, -107.5),
        c32x4!(69.625, 136.875, 127.25, -126.375, -16.375, 125.125, -17.5, -36.875),
        c32x4!(189.0, 87.25, 143.125, -46.75, 375.0, -16.125, 383.375, -60.0),
        c32x4!(16.25, 10.125, 16.0, -47.5, 11.5, -62.25, -35.75, 54.375),
        c32x4!(-147.875, -163.25, -67.0, -7.375, -116.5, -3.0, -54.25, -257.75),
        c32x4!(-318.0, 29.375, -6.25, -26.625, -42.75, 104.0, 53.0, -193.125),
        c32x4!(33.25, -105.125, 0.375, 63.125, 17.125, 127.75, -85.25, -46.75),
        c32x4!(71.5, -114.75, -167.0, 163.625, 29.125, -23.125, 53.25, -131.0),
        c32x4!(-69.875, -133.25, -67.375, 28.625, 207.5, 107.875, -54.625, -17.625),
        c32x4!(-132.25, -200.375, -12.25, -76.125, 147.75, 100.25, -162.75, -60.0),
        c32x4!(-178.125, -40.25, 12.625, -119.875, 1.375, 65.875, -45.375, -84.75),
        c32x4!(-80.75, -105.625, -53.0, 168.5, 104.125, 41.75, 218.875, -151.75),
        c32x4!(-55.25, -59.75, 78.875, -153.75, -79.625, 10.125, -86.5, -152.0),
        c32x4!(13.625, -81.0, -61.125, -137.375, 283.25, -61.5, -30.5, -54.75),
        c32x4!(-107.375, 100.75, 31.375, 41.0, -11.75, -141.25, 57.375, -155.625),
        c32x4!(-58.375, -82.75, -48.625, -0.25, -50.375, -77.125, 26.375, -77.75),
        c32x4!(-19.75, -36.625, 26.375, 51.0, -50.375, 167.625, -60.375, -64.625),
        c32x4!(26.875, 22.75, -31.625, -18.875, -6.0, 7.25, 6.875, -48.0),
        c32x4!(-45.625, 9.125, 69.75, -37.0, -17.0, 57.375, 33.625, -58.0),
        c32x4!(-43.5, -42.5, 24.0, 26.25, -25.0, 18.5, -35.125, -157.25),
        c32x4!(13.125, -26.875, -20.25, 35.0, 27.75, 26.875, 22.25, -23.625),
        c32x4!(-241.25, 66.25, 35.875, 36.125, -33.0, -10.375, -8.125, -58.875),
        c32x4!(-39.75, -72.125, 58.875, -187.25, 151.125, 51.5, 18.125, -114.125),
        c32x4!(-44.0, 14.25, 89.875, 115.125, 70.75, 146.125, -163.625, -101.75),
        c32x4!(-97.375, -62.75, 113.875, 106.875, -151.25, 47.5, -60.25, 136.625),
        c32x4!(64.5, -15.375, 42.375, -88.5, 58.0, 77.375, -30.25, -379.125),
        c32x4!(-15.0, -4.625, 34.625, 29.625, -25.875, -92.5, -145.5, -118.125),
        c32x4!(-100.5, 114.375, -109.625, -71.75, 30.75, 180.875, 15.625, -54.375),
        c32x4!(-63.5, -22.875, -66.625, 148.625, 201.5, 66.875, 137.375, -207.0),
        c32x4!(55.75, -110.375, -106.375, 134.875, 70.375, 3.125, 50.625, -372.875),
        c32x4!(-125.125, -25.5, -228.375, -138.5, 11.0, -19.5, -16.5, -66.0),
        c32x4!(-23.5, -99.625, 112.0, 20.25, -49.75, 141.125, -51.5, -62.5),
        c32x4!(-70.0, 25.625, -89.25, 97.375, 58.75, -118.0, 46.875, -12.0),
        c32x4!(-61.5, -47.0, -11.75, -58.75, 23.625, -18.75, 97.5, -168.0),
        c32x4!(360.125, -117.875, 298.25, 143.75, 240.625, -75.75, 219.25, -45.125),
        c32x4!(-101.0, 32.0, 168.0, -3.375, -68.125, -15.0, -20.0, -136.125),
        c32x4!(-32.875, -140.375, -42.375, -76.625, -28.125, -108.0, -1.375, -167.25),
        c32x4!(-139.875, 6.625, -90.125, 159.625, 6.75, -2.875, -35.375, -11.375),
        c32x4!(-132.375, 56.25, 26.875, 111.5, 46.375, -66.625, 139.125, -151.25),
        c32x4!(-71.875, -90.0, 22.5, 48.0, 159.625, 44.0, 24.625, -293.625),
        c32x4!(-178.125, -40.875, 78.25, 40.0, 62.375, 73.5, 44.375, -197.0),
        c32x4!(-122.125, -13.625, 164.5, -32.0, -40.875, -39.0, 7.125, -380.0),
        c32x4!(-141.625, 47.625, -19.125, 115.375, 13.5, -11.625, -75.125, -30.625),
        c32x4!(-43.375, 127.375, 1.0, -10.75, 154.25, 104.0, 270.125, -244.0),
        c32x4!(-221.875, -98.125, 11.375, 20.375, -72.5, 71.0, 2.0, -66.25),
        c32x4!(-103.875, -60.375, 91.25, 65.875, 31.875, 104.5, 29.125, -175.875),
        c32x4!(-6.125, -68.125, -50.75, 86.125, 97.0, 60.875, 62.625, -161.75),
        c32x4!(-69.125, 35.375, -17.125, 103.375, 15.625, -53.0, 10.5, -158.25),
        c32x4!(-28.625, -97.375, -53.625, 118.625, 2.625, -23.875, -104.25, -161.5),
        c32x4!(-54.25, -51.125, -23.25, -90.25, -7.875, 22.375, 18.5, 24.875),
        c32x4!(-77.25, -54.0, 80.0, 13.625, 28.25, -13.125, -18.625, -84.0),
        c32x4!(3.5, -70.0, -47.125, 15.125, -24.625, -8.75, 33.0, -49.75),
        c32x4!(24.0, -66.0, -17.875, -44.5, 21.875, -9.0, 77.0, -171.875),
        c32x4!(120.375, -88.875, -6.75, 4.75, 118.375, -6.5, -91.5, -7.75),
        c32x4!(-45.125, -86.625, 52.0, -12.125, 47.875, 88.25, -114.125, -161.625),
        c32x4!(-49.375, -101.0, 153.25, -64.125, -5.625, -211.375, -101.5, -177.875),
        c32x4!(-89.75, -5.625, -85.5, 12.875, 78.875, 99.875, 180.625, 27.5),
        c32x4!(-48.875, 11.125, -54.25, 37.625, -47.375, -34.125, 61.875, -52.625),
        c32x4!(-143.375, -82.875, -112.75, -2.875, 65.75, 295.75, 4.875, -141.625),
        c32x4!(5.0, -64.875, -130.375, -94.375, 42.125, 29.5, -10.875, 6.125),
        c32x4!(14.75, -42.125, -65.0, -109.125, 99.0, -29.125, 21.875, -297.0),
        c32x4!(-29.375, -58.625, 38.875, 24.125, -76.875, -3.625, -81.875, -47.0),
        c32x4!(-54.625, -46.0, -97.375, 17.875, 19.875, -130.5, 25.25, -86.875),
        c32x4!(-97.0, 11.5, 10.625, 10.0, -2.625, -29.0, 126.25, -193.625),
        c32x4!(25.25, -99.625, 67.375, 82.25, -23.75, -86.875, -145.75, -32.0),
        c32x4!(-2.625, -206.125, 40.0, 83.125, 56.25, 59.75, 98.75, -192.875),
        c32x4!(213.5, 18.75, 403.625, 24.125, 319.75, 6.25, 342.75, -78.375),
        c32x4!(-38.5, -75.625, -110.625, -8.0, 122.125, -8.125, 78.875, -172.625),
        c32x4!(-120.0, -122.25, -70.75, 15.875, 150.25, 70.25, 138.0, -61.75),
        c32x4!(-71.125, -76.75, -8.5, -85.0, -42.75, 37.0, 141.25, -188.875),
        c32x4!(-112.5, -23.625, 4.25, -34.5, 122.375, 34.625, 7.5, 26.5),
        c32x4!(-65.0, 41.0, -8.375, -256.25, 15.375, 16.625, 152.625, -116.5),
        c32x4!(78.875, 143.5, -34.5, 78.5, -100.0, -38.625, 31.75, -204.125),
        c32x4!(18.125, -4.375, 74.875, 13.0, -48.5, 62.0, -68.625, -31.625),
        c32x4!(-146.0, 8.125, -83.875, 17.5, 81.25, 18.25, -45.125, -159.5),
        c32x4!(3.375, -91.5, -30.75, -110.75, 56.375, 26.375, -3.0, -101.125),
        c32x4!(-84.0, -92.0, -54.375, 134.375, 157.25, -46.125, 79.25, 16.5),
        c32x4!(71.875, -63.5, -105.0, -31.875, 8.5, 163.875, 33.625, 49.75),
        c32x4!(-79.25, 24.375, -97.875, 33.875, -20.375, -34.375, 71.375, -104.375),
        c32x4!(-158.0, -12.875, -74.0, -11.875, 75.125, 90.625, 21.875, -102.25),
        c32x4!(-78.5, 4.75, 5.25, 56.25, 27.0, 0.5, -16.75, -74.375),
        c32x4!(-81.25, -55.5, -94.0, -37.75, -26.5, -46.875, 27.875, -90.375),
    ];

    #[test]
    fn test_write_main_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let phase_centre = RADec::new(0., -0.471238898038468967);
        let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let flags = Array::from_shape_fn((768, 4), |(c, _)| c > 546);
        let weights = Array::from_elem((768, 4), 8_f32);

        let mut main_table = Table::open(&table_path, TableOpenMode::ReadWrite).unwrap();

        main_table.add_rows(1).unwrap();

        ms_writer
            .write_main_row(
                &mut main_table,
                0,
                5077351976.000001,
                5077351976.000001,
                0,
                1,
                0,
                &vec![-54.23287132681204, -1.1046021021675756, 6.190270793715456],
                2.,
                -1,
                1,
                -1,
                &vec![1., 1., 1., 1.],
                arr2(VIS_DATA_1254670392),
                flags,
                weights,
                false,
            )
            .unwrap();
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

    #[cfg(feature = "mwalib")]
    #[test]
    fn test_write_vis_from_mwalib() {
        use ndarray::Array4;

        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let array_pos = Some(LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        });

        let context = CorrelatorContext::new(
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

        let phase_centre = RADec::from_mwalib_phase_or_pointing(&context.metafits_context);
        let mut ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, array_pos);

        let mwalib_timestep_range = 0..1_usize;
        let mwalib_coarse_chan_range = *context.provided_coarse_chan_indices.first().unwrap()
            ..(*context.provided_coarse_chan_indices.last().unwrap() + 1);
        let mwalib_baseline_idxs = vec![1_usize];

        ms_writer
            .initialize_from_mwalib(&context, &mwalib_timestep_range, &mwalib_coarse_chan_range)
            .unwrap();

        let jones_array = Array3::from_shape_fn((1, 768, 1), |(_, c, _)| {
            Jones::from([
                VIS_DATA_1254670392[c][0],
                VIS_DATA_1254670392[c][1],
                VIS_DATA_1254670392[c][2],
                VIS_DATA_1254670392[c][3],
            ])
        });

        let weight_array = Array4::from_elem((1, 768, 1, 4), 8.);
        let flag_array = Array::from_shape_fn((1, 768, 1, 4), |(_, c, _, _)| c > 546);

        ms_writer
            .write_vis_mwalib(
                jones_array.view(),
                weight_array.view(),
                flag_array.view(),
                &context,
                &mwalib_timestep_range,
                &mwalib_coarse_chan_range,
                &mwalib_baseline_idxs,
            )
            .unwrap();

        let mut main_table = Table::open(&table_path, TableOpenMode::Read).unwrap();
        let mut expected_table =
            Table::open(PATH_1254670392.join(""), TableOpenMode::Read).unwrap();

        assert_table_nrows_match!(main_table, expected_table);
        for col_name in [
            // Time is wrong in Cotter
            // "TIME",
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
            "WEIGHT_SPECTRUM",
            "FLAG",
            // Cotter doesn't write this anyway
            // "FLAG_CATEGORY",
            // Cotter is wrong, it doesn't flag a row when all of its' flags are set
            // "FLAG_ROW",
        ] {
            if col_name == "TIME_CENTROID" {
                assert_table_columns_match!(main_table, expected_table, col_name, 1e-5);
            } else {
                assert_table_columns_match!(main_table, expected_table, col_name);
            }
        }
    }
}
