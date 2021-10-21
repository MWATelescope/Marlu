use crate::{
    c32,
    io::error::MeasurementSetWriteError,
    ndarray::{Array2, Array3, Axis},
};
use flate2::read::GzDecoder;
use ndarray::Array1;
use rubbl_casatables::{
    GlueDataType, Table, TableCreateMode, TableDesc, TableDescCreateMode, TableOpenMode,
};
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

        // TODO: TableMeasDesc for REST_FREQUENCY
        // TableMeasRefDesc measRef(MFrequency::DEFAULT);
        // TableMeasValueDesc measVal(sourceTableDesc, MSSource::columnName(MSSourceEnums::REST_FREQUENCY));
        // TableMeasDesc<MFrequency> restFreqColMeas(measVal, measRef);
        // // write makes the Measure column persistent.
        // restFreqColMeas.write(sourceTableDesc);

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

        // TODO: TableMeasDesc for MWA_DATE_REQUESTED
        // casacore::Vector<Unit> unitVec(1);
        // unitVec[0] = Unit("s");
        // TableMeasRefDesc measRef(MEpoch::DEFAULT);
        // TableMeasValueDesc measVal(columnName(MWAMSEnums::MWA_DATE_REQUESTED));
        // TableMeasDesc<MEpoch> intervalColMeas(measVal, measRef, unitVec);
        // intervalColMeas.write(obsTable);
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

        let pointing_table_path = self.path.join("MWA_TILE_POINTING");
        let pointing_table = Table::new(
            pointing_table_path,
            pointing_table_desc,
            0,
            TableCreateMode::New,
        )
        .unwrap();

        // TODO: TableMeasDesc for INTERVAL
        // casacore::Vector<Unit> unitVec(1);
        // unitVec[0] = Unit("s");
        // TableMeasRefDesc measRef(MEpoch::DEFAULT);
        // TableMeasValueDesc measVal(tilePointingTableDesc, columnName(MWAMSEnums::INTERVAL));
        // TableMeasDesc<MEpoch> intervalColMeas(measVal, measRef, unitVec);
        // intervalColMeas.write(tilePointingTableDesc);

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

    /// Write a row into the SPECTRAL_WINDOW table.
    /// Return the row index.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
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
        name: &str,
        ref_freq: f64,
        chan_info: Array2<f64>,
        total_bw: f64,
        flag: bool,
    ) -> Result<u64, MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        let spw_idx = table.n_rows();

        match chan_info.shape() {
            [num_chans, 4] => {
                table.add_rows(1).unwrap();
                table
                    .put_cell("NUM_CHAN", spw_idx, &(*num_chans as i32))
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

    /// Write a row into the SPECTRAL_WINDOW table with extra mwa columns enabled.
    /// Return the row index.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
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
        name: &str,
        ref_freq: f64,
        chan_info: Array2<f64>,
        total_bw: f64,
        centre_subband_nr: i32,
        flag: bool,
    ) -> Result<u64, MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        let spw_idx = self
            .write_spectral_window_row(table, name, ref_freq, chan_info, total_bw, flag)
            .unwrap();

        table
            .put_cell("MWA_CENTRE_SUBBAND_NR", spw_idx, &centre_subband_nr)
            .unwrap();

        Ok(spw_idx)
    }

    /// Write a row into the DATA_DESCRIPTION table.
    /// Return the row index.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `spectral_window_id` - Pointer to spectralwindow table
    /// - `polarization_id` - Pointer to polarization table
    /// - `flag_row` - Flag this row
    pub fn write_data_description_row(
        &self,
        table: &mut Table,
        spectral_window_id: i32,
        polarization_id: i32,
        flag_row: bool,
    ) -> Result<u64, MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

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

    /// Write a row into the `ANTENNA` table.
    /// Return the row index.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
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
        name: &str,
        station: &str,
        ant_type: &str,
        mount: &str,
        position: &Vec<f64>,
        dish_diameter: f64,
        flag_row: bool,
    ) -> Result<u64, MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        let ant_idx = table.n_rows();
        table.add_rows(1).unwrap();

        table.put_cell("NAME", ant_idx, &name.to_string()).unwrap();
        table
            .put_cell("STATION", ant_idx, &station.to_string())
            .unwrap();
        table
            .put_cell("TYPE", ant_idx, &ant_type.to_string())
            .unwrap();
        table
            .put_cell("MOUNT", ant_idx, &mount.to_string())
            .unwrap();
        table.put_cell("POSITION", ant_idx, position).unwrap();
        table
            .put_cell("DISH_DIAMETER", ant_idx, &dish_diameter)
            .unwrap();
        table.put_cell("FLAG_ROW", ant_idx, &flag_row).unwrap();
        Ok(ant_idx)
    }

    /// Write a row into the `ANTENNA` table with extra mwa columns enabled.
    /// Return the row index.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
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
    ) -> Result<u64, MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        let ant_idx = self
            .write_antenna_row(
                table,
                name,
                station,
                ant_type,
                mount,
                position,
                dish_diameter,
                flag_row,
            )
            .unwrap();

        table.put_cell("MWA_INPUT", ant_idx, input).unwrap();
        table.put_cell("MWA_TILE_NR", ant_idx, &tile_nr).unwrap();
        table.put_cell("MWA_RECEIVER", ant_idx, &receiver).unwrap();
        table.put_cell("MWA_SLOT", ant_idx, slot).unwrap();
        table
            .put_cell("MWA_CABLE_LENGTH", ant_idx, cable_length)
            .unwrap();

        Ok(ant_idx)
    }

    /// Write a row into the `POLARIZATION` table.
    /// Return the row index.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
    /// - `corr_type` - The polarization type for each correlation product, as a Stokes enum.
    /// - `corr_product` - Indices describing receptors of feed going into correlation.
    ///     Shape should be [n, 2] where n is the length of `corr_type`
    /// - `flag_row` - Row flag
    pub fn write_polarization_row(
        &self,
        table: &mut Table,
        corr_type: &Vec<i32>,
        corr_product: &Array2<i32>,
        flag_row: bool,
    ) -> Result<u64, MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        let pol_idx = table.n_rows();

        let num_corr_type = corr_type.len();

        match corr_product.shape() {
            [num_corr, 2] if *num_corr == num_corr_type => {
                table.add_rows(1).unwrap();
                table
                    .put_cell("NUM_CORR", pol_idx, &(*num_corr as i32))
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

        table.put_cell("CORR_TYPE", pol_idx, corr_type).unwrap();
        table
            .put_cell("CORR_PRODUCT", pol_idx, corr_product)
            .unwrap();
        table.put_cell("FLAG_ROW", pol_idx, &flag_row).unwrap();
        Ok(pol_idx)
    }

    /// Write a row into the `FIELD` table.
    /// Return the row index.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
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
    pub fn write_field_row(
        &self,
        table: &mut Table,
        name: &str,
        code: &str,
        time: f64,
        dir_info: &Array3<f64>,
        source_id: i32,
        flag_row: bool,
    ) -> Result<u64, MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        let field_idx = table.n_rows();

        match dir_info.shape() {
            [3, p, 2] if *p > 0 => {
                table.add_rows(1).unwrap();
                table
                    .put_cell("NUM_POLY", field_idx, &((*p - 1) as i32))
                    .unwrap();
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

        table
            .put_cell("NAME", field_idx, &name.to_string())
            .unwrap();
        table
            .put_cell("CODE", field_idx, &code.to_string())
            .unwrap();
        table.put_cell("TIME", field_idx, &time).unwrap();

        let col_names = ["DELAY_DIR", "PHASE_DIR", "REFERENCE_DIR"];
        for (value, &col_name) in dir_info.outer_iter().zip(col_names.iter()) {
            // println!("{:?}", value.shape());
            table
                .put_cell(col_name, field_idx, &value.to_owned())
                .unwrap();
        }

        table.put_cell("SOURCE_ID", field_idx, &source_id).unwrap();
        table.put_cell("FLAG_ROW", field_idx, &flag_row).unwrap();
        Ok(field_idx)
    }

    /// Write a row into the `FIELD` table with extra mwa columns enabled.
    /// Return the row index.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
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
        name: &str,
        code: &str,
        time: f64,
        dir_info: &Array3<f64>,
        source_id: i32,
        has_calibrator: bool,
        flag_row: bool,
    ) -> Result<u64, MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        let field_idx = self
            .write_field_row(table, name, code, time, dir_info, source_id, flag_row)
            .unwrap();

        table
            .put_cell("MWA_HAS_CALIBRATOR", field_idx, &has_calibrator)
            .unwrap();

        Ok(field_idx)
    }

    /// Write a row into the `OBSERVATION` table.
    /// Return the row index.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
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
        telescope_name: &str,
        time_range: (f64, f64),
        observer: &str,
        schedule_type: &str,
        project: &str,
        release_date: f64,
        flag_row: bool,
    ) -> Result<u64, MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        let obs_idx = table.n_rows();
        table.add_rows(1).unwrap();

        table
            .put_cell("TELESCOPE_NAME", obs_idx, &telescope_name.to_string())
            .unwrap();
        let time_range = vec![time_range.0, time_range.1];
        table.put_cell("TIME_RANGE", obs_idx, &time_range).unwrap();
        table
            .put_cell("OBSERVER", obs_idx, &observer.to_string())
            .unwrap();
        table
            .put_cell("SCHEDULE_TYPE", obs_idx, &schedule_type.to_string())
            .unwrap();
        table
            .put_cell("PROJECT", obs_idx, &project.to_string())
            .unwrap();
        table
            .put_cell("RELEASE_DATE", obs_idx, &release_date)
            .unwrap();
        table.put_cell("FLAG_ROW", obs_idx, &flag_row).unwrap();
        Ok(obs_idx)
    }

    /// Write a row into the `OBSERVATION` table with extra mwa columns enabled.
    /// Return the row index.
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
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
    ) -> Result<u64, MeasurementSetWriteError> {
        // TODO: fix all these unwraps after https://github.com/pkgw/rubbl/pull/148

        let obs_idx = self
            .write_observation_row(
                table,
                telescope_name,
                time_range,
                observer,
                schedule_type,
                project,
                release_date,
                flag_row,
            )
            .unwrap();

        table.put_cell("MWA_GPS_TIME", obs_idx, &gps_time).unwrap();
        table
            .put_cell("MWA_FILENAME", obs_idx, &filename.to_string())
            .unwrap();
        table
            .put_cell(
                "MWA_OBSERVATION_MODE",
                obs_idx,
                &observation_mode.to_string(),
            )
            .unwrap();
        table
            .put_cell("MWA_FLAG_WINDOW_SIZE", obs_idx, &flag_window_size)
            .unwrap();
        table
            .put_cell("MWA_DATE_REQUESTED", obs_idx, &date_requested)
            .unwrap();
        Ok(obs_idx)
    }

    /// TODO
    ///
    /// Write a row into the `HISTORY_ITERM` table.
    /// Return the row index.
    pub fn write_history_item_row(
        &self,
        table: &mut Table,
    ) -> Result<u64, MeasurementSetWriteError> {
        Ok(0)
    }

    /// Write a row into the `MWA_TILE_POINTING` table.
    /// Return the row index.
    ///
    /// - `start` - start MJD of observation
    /// - `end` - end MJD of observation
    /// - `delays` - beamformer delays, from metafits:DELAYS
    /// - `direction_{ra|dec}` - pointing direction [Ra/Dec]
    pub fn write_mwa_tile_pointing_row(
        &self,
        table: &mut Table,
        start: f64,
        end: f64,
        delays: &Vec<i32>,
        direction_ra: f64,
        direction_dec: f64,
    ) -> Result<u64, MeasurementSetWriteError> {
        let point_idx = table.n_rows();
        table.add_rows(1).unwrap();

        table
            .put_cell("INTERVAL", point_idx, &vec![start, end])
            .unwrap();
        table.put_cell("DELAYS", point_idx, delays).unwrap();
        table
            .put_cell("DIRECTION", point_idx, &vec![direction_ra, direction_dec])
            .unwrap();

        Ok(point_idx)
    }

    /// Write a row into the `MWA_SUBBAND` table.
    /// Return the row index.
    ///
    /// - `number` - Subband (coarse channel) index
    /// - `gain` - (deprecated) - from metafits:CHANGAIN, use 0.
    /// - `flag_row` - flag this subband
    pub fn write_mwa_subband_row(
        &self,
        table: &mut Table,
        number: i32,
        gain: f64,
        flag_row: bool,
    ) -> Result<u64, MeasurementSetWriteError> {
        let point_idx = table.n_rows();
        table.add_rows(1).unwrap();

        table.put_cell("NUMBER", point_idx, &number).unwrap();
        table.put_cell("GAIN", point_idx, &gain).unwrap();
        table.put_cell("FLAG_ROW", point_idx, &flag_row).unwrap();

        Ok(point_idx)
    }

    /// Write a row into the main table.
    /// Return the row index.
    ///
    /// The main table holds measurements from a Telescope
    ///
    /// - `table` - [`rubbl_casatables::Table`] object to write to.
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
        weights: Array1<f32>,
    ) -> Result<u64, MeasurementSetWriteError> {
        let idx = table.n_rows();
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
            ([d0, d1], [f0, f1], [w0])
                if f0 == d0 && d1 == &num_pols && f1 == &num_pols && w0 == &num_pols => {}
            (dsh, fsh, wsh) => {
                return Err(MeasurementSetWriteError::BadArrayShape {
                    argument: "data|flags|weights".into(),
                    function: "write_main_row".into(),
                    expected: format!(
                        "[n, p]|[n, p]|[p] where n=num_chans, p=num_pols({})",
                        num_pols
                    )
                    .into(),
                    received: format!("{:?}|{:?}|{:?}", dsh, fsh, wsh).into(),
                })
            }
        }

        // TODO: calculate weight aggregate

        table.add_rows(1).unwrap();

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
        table.put_cell("WEIGHT", idx, &weights).unwrap();
        table.put_cell("FLAG", idx, &flags).unwrap();

        Ok(idx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    use crate::{
        approx::abs_diff_eq,
        c32,
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
                    $left, $right, $col_name, col_desc, c32
                ),
                GlueDataType::TpDComplex => assert_table_column_values_match_approx!(
                    $left, $right, $col_name, col_desc, c32
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

        // TODO: check TableMeasDesc

        let mut main_table = Table::open(&table_path.clone(), TableOpenMode::Read).unwrap();
        let main_table_keywords = main_table.table_keyword_names().unwrap();
        assert!(main_table_keywords.contains(&"SOURCE".into()));
    }

    #[test]
    fn test_add_mwa_mods() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
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
        let spw_table_path = table_path.join("SPECTRAL_WINDOW");
        let mut spw_table = Table::open(&spw_table_path, TableOpenMode::ReadWrite).unwrap();

        let result = ms_writer
            .write_spectral_window_row(
                &mut spw_table,
                "MWA_BAND_182.4",
                182395000.,
                chan_info,
                30720000.,
                false,
            )
            .unwrap();
        drop(ms_writer);

        assert_eq!(result, 0);

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
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
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

        let result = ms_writer
            .write_spectral_window_row_mwa(
                &mut spw_table,
                "MWA_BAND_182.4",
                182395000.,
                chan_info,
                30720000.,
                143,
                false,
            )
            .unwrap();
        drop(ms_writer);

        assert_eq!(result, 0);

        let mut spw_table = Table::open(&spw_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("SPECTRAL_WINDOW"), TableOpenMode::Read).unwrap();

        assert_tables_match!(spw_table, expected_table);
    }

    #[test]
    fn handle_bad_spw_chan_info() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let chan_info = Array2::from_shape_fn((768, 3), |(_, _)| 40000.);

        let spw_table_path = table_path.join("SPECTRAL_WINDOW");
        let mut spw_table = Table::open(&spw_table_path, TableOpenMode::ReadWrite).unwrap();

        let result = ms_writer.write_spectral_window_row(
            &mut spw_table,
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
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let ddesc_table_path = table_path.join("DATA_DESCRIPTION");
        let mut ddesc_table = Table::open(&ddesc_table_path, TableOpenMode::ReadWrite).unwrap();

        let result = ms_writer
            .write_data_description_row(&mut ddesc_table, 0, 0, false)
            .unwrap();
        drop(ms_writer);

        assert_eq!(result, 0);

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
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let ant_table_path = ms_writer.path.join("ANTENNA");
        let mut ant_table = Table::open(ant_table_path, TableOpenMode::ReadWrite).unwrap();

        for (idx, (name, position)) in izip!(ANT_NAMES, ANT_POSITIONS).enumerate() {
            let position = position.iter().cloned().collect();

            let result = ms_writer
                .write_antenna_row(
                    &mut ant_table,
                    name,
                    "MWA",
                    "GROUND-BASED",
                    "ALT-AZ",
                    &position,
                    4.0,
                    false,
                )
                .unwrap();
            assert_eq!(result, idx as u64);
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
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);
        ms_writer.add_mwa_mods();

        let ant_table_path = ms_writer.path.join("ANTENNA");
        let mut ant_table = Table::open(ant_table_path.clone(), TableOpenMode::ReadWrite).unwrap();

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

            let result = ms_writer
                .write_antenna_row_mwa(
                    &mut ant_table,
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
            assert_eq!(result, idx as u64);
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
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let corr_product = array![[0, 0], [0, 1], [1, 0], [1, 1]];
        let corr_type = vec![9, 10, 11, 12];

        let pol_table_path = table_path.join("POLARIZATION");
        let mut pol_table = Table::open(&pol_table_path, TableOpenMode::ReadWrite).unwrap();

        let result = ms_writer
            .write_polarization_row(&mut pol_table, &corr_type, &corr_product, false)
            .unwrap();
        drop(ms_writer);

        assert_eq!(result, 0);

        let mut pol_table = Table::open(&pol_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("POLARIZATION"), TableOpenMode::Read).unwrap();

        assert_tables_match!(pol_table, expected_table);
    }

    #[test]
    fn handle_bad_pol_small_corr_type() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let corr_product = array![[0, 0], [0, 1], [1, 0], [1, 1]];
        let corr_type = vec![9, 10, 11];

        let pol_table_path = table_path.join("POLARIZATION");
        let mut pol_table = Table::open(&pol_table_path, TableOpenMode::ReadWrite).unwrap();

        let result =
            ms_writer.write_polarization_row(&mut pol_table, &corr_type, &corr_product, false);

        assert!(matches!(
            result,
            Err(MeasurementSetWriteError::BadArrayShape { .. })
        ))
    }

    #[test]
    fn handle_bad_pol_big_corr_product() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let corr_product = array![[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]];
        let corr_type = vec![9, 10, 11, 12];

        let pol_table_path = table_path.join("POLARIZATION");
        let mut pol_table = Table::open(&pol_table_path, TableOpenMode::ReadWrite).unwrap();

        let result =
            ms_writer.write_polarization_row(&mut pol_table, &corr_type, &corr_product, false);

        assert!(matches!(
            result,
            Err(MeasurementSetWriteError::BadArrayShape { .. })
        ))
    }

    #[test]
    fn test_write_field_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let field_table_path = table_path.join("FIELD");
        let mut field_table = Table::open(&field_table_path, TableOpenMode::ReadWrite).unwrap();

        let dir_info = array![
            [[0., -0.471238898038468967]],
            [[0., -0.471238898038468967]],
            [[0., -0.471238898038468967]]
        ];
        let result = ms_writer
            .write_field_row(
                &mut field_table,
                "high_2019B_2458765_EOR0_RADec0.0,-27.0",
                "",
                5077351974.,
                &dir_info,
                -1,
                false,
            )
            .unwrap();
        drop(ms_writer);

        assert_eq!(result, 0);

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
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);
        ms_writer.add_mwa_mods();

        let field_table_path = table_path.join("FIELD");
        let mut field_table = Table::open(&field_table_path, TableOpenMode::ReadWrite).unwrap();

        let dir_info = array![
            [[0., -0.471238898038468967]],
            [[0., -0.471238898038468967]],
            [[0., -0.471238898038468967]]
        ];
        let result = ms_writer
            .write_field_row_mwa(
                &mut field_table,
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

        assert_eq!(result, 0);

        let mut field_table = Table::open(&field_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("FIELD"), TableOpenMode::Read).unwrap();

        assert_tables_match!(field_table, expected_table);
    }

    #[test]
    fn handle_bad_field_shape() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let field_table_path = table_path.join("FIELD");
        let mut field_table = Table::open(&field_table_path, TableOpenMode::ReadWrite).unwrap();

        let dir_info = array![
            [[0., -0.471238898038468967]],
            [[0., -0.471238898038468967]],
            [[0., -0.471238898038468967]],
            [[0., 1.]]
        ];
        let result = ms_writer.write_field_row(
            &mut field_table,
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
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let obs_table_path = table_path.join("OBSERVATION");
        let mut obs_table = Table::open(&obs_table_path, TableOpenMode::ReadWrite).unwrap();

        let result = ms_writer
            .write_observation_row(
                &mut obs_table,
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

        assert_eq!(result, 0);

        let mut obs_table = Table::open(&obs_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("OBSERVATION"), TableOpenMode::Read).unwrap();

        assert_table_nrows_match!(obs_table, expected_table);
        for col_name in [
            "TIME_RANGE",
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
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);
        ms_writer.add_mwa_mods();

        let obs_table_path = table_path.join("OBSERVATION");
        let mut obs_table = Table::open(&obs_table_path, TableOpenMode::ReadWrite).unwrap();

        let result = ms_writer
            .write_observation_row_mwa(
                &mut obs_table,
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

        assert_eq!(result, 0);

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
    fn test_write_mwa_tile_pointing_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);
        ms_writer.add_mwa_mods();

        let point_table_path = table_path.join("MWA_TILE_POINTING");
        let mut point_table = Table::open(&point_table_path, TableOpenMode::ReadWrite).unwrap();

        let result = ms_writer
            .write_mwa_tile_pointing_row(
                &mut point_table,
                5077351975.,
                5077351984.,
                &vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0],
                6.28310690918895887,
                -0.464403366228935188,
            )
            .unwrap();
        drop(ms_writer);

        assert_eq!(result, 0);

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
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);
        ms_writer.add_mwa_mods();

        let subband_table_path = table_path.join("MWA_SUBBAND");
        let mut subband_table = Table::open(&subband_table_path, TableOpenMode::ReadWrite).unwrap();

        for idx in 0..24 {
            let result = ms_writer
                .write_mwa_subband_row(&mut subband_table, idx, 0., false)
                .unwrap();
            assert_eq!(result, idx as u64);
        }

        drop(ms_writer);

        let mut subband_table = Table::open(&subband_table_path, TableOpenMode::Read).unwrap();

        let mut expected_table =
            Table::open(PATH_1254670392.join("MWA_SUBBAND"), TableOpenMode::Read).unwrap();

        assert_tables_match!(subband_table, expected_table);
    }

    #[test]
    fn test_write_main_row() {
        let temp_dir = tempdir().unwrap();
        let table_path = temp_dir.path().join("test.ms");
        let ms_writer = MeasurementSetWriter::new(table_path.clone());
        ms_writer.decompress_default_tables().unwrap();
        ms_writer.decompress_source_table().unwrap();
        ms_writer.add_cotter_mods(768);

        let data = arr2(&[
            c32x4!(26134.125, 0., 0., 0., -867.5, 154., 26266.875, 0.),
            c32x4!(28940.875, 0., 0., 0., -944.875, 157.625, 29159.75, 0.),
            c32x4!(34085., 0., 0., 0., -1165., -134.625, 34177.75, 0.),
            c32x4!(39534., 0., 0., 0., -1411.75, 37.625, 39789.5, 0.),
            c32x4!(44169.625, 0., 0., 0., -1381.5, -61.375, 44585., 0.),
            c32x4!(47219., 0., 0., 0., -1760.875, 67.875, 47436.375, 0.),
            c32x4!(48642.375, 0., 0., 0., -1814.5, 58.75, 49032.625, 0.),
            c32x4!(49047.25, 0., 0., 0., -1800., 41.125, 49140., 0.),
            c32x4!(48757.25, 0., 0., 0., -1660.75, -50.125, 49267., 0.),
            c32x4!(48756., 0., 0., 0., -1854.125, 9.875, 48875.75, 0.),
            c32x4!(48865.25, 0., 0., 0., -1648., 115.375, 48852.875, 0.),
            c32x4!(48828.25, 0., 0., 0., -1946., 169.5, 48984.5, 0.),
            c32x4!(49121.375, 0., 0., 0., -1682., 62.5, 49216.625, 0.),
            c32x4!(48863.75, 0., 0., 0., -1808.375, 62.375, 49011.375, 0.),
            c32x4!(49068.625, 0., 0., 0., -1867.625, 58.375, 49006.75, 0.),
            c32x4!(48936.875, 0., 0., 0., -1866.25, -113.5, 49516., 0.),
            c32x4!(49001.625, 0., 0., 0., -1315.375, 244.125, 49147.125, 0.),
            c32x4!(49034.5, 0., 0., 0., -1850.5, 71.625, 49491.375, 0.),
            c32x4!(49019.5, 0., 0., 0., -2073.75, 227.25, 49396.625, 0.),
            c32x4!(49248., 0., 0., 0., -1817.125, -137.25, 49431.5, 0.),
            c32x4!(49029.375, 0., 0., 0., -1597.875, 102., 49570.125, 0.),
            c32x4!(48731.25, 0., 0., 0., -1552.125, 85.375, 49542.75, 0.),
            c32x4!(49431.375, 0., 0., 0., -1600.375, -7.375, 49484.25, 0.),
            c32x4!(49294., 0., 0., 0., -1979.875, 201.375, 49633.875, 0.),
            c32x4!(49568.375, 0., 0., 0., -1624.125, 160.625, 49627.875, 0.),
            c32x4!(49000.375, 0., 0., 0., -1609.5, -108.125, 49340.625, 0.),
            c32x4!(47914.875, 0., 0., 0., -1622.625, 19.625, 48311.375, 0.),
            c32x4!(44976.5, 0., 0., 0., -1528.125, 12., 45825.75, 0.),
            c32x4!(40752.125, 0., 0., 0., -1263., -39.875, 41233.75, 0.),
            c32x4!(35286.125, 0., 0., 0., -1276.5, -88.375, 35738.375, 0.),
            c32x4!(29878.375, 0., 0., 0., -979.75, 63.125, 30324.125, 0.),
            c32x4!(26534.75, 0., 0., 0., -931.75, -37.875, 26891.5, 0.),
            c32x4!(26492., 0., 0., 0., -927.25, 93.75, 26953., 0.),
            c32x4!(29331.375, 0., 0., 0., -1007.5, 41., 29891.625, 0.),
            c32x4!(34517.125, 0., 0., 0., -1166.25, 69.5, 35022.125, 0.),
            c32x4!(40148.375, 0., 0., 0., -1470.25, 35.625, 40891.25, 0.),
            c32x4!(44521.25, 0., 0., 0., -1715., 56., 45387.5, 0.),
            c32x4!(47516., 0., 0., 0., -1705.625, 98.375, 48694.25, 0.),
            c32x4!(49433.375, 0., 0., 0., -1663.375, -26.25, 49870., 0.),
            c32x4!(49397.875, 0., 0., 0., -1672.75, 243.25, 50122.5, 0.),
            c32x4!(49684.125, 0., 0., 0., -1947.75, 118.5, 49875.375, 0.),
            c32x4!(49532.75, 0., 0., 0., -1771.25, 126.625, 49996.375, 0.),
            c32x4!(49629.625, 0., 0., 0., -1689.875, -68.25, 49910.75, 0.),
            c32x4!(49366.125, 0., 0., 0., -1638.875, 235.875, 50207.625, 0.),
            c32x4!(49502.75, 0., 0., 0., -1613.5, -182., 50043.125, 0.),
            c32x4!(49425.75, 0., 0., 0., -1772.125, 255.625, 50069.625, 0.),
            c32x4!(49572., 0., 0., 0., -1954., 13.75, 50249.625, 0.),
            c32x4!(49597.25, 0., 0., 0., -1687.125, 236.375, 50100.5, 0.),
            c32x4!(49914.25, 0., 0., 0., -1471.875, 246.625, 50550.5, 0.),
            c32x4!(49339.375, 0., 0., 0., -1629., 6.75, 50034.5, 0.),
            c32x4!(49689.375, 0., 0., 0., -1914.25, 316.75, 50231.875, 0.),
            c32x4!(49721.75, 0., 0., 0., -1923.75, 183.625, 50414.375, 0.),
            c32x4!(49691., 0., 0., 0., -1711.875, 139.375, 50329.875, 0.),
            c32x4!(49842.125, 0., 0., 0., -1680.75, 120.125, 50100.875, 0.),
            c32x4!(50036.125, 0., 0., 0., -1806.25, 48.375, 50454.625, 0.),
            c32x4!(50072.5, 0., 0., 0., -1630.25, 215., 50734.25, 0.),
            c32x4!(50078.75, 0., 0., 0., -1847., -20.25, 50360.125, 0.),
            c32x4!(49793.25, 0., 0., 0., -1619.125, 134.75, 50273.25, 0.),
            c32x4!(48479.5, 0., 0., 0., -1740.375, -7.375, 49218.25, 0.),
            c32x4!(45842.125, 0., 0., 0., -1764., 2.625, 46629.625, 0.),
            c32x4!(41494.625, 0., 0., 0., -1476.625, 84.625, 41988.875, 0.),
            c32x4!(35740.125, 0., 0., 0., -1377.625, -82.875, 36524.25, 0.),
            c32x4!(30237.875, 0., 0., 0., -984.5, 161.25, 30872.25, 0.),
            c32x4!(26857., 0., 0., 0., -1028.125, 42.125, 27267.875, 0.),
            c32x4!(25889.625, 0., 0., 0., -881.25, 136.5, 26453.625, 0.),
            c32x4!(28979.75, 0., 0., 0., -1149., 88.875, 29474.875, 0.),
            c32x4!(33682.625, 0., 0., 0., -1279., 62.375, 34679.375, 0.),
            c32x4!(39557.75, 0., 0., 0., -1512.625, 130.625, 40141.75, 0.),
            c32x4!(44098.125, 0., 0., 0., -1462., 62.25, 44656.625, 0.),
            c32x4!(47034., 0., 0., 0., -1750.875, 242.5, 47682.5, 0.),
            c32x4!(48190., 0., 0., 0., -1677.125, 321.625, 49401.5, 0.),
            c32x4!(48665., 0., 0., 0., -1734., 101.125, 49856.5, 0.),
            c32x4!(48650.75, 0., 0., 0., -1740.875, -220.5, 49570.25, 0.),
            c32x4!(48447.75, 0., 0., 0., -1755.25, -1.75, 49415.125, 0.),
            c32x4!(48434.875, 0., 0., 0., -1753.375, 110.5, 49349.625, 0.),
            c32x4!(48512.375, 0., 0., 0., -1638.375, 241.625, 49595.25, 0.),
            c32x4!(48571.25, 0., 0., 0., -1533.375, 105., 49366.375, 0.),
            c32x4!(48160.5, 0., 0., 0., -1919.75, 138.75, 49607.75, 0.),
            c32x4!(48049.75, 0., 0., 0., -1764.25, 79.5, 49088.625, 0.),
            c32x4!(48494.75, 0., 0., 0., -1668.875, 171.875, 49501.875, 0.),
            c32x4!(48982.75, 0., 0., 0., -1349., 394.375, 49739.875, 0.),
            c32x4!(48512.5, 0., 0., 0., -1825.375, 289.75, 49085.125, 0.),
            c32x4!(48535.875, 0., 0., 0., -1883.875, 141., 49219.875, 0.),
            c32x4!(48717.75, 0., 0., 0., -1512., 273.625, 49441.625, 0.),
            c32x4!(48586.25, 0., 0., 0., -1526.75, 244.125, 49590.25, 0.),
            c32x4!(48413.5, 0., 0., 0., -1722.375, 391., 49733.625, 0.),
            c32x4!(48992.25, 0., 0., 0., -1738.25, 174.5, 49313.125, 0.),
            c32x4!(48929.875, 0., 0., 0., -1757.25, -153.25, 49494.625, 0.),
            c32x4!(49031.375, 0., 0., 0., -1637.75, 1.625, 49893.125, 0.),
            c32x4!(48798.25, 0., 0., 0., -1674.625, 59.125, 49291.875, 0.),
            c32x4!(47650., 0., 0., 0., -1842.125, 108.375, 48311.875, 0.),
            c32x4!(45108.125, 0., 0., 0., -1808., 126.875, 45435.875, 0.),
            c32x4!(40429.875, 0., 0., 0., -1365.375, 95.5, 41316.875, 0.),
            c32x4!(35334.75, 0., 0., 0., -1320.875, 130.5, 36063.875, 0.),
            c32x4!(29831.375, 0., 0., 0., -954.5, 215.75, 30397.125, 0.),
            c32x4!(26414.5, 0., 0., 0., -908.75, 102.25, 26821., 0.),
            c32x4!(24709.125, 0., 0., 0., -842., 116.625, 25140.875, 0.),
            c32x4!(27489.25, 0., 0., 0., -1059.5, 171.25, 27891.75, 0.),
            c32x4!(32120.375, 0., 0., 0., -1164.75, 244.625, 32858.125, 0.),
            c32x4!(37593.25, 0., 0., 0., -1221.25, -31.25, 38179.25, 0.),
            c32x4!(41812.5, 0., 0., 0., -1746.625, 95., 42359., 0.),
            c32x4!(44467.75, 0., 0., 0., -1641.25, 177.125, 45359.125, 0.),
            c32x4!(46023.125, 0., 0., 0., -1688.875, 213.625, 46648.375, 0.),
            c32x4!(45897.25, 0., 0., 0., -1782., 79.75, 46845.875, 0.),
            c32x4!(46162.75, 0., 0., 0., -1424.75, 198., 47028.375, 0.),
            c32x4!(46203.625, 0., 0., 0., -1651.5, 217.875, 46533.25, 0.),
            c32x4!(45645.125, 0., 0., 0., -1867.875, 300.75, 46546.375, 0.),
            c32x4!(46094.5, 0., 0., 0., -1557.625, 155., 46453.75, 0.),
            c32x4!(45814., 0., 0., 0., -1654.125, 209., 46753.875, 0.),
            c32x4!(45987.375, 0., 0., 0., -1808.125, -72., 46312.5, 0.),
            c32x4!(45699.875, 0., 0., 0., -1645.125, 276., 46211.875, 0.),
            c32x4!(45921., 0., 0., 0., -1703.75, 129.875, 46500.5, 0.),
            c32x4!(46374.375, 0., 0., 0., -1180.625, 171.75, 46649.375, 0.),
            c32x4!(45584.375, 0., 0., 0., -1627.625, 210.875, 46245.125, 0.),
            c32x4!(46114.125, 0., 0., 0., -1620.5, 196.375, 46619.5, 0.),
            c32x4!(46004.125, 0., 0., 0., -1814., 151.5, 46420.25, 0.),
            c32x4!(46107.875, 0., 0., 0., -1585., 43.875, 46505.375, 0.),
            c32x4!(45852.75, 0., 0., 0., -1755.375, 213.5, 46341.75, 0.),
            c32x4!(46173.375, 0., 0., 0., -1619.25, 282.75, 46565., 0.),
            c32x4!(46114.75, 0., 0., 0., -1675., 251.25, 46829.5, 0.),
            c32x4!(46478., 0., 0., 0., -1726.375, 266., 46916.25, 0.),
            c32x4!(46178.25, 0., 0., 0., -1619.75, 213.875, 46873.75, 0.),
            c32x4!(44984.25, 0., 0., 0., -1562.375, 352.75, 45543.625, 0.),
            c32x4!(42465.625, 0., 0., 0., -1596., 123.5, 43310.875, 0.),
            c32x4!(38271.75, 0., 0., 0., -1429.75, 111.25, 38842.625, 0.),
            c32x4!(33223.25, 0., 0., 0., -1147.5, 77.125, 33704.875, 0.),
            c32x4!(28216., 0., 0., 0., -1012.875, 331., 28525.375, 0.),
            c32x4!(25014.25, 0., 0., 0., -903., 73.125, 25541.375, 0.),
            c32x4!(24015.25, 0., 0., 0., -838.125, 155.375, 24546.375, 0.),
            c32x4!(26790.125, 0., 0., 0., -856., 111.5, 27310.25, 0.),
            c32x4!(31208.125, 0., 0., 0., -1109., 130.75, 31861.875, 0.),
            c32x4!(36162.125, 0., 0., 0., -1225.625, 176.75, 37212.25, 0.),
            c32x4!(40178.375, 0., 0., 0., -1261.625, 107.875, 41282., 0.),
            c32x4!(43125.25, 0., 0., 0., -1598.5, -70.25, 44425.5, 0.),
            c32x4!(44448.5, 0., 0., 0., -1627.125, 349.75, 45540.375, 0.),
            c32x4!(44373.75, 0., 0., 0., -1829.375, 192.625, 45763.375, 0.),
            c32x4!(44222.375, 0., 0., 0., -1787.375, 163.125, 45916.125, 0.),
            c32x4!(44524.75, 0., 0., 0., -1859., 163.125, 45272.125, 0.),
            c32x4!(44134.375, 0., 0., 0., -1659.875, 134.25, 45397., 0.),
            c32x4!(44510.375, 0., 0., 0., -1538.75, -55., 45363.375, 0.),
            c32x4!(44252.25, 0., 0., 0., -1646., 181.375, 44829., 0.),
            c32x4!(44123.625, 0., 0., 0., -1668.75, 178.875, 45187.125, 0.),
            c32x4!(43978.75, 0., 0., 0., -1593.75, 242.5, 45242.875, 0.),
            c32x4!(44157.25, 0., 0., 0., -1729.75, 281., 45092.5, 0.),
            c32x4!(44536.875, 0., 0., 0., -1218.875, 305.625, 45192., 0.),
            c32x4!(44240.625, 0., 0., 0., -1548., 89., 44976.375, 0.),
            c32x4!(44447.5, 0., 0., 0., -1428.25, 267.25, 44753.125, 0.),
            c32x4!(44128.25, 0., 0., 0., -1615.375, 312., 45518.25, 0.),
            c32x4!(44348.75, 0., 0., 0., -1587.75, 328.625, 45194., 0.),
            c32x4!(44406.625, 0., 0., 0., -1659.125, 234.125, 45239.5, 0.),
            c32x4!(44283.375, 0., 0., 0., -1679.75, 333.75, 45244., 0.),
            c32x4!(44559.125, 0., 0., 0., -1720.875, 218.75, 45093.875, 0.),
            c32x4!(44739.75, 0., 0., 0., -1466.375, 260.625, 45420.125, 0.),
            c32x4!(44435.875, 0., 0., 0., -1580.75, 124.5, 45199.125, 0.),
            c32x4!(43523.375, 0., 0., 0., -1485.625, -1.75, 44271.5, 0.),
            c32x4!(40903.625, 0., 0., 0., -1284.625, 220.75, 41723.625, 0.),
            c32x4!(36990.375, 0., 0., 0., -1457.625, 421.25, 37751.875, 0.),
            c32x4!(32147.5, 0., 0., 0., -1041.125, 296.125, 32802.875, 0.),
            c32x4!(27394.625, 0., 0., 0., -1063.375, 251.625, 27753.875, 0.),
            c32x4!(24359.125, 0., 0., 0., -798.75, 92.125, 24723.875, 0.),
            c32x4!(23407.375, 0., 0., 0., -744.125, 129.5, 23769., 0.),
            c32x4!(25960.75, 0., 0., 0., -962., 80.25, 26264.875, 0.),
            c32x4!(30354.75, 0., 0., 0., -986., 169.25, 31005.125, 0.),
            c32x4!(35424., 0., 0., 0., -1254.5, 222.625, 35849.125, 0.),
            c32x4!(39413.375, 0., 0., 0., -1338.125, 401.75, 40161.5, 0.),
            c32x4!(41939.625, 0., 0., 0., -1544.75, 136., 43051.875, 0.),
            c32x4!(43172.375, 0., 0., 0., -1611.75, 10.375, 44115.75, 0.),
            c32x4!(43200., 0., 0., 0., -1521., 329.875, 44125.75, 0.),
            c32x4!(43182.75, 0., 0., 0., -1530.375, 89., 44050., 0.),
            c32x4!(43037., 0., 0., 0., -1376.375, 135.75, 43865.625, 0.),
            c32x4!(42772.125, 0., 0., 0., -1441., 137.5, 43939.375, 0.),
            c32x4!(42875.875, 0., 0., 0., -1457.25, 206.5, 43635.75, 0.),
            c32x4!(42945.25, 0., 0., 0., -1445.125, 153., 43423.875, 0.),
            c32x4!(43015., 0., 0., 0., -1439.5, 261.75, 43584.375, 0.),
            c32x4!(42764.625, 0., 0., 0., -1828.25, 360.25, 43191.5, 0.),
            c32x4!(42949.625, 0., 0., 0., -1736.25, 138.75, 43342.75, 0.),
            c32x4!(43141.625, 0., 0., 0., -1177.5, 343.125, 43430., 0.),
            c32x4!(42771.75, 0., 0., 0., -1530.5, 217.625, 43383.875, 0.),
            c32x4!(42705.875, 0., 0., 0., -1512.125, 255.875, 43337., 0.),
            c32x4!(42822.375, 0., 0., 0., -1252.625, 165., 43225.625, 0.),
            c32x4!(43057.75, 0., 0., 0., -1320.375, 308.5, 43377.375, 0.),
            c32x4!(42983., 0., 0., 0., -1542.125, 332., 43092.125, 0.),
            c32x4!(42990.5, 0., 0., 0., -1448.5, 167.5, 43340.625, 0.),
            c32x4!(42881.25, 0., 0., 0., -1522.875, 276.75, 43344.75, 0.),
            c32x4!(43268.875, 0., 0., 0., -1524.875, 85.75, 43745.75, 0.),
            c32x4!(43076.375, 0., 0., 0., -1603.75, 241.625, 43496.875, 0.),
            c32x4!(41966.25, 0., 0., 0., -1551.5, 13.375, 42387.75, 0.),
            c32x4!(39854., 0., 0., 0., -1449.5, 53.625, 39905.5, 0.),
            c32x4!(35783.25, 0., 0., 0., -1431.875, 138.375, 36383.125, 0.),
            c32x4!(31189.375, 0., 0., 0., -1115.125, 151.5, 31673.25, 0.),
            c32x4!(26534.25, 0., 0., 0., -977.625, 161., 26994.125, 0.),
            c32x4!(23632.75, 0., 0., 0., -766.75, 166.25, 23903.625, 0.),
            c32x4!(23300.75, 0., 0., 0., -840.25, 112.75, 23622.375, 0.),
            c32x4!(25865.875, 0., 0., 0., -972.125, 158.875, 26461.125, 0.),
            c32x4!(30336.5, 0., 0., 0., -1082., 209.125, 31190., 0.),
            c32x4!(35487.25, 0., 0., 0., -1272.125, 127.25, 35949.625, 0.),
            c32x4!(39218., 0., 0., 0., -1435.875, 331.625, 40281.125, 0.),
            c32x4!(42116.625, 0., 0., 0., -1515.375, 244.625, 42979., 0.),
            c32x4!(43435.375, 0., 0., 0., -1556.75, 378.875, 44259.25, 0.),
            c32x4!(43370.5, 0., 0., 0., -1703.375, 114., 44269.75, 0.),
            c32x4!(43299.5, 0., 0., 0., -1630.125, 166.625, 44090.375, 0.),
            c32x4!(43097.125, 0., 0., 0., -1369.375, 330., 43802.25, 0.),
            c32x4!(42906., 0., 0., 0., -1518., 162.125, 43653.375, 0.),
            c32x4!(42599.75, 0., 0., 0., -1509.125, 344.625, 43985.375, 0.),
            c32x4!(42796.375, 0., 0., 0., -1465.75, 123.125, 43842.75, 0.),
            c32x4!(42599.75, 0., 0., 0., -1535., 204.375, 43346.75, 0.),
            c32x4!(42539.875, 0., 0., 0., -1526.625, 270.5, 43585.75, 0.),
            c32x4!(42463.75, 0., 0., 0., -1640.75, 263.25, 43625., 0.),
            c32x4!(42899., 0., 0., 0., -984.5, -10.75, 43694.125, 0.),
            c32x4!(42822.5, 0., 0., 0., -1547.625, -0.875, 43066.125, 0.),
            c32x4!(42941.75, 0., 0., 0., -1415.875, 117., 43469.625, 0.),
            c32x4!(42651.625, 0., 0., 0., -1618.875, 356.75, 42897.375, 0.),
            c32x4!(42760.75, 0., 0., 0., -1558.75, 261.375, 43153.625, 0.),
            c32x4!(42822., 0., 0., 0., -1442.5, -167.25, 43397.125, 0.),
            c32x4!(42799.875, 0., 0., 0., -1592.25, 244.75, 43172.875, 0.),
            c32x4!(42949.375, 0., 0., 0., -1293.625, 160.375, 43386.125, 0.),
            c32x4!(43213.625, 0., 0., 0., -1539.25, 153.375, 43674.5, 0.),
            c32x4!(43081.75, 0., 0., 0., -1556.875, 144.5, 43494., 0.),
            c32x4!(41881.75, 0., 0., 0., -1488.375, 66.625, 42347.25, 0.),
            c32x4!(39886.125, 0., 0., 0., -1524.25, 253.5, 40027.75, 0.),
            c32x4!(36206.125, 0., 0., 0., -1154.375, 225.75, 36067.5, 0.),
            c32x4!(31212.75, 0., 0., 0., -1061.625, 65.75, 31218.125, 0.),
            c32x4!(26710.625, 0., 0., 0., -922.625, -16.375, 26826.375, 0.),
            c32x4!(23722.625, 0., 0., 0., -762.25, 108.625, 23897.75, 0.),
            c32x4!(23227.5, 0., 0., 0., -854.25, 33.75, 23544.5, 0.),
            c32x4!(26091.25, 0., 0., 0., -836.25, 127.25, 26225.5, 0.),
            c32x4!(30565.625, 0., 0., 0., -1053.125, 338.125, 30880.25, 0.),
            c32x4!(35313.875, 0., 0., 0., -1431.5, 101.375, 35798.25, 0.),
            c32x4!(39391., 0., 0., 0., -1325.625, 68.125, 40056.5, 0.),
            c32x4!(42078.625, 0., 0., 0., -1181.5, 209.5, 42596.375, 0.),
            c32x4!(43459.875, 0., 0., 0., -1576.25, 350.5, 43665.125, 0.),
            c32x4!(43568.25, 0., 0., 0., -1686.75, 152.125, 44535.625, 0.),
            c32x4!(43129., 0., 0., 0., -1571.25, 373., 44078., 0.),
            c32x4!(43025.75, 0., 0., 0., -1588.75, 135.625, 43750.5, 0.),
            c32x4!(42759.625, 0., 0., 0., -1778., 93.875, 43265.75, 0.),
            c32x4!(42894.25, 0., 0., 0., -1402.625, 249.75, 43356.125, 0.),
            c32x4!(42682.75, 0., 0., 0., -1569.125, 207.375, 43599.375, 0.),
            c32x4!(42809.125, 0., 0., 0., -1471.25, 300.875, 43304.625, 0.),
            c32x4!(42481.5, 0., 0., 0., -1591.125, 139.625, 43196., 0.),
            c32x4!(42160.5, 0., 0., 0., -1334.375, 204.625, 43095.75, 0.),
            c32x4!(42705.625, 0., 0., 0., -1232.375, 181., 43330.75, 0.),
            c32x4!(42367.625, 0., 0., 0., -1528.875, 110.75, 42924.625, 0.),
            c32x4!(42600.875, 0., 0., 0., -1602.75, 319., 43161.375, 0.),
            c32x4!(42445., 0., 0., 0., -1780.375, 316.125, 43078., 0.),
            c32x4!(42723.375, 0., 0., 0., -1478., -2.75, 43141.125, 0.),
            c32x4!(42235.125, 0., 0., 0., -1357.125, 67.875, 42905., 0.),
            c32x4!(42604.625, 0., 0., 0., -1622.25, 343.5, 43107., 0.),
            c32x4!(42851.125, 0., 0., 0., -1451.75, 299.25, 43425.125, 0.),
            c32x4!(42995.5, 0., 0., 0., -1619.875, 161.5, 43211.375, 0.),
            c32x4!(42648.875, 0., 0., 0., -1503., 57.375, 43013.875, 0.),
            c32x4!(41880.5, 0., 0., 0., -1632., 54., 41991.75, 0.),
            c32x4!(39693., 0., 0., 0., -1370.25, 311., 39807.5, 0.),
            c32x4!(35883.125, 0., 0., 0., -1217., 101.125, 36258.75, 0.),
            c32x4!(31513.75, 0., 0., 0., -975.625, 127.25, 31427.875, 0.),
            c32x4!(26745.25, 0., 0., 0., -890.75, 117.25, 26644.875, 0.),
            c32x4!(23648.125, 0., 0., 0., -830., 146., 23740.625, 0.),
            c32x4!(23979.25, 0., 0., 0., -841.875, 86.5, 24315.875, 0.),
            c32x4!(26653.25, 0., 0., 0., -917.25, 87.5, 26690.75, 0.),
            c32x4!(31401.125, 0., 0., 0., -1101., 363.625, 31534.375, 0.),
            c32x4!(36328.25, 0., 0., 0., -1347.25, 194.125, 36877.25, 0.),
            c32x4!(40466.125, 0., 0., 0., -1431.125, 131.375, 41179.625, 0.),
            c32x4!(42987.375, 0., 0., 0., -1501.75, 299.875, 43610.5, 0.),
            c32x4!(44200.5, 0., 0., 0., -1255.25, 274.875, 45191.625, 0.),
            c32x4!(44564.5, 0., 0., 0., -1613.75, 169.375, 45663.25, 0.),
            c32x4!(44709.5, 0., 0., 0., -1544.375, 279.75, 45122.125, 0.),
            c32x4!(44202.375, 0., 0., 0., -1836., 406.875, 44793.625, 0.),
            c32x4!(43981.625, 0., 0., 0., -1572.125, -60.25, 44537.75, 0.),
            c32x4!(44146.875, 0., 0., 0., -1475.125, 193., 44596.125, 0.),
            c32x4!(43959.625, 0., 0., 0., -1474.375, 256.625, 44612.125, 0.),
            c32x4!(43836.25, 0., 0., 0., -1481.125, 115., 44271.875, 0.),
            c32x4!(43612.75, 0., 0., 0., -1382.125, 285.625, 44339., 0.),
            c32x4!(43425.5, 0., 0., 0., -1345.25, 247.625, 44048.25, 0.),
            c32x4!(43908.5, 0., 0., 0., -1240.5, 256.5, 44092., 0.),
            c32x4!(43544.125, 0., 0., 0., -1384.75, 225.5, 43938.25, 0.),
            c32x4!(43806.875, 0., 0., 0., -1453., 292.75, 44023.75, 0.),
            c32x4!(43250., 0., 0., 0., -1418.625, 153.5, 44362.75, 0.),
            c32x4!(43576.625, 0., 0., 0., -1690.125, 367.75, 44258.125, 0.),
            c32x4!(43509.5, 0., 0., 0., -1598.75, 237.375, 43828.5, 0.),
            c32x4!(43606.125, 0., 0., 0., -1734.625, 246.25, 43822.625, 0.),
            c32x4!(43757.875, 0., 0., 0., -1568.125, 300.125, 43920., 0.),
            c32x4!(43968.25, 0., 0., 0., -1376.75, -104.375, 44367.125, 0.),
            c32x4!(43827., 0., 0., 0., -1372., 105.375, 43953.875, 0.),
            c32x4!(42936.375, 0., 0., 0., -1355.75, 427.875, 43205.75, 0.),
            c32x4!(40598.25, 0., 0., 0., -1500.75, 209.625, 41030.875, 0.),
            c32x4!(36670.5, 0., 0., 0., -1285.875, 96.25, 37062.625, 0.),
            c32x4!(31854.375, 0., 0., 0., -999.875, 103.75, 32281.125, 0.),
            c32x4!(27261.625, 0., 0., 0., -938.625, 122.875, 27275.375, 0.),
            c32x4!(24204.625, 0., 0., 0., -784.625, 96.25, 24436., 0.),
            c32x4!(23706.5, 0., 0., 0., -806., 76., 24002.75, 0.),
            c32x4!(26428.375, 0., 0., 0., -835., 79.5, 26757.5, 0.),
            c32x4!(31140.125, 0., 0., 0., -1082.125, 254.625, 31235.125, 0.),
            c32x4!(36282.125, 0., 0., 0., -1355.75, 176.875, 36557.125, 0.),
            c32x4!(40156.875, 0., 0., 0., -1184.75, 218.625, 40971.875, 0.),
            c32x4!(42798., 0., 0., 0., -1376.625, 270.875, 43805.875, 0.),
            c32x4!(44334., 0., 0., 0., -1395.375, 251.625, 44623.875, 0.),
            c32x4!(44419.125, 0., 0., 0., -1447.625, 243.125, 45257.25, 0.),
            c32x4!(44339.75, 0., 0., 0., -1656.75, 47.375, 44972., 0.),
            c32x4!(44307.125, 0., 0., 0., -1510.25, 255., 44800.125, 0.),
            c32x4!(43364., 0., 0., 0., -1466.25, 135., 44389.125, 0.),
            c32x4!(43638.25, 0., 0., 0., -1360.25, 207., 44658., 0.),
            c32x4!(43315.125, 0., 0., 0., -1598., 259.5, 44341.75, 0.),
            c32x4!(43427.25, 0., 0., 0., -1543.25, 301.125, 44448.5, 0.),
            c32x4!(43068.875, 0., 0., 0., -1512.25, -31., 44005.5, 0.),
            c32x4!(43072.25, 0., 0., 0., -1540.875, 85.625, 43911.25, 0.),
            c32x4!(43500.875, 0., 0., 0., -1067.75, 112.125, 44462.5, 0.),
            c32x4!(42851.625, 0., 0., 0., -1531.125, 299.625, 43651.375, 0.),
            c32x4!(43037.125, 0., 0., 0., -1387.375, 221.5, 43620.75, 0.),
            c32x4!(43112.75, 0., 0., 0., -1527.875, 208.125, 43689.125, 0.),
            c32x4!(43137.25, 0., 0., 0., -1676.125, 458.5, 43792.5, 0.),
            c32x4!(43058.875, 0., 0., 0., -1660.125, 235.375, 43751.75, 0.),
            c32x4!(42807.375, 0., 0., 0., -1390.875, 175.375, 43488.375, 0.),
            c32x4!(43432.25, 0., 0., 0., -1599.375, 400.5, 43682.25, 0.),
            c32x4!(43334.375, 0., 0., 0., -1585.125, 64.5, 44028.375, 0.),
            c32x4!(43280.125, 0., 0., 0., -1593.125, 221.375, 44009.5, 0.),
            c32x4!(42820.25, 0., 0., 0., -1447.75, 184., 42808.125, 0.),
            c32x4!(40309.375, 0., 0., 0., -1236., 69.875, 40615.25, 0.),
            c32x4!(36223.125, 0., 0., 0., -1314.5, 96.25, 37031.25, 0.),
            c32x4!(31608.5, 0., 0., 0., -1097.625, 169.25, 32154.875, 0.),
            c32x4!(26928.625, 0., 0., 0., -988.875, 78.25, 27307.875, 0.),
            c32x4!(23969.75, 0., 0., 0., -695.875, 160.125, 24360.125, 0.),
            c32x4!(24009.75, 0., 0., 0., -888.875, 99.875, 24232.875, 0.),
            c32x4!(26611.125, 0., 0., 0., -818.25, 154.625, 27060.25, 0.),
            c32x4!(31395.875, 0., 0., 0., -1015.875, 158.875, 32017.625, 0.),
            c32x4!(36504.375, 0., 0., 0., -1311.75, 128.375, 37257.75, 0.),
            c32x4!(40709., 0., 0., 0., -1480.625, 331.5, 41512.5, 0.),
            c32x4!(43526.25, 0., 0., 0., -1576.125, -3.25, 44367.125, 0.),
            c32x4!(44777., 0., 0., 0., -1539.125, 197.25, 45485.375, 0.),
            c32x4!(44802.375, 0., 0., 0., -1795., 135.625, 45818.5, 0.),
            c32x4!(44747.5, 0., 0., 0., -1597.875, 397.625, 45427.375, 0.),
            c32x4!(44310.375, 0., 0., 0., -1857.625, 394.875, 45284.125, 0.),
            c32x4!(44187.5, 0., 0., 0., -1637.5, 215.25, 45420.875, 0.),
            c32x4!(43980.25, 0., 0., 0., -1787.375, 251.625, 45162.375, 0.),
            c32x4!(44245.375, 0., 0., 0., -1746.625, 60.625, 45235.125, 0.),
            c32x4!(43790., 0., 0., 0., -1801.625, 151.5, 44991.5, 0.),
            c32x4!(44015.25, 0., 0., 0., -1654.25, 80.75, 44864.25, 0.),
            c32x4!(43380.375, 0., 0., 0., -1722.75, 165.375, 44753.625, 0.),
            c32x4!(43900.625, 0., 0., 0., -1399., 276.375, 44645.25, 0.),
            c32x4!(43075.75, 0., 0., 0., -1604.625, 126.625, 44483.125, 0.),
            c32x4!(43436.125, 0., 0., 0., -1624.125, 200.125, 44344., 0.),
            c32x4!(43392.875, 0., 0., 0., -1481.75, 93., 43947.75, 0.),
            c32x4!(43201.125, 0., 0., 0., -1714.75, 309., 44369.25, 0.),
            c32x4!(43476.75, 0., 0., 0., -1603.75, 114.625, 43938.75, 0.),
            c32x4!(43283.625, 0., 0., 0., -1544.125, 228., 44177.25, 0.),
            c32x4!(43598.375, 0., 0., 0., -1631.375, 187.75, 44300.5, 0.),
            c32x4!(43814.625, 0., 0., 0., -1625.375, 101.875, 44257.5, 0.),
            c32x4!(43641.625, 0., 0., 0., -1641.625, 264.625, 44362.875, 0.),
            c32x4!(42325.125, 0., 0., 0., -1620.75, 288.625, 43360.25, 0.),
            c32x4!(40252., 0., 0., 0., -1668., 67.875, 40949.875, 0.),
            c32x4!(36817.625, 0., 0., 0., -1421.75, 81., 37377.375, 0.),
            c32x4!(31696.75, 0., 0., 0., -1058.5, 102.625, 32489.625, 0.),
            c32x4!(27255.75, 0., 0., 0., -1012.375, 135.625, 27538.875, 0.),
            c32x4!(24108.75, 0., 0., 0., -843.625, 218.5, 24540.125, 0.),
            c32x4!(23717.375, 0., 0., 0., -860.375, 49.125, 24154.375, 0.),
            c32x4!(26591.25, 0., 0., 0., -997.5, 98.75, 27122.125, 0.),
            c32x4!(31300.625, 0., 0., 0., -1056.75, 230.875, 31959.375, 0.),
            c32x4!(36454.75, 0., 0., 0., -1149.5, 282.375, 37237.125, 0.),
            c32x4!(40559.75, 0., 0., 0., -1547.25, 85.625, 41808.125, 0.),
            c32x4!(43188.75, 0., 0., 0., -1511.75, 326.875, 44638., 0.),
            c32x4!(44551.25, 0., 0., 0., -1830.75, 185.75, 45960.625, 0.),
            c32x4!(44759.5, 0., 0., 0., -1796.125, 141.75, 46108.25, 0.),
            c32x4!(44619.5, 0., 0., 0., -1607.875, 327.625, 46041., 0.),
            c32x4!(44430.375, 0., 0., 0., -1632.25, 237.75, 45368.375, 0.),
            c32x4!(44087., 0., 0., 0., -1662.625, 349.875, 45245.5, 0.),
            c32x4!(44008.625, 0., 0., 0., -1713.125, -21.375, 44952.25, 0.),
            c32x4!(43658.125, 0., 0., 0., -1742.625, 185.375, 45330.625, 0.),
            c32x4!(43733.25, 0., 0., 0., -1449.625, 163.75, 45009.5, 0.),
            c32x4!(43407.25, 0., 0., 0., -1561.5, 247.75, 44739.625, 0.),
            c32x4!(43203.625, 0., 0., 0., -1768.75, 155.75, 44657.25, 0.),
            c32x4!(43148.75, 0., 0., 0., -1259.375, 271.875, 44534.375, 0.),
            c32x4!(43094.375, 0., 0., 0., -1518.375, 114.25, 44105.125, 0.),
            c32x4!(42903.25, 0., 0., 0., -1737., 0.75, 44320.5, 0.),
            c32x4!(42856.25, 0., 0., 0., -1665.875, 116.875, 44268.75, 0.),
            c32x4!(42910.75, 0., 0., 0., -1611.25, 168.125, 44305.875, 0.),
            c32x4!(42831.5, 0., 0., 0., -1449.375, 167., 43936.375, 0.),
            c32x4!(42584.625, 0., 0., 0., -1565.5, 156.625, 44086.25, 0.),
            c32x4!(43047.625, 0., 0., 0., -1585.125, 178., 44218.125, 0.),
            c32x4!(43117.25, 0., 0., 0., -1705.375, 60.125, 44423.5, 0.),
            c32x4!(42847.125, 0., 0., 0., -1437.125, 84.5, 43913.75, 0.),
            c32x4!(42303.125, 0., 0., 0., -1669.25, 226.375, 43010.625, 0.),
            c32x4!(39904.75, 0., 0., 0., -1515.625, 236.75, 40756.5, 0.),
            c32x4!(36179.75, 0., 0., 0., -1176., 87., 37136.375, 0.),
            c32x4!(31511.875, 0., 0., 0., -1281.25, 129.125, 32117.875, 0.),
            c32x4!(27015.875, 0., 0., 0., -904.625, 136.75, 27433.25, 0.),
            c32x4!(23923.75, 0., 0., 0., -847., 44., 24569.125, 0.),
            c32x4!(23081.875, 0., 0., 0., -924.75, 11.125, 23415.375, 0.),
            c32x4!(25732.5, 0., 0., 0., -904.25, -13.625, 26165.5, 0.),
            c32x4!(30170.625, 0., 0., 0., -1067.5, 165.75, 31094.75, 0.),
            c32x4!(35517.75, 0., 0., 0., -1351.5, 335.625, 36335.625, 0.),
            c32x4!(39351.875, 0., 0., 0., -1490., 234.125, 40515.5, 0.),
            c32x4!(42243., 0., 0., 0., -1683.5, 147.75, 43276.875, 0.),
            c32x4!(43381.625, 0., 0., 0., -1853., 131.625, 44467.25, 0.),
            c32x4!(43504., 0., 0., 0., -1708.375, 261.625, 45027.5, 0.),
            c32x4!(43109.75, 0., 0., 0., -1561.75, 123.75, 44722.125, 0.),
            c32x4!(42955.75, 0., 0., 0., -1543.5, 87.75, 44499.875, 0.),
            c32x4!(42604.625, 0., 0., 0., -1540.5, -15.125, 44156.875, 0.),
            c32x4!(42924.875, 0., 0., 0., -1433.875, 94.625, 44401.25, 0.),
            c32x4!(42293.125, 0., 0., 0., -1573.375, 23.75, 43989.625, 0.),
            c32x4!(42248.25, 0., 0., 0., -1435.125, -26.75, 43791.875, 0.),
            c32x4!(41936.875, 0., 0., 0., -1729.625, 123.75, 43332.75, 0.),
            c32x4!(41952.25, 0., 0., 0., -1602.375, 97.125, 43499.5, 0.),
            c32x4!(41999.25, 0., 0., 0., -1287., 102., 43615.25, 0.),
            c32x4!(41574.125, 0., 0., 0., -1646.25, 206.125, 43077.625, 0.),
            c32x4!(41704.375, 0., 0., 0., -1566.5, 23., 42987.25, 0.),
            c32x4!(41539.375, 0., 0., 0., -1653.625, 94.125, 42987.25, 0.),
            c32x4!(41168.75, 0., 0., 0., -1391., 28.25, 42746.375, 0.),
            c32x4!(41150.375, 0., 0., 0., -1423.375, -20.625, 42987.375, 0.),
            c32x4!(41362.625, 0., 0., 0., -1486.875, 140.75, 42386.75, 0.),
            c32x4!(41491., 0., 0., 0., -1518.5, 144.75, 42868.125, 0.),
            c32x4!(41405.375, 0., 0., 0., -1575.375, 25.875, 42750.625, 0.),
            c32x4!(41304.125, 0., 0., 0., -1592.75, 49.375, 42077.125, 0.),
            c32x4!(40587.125, 0., 0., 0., -1533.625, 129.125, 41724.75, 0.),
            c32x4!(38099.375, 0., 0., 0., -1469.25, 63.75, 39533.625, 0.),
            c32x4!(34564.5, 0., 0., 0., -1231.875, 273.5, 35629.25, 0.),
            c32x4!(30324.25, 0., 0., 0., -1113.75, 97.125, 31122.25, 0.),
            c32x4!(26008.25, 0., 0., 0., -891.75, 130.125, 26484.125, 0.),
            c32x4!(23100.625, 0., 0., 0., -807., -30.75, 23645.875, 0.),
            c32x4!(22210.5, 0., 0., 0., -701.5, -32.75, 22701.5, 0.),
            c32x4!(24793.125, 0., 0., 0., -900.375, 133.625, 25392.75, 0.),
            c32x4!(29118.125, 0., 0., 0., -957., -14., 30242.875, 0.),
            c32x4!(34204., 0., 0., 0., -1139.875, 245.625, 35363.625, 0.),
            c32x4!(37945.125, 0., 0., 0., -1416.75, 314.125, 39357., 0.),
            c32x4!(40300., 0., 0., 0., -1274., 212.5, 41993.125, 0.),
            c32x4!(41846.125, 0., 0., 0., -1549.75, 36.375, 43237.125, 0.),
            c32x4!(41766.75, 0., 0., 0., -1603.875, 158.625, 43618.375, 0.),
            c32x4!(41964.25, 0., 0., 0., -1451.75, 254.375, 43373., 0.),
            c32x4!(41482., 0., 0., 0., -1552.875, -39., 42999.875, 0.),
            c32x4!(41182.5, 0., 0., 0., -1611.25, 73.5, 42735.125, 0.),
            c32x4!(40938.875, 0., 0., 0., -1566.5, 239.375, 42691.625, 0.),
            c32x4!(40893.125, 0., 0., 0., -1616.625, 165.375, 42996.25, 0.),
            c32x4!(40604.375, 0., 0., 0., -1464.625, 262.375, 42551.125, 0.),
            c32x4!(40735.875, 0., 0., 0., -1513.125, 48.625, 42125.125, 0.),
            c32x4!(40393.5, 0., 0., 0., -1553.875, 146.625, 42179., 0.),
            c32x4!(40659.25, 0., 0., 0., -972.625, 304.75, 42357., 0.),
            c32x4!(40029.375, 0., 0., 0., -1470.125, 208.5, 41671.375, 0.),
            c32x4!(39957.25, 0., 0., 0., -1555.25, 276.5, 41839.375, 0.),
            c32x4!(39912.125, 0., 0., 0., -1564., 146.75, 41519.25, 0.),
            c32x4!(39732.125, 0., 0., 0., -1175.125, 158.5, 41238.25, 0.),
            c32x4!(39464.125, 0., 0., 0., -1433.625, 146.625, 41413.25, 0.),
            c32x4!(39208.875, 0., 0., 0., -1730.25, 185.5, 41068., 0.),
            c32x4!(39442.625, 0., 0., 0., -1563.625, 141.875, 41151.125, 0.),
            c32x4!(39770.75, 0., 0., 0., -1476.25, 71.25, 41318.5, 0.),
            c32x4!(39560.75, 0., 0., 0., -1498.5, 169.5, 41064.875, 0.),
            c32x4!(38732.875, 0., 0., 0., -1327.375, 68.375, 40003.375, 0.),
            c32x4!(36666., 0., 0., 0., -1264.25, 113.75, 37870.125, 0.),
            c32x4!(33244.5, 0., 0., 0., -1212.25, 271.875, 34313.25, 0.),
            c32x4!(28899.125, 0., 0., 0., -1015.625, 272.25, 29691.125, 0.),
            c32x4!(24948.25, 0., 0., 0., -894.25, 77.375, 25617.125, 0.),
            c32x4!(22182.625, 0., 0., 0., -682., 41.75, 22783., 0.),
            c32x4!(21292.5, 0., 0., 0., -648.75, 163.375, 21898.125, 0.),
            c32x4!(23657.25, 0., 0., 0., -802.875, 240.5, 24631.75, 0.),
            c32x4!(27857.875, 0., 0., 0., -809.25, 108.625, 28813.75, 0.),
            c32x4!(32561.25, 0., 0., 0., -1182., 158.25, 33586., 0.),
            c32x4!(36462.5, 0., 0., 0., -1322.125, 38., 37947.875, 0.),
            c32x4!(38814., 0., 0., 0., -1382.125, 237.875, 40232., 0.),
            c32x4!(40195.125, 0., 0., 0., -1591.5, 40.625, 41873.125, 0.),
            c32x4!(40165.625, 0., 0., 0., -1449.125, 121.625, 42131.625, 0.),
            c32x4!(40429., 0., 0., 0., -1545.75, 69.625, 42036.875, 0.),
            c32x4!(40089.875, 0., 0., 0., -1355.875, 130., 42130., 0.),
            c32x4!(39623.75, 0., 0., 0., -1594.75, 239.875, 41890.875, 0.),
            c32x4!(39645.875, 0., 0., 0., -1401.5, 211.25, 41870.375, 0.),
            c32x4!(39551.875, 0., 0., 0., -1609.25, 195.625, 41662.5, 0.),
            c32x4!(39452.25, 0., 0., 0., -1524.125, 156.125, 41595.375, 0.),
            c32x4!(39052.75, 0., 0., 0., -1322.5, 326.875, 41083.75, 0.),
            c32x4!(38766.125, 0., 0., 0., -1401.375, 273.625, 41159., 0.),
            c32x4!(39327.375, 0., 0., 0., -868.75, 203.25, 41122.625, 0.),
            c32x4!(38854.75, 0., 0., 0., -1063.25, 478.625, 40695.625, 0.),
            c32x4!(38504.125, 0., 0., 0., -1372.375, 285., 40731.875, 0.),
            c32x4!(38546.875, 0., 0., 0., -1344.5, 364.25, 40489.375, 0.),
            c32x4!(38382.875, 0., 0., 0., -1305.625, 322.875, 40266.875, 0.),
            c32x4!(38266.375, 0., 0., 0., -1300.625, 219.625, 39999.375, 0.),
            c32x4!(37982.5, 0., 0., 0., -1544.375, 112.875, 40020.125, 0.),
            c32x4!(38105.875, 0., 0., 0., -1221.625, 160.25, 39939., 0.),
            c32x4!(37993., 0., 0., 0., -1164.625, 248., 40425.375, 0.),
            c32x4!(37946.875, 0., 0., 0., -1296., 202.625, 40023.375, 0.),
            c32x4!(37149.375, 0., 0., 0., -1263.5, 167.625, 38929.625, 0.),
            c32x4!(35211.5, 0., 0., 0., -1235.75, 56.75, 36832.125, 0.),
            c32x4!(31929.125, 0., 0., 0., -1095.5, 217.75, 33279.625, 0.),
            c32x4!(27846.375, 0., 0., 0., -941.125, 119., 29100.125, 0.),
            c32x4!(23737.5, 0., 0., 0., -765.375, 82.375, 24896.75, 0.),
            c32x4!(21336.875, 0., 0., 0., -648.25, 82.125, 22198.25, 0.),
            c32x4!(20872.125, 0., 0., 0., -706.75, 156.75, 21526.875, 0.),
            c32x4!(23173.875, 0., 0., 0., -704.25, 104.875, 24282.25, 0.),
            c32x4!(27350.75, 0., 0., 0., -894.375, 243.125, 28443.125, 0.),
            c32x4!(31796., 0., 0., 0., -1020.375, 165.125, 33168.25, 0.),
            c32x4!(35704.625, 0., 0., 0., -1182.625, 261.75, 37239., 0.),
            c32x4!(38242.875, 0., 0., 0., -1076.125, 195.125, 39784.625, 0.),
            c32x4!(39192., 0., 0., 0., -1195.125, 164., 41089.125, 0.),
            c32x4!(39820.25, 0., 0., 0., -1352.875, 330.25, 41333.125, 0.),
            c32x4!(39285.5, 0., 0., 0., -1280.75, 230.375, 41288.75, 0.),
            c32x4!(39415.625, 0., 0., 0., -1298.625, 450., 41028.25, 0.),
            c32x4!(39414.5, 0., 0., 0., -1228.875, 168.375, 40890.5, 0.),
            c32x4!(39172.375, 0., 0., 0., -1374.5, 239.625, 40976.75, 0.),
            c32x4!(38787.875, 0., 0., 0., -1253.875, 186.25, 40868.625, 0.),
            c32x4!(38623.625, 0., 0., 0., -1475.25, 174., 40627.125, 0.),
            c32x4!(38692.125, 0., 0., 0., -1369.625, 190., 40524.625, 0.),
            c32x4!(38336.375, 0., 0., 0., -1333.875, 229.875, 39975., 0.),
            c32x4!(38177.875, 0., 0., 0., -817.75, 283.75, 40330., 0.),
            c32x4!(38015.5, 0., 0., 0., -1230.75, 238.125, 40049.375, 0.),
            c32x4!(38010., 0., 0., 0., -1286.75, 290., 39484.5, 0.),
            c32x4!(37633.125, 0., 0., 0., -1228.125, 325.75, 39605.875, 0.),
            c32x4!(37439.875, 0., 0., 0., -1080.5, 232.75, 39035.75, 0.),
            c32x4!(37266., 0., 0., 0., -1052.125, 351.625, 39196.5, 0.),
            c32x4!(37187.125, 0., 0., 0., -1175.625, 112.875, 39016., 0.),
            c32x4!(37340.625, 0., 0., 0., -1274.25, 243.5, 38876.125, 0.),
            c32x4!(37214.25, 0., 0., 0., -1210.875, 454.375, 38865.375, 0.),
            c32x4!(37008.5, 0., 0., 0., -1150.125, 371.75, 38805.375, 0.),
            c32x4!(36274., 0., 0., 0., -1191.5, 129.75, 37674.25, 0.),
            c32x4!(34554., 0., 0., 0., -971.375, 240.25, 35718.75, 0.),
            c32x4!(30987.75, 0., 0., 0., -1029.875, 134.875, 32503.625, 0.),
            c32x4!(27209.875, 0., 0., 0., -727.125, 132.875, 28130.25, 0.),
            c32x4!(23303.875, 0., 0., 0., -794.25, 55.25, 24154.125, 0.),
            c32x4!(20760.625, 0., 0., 0., -578.875, 143.875, 21663.625, 0.),
            c32x4!(20476.75, 0., 0., 0., -571.625, 78., 21166., 0.),
            c32x4!(22929.75, 0., 0., 0., -809.25, 141.5, 23623.625, 0.),
            c32x4!(27122., 0., 0., 0., -916.625, 227., 27869., 0.),
            c32x4!(31445.375, 0., 0., 0., -973.125, 202.25, 32687.375, 0.),
            c32x4!(35526.125, 0., 0., 0., -1110.125, 144.75, 36582.625, 0.),
            c32x4!(37845.5, 0., 0., 0., -989.25, 133.625, 39113.625, 0.),
            c32x4!(39025., 0., 0., 0., -1376.625, 91.875, 40224.25, 0.),
            c32x4!(39145.75, 0., 0., 0., -1322.875, 161.625, 40969.875, 0.),
            c32x4!(39191.75, 0., 0., 0., -1263.125, 445.5, 40694., 0.),
            c32x4!(39193., 0., 0., 0., -1236.375, 194.625, 40787.125, 0.),
            c32x4!(38913.375, 0., 0., 0., -1249.625, 331.5, 40460.375, 0.),
            c32x4!(38984.25, 0., 0., 0., -1343.125, 264.375, 40899., 0.),
            c32x4!(38878.5, 0., 0., 0., -1232.875, 217.625, 40415.875, 0.),
            c32x4!(38567., 0., 0., 0., -1052., 192.25, 40160.875, 0.),
            c32x4!(38675.5, 0., 0., 0., -1273.75, 317.5, 40059.625, 0.),
            c32x4!(38170.125, 0., 0., 0., -1201.125, 131.875, 40005.25, 0.),
            c32x4!(38276.375, 0., 0., 0., -832.125, 189.375, 40195.125, 0.),
            c32x4!(38185.75, 0., 0., 0., -1170., 228.75, 39551., 0.),
            c32x4!(37809.5, 0., 0., 0., -1151.75, 355.875, 39623., 0.),
            c32x4!(37579.625, 0., 0., 0., -1122.125, 262.625, 39217., 0.),
            c32x4!(37580.75, 0., 0., 0., -1238., 332., 39035.625, 0.),
            c32x4!(37303.75, 0., 0., 0., -1089., 241.625, 38993.125, 0.),
            c32x4!(36796.25, 0., 0., 0., -1176.375, 152.875, 38595.25, 0.),
            c32x4!(36940.625, 0., 0., 0., -1070.125, 287.875, 38571.125, 0.),
            c32x4!(37165.25, 0., 0., 0., -1241.75, 255.625, 38489.875, 0.),
            c32x4!(36826.25, 0., 0., 0., -1140.375, 249.125, 38594.5, 0.),
            c32x4!(36318.25, 0., 0., 0., -970.75, 188.25, 37410.75, 0.),
            c32x4!(34030.875, 0., 0., 0., -1115.625, 29.625, 35487.625, 0.),
            c32x4!(30933.625, 0., 0., 0., -1038.125, 221.875, 31930.875, 0.),
            c32x4!(26962., 0., 0., 0., -761.375, 95.5, 27880.125, 0.),
            c32x4!(23115.25, 0., 0., 0., -715.75, 29.5, 23726.75, 0.),
            c32x4!(20538., 0., 0., 0., -596.375, 171.375, 21357., 0.),
            c32x4!(20901.625, 0., 0., 0., -515.125, 96.625, 21474., 0.),
            c32x4!(23382.25, 0., 0., 0., -709.5, 251., 23938.25, 0.),
            c32x4!(27672.125, 0., 0., 0., -850.375, 92.5, 28526.75, 0.),
            c32x4!(32074.875, 0., 0., 0., -987.625, 128., 33040.875, 0.),
            c32x4!(36290.625, 0., 0., 0., -949.875, 206.625, 37058.75, 0.),
            c32x4!(38549.625, 0., 0., 0., -1226.125, 240.375, 39745.5, 0.),
            c32x4!(40019.875, 0., 0., 0., -1014.875, 254.75, 41100.625, 0.),
            c32x4!(40133.875, 0., 0., 0., -1064.875, 330.125, 41034.75, 0.),
            c32x4!(39971.375, 0., 0., 0., -1391.5, 252.75, 41480.75, 0.),
            c32x4!(40207.875, 0., 0., 0., -1321., 292.625, 41134.25, 0.),
            c32x4!(40195.625, 0., 0., 0., -1322.75, 226.375, 40976.5, 0.),
            c32x4!(39652.25, 0., 0., 0., -1120., 363.625, 41256.875, 0.),
            c32x4!(39865.5, 0., 0., 0., -1144.125, 276.125, 40875.125, 0.),
            c32x4!(39566.25, 0., 0., 0., -1165.75, 169.5, 40716.375, 0.),
            c32x4!(39814.125, 0., 0., 0., -1002.625, 81.875, 40604.25, 0.),
            c32x4!(39341.625, 0., 0., 0., -1238., 472.5, 40509.75, 0.),
            c32x4!(39119.125, 0., 0., 0., -1060.75, 280.125, 40643.375, 0.),
            c32x4!(38714.625, 0., 0., 0., -1090.25, 260.375, 40264.125, 0.),
            c32x4!(38593.375, 0., 0., 0., -1139.375, 416.875, 40045.125, 0.),
            c32x4!(38339.5, 0., 0., 0., -1247.625, 392., 39717.875, 0.),
            c32x4!(38073.75, 0., 0., 0., -1236.375, 304.375, 39718.75, 0.),
            c32x4!(37963.375, 0., 0., 0., -1162.5, 339.625, 39591.25, 0.),
            c32x4!(37854.375, 0., 0., 0., -1211.625, 194.75, 39345.5, 0.),
            c32x4!(37863.5, 0., 0., 0., -1185.125, 483.625, 39073.375, 0.),
            c32x4!(37713.25, 0., 0., 0., -995., 218.5, 39286.75, 0.),
            c32x4!(37537.875, 0., 0., 0., -1099.875, 235.5, 38776.875, 0.),
            c32x4!(36682.875, 0., 0., 0., -1121.75, 237.5, 37668.875, 0.),
            c32x4!(34755.625, 0., 0., 0., -1070., 201.125, 35748.25, 0.),
            c32x4!(31552.5, 0., 0., 0., -848.75, 157.5, 32339., 0.),
            c32x4!(27156.625, 0., 0., 0., -785.125, 236.625, 28056.5, 0.),
            c32x4!(23478.25, 0., 0., 0., -680.75, 139.125, 23972.375, 0.),
            c32x4!(20946.875, 0., 0., 0., -612.25, 207.5, 21621.375, 0.),
            c32x4!(20486.75, 0., 0., 0., -511.75, 125.375, 20955.25, 0.),
            c32x4!(23034.125, 0., 0., 0., -569.5, 193.75, 23402.75, 0.),
            c32x4!(26946.875, 0., 0., 0., -805., 211.625, 27714.625, 0.),
            c32x4!(31502.125, 0., 0., 0., -757.75, 342.125, 32461., 0.),
            c32x4!(35248.875, 0., 0., 0., -919.25, 322.375, 36263.375, 0.),
            c32x4!(37940.375, 0., 0., 0., -1281.875, 233., 38731.25, 0.),
            c32x4!(39152.75, 0., 0., 0., -1276.625, 382.5, 40115.875, 0.),
            c32x4!(39505., 0., 0., 0., -1291.625, 196.125, 40566.875, 0.),
            c32x4!(39156.5, 0., 0., 0., -1312.5, 214., 40381.625, 0.),
            c32x4!(39178.625, 0., 0., 0., -1348., 389.375, 40314.125, 0.),
            c32x4!(39445.875, 0., 0., 0., -1095., 414.625, 40342., 0.),
            c32x4!(39051.125, 0., 0., 0., -1053.75, 242.75, 40364.25, 0.),
            c32x4!(39136.625, 0., 0., 0., -1247.625, 74.25, 40320.25, 0.),
            c32x4!(38888.875, 0., 0., 0., -1114., 187., 40140.625, 0.),
            c32x4!(38760.625, 0., 0., 0., -1308.625, 268., 39711.375, 0.),
            c32x4!(38410.5, 0., 0., 0., -1210.25, 315.125, 39748.25, 0.),
            c32x4!(38778.125, 0., 0., 0., -933.625, 325.375, 39879.625, 0.),
            c32x4!(37947.75, 0., 0., 0., -1097.875, 240.125, 39385.875, 0.),
            c32x4!(38196.25, 0., 0., 0., -1092.625, 332.25, 39152.875, 0.),
            c32x4!(37802.25, 0., 0., 0., -946.25, 128.25, 39071.5, 0.),
            c32x4!(37673.75, 0., 0., 0., -1169.5, 192.75, 39047.375, 0.),
            c32x4!(37350.25, 0., 0., 0., -1023.25, 218.5, 38499.5, 0.),
            c32x4!(37144.25, 0., 0., 0., -1025.125, 273.625, 38320.625, 0.),
            c32x4!(36945.5, 0., 0., 0., -1121.375, 289.25, 38243., 0.),
            c32x4!(37061., 0., 0., 0., -962., 401.125, 38187.125, 0.),
            c32x4!(36784.125, 0., 0., 0., -1039.75, 234.375, 37829.75, 0.),
            c32x4!(35775.375, 0., 0., 0., -1064.375, 294.875, 37048.625, 0.),
            c32x4!(33761.625, 0., 0., 0., -937., 215.25, 35000.625, 0.),
            c32x4!(30812.5, 0., 0., 0., -869., 89.875, 31602.875, 0.),
            c32x4!(26727.25, 0., 0., 0., -792.625, 224.25, 27404.875, 0.),
            c32x4!(22829.125, 0., 0., 0., -582.875, 246.5, 23520.125, 0.),
            c32x4!(20530.5, 0., 0., 0., -508.625, 67.625, 20960.125, 0.),
            c32x4!(20425.875, 0., 0., 0., -532.875, 183.875, 20745.375, 0.),
            c32x4!(22724., 0., 0., 0., -554.875, 91.25, 23449., 0.),
            c32x4!(27107.75, 0., 0., 0., -839.5, 137.5, 27403.625, 0.),
            c32x4!(31286.5, 0., 0., 0., -770.5, 238.875, 31968.875, 0.),
            c32x4!(35360.375, 0., 0., 0., -1203.625, 304.5, 36087.375, 0.),
            c32x4!(37582.875, 0., 0., 0., -1188.125, 366.5, 38424.125, 0.),
            c32x4!(39162.125, 0., 0., 0., -1154.875, 103.375, 39894.25, 0.),
            c32x4!(39580.125, 0., 0., 0., -1270.5, 331.875, 39952.375, 0.),
            c32x4!(39322.5, 0., 0., 0., -1264.75, 104.75, 40279., 0.),
            c32x4!(39065.25, 0., 0., 0., -1179.75, 195.625, 40267.625, 0.),
            c32x4!(39436.25, 0., 0., 0., -1173.25, 208.5, 40401.875, 0.),
            c32x4!(39104.375, 0., 0., 0., -1101.5, 239.25, 40385., 0.),
            c32x4!(39156.375, 0., 0., 0., -1113.75, 326.25, 40516.625, 0.),
            c32x4!(38990.25, 0., 0., 0., -1016., 273.75, 40125.125, 0.),
            c32x4!(39029.625, 0., 0., 0., -1018.75, 220.5, 40121.75, 0.),
            c32x4!(38607.75, 0., 0., 0., -1068.5, 268.5, 40223., 0.),
            c32x4!(38894.625, 0., 0., 0., -792.5, 233.125, 40477.125, 0.),
            c32x4!(38466.875, 0., 0., 0., -1124.875, 189.25, 39694.5, 0.),
            c32x4!(38330.75, 0., 0., 0., -1259.75, 60.375, 39960., 0.),
            c32x4!(38145.875, 0., 0., 0., -1160.5, 206., 39353., 0.),
            c32x4!(37785.5, 0., 0., 0., -960.125, 212.875, 39277.75, 0.),
            c32x4!(37801.125, 0., 0., 0., -1256.5, 210., 39082.875, 0.),
            c32x4!(37463.125, 0., 0., 0., -1062.875, 185., 38938.25, 0.),
            c32x4!(37355.375, 0., 0., 0., -1012.5, 265.125, 38750.625, 0.),
            c32x4!(37237., 0., 0., 0., -1015.5, 302.625, 38729.875, 0.),
            c32x4!(37015.5, 0., 0., 0., -1045.875, 173.125, 38440.75, 0.),
            c32x4!(36215.25, 0., 0., 0., -879.5, 130.875, 37396.25, 0.),
            c32x4!(33947.75, 0., 0., 0., -1000.25, 262.625, 35013.75, 0.),
            c32x4!(31005.75, 0., 0., 0., -845.875, 142.5, 31924.125, 0.),
            c32x4!(26871.375, 0., 0., 0., -662.125, 160.875, 27870.875, 0.),
            c32x4!(22903.75, 0., 0., 0., -622.875, 144.375, 23613., 0.),
            c32x4!(20521.125, 0., 0., 0., -592.5, 75., 21068.75, 0.),
            c32x4!(19873.75, 0., 0., 0., -582.125, 149.25, 20620.25, 0.),
            c32x4!(22330.625, 0., 0., 0., -702., 111.875, 23059.875, 0.),
            c32x4!(26352.5, 0., 0., 0., -738.75, 44.375, 27162.5, 0.),
            c32x4!(30681.25, 0., 0., 0., -780.125, -16., 31598.25, 0.),
            c32x4!(34520.75, 0., 0., 0., -915., 407.125, 35515., 0.),
            c32x4!(36988.875, 0., 0., 0., -1053.875, 334., 37891.75, 0.),
            c32x4!(38021.875, 0., 0., 0., -1202.25, -32.125, 39332.375, 0.),
            c32x4!(38465.5, 0., 0., 0., -1145.25, 282., 39423.875, 0.),
            c32x4!(38619.25, 0., 0., 0., -1076.5, 178.75, 39824.125, 0.),
            c32x4!(38447.25, 0., 0., 0., -1086.375, 235.25, 39635., 0.),
            c32x4!(38553.625, 0., 0., 0., -1080.875, 290.875, 39434.5, 0.),
            c32x4!(38195.75, 0., 0., 0., -977.125, 252.5, 39704.875, 0.),
            c32x4!(38304.125, 0., 0., 0., -1068.875, 12., 39712.625, 0.),
            c32x4!(38358.875, 0., 0., 0., -1327.625, 410.75, 39837., 0.),
            c32x4!(38366.625, 0., 0., 0., -1039.25, 265., 39414.875, 0.),
            c32x4!(37972.375, 0., 0., 0., -1131.25, 249.875, 39605., 0.),
            c32x4!(38281.125, 0., 0., 0., -740.5, -63.75, 39617.5, 0.),
            c32x4!(37740.375, 0., 0., 0., -1013.5, 110.25, 39177.75, 0.),
            c32x4!(37487.875, 0., 0., 0., -1103.875, 348.125, 39241.25, 0.),
            c32x4!(37179., 0., 0., 0., -1034.375, 168.375, 38975.5, 0.),
            c32x4!(37192.875, 0., 0., 0., -950.125, 16., 38723.875, 0.),
            c32x4!(36680.625, 0., 0., 0., -953.75, 124.75, 38589.375, 0.),
            c32x4!(36685.5, 0., 0., 0., -907.375, 169.875, 38232.25, 0.),
            c32x4!(36529.625, 0., 0., 0., -992.625, 160.125, 38497.25, 0.),
            c32x4!(36323.75, 0., 0., 0., -1259.5, -11.75, 38172.625, 0.),
            c32x4!(35959.875, 0., 0., 0., -865.625, 227.5, 37912.125, 0.),
            c32x4!(35222.625, 0., 0., 0., -1068.875, 251.375, 36980.125, 0.),
            c32x4!(33109.625, 0., 0., 0., -843.25, 185.125, 34959.75, 0.),
            c32x4!(30121., 0., 0., 0., -867., 179.75, 31412.125, 0.),
            c32x4!(26028.5, 0., 0., 0., -693.875, 161.625, 27213., 0.),
            c32x4!(22304.875, 0., 0., 0., -546.625, 142.75, 23279.5, 0.),
            c32x4!(20020.625, 0., 0., 0., -495.625, 95.375, 20724., 0.),
            c32x4!(19609.375, 0., 0., 0., -464.75, 43., 20349.875, 0.),
            c32x4!(21877.375, 0., 0., 0., -473.25, 99.375, 22800.75, 0.),
            c32x4!(25727.75, 0., 0., 0., -786.625, 165.875, 26834.625, 0.),
            c32x4!(30171.5, 0., 0., 0., -800.75, 321.625, 31180.625, 0.),
            c32x4!(33535.5, 0., 0., 0., -1101.625, 222.625, 34919.875, 0.),
            c32x4!(36049.75, 0., 0., 0., -935.75, 243.125, 37195.5, 0.),
            c32x4!(37357.75, 0., 0., 0., -1111.125, 194.625, 38766.75, 0.),
            c32x4!(37662.625, 0., 0., 0., -1092.75, 331.5, 39024.5, 0.),
            c32x4!(37647.75, 0., 0., 0., -1086.375, 222.75, 39228., 0.),
            c32x4!(37814.125, 0., 0., 0., -1130.5, 305.5, 39100.875, 0.),
            c32x4!(37848.75, 0., 0., 0., -1038.875, 62.5, 39073.125, 0.),
            c32x4!(37892.75, 0., 0., 0., -1153.75, 85.125, 39194.625, 0.),
            c32x4!(38280.875, 0., 0., 0., -1033.875, 290.875, 39222.625, 0.),
            c32x4!(37939.125, 0., 0., 0., -945., 161., 39357.375, 0.),
            c32x4!(37728.875, 0., 0., 0., -1101.5, 107.5, 39338.625, 0.),
            c32x4!(37554., 0., 0., 0., -1205.875, -58., 39034.125, 0.),
            c32x4!(37910.625, 0., 0., 0., -830.625, 395.625, 39308.375, 0.),
            c32x4!(37290.5, 0., 0., 0., -994.625, 27.75, 39024.25, 0.),
            c32x4!(37251.625, 0., 0., 0., -950.125, 226., 38815.875, 0.),
            c32x4!(37189.625, 0., 0., 0., -1058.625, 253.875, 38753.25, 0.),
            c32x4!(37327., 0., 0., 0., -1132.125, 373.375, 38510., 0.),
            c32x4!(36553.375, 0., 0., 0., -742.625, 238.25, 38211.5, 0.),
            c32x4!(36380.25, 0., 0., 0., -1003.625, 376.625, 37757.25, 0.),
            c32x4!(36311.875, 0., 0., 0., -929.875, 273.875, 37697.75, 0.),
            c32x4!(35952.25, 0., 0., 0., -990.625, 78.625, 37993.25, 0.),
            c32x4!(35701., 0., 0., 0., -1025.25, 192.25, 37368.25, 0.),
            c32x4!(34809., 0., 0., 0., -1057., 141.375, 36447.125, 0.),
            c32x4!(32742.75, 0., 0., 0., -880.5, 230.25, 34409.375, 0.),
            c32x4!(29643.25, 0., 0., 0., -582.375, 218.25, 30782.5, 0.),
            c32x4!(25734.375, 0., 0., 0., -746.75, 85.25, 26752.25, 0.),
            c32x4!(22020.875, 0., 0., 0., -647.375, 53., 22854.875, 0.),
            c32x4!(19723.5, 0., 0., 0., -447.25, 138.375, 20598.75, 0.),
            c32x4!(18779.375, 0., 0., 0., -412.75, 58.75, 19629.125, 0.),
            c32x4!(20881.75, 0., 0., 0., -527.125, 123.5, 21877.5, 0.),
            c32x4!(24721., 0., 0., 0., -680.125, 155.875, 25662.625, 0.),
            c32x4!(28741.875, 0., 0., 0., -703.125, 140.875, 29785.25, 0.),
            c32x4!(32157.875, 0., 0., 0., -843., 96.875, 33580.5, 0.),
            c32x4!(34379.875, 0., 0., 0., -792.875, 308.5, 35707.75, 0.),
            c32x4!(35601.875, 0., 0., 0., -894.75, 40.5, 37125., 0.),
            c32x4!(35875.5, 0., 0., 0., -965.25, 283., 37364.375, 0.),
            c32x4!(36192.875, 0., 0., 0., -860.25, 178.625, 37416.25, 0.),
            c32x4!(36028.375, 0., 0., 0., -868., 120.875, 37331., 0.),
            c32x4!(36209.75, 0., 0., 0., -983.125, 284.5, 37488.375, 0.),
            c32x4!(36209.25, 0., 0., 0., -758.25, 239.625, 37368.5, 0.),
            c32x4!(36492.75, 0., 0., 0., -1011.5, 120.25, 37770.25, 0.),
            c32x4!(36316.875, 0., 0., 0., -946.875, 124.625, 37745.875, 0.),
            c32x4!(36146.25, 0., 0., 0., -973.375, 315.25, 37622., 0.),
            c32x4!(35863.25, 0., 0., 0., -910.875, 82.625, 37283.25, 0.),
            c32x4!(36186.625, 0., 0., 0., -617.75, 248.5, 37650.5, 0.),
            c32x4!(35771.625, 0., 0., 0., -842.75, 124.375, 37294.625, 0.),
            c32x4!(35511.625, 0., 0., 0., -810.75, 258., 37246.375, 0.),
            c32x4!(35367.125, 0., 0., 0., -983.125, 266.625, 37181.625, 0.),
            c32x4!(35178.25, 0., 0., 0., -1000.125, 261.5, 37217., 0.),
            c32x4!(34815.25, 0., 0., 0., -901.375, 157.375, 36909., 0.),
            c32x4!(34829.5, 0., 0., 0., -867.125, 135., 36729.75, 0.),
            c32x4!(34579.375, 0., 0., 0., -909.25, 124.375, 36442.875, 0.),
            c32x4!(34327.875, 0., 0., 0., -748., 147.25, 36361.75, 0.),
            c32x4!(34155.375, 0., 0., 0., -909.375, 164., 35979., 0.),
            c32x4!(33196.875, 0., 0., 0., -826.5, 242.5, 34978.25, 0.),
            c32x4!(31463.5, 0., 0., 0., -721., 53.5, 33162.5, 0.),
            c32x4!(28267., 0., 0., 0., -628.375, 185.75, 29799.375, 0.),
            c32x4!(24576.75, 0., 0., 0., -530.875, 195.5, 25937.5, 0.),
            c32x4!(20988.75, 0., 0., 0., -393.625, 87.125, 22245.375, 0.),
            c32x4!(18847.125, 0., 0., 0., -382.25, 109.375, 19831.75, 0.),
            c32x4!(18083.5, 0., 0., 0., -403.875, -52.625, 18720.875, 0.),
            c32x4!(20019.875, 0., 0., 0., -444.5, 81.75, 20919.875, 0.),
            c32x4!(23688., 0., 0., 0., -536.875, 88.875, 24603.625, 0.),
            c32x4!(27465., 0., 0., 0., -629.125, 189.25, 28677.5, 0.),
            c32x4!(30625.625, 0., 0., 0., -678.25, 171.125, 31969., 0.),
            c32x4!(32780.75, 0., 0., 0., -639.625, 182.25, 33969.875, 0.),
            c32x4!(33722.875, 0., 0., 0., -805.875, 125.5, 35451.5, 0.),
            c32x4!(34173.875, 0., 0., 0., -808.5, 96., 35879.75, 0.),
            c32x4!(34272.25, 0., 0., 0., -941.25, 108.25, 35837.625, 0.),
            c32x4!(34365.5, 0., 0., 0., -700.375, 154.375, 35810.625, 0.),
            c32x4!(34220.75, 0., 0., 0., -882., 199.375, 35714.625, 0.),
            c32x4!(34378.5, 0., 0., 0., -706.875, 152.375, 35887.75, 0.),
            c32x4!(34536., 0., 0., 0., -852.625, 402.75, 36140.125, 0.),
            c32x4!(34348., 0., 0., 0., -864.125, 153., 36012.25, 0.),
            c32x4!(34319.625, 0., 0., 0., -736.875, 338.25, 35969.25, 0.),
            c32x4!(34384.625, 0., 0., 0., -594.5, 149.375, 35886.375, 0.),
            c32x4!(34625.25, 0., 0., 0., -365.375, 208.125, 36315.25, 0.),
            c32x4!(34305.75, 0., 0., 0., -741.125, 236.125, 35855.25, 0.),
            c32x4!(34246.875, 0., 0., 0., -694.625, 109.25, 35866.25, 0.),
            c32x4!(33968.125, 0., 0., 0., -836.25, 148., 35625.5, 0.),
            c32x4!(33876.375, 0., 0., 0., -762.125, 119., 35671., 0.),
            c32x4!(33568., 0., 0., 0., -666., 200.625, 35172.125, 0.),
            c32x4!(33211.75, 0., 0., 0., -603.5, 99.625, 35135.25, 0.),
            c32x4!(33290.75, 0., 0., 0., -690.5, 132.5, 34923.25, 0.),
            c32x4!(33162.75, 0., 0., 0., -808.75, 170.875, 35102.625, 0.),
            c32x4!(32798.375, 0., 0., 0., -678., 300.5, 34476.75, 0.),
            c32x4!(31804.875, 0., 0., 0., -623.875, 39.25, 33628.875, 0.),
            c32x4!(30186.625, 0., 0., 0., -652.25, 175.75, 31692.625, 0.),
            c32x4!(27166.875, 0., 0., 0., -554.375, 222., 28661., 0.),
            c32x4!(23605.125, 0., 0., 0., -534.875, 68.25, 24807.25, 0.),
            c32x4!(20214.625, 0., 0., 0., -464.625, 101.375, 21283.25, 0.),
            c32x4!(18031.625, 0., 0., 0., -327.125, 44., 18923.75, 0.),
        ]);

        let flags = Array::from_elem((768, 4), true);
        let weights = Array::from_elem((4,), 6144_f32);

        let mut main_table = Table::open(&table_path, TableOpenMode::ReadWrite).unwrap();

        ms_writer
            .write_main_row(
                &mut main_table,
                5077351976.000001,
                5077351976.000001,
                0,
                0,
                0,
                &vec![0., 0., 0.],
                2.,
                -1,
                1,
                -1,
                &vec![1., 1., 1., 1.],
                data,
                flags,
                weights,
            )
            .unwrap();
        drop(ms_writer);

        let mut main_table = Table::open(&table_path, TableOpenMode::Read).unwrap();
        let mut expected_table =
            Table::open(PATH_1254670392.join(""), TableOpenMode::Read).unwrap();

        assert_table_nrows_match!(main_table, expected_table);
        for col_name in [
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
        ] {
            assert_table_columns_match!(main_table, expected_table, col_name);
        }
    }
}
