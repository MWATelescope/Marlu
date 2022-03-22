use hifitime::{Duration, Epoch, TimeSeries};
use ndarray::Array2;

use crate::{LatLngHeight, RADec, XyzGeocentric, XyzGeodetic, ENH};

cfg_if::cfg_if! {
    if #[cfg(feature = "mwalib")] {
        use std::ops::Range;
        use mwalib::{CorrelatorContext, MetafitsContext};
        use hifitime::Unit::Millisecond;
        use itertools::izip;
        use ndarray::array;
    }
}

/// A container for observation metadata common across most file types
pub struct ObsContext {
    /// Scheduled start time
    pub sched_start_timestamp: Epoch,

    /// Scheduled duration
    pub sched_duration: Duration,

    /// Observation name
    pub name: Option<String>,

    /// Name of field being observed
    pub field_name: Option<String>,

    /// The project ID of the observation
    pub project_id: Option<String>,

    /// The observer or creator of the observation
    pub observer: Option<String>,

    /// The phase centre.
    pub phase_centre: RADec,

    /// The pointing centre.
    pub pointing_centre: Option<RADec>,

    /// The Earth position of the instrumental array
    pub array_pos: LatLngHeight,

    /// TODO: store in ENH or geodetic?
    /// The geodetic position of each antenna.
    // pub tiles_xyz_geod: Vec<XyzGeodetic>,
    /// The east-north-height position of each antenna
    pub ant_positions_enh: Vec<ENH>,

    /// The name of each antenna / tile.
    pub ant_names: Vec<String>,
}

impl ObsContext {
    #[cfg(feature = "mwalib")]
    pub fn from_mwalib(meta_ctx: &MetafitsContext) -> Self {
        let obs_name = meta_ctx.obs_name.clone();
        let field_name: String = obs_name
            .rsplit_once('_')
            .unwrap_or((obs_name.as_str(), ""))
            .0
            .into();

        let ants = &meta_ctx.antennas;
        let mut ant_positions_enh = Vec::<ENH>::with_capacity(ants.len());
        let mut ant_names = Vec::<String>::with_capacity(ants.len());
        for ant in ants {
            ant_positions_enh.push(ENH {
                e: ant.east_m,
                n: ant.north_m,
                h: ant.height_m,
            });
            ant_names.push(ant.tile_name.clone());
        }

        Self {
            sched_start_timestamp: Epoch::from_gpst_seconds(
                meta_ctx.sched_start_gps_time_ms as f64 / 1e3,
            ),
            sched_duration: Duration::from_f64(meta_ctx.sched_duration_ms as f64, Millisecond),
            name: Some(obs_name),
            field_name: Some(field_name),
            project_id: Some(meta_ctx.project_id.clone()),
            observer: Some(meta_ctx.creator.clone()),
            phase_centre: RADec::from_mwalib_phase_or_pointing(meta_ctx),
            pointing_centre: Some(RADec::from_mwalib_tile_pointing(meta_ctx)),
            array_pos: LatLngHeight::new_mwa(),
            ant_positions_enh,
            ant_names,
        }
    }

    pub fn ant_positions_geodetic(&self) -> Vec<XyzGeodetic> {
        self.ant_positions_enh
            .iter()
            .map(|enh| enh.to_xyz(self.array_pos.latitude_rad))
            .collect()
    }

    pub fn ant_positions_geocentric(&self) -> Vec<XyzGeocentric> {
        self.ant_positions_enh
            .iter()
            .map(|enh| {
                enh.to_xyz(self.array_pos.latitude_rad)
                    .to_geocentric(self.array_pos)
                    .unwrap()
            })
            .collect()
    }

    pub fn num_ants(&self) -> usize {
        self.ant_positions_enh.len()
    }
}

/// An extension of [`ObsContext`] that for MWA-specific metadata that is not
/// present in some file types like uvfits.
pub struct MwaObsContext {
    /// Antenna input numbers. [ant_idx][pol]
    pub ant_inputs: Array2<usize>,

    /// Antenna tile numbers
    pub ant_numbers: Vec<usize>,

    /// Antenna receiver numbers
    pub ant_receivers: Vec<usize>,

    /// Antenna slot numbers. [ant_idx][pol]
    pub ant_slots: Array2<usize>,

    /// Antenna slot numbers. [ant_idx][pol]
    pub ant_cable_lengths: Array2<f64>,

    /// Coarse Channel Receiver Numbers
    pub coarse_chan_recs: Vec<usize>,

    /// Whether the observation has a calibrator
    pub has_calibrator: bool,

    /// MWA Observation Mode. See: [`mwalib::MetafitsContext::obs_mode`]
    pub mode: String,

    /// Tile pointing delays
    pub delays: Vec<u32>,
}

impl MwaObsContext {
    #[cfg(feature = "mwalib")]
    pub fn from_mwalib(meta_ctx: &MetafitsContext) -> Self {
        let ants = &meta_ctx.antennas;

        let mut result = Self {
            ant_inputs: Array2::zeros((ants.len(), 2)),
            ant_numbers: vec![0; ants.len()],
            ant_receivers: vec![0; ants.len()],
            ant_slots: Array2::zeros((ants.len(), 2)),
            ant_cable_lengths: Array2::zeros((ants.len(), 2)),
            coarse_chan_recs: meta_ctx
                .metafits_coarse_chans
                .iter()
                .map(|c| c.rec_chan_number)
                .collect(),
            has_calibrator: meta_ctx.calibrator,
            mode: meta_ctx.mode.to_string(),
            delays: meta_ctx.delays.clone(),
        };

        for (ant, mut input, number, receiver, mut slot, mut length) in izip!(
            ants,
            result.ant_inputs.outer_iter_mut(),
            result.ant_numbers.iter_mut(),
            result.ant_receivers.iter_mut(),
            result.ant_slots.outer_iter_mut(),
            result.ant_cable_lengths.outer_iter_mut(),
        ) {
            let (rf_x, rf_y) = (ant.rfinput_x.clone(), ant.rfinput_y.clone());
            input.assign(&array![rf_x.input as usize, rf_y.input as _]);
            *number = ant.tile_id as _;
            *receiver = rf_x.rec_number as _;
            slot.assign(&array![
                rf_x.rec_slot_number as usize,
                rf_y.rec_slot_number as _
            ]);
            length.assign(&array![
                rf_x.electrical_length_m as f64,
                rf_y.electrical_length_m as _
            ]);
        }

        result
    }
}

/// A lightweight container for correlator visibility metadata used in Marlu operations.
///
/// This is intended to describe an accompanying visibility and weight ndarray.
///
/// This is similar to, and can be derived from [`mwalib::CorrelatorContext`],
/// but applies to a specific selection of baselines from a contiguous range of
/// timesteps and frequencies with a given configuration of averaging, and other
/// pre-processing settings.
///
/// A `VisContext` is oblivious to mwalib concepts like coarse channels.
pub struct VisContext {
    /// The number of selected timesteps (Axis 0) in the accompanying visibility and weight ndarrays.
    pub num_sel_timesteps: usize,
    /// The timestamp at the start of the first selected pre-averaging timestep
    pub start_timestamp: Epoch,
    /// Duration between each pre-averaging timestep [milliseconds]
    pub int_time: Duration,
    /// The number of selected channels (Axis 1) in the accompanying visibility and weight ndarrays.
    pub num_sel_chans: usize,
    /// The centre frequency of the first selected pre-averaging channel [Hz]
    pub start_freq_hz: f64,
    /// The bandwidth between each pre-averaging channel [Hz]
    pub freq_resolution_hz: f64,
    /// The tile index pairs for each selected baseline
    pub sel_baselines: Vec<(usize, usize)>,
    /// Time averaging factor
    pub avg_time: usize,
    /// Frequency averaging factor
    pub avg_freq: usize,
    /// Number of polarisation combinations in the visibilities e.g. XX,XY,YX,YY == 4
    pub num_vis_pols: usize,
}

impl VisContext {
    #[cfg(feature = "mwalib")]
    pub fn from_mwalib(
        corr_ctx: &CorrelatorContext,
        timestep_range: &Range<usize>,
        coarse_chan_range: &Range<usize>,
        baseline_idxs: &[usize],
        avg_time: usize,
        avg_freq: usize,
    ) -> Self {
        // Time axis
        let num_sel_timesteps = timestep_range.len();

        let start_timestamp = Epoch::from_gpst_seconds(
            corr_ctx.timesteps[timestep_range.start].gps_time_ms as f64 / 1e3,
        );

        let int_time =
            Duration::from_f64(corr_ctx.metafits_context.corr_int_time_ms as _, Millisecond);

        // Frequency axis
        let num_sel_coarse_chans = coarse_chan_range.len();
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let num_sel_chans = fine_chans_per_coarse * num_sel_coarse_chans;
        let start_freq_hz = corr_ctx.metafits_context.metafits_fine_chan_freqs_hz
            [coarse_chan_range.start * fine_chans_per_coarse];
        let freq_resolution_hz = corr_ctx.metafits_context.corr_fine_chan_width_hz as f64;

        // baseline axis
        let sel_baselines = baseline_idxs
            .iter()
            .map(|&idx| {
                let baseline = &corr_ctx.metafits_context.baselines[idx];
                (baseline.ant1_index, baseline.ant2_index)
            })
            .collect::<Vec<_>>();

        let num_vis_pols = corr_ctx.metafits_context.num_visibility_pols;

        VisContext {
            num_sel_timesteps,
            start_timestamp,
            int_time,
            num_sel_chans,
            start_freq_hz,
            freq_resolution_hz,
            sel_baselines,
            avg_time,
            avg_freq,
            num_vis_pols,
        }
    }

    /// The expected dimensions of the visibility and weight ndarray selection.
    pub fn sel_dims(&self) -> (usize, usize, usize) {
        (
            self.num_sel_timesteps,
            self.num_sel_chans,
            self.sel_baselines.len(),
        )
    }

    /// The expected dimensions of the visibility and weight ndarray selection.
    pub fn avg_dims(&self) -> (usize, usize, usize) {
        (
            self.num_avg_timesteps(),
            self.num_avg_chans(),
            self.sel_baselines.len(),
        )
    }

    /// Whether this context corresponds with the the case where the time and
    /// frequency averaging factors are both 1, and no averaging is performed.
    pub fn trivial_averaging(&self) -> bool {
        self.avg_time == 1 && self.avg_freq == 1
    }

    /// The number of timesteps in the post-averaging time dimension
    pub fn num_avg_timesteps(&self) -> usize {
        (self.num_sel_timesteps as f64 / self.avg_time as f64).ceil() as usize
    }

    /// The integration time of the post-averaging data.
    pub fn avg_int_time(&self) -> Duration {
        self.int_time * (self.avg_time as i64)
    }

    pub fn timeseries(&self, averaging: bool, centroid: bool) -> TimeSeries {
        let (num_timesteps, int_time) = if averaging {
            (self.num_avg_timesteps(), self.avg_int_time())
        } else {
            (self.num_sel_timesteps, self.int_time)
        };
        let offset = if centroid { 0.5 } else { 0.0 };
        let start_timestamp = self.start_timestamp + offset * int_time;
        TimeSeries::inclusive(
            start_timestamp,
            start_timestamp + (num_timesteps as f64) * int_time,
            int_time,
        )
    }

    /// The number of channels in the post-averaging frequency dimension
    pub fn num_avg_chans(&self) -> usize {
        (self.num_sel_chans as f64 / self.avg_freq as f64).ceil() as usize
    }

    /// The frequency resolution after averaging
    pub fn avg_freq_resolution_hz(&self) -> f64 {
        self.freq_resolution_hz * self.avg_freq as f64
    }

    /// An iterator over all selected frequencies
    ///
    /// TODO: iterator return type?
    pub fn frequencies_hz(&self) -> Vec<f64> {
        (0..self.num_sel_chans)
            .map(|i| self.start_freq_hz + i as f64 * self.freq_resolution_hz)
            .collect()
    }

    /// An iterator over averaged frequencies
    ///
    /// TODO: iterator return type?
    pub fn avg_frequencies_hz(&self) -> Vec<f64> {
        self.frequencies_hz()
            .chunks(self.avg_freq)
            .map(|chunk| chunk.iter().sum::<f64>() / chunk.len() as f64)
            .collect()
    }

    /// Get the weight factor: a measure of the resolution relative to the base
    /// resolution of the legacy MWA correlator (1s / 10kHZ).
    ///
    /// This is a conceptfrom Cotter, and the legacy MWA correlator where the value
    /// is a multiple of the frequency resolution (relative to 10kHz), and the
    /// time averaging factor (relative to 1s).
    pub fn weight_factor(&self) -> f64 {
        self.freq_resolution_hz * self.int_time.in_seconds() / 10000.0
    }
}
