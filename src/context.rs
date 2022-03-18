use hifitime::{Duration, Epoch, TimeSeries};

use crate::XyzGeodetic;

cfg_if::cfg_if! {
    if #[cfg(feature = "mwalib")] {
        use std::ops::Range;
        use mwalib::CorrelatorContext;
        use hifitime::Unit::Millisecond;
    }
}

// use crate::{LatLngHeight, RADec};
// /// A lightweight container for observation metadata used in Marlu operations.
// pub struct MarluObsContext {
//     /// The observation ID, which is also the observation's scheduled start GPS
//     /// time.
//     pub(crate) obsid: Option<u32>,

//     array_pos: LatLngHeight,
//     phase_centre: RADec,
//     ...
// }

/// A lightweight container for correlator visibility metadata used in Marlu operations.
///
/// This is intended to describe an accompanying visibility and weight ndarray.
///
/// This is similar to, and can be derived from [`mwalib::CorrelatorContext`],
/// but applies to a specific selection of baselines from a contiguous range of
/// timesteps and frequencies with a given configuration of averaging, and other
/// pre-processing settings.
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
    /// The geodetic position of each antenna.
    // TODO: maybe this doesn't belong here.
    pub tiles_xyz_geod: Vec<XyzGeodetic>,
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
        let tiles_xyz_geod = XyzGeodetic::get_tiles_mwa(&corr_ctx.metafits_context);
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
            tiles_xyz_geod,
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
