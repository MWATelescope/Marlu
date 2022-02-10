use hifitime::{Duration, Epoch};

use crate::XyzGeodetic;

use std::ops::Range;

#[cfg(feature = "mwalib")]
use mwalib::CorrelatorContext;

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
pub struct MarluVisContext {
    /// The number of selected timesteps (Axis 0) in the accompanying visibility and weight ndarrays.
    pub num_sel_timesteps: usize,
    /// The timestamp at the start of the first selected pre-averaging timestep
    pub start_timestamp: Epoch,
    /// Duration between each pre-averaging timestep [milliseconds]
    pub int_time: Duration,
    /// The number of selected channels (Axis 1) in the accompanying visibility and weight ndarrays.
    pub num_sel_chans: usize,
    /// The centre frequency of the first selected pre-averaging channel [Hz]
    // pub start_freq_hz: f64,
    /// The bandwidth between each pre-averaging channel [Hz]
    // pub freq_resolution_hz: f64,
    /// The tile index pairs for each selected baseline
    pub sel_baselines: Vec<(usize, usize)>,
    /// The geodetic position of each antenna.
    pub tiles_xyz_geod: Vec<XyzGeodetic>,
    /// Time averaging factor
    pub avg_time: usize,
    /// Frequency averaging factor
    pub avg_freq: usize,
    /// Number of polarisation combinations in the visibilities e.g. XX,XY,YX,YY == 4
    pub num_vis_pols: usize,
}

impl MarluVisContext {
    #[cfg(feature = "mwalib")]
    pub fn from_mwalib(
        context: &CorrelatorContext,
        timestep_range: &Range<usize>,
        coarse_chan_range: &Range<usize>,
        baseline_idxs: &[usize],
        avg_time: usize,
        avg_freq: usize,
    ) -> Self {
        let num_sel_timesteps = timestep_range.len();
        let num_sel_coarse_chans = coarse_chan_range.len();
        let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
        let num_sel_chans = fine_chans_per_coarse * num_sel_coarse_chans;

        let start_timestamp = Epoch::from_gpst_seconds(
            context.timesteps[timestep_range.start].gps_time_ms as f64 / 1e3,
        );

        let int_time =
            Duration::from_milliseconds(context.metafits_context.corr_int_time_ms.into());

        // baselines
        let tiles_xyz_geod = XyzGeodetic::get_tiles_mwa(&context.metafits_context);
        let sel_baselines = baseline_idxs
            .iter()
            .map(|&idx| {
                let baseline = &context.metafits_context.baselines[idx];
                (baseline.ant1_index, baseline.ant2_index)
            })
            .collect::<Vec<_>>();

        let num_vis_pols = context.metafits_context.num_visibility_pols;

        MarluVisContext {
            num_sel_timesteps,
            start_timestamp,
            int_time,
            num_sel_chans,
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

    /// Whether this context corresponds with the the case where the time and
    /// frequency averaging factors are both 1, and no averaging is performed.
    pub fn trivial_averaging(&self) -> bool {
        self.avg_time == 1 && self.avg_freq == 1
    }

    /// The number of timesteps in the post-averaging time dimension
    pub fn num_avg_timesteps(&self) -> usize {
        (self.num_sel_timesteps as f64 / self.avg_time as f64).ceil() as usize
    }

    /// The number of channels in the post-averaging frequency dimension
    pub fn num_avg_chans(&self) -> usize {
        (self.num_sel_chans as f64 / self.avg_freq as f64).ceil() as usize
    }

    /// The expected dimensions of the visibility and weight ndarray selection.
    pub fn avg_dims(&self) -> (usize, usize, usize) {
        (
            self.num_avg_timesteps(),
            self.num_avg_chans(),
            self.sel_baselines.len(),
        )
    }

    /// The integration time of the post-averaging data.
    pub fn avg_int_time(&self) -> Duration {
        self.int_time * self.avg_time as u64
    }

    /// A vector of centroid timesteps for each post-averaging timestep.
    pub fn avg_centroid_timestamps(&self) -> Vec<Epoch> {
        let avg_int_time = self.avg_int_time();
        (0..self.num_avg_timesteps())
            .step_by(self.avg_time)
            .map(|timestep_idx| self.start_timestamp + avg_int_time * (timestep_idx as f64 + 0.5))
            .collect()
    }
}
