use crate::{context::VisContext, XyzGeodetic};

pub mod error;
pub mod ms;
pub mod uvfits;

// re-exports
pub use ms::MeasurementSetWriter;
pub use uvfits::UvfitsWriter;

cfg_if::cfg_if! {
    if #[cfg(feature = "mwalib")] {
        use std::ops::Range;

        use crate::{mwalib::CorrelatorContext, Jones};
        use ndarray::{ArrayView3, ArrayView4, ArrayViewMut3, Axis};
        use self::error::IOError;
        use approx::abs_diff_eq;
        use log::trace;
    }
}

/// The container has visibilities which can be read by passing in a mwalib
/// context and the range of values to read.
pub trait VisReadable: Sync + Send {
    /// Read the visibilities and weights for the selected timesteps, coarse
    /// channels and baselines into the provided arrays.
    ///
    /// # Errors
    ///
    /// Can throw `IOError` if there is an issue reading.
    ///
    /// TODO: reduce number of arguments.
    #[allow(clippy::too_many_arguments)]
    #[cfg(feature = "mwalib")]
    fn read_vis_mwalib(
        &self,
        jones_array: ArrayViewMut3<Jones<f32>>,
        weight_array: ArrayViewMut3<f32>,
        context: &CorrelatorContext,
        timestep_range: &Range<usize>,
        coarse_chan_range: &Range<usize>,
        baseline_idxs: &[usize],
    ) -> Result<(), IOError>;
}

/// The container can accept a chunk of visibilities to be written.
pub trait VisWritable: Sync + Send {
    /// Write a chunk of visibilities, contextualised with a [`VisContext`].
    ///
    /// This makes use of the heuristic that the weights of all pols of a visibility
    /// should be equal, and thus, only three dimensions of weights / flags are used.
    ///
    /// `vis` - a three dimensional array of jones matrix visibilities.
    ///     The dimensions of the array are `[timestep][channel][baseline]`
    ///
    /// `weights` - a three dimensional array of visibility weights.
    ///     The dimensions of the array are `[timestep][channel][baseline]`
    ///
    /// `flags` - a three dimensional array of visibility flags.
    ///     The dimensions of the array are `[timestep][channel][baseline]`
    ///
    /// `vis_ctx` - a [`VisContext`] which contextualises each axis of the visibilities.
    ///
    /// `draw_progress` - whether or not to draw a progress bar.
    ///
    /// Hyperdrive:
    /// - shim for baseline
    /// - shim for timestep indices
    fn write_vis_marlu(
        &mut self,
        vis: ArrayView3<Jones<f32>>,
        weights: ArrayView3<f32>,
        flags: ArrayView3<bool>,
        vis_ctx: &VisContext,
        tiles_xyz_geod: &[XyzGeodetic],
        draw_progress: bool,
    ) -> Result<(), IOError>;

    /// Contract the weight and flag arrays along the pol axis, and write using
    /// [`VisWritable::write_vis_marlu`]
    ///
    /// `jones_array` - a three dimensional array of jones matrix visibilities.
    ///     The dimensions of the array are `[timestep][channel][baseline]`
    ///
    /// `weight_array` - a four dimensional array of visibility weights.
    ///     The dimensions of the array are `[timestep][channel][baseline][pol]`
    ///
    /// `flag_array` - a four dimensional array of visibility flags.
    ///     The dimensions of the array are `[timestep][channel][baseline][pol]`
    ///
    /// `context` - a [`mwalib::CorrelatorContext`] object containing information
    ///     about the timesteps, channels and baselines referred to in the jones
    ///     array.
    ///
    /// `timestep_range` - the range of indices into `CorrelatorContext.timesteps`
    ///     corresponding to the first dimension of the jones array.
    ///
    /// `coarse_chan_range` - the range of indices into `CorrelatorContext.coarse_chans`
    ///     corresponding to the second dimension of the jones array.
    ///     Note: this is a range of coarse channels, where the jones array is a range
    ///     of fine channels
    ///
    /// `baseline_idxs` - the range of indices into `CorrelatorContext.metafits_context.baselines`
    ///     corresponding to the third dimension of the jones array.
    ///
    /// `avg_time` - the number of timesteps to average together.
    ///
    /// `avg_freq` - the number of channels to average together.
    ///
    /// `draw_progress` - whether or not to draw a progress bar.
    ///
    /// # Errors
    ///
    /// Can throw [`IOError`] if there is an issue writing to the file, or the indices
    /// into `context` are invalid.
    ///
    // TODO: deprecate this entirely
    #[allow(clippy::too_many_arguments)]
    #[cfg(feature = "mwalib")]
    fn write_vis_mwalib(
        &mut self,
        jones_array: ArrayView3<Jones<f32>>,
        weight_array: ArrayView4<f32>,
        flag_array: ArrayView4<bool>,
        corr_ctx: &CorrelatorContext,
        timestep_range: &Range<usize>,
        coarse_chan_range: &Range<usize>,
        baseline_idxs: &[usize],
        avg_time: usize,
        avg_freq: usize,
        draw_progress: bool,
    ) -> Result<(), IOError> {
        trace!("write_vis_mwalib");

        let vis_ctx = VisContext::from_mwalib(
            corr_ctx,
            timestep_range,
            coarse_chan_range,
            baseline_idxs,
            avg_time,
            avg_freq,
        );

        let weight_array = weight_array.map_axis(Axis(3), |weights| {
            assert!(weights.iter().all(|&w| abs_diff_eq!(weights[0], w)));
            weights[0]
        });

        let flag_array = flag_array.map_axis(Axis(3), |flags| {
            assert!(flags.iter().all(|&w| flags[0] == w));
            flags[0]
        });

        self.write_vis_marlu(
            jones_array,
            weight_array.view(),
            flag_array.view(),
            &vis_ctx,
            &XyzGeodetic::get_tiles_mwa(&corr_ctx.metafits_context),
            draw_progress,
        )
    }
}
