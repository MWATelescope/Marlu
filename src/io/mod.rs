mod error;
mod ms;

use std::ops::Range;

use crate::{Jones, mwalib::CorrelatorContext};
use ndarray::{ArrayView3, ArrayViewMut3};

use self::error::IOError;

/// The container has visibilities which can be read by passing in a mwalib
/// context and the range of values to read.
pub trait VisReadable: Sync + Send {
    /// Read the visibilities and weights for the selected timesteps, coarse
    /// channels and baselines into the provided arrays.
    ///
    /// # Errors
    ///
    /// Can throw IOError if there is an issue reading.
    ///
    /// TODO: reduce number of arguments.
    #[allow(clippy::too_many_arguments)]
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

/// The container can accept visibilities by passing in the range of mwalib
/// indices corresponding to the visibilities being written.
pub trait VisWritable: Sync + Send {
    /// Write visibilities and weights from the arrays. Timestep, coarse channel
    /// and baseline indices are needed for labelling the visibility array
    ///
    /// `jones_array` - a three dimensional array of jones matrix visibilities. 
    ///     The dimensions of the array are `[timestep][channel][baseline]`
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
    /// # Errors
    ///
    /// Can throw IOError if there is an issue writing to the file, or the indices
    /// into `context` are invalid.
    #[allow(clippy::too_many_arguments)]
    fn write_vis_mwalib(
        &mut self,
        jones_array: ArrayView3<Jones<f32>>,
        weight_array: ArrayView3<f32>,
        context: &CorrelatorContext,
        timestep_range: &Range<usize>,
        coarse_chan_range: &Range<usize>,
        baseline_idxs: &[usize],
    ) -> Result<(), IOError>;
}
