// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

pub mod error;
use ndarray::prelude::*;

use crate::{context::VisContext, Jones};
use error::IOError;

cfg_if::cfg_if! {
    if #[cfg(feature = "cfitsio")] {
        pub mod uvfits;

        pub use error::UvfitsWriteError;
        pub use uvfits::UvfitsWriter;
    }
}

cfg_if::cfg_if! {
    if #[cfg(feature = "ms")] {
        pub mod ms;

        pub use error::MeasurementSetWriteError;
        pub use ms::MeasurementSetWriter;
    }
}

/// The container can accept a chunk of visibilities to be written.
pub trait VisWrite {
    /// Write a chunk of visibilities, contextualised with a [`VisContext`].
    ///
    /// This makes use of the heuristic that the weights of all pols of a visibility
    /// should be equal, and thus, only three dimensions of weights / flags are used.
    ///
    /// `vis` - a three dimensional array of jones matrix visibilities.
    ///     The dimensions of the array are `[timestep][channel][baseline]`
    ///
    /// `weights` - a three dimensional array of visibility weights, where the sign
    ///     of the double is the flag. Assumes all weights and flags are the same
    ///     for a given pol. The dimensions of the array are `[timestep][channel][baseline]`
    ///
    /// `vis_ctx` - a [`VisContext`] which contextualises each axis of the visibilities.
    ///
    fn write_vis(
        &mut self,
        vis: ArrayView3<Jones<f32>>,
        weights: ArrayView3<f32>,
        vis_ctx: &VisContext,
    ) -> Result<(), IOError>;

    /// When all visibilities have been given to this [`VisWrite`] implementor,
    /// calling this function will perform any remaining tasks before the writer
    /// can be dropped.
    fn finalise(&mut self) -> Result<(), IOError>;
}
