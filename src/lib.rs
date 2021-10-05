// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Core code to describe coordinate transformations, Jones matrices, etc.

pub mod constants;
pub mod jones;
pub mod math;
pub mod pos;
pub mod sexagesimal;
pub mod time;
pub mod io;

// Re-exports.
pub use jones::Jones;
pub use pos::{
    azel::AzEl,
    earth::{Ellipsoid, LatLngHeight},
    enh::ENH,
    hadec::HADec,
    lmn::LMN,
    pal, precession,
    radec::RADec,
    uvw::UVW,
    xyz::{XyzGeocentric, XyzGeodetic},
};

pub use num_complex::{Complex, Complex32 as c32, Complex64 as c64};

pub use erfa_sys;
pub use hifitime;

// If "mwalib" is enabled, re-export the crate here, as well its re-exported
// crates.
#[cfg(feature = "mwalib")]
pub use mwalib;
#[cfg(feature = "mwalib")]
pub use mwalib::{fitsio, fitsio_sys};
