// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Core code to describe coordinate transformations, Jones matrices, etc.

pub mod constants;
pub mod jones;
pub mod math;
pub mod pos;
pub mod precession;
pub mod sexagesimal;
pub mod time;

// Re-exports.
pub use jones::Jones;
pub use pos::{
    azel::AzEl,
    earth::{Ellipsoid, LatLngHeight},
    enh::ENH,
    hadec::HADec,
    lmn::LMN,
    pal,
    radec::RADec,
    uvw::UVW,
    xyz::{XyzGeocentric, XyzGeodetic},
};

pub use num_complex::{Complex32 as c32, Complex64 as c64};

pub use erfa_sys;
#[cfg(feature = "mwalib")]
pub use mwalib;
