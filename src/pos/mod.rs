// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Super module for all positional code.

pub mod azel;
pub mod earth;
pub mod enh;
pub mod hadec;
pub mod lmn;
pub mod pal;
pub mod radec;
pub mod uvw;
pub mod xyz;

use thiserror::Error;

#[derive(Error, Debug)]
#[error(
    "{source_file}:{source_line} Call to ERFA function {function} returned status code {status}"
)]
pub struct ErfaError {
    source_file: &'static str,
    source_line: u32,
    status: i32,
    function: &'static str,
}
