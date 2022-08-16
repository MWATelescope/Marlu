// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

use rubbl_casatables::CasacoreError;
use thiserror::Error;

#[derive(Error, Debug)]
#[error("bad array shape supplied to argument {argument} of function {function}. expected {expected}, received {received}")]
pub struct BadArrayShape {
    pub argument: &'static str,
    pub function: &'static str,
    pub expected: String,
    pub received: String,
}

// TODO: there are plenty of panics in ms that need enums
#[derive(Error, Debug)]
pub enum MeasurementSetWriteError {
    /// An error when trying to write to an unexpected row.
    #[error("Tried to write {rows_attempted} rows, but only {rows_remaining} rows are remaining out of {rows_total}")]
    MeasurementSetFull {
        /// The row number (0-indexed)
        rows_attempted: usize,
        /// Rows remaining.
        rows_remaining: usize,
        /// Total capacity of measurement set
        rows_total: usize,
    },

    // TODO: https://github.com/pkgw/rubbl/pull/148
    /// From Rubbl Casacore
    #[error("Rubbl CASACore error {inner:?}")]
    CasacoreError { inner: CasacoreError },

    /// From Rubbl
    #[error("Rubbl error {inner:?}")]
    RubblError { inner: failure::Error },

    /// Tried to create a directory where a file already exists
    #[error("cannot create directory, path={path} already exists and is not a directory")]
    NotADirectory { path: String },

    #[error(transparent)]
    BadArrayShape(#[from] BadArrayShape),

    /// An IO error.
    #[error(transparent)]
    StdIo(#[from] std::io::Error),

    #[error(transparent)]
    SystemTimeError(#[from] std::time::SystemTimeError),
}

impl From<failure::Error> for MeasurementSetWriteError {
    fn from(inner: failure::Error) -> Self {
        Self::RubblError { inner }
    }
}

impl From<CasacoreError> for MeasurementSetWriteError {
    fn from(inner: CasacoreError) -> Self {
        Self::CasacoreError { inner }
    }
}

#[derive(Error, Debug)]
pub enum UvfitsWriteError {
    /// An error when trying to write to an unexpected row.
    #[error("Tried to write to row number {row_num}, but only {num_rows} rows are expected")]
    BadRowNum {
        /// The row number (0-indexed)
        row_num: usize,
        /// Total number of rows expected.
        num_rows: usize,
    },

    /// An error when less rows were written to an HDU than expected.
    #[error("Expected {total} uvfits rows to be written, but only {current} were written")]
    NotEnoughRowsWritten {
        /// Number of rows written
        current: usize,
        /// Total number of rows expected.
        total: usize,
    },

    /// An error associated with ERFA.
    #[error(transparent)]
    Erfa(#[from] crate::pos::ErfaError),

    /// An error associated with fitsio.
    #[error(transparent)]
    Fitsio(#[from] mwalib::fitsio::errors::Error),

    /// An error when converting a Rust string to a C string.
    #[error(transparent)]
    BadString(#[from] std::ffi::NulError),

    /// An IO error.
    #[error(transparent)]
    StdIo(#[from] std::io::Error),
}

impl From<crate::io::uvfits::FitsioOrCStringError> for UvfitsWriteError {
    fn from(e: crate::io::uvfits::FitsioOrCStringError) -> Self {
        match e {
            super::uvfits::FitsioOrCStringError::Fitsio(e) => Self::Fitsio(e),
            super::uvfits::FitsioOrCStringError::Nul(e) => Self::BadString(e),
        }
    }
}

#[derive(Error, Debug)]
#[allow(clippy::upper_case_acronyms)]
/// All the errors that can occur in file io operations
pub enum IOError {
    #[error(transparent)]
    /// Error derived from [`io::errors::MeasurementSetWriteError`]
    MeasurementSetWriteError(#[from] MeasurementSetWriteError),

    #[error(transparent)]
    /// Error derived from [`marlu::mwalib::FitsError`]
    FitsError(#[from] mwalib::FitsError),

    #[error(transparent)]
    /// Error derived from [`fitsio::errors::Error`]
    FitsioError(#[from] mwalib::fitsio::errors::Error),

    #[error(transparent)]
    /// Error derived from [`io::errors::UvfitsWriteError`]
    UvfitsWriteError(#[from] UvfitsWriteError),

    #[error(transparent)]
    BadArrayShape(#[from] BadArrayShape),

    /// From Rubbl
    #[error("Rubbl error {inner:?}")]
    RubblError { inner: failure::Error },
}

impl From<failure::Error> for IOError {
    fn from(inner: failure::Error) -> Self {
        Self::RubblError { inner }
    }
}
