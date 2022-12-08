// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
#[cfg(feature = "ms")]
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

    /// An error from within rubbl's casacore
    #[error(transparent)]
    Casacore(#[from] rubbl_casatables::CasacoreError),

    /// An error from a CASA table
    #[error(transparent)]
    Table(#[from] rubbl_casatables::TableError),

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

#[derive(Error, Debug)]
#[cfg(feature = "cfitsio")]
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

    /// An error associated with fitsio.
    #[error(transparent)]
    Fitsio(#[from] fitsio::errors::Error),

    /// An error when converting a Rust string to a C string.
    #[error(transparent)]
    BadString(#[from] std::ffi::NulError),

    /// An IO error.
    #[error(transparent)]
    StdIo(#[from] std::io::Error),
}

#[cfg(feature = "cfitsio")]
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
    #[cfg(feature = "ms")]
    /// Error derived from [`io::errors::MeasurementSetWriteError`]
    MeasurementSetWriteError(#[from] MeasurementSetWriteError),

    #[cfg(feature = "mwalib")]
    #[error(transparent)]
    /// Error derived from [`mwalib::FitsError`]
    FitsError(#[from] mwalib::FitsError),

    #[error(transparent)]
    #[cfg(feature = "cfitsio")]
    /// Error derived from [`fitsio::errors::Error`]
    FitsioError(#[from] fitsio::errors::Error),

    #[error(transparent)]
    #[cfg(feature = "cfitsio")]
    /// Error derived from [`io::errors::UvfitsWriteError`]
    UvfitsWriteError(#[from] UvfitsWriteError),

    #[error(transparent)]
    BadArrayShape(#[from] BadArrayShape),

    #[cfg(feature = "ms")]
    #[error(transparent)]
    Table(#[from] rubbl_casatables::TableError),
}
