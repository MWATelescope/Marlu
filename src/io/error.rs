use thiserror::Error;

// TODO: there are plenty of panics in ms that need enums
#[derive(Error, Debug)]
pub enum MeasurementSetWriteError {
    // TODO: https://github.com/pkgw/rubbl/pull/148
    // #[error("{0}")]
    // RubblError(#[from] CasacoreError)
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
    IO(#[from] std::io::Error),
}

#[derive(Error, Debug)]
#[allow(clippy::upper_case_acronyms)]
/// All the errors that can occur in file io operations
pub enum IOError {
    #[error("bad array shape supplied to argument {argument} of function {function}. expected {expected}, received {received}")]
    BadArrayShape {
        argument: String,
        function: String,
        expected: String,
        received: String,
    },

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
}
