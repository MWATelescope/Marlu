use thiserror::Error;

// use rubbl_casatables::CasacoreError;

#[derive(Error, Debug)]
pub enum MeasurementSetWriteError {
    #[error("bad channel info shape supplied to argument {argument} of function {function}. expected {expected}, received {received}")]
    BadArrayShape {
        argument: String,
        function: String,
        expected: String,
        received: String
    },

    // TODO: https://github.com/pkgw/rubbl/pull/148
    // #[error("{0}")]
    // RubblError(#[from] CasacoreError)
        
}
