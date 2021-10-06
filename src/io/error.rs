use thiserror::Error;

// use rubbl_casatables::CasacoreError;

#[derive(Error, Debug)]
pub enum MeasurementSetWriteError {
    #[error("bad channel info shape supplied. expected {expected}, received {received}")]
    BadSpwChannelInfoShape {
        expected: String,
        received: String
    },

    // TODO: https://github.com/pkgw/rubbl/pull/148
    // #[error("{0}")]
    // RubblError(#[from] CasacoreError)
        
}
