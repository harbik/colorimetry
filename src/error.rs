#[derive(thiserror::Error, Debug, PartialEq)]
#[non_exhaustive]
pub enum Error {
    #[error("{name} should be within in range from {low} to {high}")]
    OutOfRange { name: String, low: f64, high: f64 },
    #[error("Error: {0}")]
    ErrorString(String),
    #[error("CCT: Distance to Planck locus >0.05")]
    CCTDuvHighError,
    #[error("CCT: Distance to Blackbody locus <-0.05")]
    CCTDuvLowError,
    #[error("CCT: Temperature too high")]
    CCTTemperatureTooHigh,
    #[error("CCT: Temperature too low")]
    CCTTemperatureTooLow,
    #[error("Data size Error: need exactly 401 data values")]
    DataSize401Error,
    #[error("Linear Interpolate: Incorrect wavelength data")]
    InterpolateWavelengthError,
    #[error("This method requires distinct points")]
    RequiresDistinctPoints,
    #[error("Arguments require the identical Standard Observer")]
    RequireSameObserver,
    #[error("No Reference White values allowed")]
    NoReferenceWhiteAllowed,
    #[error("Lines do not intersect")]
    NoIntersection,
    #[error("Wavelength out of range")]
    WavelengthOutOfRange,
    #[error("Invalid Chromaticity Values")]
    InvalidChromaticityValues,
    #[error("This Method Requires CIE 1931-based XYZ values")]
    RequiresCIE1931XYZ,
    #[error("Colorant is required here")]
    NoColorant,
    #[error("RequiresSameIlluminant")]
    RequiresSameIlluminant,
    #[error("Spectrum {0} not found in Collection")]
    SpectrumNotFound(String),
    #[error("Provide at least {0} values")]
    ProvideAtLeastNValues(usize),
    #[error("Invalid RGB value")]
    InvalidRgbValue,
}

impl From<&str> for Error {
    fn from(s: &str) -> Self {
        Error::ErrorString(s.to_string())
    }
}

#[cfg(target_arch = "wasm32")]
impl From<wasm_bindgen::JsValue> for Error {
    fn from(s: wasm_bindgen::JsValue) -> Self {
        Error::ErrorString(s.as_string().expect("Sorry, Unknown Error Encountered"))
    }
}

#[cfg(target_arch = "wasm32")]
impl From<Error> for wasm_bindgen::JsValue {
    fn from(value: Error) -> Self {
        value.to_string().into()
    }
}
