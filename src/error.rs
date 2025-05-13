use wasm_bindgen::JsValue;

#[derive(thiserror::Error, Debug, PartialEq)]
pub enum CmtError {
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
    #[error("Allowed wavelength range for this function is {0} to {1} nanometer")]
    NoUniqueSpectralLocus(usize, usize),
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

impl From<&str> for CmtError {
    fn from(s: &str) -> Self {
        CmtError::ErrorString(s.to_string())
    }
}

impl From<JsValue> for CmtError {
    fn from(s: JsValue) -> Self {
        CmtError::ErrorString(s.as_string().expect("Sorry, Unknown Error Encountered"))
    }
}

impl From<CmtError> for JsValue {
    fn from(value: CmtError) -> Self {
        value.to_string().into()
    }
}
