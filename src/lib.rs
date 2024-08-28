
#![allow(dead_code, unused_variables, unused_imports)]
#![doc = include_str!("../README.md")]


pub use spc::*;
pub use xyz::*;
pub use obs::*;
pub use data::*;
pub use rgb::*;
pub use geometry::*;
use wasm_bindgen::JsValue;


pub mod obs;
pub mod spc;
pub mod xyz;
pub mod lab;
pub mod rgb;
pub mod gamma;
pub mod physics;
pub mod data;
pub mod cri;
pub mod geometry;


#[derive(thiserror::Error, Debug)]
pub enum CmError {
    #[error("please provide {0} arguments only")]
    ProvideOnlyArguments(String),
    #[error("please provide a single {0} argument only")]
    ProvideArgumentOnly(String),
    #[error("{name} should be within in range from {low} to {high}")]
    OutOfRange{name: String, low: f64,  high: f64},
    #[error("Currently only support for 3 primaries")]
    ThreePrimariesOnly,
    #[error("Out of Gamut")]
    OutOfGamut,
    #[error("Not yet implemented")]
    NotYetImplemented,
    #[error("Provide at least one {0} argument")]
    AtLeastOne(String),
    #[error("CCT Robertson Slope Error")]
    CCTRobertsonSlopeError,
    #[error("Error: {0}")]
    ErrorString(String),
    #[error("CCT: Distance to Planck locus >0.05")]
    CCTDuvHighError,
    #[error("CCT: Distance to Blackbody locus <-0.05")]
    CCTDuvLowError,
    #[error("CCT: Temperature too high")]
    CCTTemperatureTooHigh,
    #[error("CCT: Temperature too low")]
    CCTTemperatureTooLow  ,
    #[error("CCT: Maximum evaluation depth for CCT calculation in this implementation is {0}")]
    CCTEvaluationDept(u32),
    #[error("MacAdam: Point outside triangle, no interpolation possible")]
    MacAdamInterpolateError,
    #[error("RGB Display: Out of gamut")]
    RgbDisplayOutOfGamutError,
    #[error("CRI: No illuminant, use 'illuminant(light)' first")]
    CriNoSourceError,
    #[error("Filter: No match found")]
    FilterNoMatch,
    #[error("Primaries: Same White Point needed for RGB-transform between Observers")]
    PrimariesRgbTransformWhiteMismatch,
    #[error("Color Matching Function {0} not found")]
    ColorMatchingFunctionNotfound(String),
    #[error("Please provide the primaries in red, green, and blue order")]
    RgbTransformNotInRgbOrder,
    #[error("White point seems to be out of gamut")]
    NoWhiteBalance,
    #[error("Could not invert the RGB Matrix")]
    CouldNotInvertRGBMatrix,
    #[error("RgbTransforms: Please provide exactly three primaries only")]
    ProvideThreePrimariesOnly,
    #[error("RgbDisplay: Provide Observer Names only")]
    ProvideObserverNamesOnly,
    #[error("RgbDisplay: Provide Observer Names only")]
    ProvideOnlyTwoObservers,
    #[error("Filter: Please provide a valid Json Filter file")]
    ProvideValidJsonFilter,
    #[error("BoxcarIntegrator: please decrease domain resolution")]
    DomainStepError,
    #[error("Data size Error: need exactly 401 data values")]
    DataSize401Error,
    #[error("Linear Interpolate: Incorrect wavelength data")]
    InterpolateWavelengthError,
    #[error("This method requires distinct points")]
    RequiresDistinctPoints,
    #[error("Arguments require the identical Standard Observer")]
    RequireSameObserver,
    #[error("Lines do not intersect")]
    NoIntersection,
    #[error("Wavelength out of range")]
    WavelengthOutOfRange,
    #[error("Allowed wavelength range for this function is {0} to {1} nanometer")]
    NoUniqueSpectralLocus(usize,usize),
    #[error("Invalid Chromaticity Values")]
    InvalidChromaticityValues,
}

impl From<&str> for CmError {
    fn from(s: &str) -> Self {
        CmError::ErrorString(s.to_string())
    }
}

impl From<JsValue> for CmError {
    fn from(s: JsValue) -> Self {
        CmError::ErrorString(s.as_string().expect("Sorry, Unknown Error Encountered"))
    }
}

impl From<CmError> for JsValue {
    fn from(value: CmError) -> Self {
        value.to_string().into()
    }
}

 