
#![allow(dead_code, unused_variables, unused_imports, )]
#![doc = include_str!("../README.md")]


#[cfg(feature="cct")]
pub use cct::*;
pub use colorant::*;
#[cfg(feature="cri")]
pub use cri::*;
pub use data::*;
pub use geometry::*;
pub use illuminant::*;
pub use observer::*;
pub use physics::*;
pub use rgb::*;
pub use rgbspace::*;
pub use spectrum::*;
pub use std_illuminants::*;
pub use stimulus::*;
pub use lab::*;
pub use xyz::*;
use wasm_bindgen::JsValue;


#[cfg(feature="cct")]
pub mod cct;
pub mod colorant;
#[cfg(feature="cri")]
pub mod cri;
pub mod data;
pub mod gamma;
pub mod geometry;
pub mod illuminant;
pub mod lab;
pub mod observer;
pub mod physics;
pub mod rgb;
pub mod rgbspace;
pub mod spectrum;
pub mod std_illuminants;
pub mod stimulus;
pub mod xyz;


#[derive(thiserror::Error, Debug, PartialEq)]
pub enum CmtError {
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
    #[error("Colorant: No match found")]
    ColorantNoMatch,
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
    #[error("No Reference White values allowed")]
    NoReferenceWhiteAllowed,
    #[error("Lines do not intersect")]
    NoIntersection,
    #[error("Wavelength out of range")]
    WavelengthOutOfRange,
    #[error("Allowed wavelength range for this function is {0} to {1} nanometer")]
    NoUniqueSpectralLocus(usize,usize),
    #[error("Invalid Chromaticity Values")]
    InvalidChromaticityValues,
    #[error("This Method Requires CIE 1931-based XYZ values")]
    RequiresCIE1931XYZ,
    #[error("Please provide a Reference Illuminant")]
    NoReferenceIlluminantProvided,
    #[error("RequiresSameIlluminant")]
    RequiresSameIlluminant,
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

 
 