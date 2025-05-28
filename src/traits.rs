use std::borrow::Cow;

use crate::{
    illuminant::Illuminant, observer::Observer, spectrum::Spectrum, xyz::XYZ,
};

/**
Spectral representation of Lights, typically in form of (standard) Illuminants.

Also allows to use lookup tristimulus values, such as the very common
[`D65`](crate::illuminant::CieIlluminant::D65) illuminant
(see [`CieIlluminant`](crate::illuminant::CieIlluminant) implemention).
Calculating them from a spectrum is the default implementation.
*/
pub trait Light {
    /// Calculates the tristimulus values of the light source, using the
    /// provided observer's color matching data, and an optional illuminance value.
    ///
    /// The illuminance value is optional, and if not provided, the actual luminous values in the
    /// spectrum are used.
    fn xyzn(&self, observer: Observer, y: Option<f64>) -> XYZ {
        let xyz = observer.data().xyz_from_spectrum(&self.spectrum());
        if let Some(illuminance) = y {
            xyz.set_illuminance(illuminance)
        } else {
            xyz
        }
    }

    fn spectrum(&self) -> Cow<Spectrum>;
}

impl From<&dyn Light> for Illuminant {
    fn from(light: &dyn Light) -> Self {
        Illuminant(light.spectrum().into_owned())
    }
}

pub trait Filter {
    fn spectrum(&self) -> Cow<Spectrum>;
}
