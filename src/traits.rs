
use std::borrow::Cow;

use crate::{Colorant, Illuminant, Spectrum};

/**
Spectral representation of Lights, typically in form of (standard) Illuminants.

Also allows to use lookup tristimulus values, such as the very common [D65] illuminant (see [crate::StdIlluminant] implemention).
Calculating them from a spectrum is the default implementation.  For standard illuminants,
especially D65, which are used so frequently, their values are obtained from buffered entries.
*/
pub trait Light {
    // Default implementation takes the illuminant spectrum, and calculates tristimulus values using the
    // provided observer's color matching data.
    fn xyzn(&self, observer: crate::Observer, y: Option<f64>) -> crate::XYZ {
        let xyz = observer.data().xyz_from_spectrum(&self.spectrum(), None);
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
