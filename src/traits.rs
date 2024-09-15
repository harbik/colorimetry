
use crate::{Colorant, Spectrum};

/**
A trait that overrides the way a light’s tristimulus values are obtained.

Calculating them from a spectrum is the default implementation.  For standard illuminants,
especially D65, which are used so frequently, their values are obtained from buffered entries.
*/
pub trait Light {
    // Default implementation takes the illuminant spectrum, and calculates tristimulus values using the
    // provided observer's color matching data.
    fn xyzn(&self, observer: crate::Observer, y: Option<f64>) -> crate::XYZ {
        let xyz = observer.data().xyz_raw(&self.spectrum(), None);
        if let Some(illuminance) = y {
            xyz.set_illuminance(illuminance)
        } else {
            xyz
        }
    }
    fn spectrum(&self) -> &Spectrum;
}

pub trait Filter {

    fn spectrum(&self) -> &Spectrum;

}
