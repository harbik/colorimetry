#![allow(unused)]
use std::{collections::HashMap, hash::Hash};

use colorimetry::{
    illuminant::D65,
    observer::{self, Observer},
    prelude::{CieIlluminant, Spectrum, Stimulus},
    rgb::{gaussian_filtered_primaries, GammaCurve},
    spectrum::NS,
    traits::Light,
    xyz::XYZ,
};
use serde::Serialize;
use strum::IntoEnumIterator;

#[derive(Serialize, Debug, Clone)]
pub struct GenRgbSpaceData {
    primaries: [Vec<f64>; 3],
    white: CieIlluminant,
    gamma: Vec<f64>,
    name: String,
    identifier: String,
    chromaticities: HashMap<Observer, [[f64; 2]; 3]>,
    white_points: HashMap<Observer, XYZ>, // XYZ coordinates for each observer
}

impl GenRgbSpaceData {
    /// The sRGB color space, created by HP and Microsoft in 1996.  It is the
    /// default color space used in an image, or a Web-page, if no color space
    /// is specified in a tag or in a color profile.
    /// In this instance, the primaries are composed of Daylight D65 filtered Gaussian filters.
    /// For the green and blue these are single Gaussians.
    /// For red primary a combination of a blue and red Gaussian is used, with
    /// the blue being equal to the blue primary.
    pub fn srgb() -> Self {
        /*
        Red Gaussian center wavelength and width, and blue lumen fraction
        added to match the red primary sRGB color space specifiation.
        See the `examples/primaries/main.rs` how these were calculated.
        */
        const RED: [f64; 3] = [627.2101041540204, 23.38636113607498, 0.006839349789397155];
        const GREEN: [f64; 2] = [541.2355906535001, 33.66683554608116];
        const BLUE: [f64; 2] = [398.0273721579992, 52.55338039394701];

        let primaries = gaussian_filtered_primaries(D65.as_ref(), RED, GREEN, BLUE);
        let white = CieIlluminant::D65;
        let gamma = vec![2.4, 1.0 / 1.055, 0.055 / 1.055, 1.0 / 12.92, 0.04045];
        Self {
            name: "sRGB".to_string(),
            identifier: "SRGB".to_string(),
            chromaticities: Self::chromaticities(&primaries),
            primaries: Self::stimuli_to_vecs(&primaries),
            white,
            white_points: GenRgbSpaceData::white_points(white),
            gamma,
        }
    }

    /// The Adobe RGB color space.
    ///
    /// The primaries are composed of Daylight D65 filtered Gaussian filters.
    /// For the green and blue these are single Gaussians.
    /// For red primary a combination of a blue and red Gaussian is used, with
    /// the blue being equal to the blue primary.
    pub fn adobe_rgb() -> Self {
        // Red Gaussian center wavelength and width, and blue lumen fraction
        // added to match the red primary sRGB color space specifiation.
        // See the `examples/primaries/main.rs` how these were calculated.
        const RED: [f64; 3] = [627.2101041540204, 23.38636113607498, 0.006839349789397155];
        const GREEN: [f64; 2] = [531.1505632157933, 20.61013919689458];
        const BLUE: [f64; 2] = [398.0273721579992, 52.55338039394701];

        let primaries = gaussian_filtered_primaries(D65.as_ref(), RED, GREEN, BLUE);
        let white = CieIlluminant::D65;
        let gamma = vec![563.0 / 256.0];
        // See https://en.wikipedia.org/wiki/Adobe_RGB_color_space#ICC_PCS_color_image_encoding
        Self {
            name: "Adobe RGB".to_string(),
            identifier: "ADOBE_RGB".to_string(),
            chromaticities: Self::chromaticities(&primaries),
            primaries: Self::stimuli_to_vecs(&primaries),
            white,
            white_points: GenRgbSpaceData::white_points(white),
            gamma,
        }
    }

    /// The Display P3 color space.
    ///
    /// Used in Apple devices, such as iPhones and iPads.
    /// The Display P3 color space is a wide-gamut RGB color space that is
    /// designed to be compatible with the DCI-P3 color space used in digital cinema.
    /// It is based on the same primaries as DCI-P3, but with a different white point,
    /// which is D65 instead of DCI-P3's DCI white point.
    ///
    /// The primaries are composed of Daylight D65 filtered Gaussian filters.
    /// For the green and blue these are single Gaussians.
    pub fn display_p3() -> Self {
        // Red Gaussian center wavelength and width, and blue lumen fraction
        // added to match the red primary sRGB color space specifiation.
        // See the `examples/primaries/main.rs` how these were calculated.
        const RED: [f64; 3] = [637.9658554073235, 23.274151215539906, 1.4261447449730065e-6];

        // Blue and Green Gaussian filters.
        const GREEN: [f64; 2] = [539.8416064376327, 21.411199475289777];
        const BLUE: [f64; 2] = [398.0273721579992, 52.55338039394701];

        let primaries = gaussian_filtered_primaries(D65.as_ref(), RED, GREEN, BLUE);
        let white = CieIlluminant::D65;
        let gamma = vec![2.4, 1.0 / 1.055, 0.055 / 1.055, 1.0 / 12.92, 0.04045];
        Self {
            name: "Display P3".to_string(),
            identifier: "DISPLAY_P3".to_string(),
            chromaticities: Self::chromaticities(&primaries),
            primaries: Self::stimuli_to_vecs(&primaries),
            white,
            white_points: GenRgbSpaceData::white_points(white),
            gamma,
        }
    }

    /// The CIE RGB color space.
    /// This color space is defined by the CIE and is based on the
    /// CIE 1931 standard observer. It is not commonly used in practice, but
    /// it is useful for color management and color science applications.
    ///
    /// The primaries are monochromatic stimuli at wavelength 700 nm (red),
    /// 546.1 nm (green), and 435.8 nm (blue), with the white point being
    /// CIE Illuminant E, which is a uniform white light source.
    ///
    /// The gamma curve is not used in CIE RGB, so it is an empty vector.
    pub fn cie_rgb() -> Self {
        let primaries = [
            Stimulus::new(Spectrum::from_wavelength_map(&[(700, 1.0)])), // 700.0 nm
            Stimulus::new(Spectrum::from_wavelength_map(&[(546, 0.9), (547, 0.1)])), // 546.1 nm
            Stimulus::new(Spectrum::from_wavelength_map(&[(435, 0.2), (436, 0.8)])), // 435.8 nm
        ];
        let white = CieIlluminant::E;
        let white_points = GenRgbSpaceData::white_points(white);
        let gamma = vec![];
        // CIE RGB does not use a gamma curve, so we use an empty vector.
        // Use with float values only!
        Self {
            name: "CIE RGB".to_string(),
            identifier: "CIE_RGB".to_string(),
            chromaticities: Self::chromaticities(&primaries),
            primaries: Self::stimuli_to_vecs(&primaries),
            white,
            white_points,
            gamma,
        }
    }

    pub fn dci_p3() -> Self {
        todo!()
    }

    pub fn rec709() -> Self {
        todo!()
    }

    pub fn rec2020() -> Self {
        todo!()
    }

    pub fn prophoto_rgb() -> Self {
        todo!()
    }

    pub fn aces() -> Self {
        todo!()
    }

    pub fn acescg() -> Self {
        todo!()
    }

    pub fn white_points(white: CieIlluminant) -> HashMap<Observer, XYZ> {
        let mut map = HashMap::new();
        for observer in Observer::iter() {
            let white_point = white.white_point(observer);
            map.insert(observer, white_point);
        }
        map
    }

    pub fn chromaticities(stimuli: &[Stimulus; 3]) -> HashMap<Observer, [[f64; 2]; 3]> {
        let mut chromaticities = HashMap::new();
        for observer in Observer::iter() {
            let mut values = [[0.0; 2]; 3];
            for (i, stimulus) in stimuli.iter().enumerate() {
                let xy = observer
                    .xyz_from_spectrum(stimulus.as_ref())
                    .chromaticity()
                    .to_array();
                values[i][0] = xy[0];
                values[i][1] = xy[1];
            }
            chromaticities.insert(observer, values);
        }
        chromaticities
    }

    pub fn stimuli_to_vecs(stimuli: &[Stimulus; 3]) -> [Vec<f64>; 3] {
        [
            stimuli[0].as_ref().values().to_vec(),
            stimuli[1].as_ref().values().to_vec(),
            stimuli[2].as_ref().values().to_vec(),
        ]
    }
    pub fn write(&self) {
        // This function would write the RGB space data to a file or output.
        // Implementation is omitted for brevity.
        println!("Writing RGB space data for: {}", self.name);
        // Example output:
        println!("Identifier: {}", self.identifier);
        println!("Primaries: {:?}", self.primaries);
        println!("White: {:?}", self.white);
        println!("White Points: {:?}", self.white_points);
        println!("Gamma curve: {:?}", self.gamma);
        println!("Chromaticities: {:?}", self.chromaticities);
    }
}

pub fn main() {
    println!("Generating RGB space data...");

    let rgbspace = GenRgbSpaceData::srgb();
    println!("Writing RGB space data for: {:?}", rgbspace.name);
    rgbspace.write();
    // for each observer
    for observer in Observer::iter() {
        println!("{}", observer);
    }
}
