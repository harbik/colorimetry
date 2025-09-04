#![allow(unused)]
use std::{collections::HashMap, hash::Hash};

use chrono::Utc;
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

/// Defines a generic RGB color space data structure
/// that is used to create a RgbSpaceData instance in the
/// colorimetry library
#[derive(Serialize, Debug, Clone)]
pub struct TemplateContext {
    date: String,
    name: String,
    identifier: String,
    red_prim: Stimulus,
    green_prim: Stimulus,
    blue_prim: Stimulus,
    white: CieIlluminant,
    gamma: String,
    chromaticities: Vec<[[f64; 2]; 3]>,
    target_chromaticities_cie1931: [[f64; 2]; 3],
    white_points: Vec<[f64; 3]>,
}

impl TemplateContext {
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

        let [red_prim, green_prim, blue_prim] =
            gaussian_filtered_primaries(D65.as_ref(), RED, GREEN, BLUE);
        let white = CieIlluminant::D65;
        let gamma = GammaCurve::new(
            5,
            [
                2.4,
                1.0 / 1.055,
                0.055 / 1.055,
                1.0 / 12.92,
                0.04045,
                0.0,
                0.0,
            ],
        );
        Self {
            date: Utc::now().format("%Y-%m-%d").to_string(),
            name: "sRGB".to_string(),
            identifier: "SRGB".to_string(),
            chromaticities: Self::chromaticities([&red_prim, &green_prim, &blue_prim]),
            red_prim,
            green_prim,
            blue_prim,
            white,
            white_points: TemplateContext::white_points(white),
            gamma: format!("{gamma:?}"),
            target_chromaticities_cie1931: [
                [0.6400, 0.3300], // Red
                [0.3000, 0.6000], // Green
                [0.1500, 0.0600], // Blue
            ],
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

        let [red_prim, green_prim, blue_prim] =
            gaussian_filtered_primaries(D65.as_ref(), RED, GREEN, BLUE);
        let white = CieIlluminant::D65;
        let gamma = GammaCurve::new(1, [563.0 / 256.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        // See https://en.wikipedia.org/wiki/Adobe_RGB_color_space#ICC_PCS_color_image_encoding
        Self {
            date: Utc::now().format("%Y-%m-%d").to_string(),
            name: "Adobe RGB".to_string(),
            identifier: "ADOBE_RGB".to_string(),
            chromaticities: Self::chromaticities([&red_prim, &green_prim, &blue_prim]),
            red_prim,
            green_prim,
            blue_prim,
            white,
            white_points: Self::white_points(white),
            gamma: format!("{gamma:?}"),
            target_chromaticities_cie1931: [
                [0.6400, 0.3300], // Red
                [0.2100, 0.7100], // Green
                [0.1500, 0.0600], // Blue
            ],
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
        //
        // The display P3 red coordinate is outside the CIE 1931 gamut using the CIE 1931 1 nanometer
        // dataset as provided by the CIE.  To still match it, we need to desature it by adding white. It also mixes
        // in a bit of blue (1.4E-6) to get to to the target.
        // const D: f64 = 0.0005620; // Desaturation ratio

        const RED: [f64; 3] = [637.9658554073235, 23.274151215539906, 1.4261447449730065e-6];

        // Blue and Green Gaussian filters.
        const GREEN: [f64; 2] = [539.8416064376327, 21.411199475289777];
        const BLUE: [f64; 2] = [398.0273721579992, 52.55338039394701];

        let [red_prim, green_prim, blue_prim] =
            gaussian_filtered_primaries(D65.as_ref(), RED, GREEN, BLUE);
        let white = CieIlluminant::D65;
        let gamma = GammaCurve::new(
            5,
            [
                2.4,
                1.0 / 1.055,
                0.055 / 1.055,
                1.0 / 12.92,
                0.04045,
                0.0,
                0.0,
            ],
        );
        Self {
            date: Utc::now().format("%Y-%m-%d").to_string(),
            name: "Display P3".to_string(),
            identifier: "DISPLAY_P3".to_string(),
            chromaticities: Self::chromaticities([&red_prim, &green_prim, &blue_prim]),
            red_prim,
            green_prim,
            blue_prim,
            white,
            white_points: TemplateContext::white_points(white),
            gamma: format!("{gamma:?}"),
            target_chromaticities_cie1931: [
                [0.6800, 0.3200], // Red
                [0.2650, 0.6900], // Green
                [0.1500, 0.0600], // Blue
            ],
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
    pub fn cie_rgb() -> Self {
        let red_prim = Stimulus::new(Spectrum::from_wavelength_map(&[(700, 1.0)])); // 700.0 nm
        let green_prim = Stimulus::new(Spectrum::from_wavelength_map(&[(546, 0.9), (547, 0.1)])); // 546.1 nm
        let blue_prim = Stimulus::new(Spectrum::from_wavelength_map(&[(435, 0.2), (436, 0.8)])); // 435.8 nm
        let white = CieIlluminant::E;
        let white_points = TemplateContext::white_points(white);
        let gamma = GammaCurve::new(1, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        Self {
            date: Utc::now().format("%Y-%m-%d").to_string(),
            name: "CIE RGB".to_string(),
            identifier: "CIE_RGB".to_string(),
            chromaticities: Self::chromaticities([&red_prim, &green_prim, &blue_prim]),
            white,
            red_prim,
            green_prim,
            blue_prim,
            white_points,
            gamma: format!("{gamma:?}"),
            target_chromaticities_cie1931: [
                [0.7347, 0.2653],   // red
                [0.27368, 0.71741], // green
                [0.16653, 0.0088],  // blue
            ],
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

    pub fn white_points(white: CieIlluminant) -> Vec<[f64; 3]> {
        let mut out = Vec::new();
        for observer in Observer::iter() {
            let white_point = white.white_point(observer);
            out.push(white_point.values());
        }
        out
    }

    pub fn chromaticities(stimuli: [&Stimulus; 3]) -> Vec<[[f64; 2]; 3]> {
        let mut out = Vec::new();
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
            out.push(values);
        }
        out
    }

    pub fn stimuli_to_vecs(stimuli: &[Stimulus; 3]) -> [Vec<f64>; 3] {
        [
            stimuli[0].as_ref().values().to_vec(),
            stimuli[1].as_ref().values().to_vec(),
            stimuli[2].as_ref().values().to_vec(),
        ]
    }
    pub fn write(&self) -> Result<(), Box<dyn std::error::Error>> {
        let mut handlebars = handlebars::Handlebars::new();
        handlebars.set_strict_mode(true);
        handlebars
            .register_template_file("rgbspace", "xtask/templates/rgbspace.rs.hbs")
            .expect("Failed to register template");

        let contents = handlebars.render("rgbspace", self)?;
        //   print!("{}", contents);
        let path = "src/rgb/rgbspace";
        let file_path = format!("{path}/{}.rs", self.identifier.to_lowercase());
        std::fs::create_dir_all(path);
        std::fs::write(file_path, contents)?;
        Ok(())
    }
}

pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    TemplateContext::srgb().write()?;
    TemplateContext::adobe_rgb().write()?;
    TemplateContext::display_p3().write()?;
    TemplateContext::cie_rgb().write()?;
    println!("Generating RGB space data...");

    Ok(())
}
