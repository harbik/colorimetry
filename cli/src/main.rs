//! A command-line tool for use with the `colorimetry` library.
//!
//! # Installation
//!
//! You can install the `colorimetry` CLI tool like this:
//!
//! ```bash
//! cargo install colorimetry-cli
//! ```
//!
//! # Usage
//!
//! You can run using:
//!
//! ```bash
//! colorimetry --help
//! ```
use colorimetry::lab::CieLab;
use colorimetry::{
    colorant::Colorant, illuminant::D65, observer::Observer, rgb::RgbSpace, spectrum::Spectrum,
};

use clap::{Parser, Subcommand};
use csv::ReaderBuilder;
use std::error::Error;
use std::fs::File;

#[derive(Parser)]
#[command(name = "colorimetry")]
#[command(about = "A CLI tool for colorimetry operations", long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Compute various color matrics from measured reflectivity data.
    Colorpatch {
        /// Input CSV file
        file: String,
    },
}

/// The input file should be a CSV file with two columns:
///
/// - Wavelength (in nm)
/// - Reflectivity (as a fraction between 0 and 1)
///
/// The output will be the XYZ, CIE Lab, and sRGB values computed from the reflectivity data.
/// The CSV file can have comments starting with `#`.
pub fn colorpatch(file: String) -> Result<(), Box<dyn Error>> {
    // Open the file relative to the current working directory
    // use cargo run colorpatch -- data/test.csv as an example
    let file = File::open(std::path::Path::new(&file))?;
    let mut rdr = ReaderBuilder::new().comment(Some(b'#')).from_reader(file);

    let mut wavelengths = Vec::new();
    let mut reflectivities = Vec::new();

    for result in rdr.records() {
        let record = result?;
        let wl: f64 = record[0].trim().parse()?;
        let val: f64 = record[1].trim().parse()?;
        wavelengths.push(wl);
        reflectivities.push(val);
    }

    let spectrum = Spectrum::linear_interpolate(&wavelengths, &reflectivities)?;
    let colorant = Colorant::new(spectrum)?;
    let xyz = Observer::Cie1931.rel_xyz(&D65, &colorant);
    let lab = CieLab::from_rxyz(xyz);
    let rgb: [u8; 3] = xyz.xyz().rgb(RgbSpace::SRGB).compress().into();

    let [x, y, z] = xyz.xyz().values();
    let xyz = format!("\t{x:.2} {y:.2} {z:.2}");
    let lab = format!("{:.2} {:.2} {:.2}", lab.l(), lab.a(), lab.b());
    let rgb = format!("\t{:3} {:3} {:3}", rgb[0], rgb[1], rgb[2]);
    println!("XYZ: {xyz}\nCieLab: {lab}\nsRGB: {rgb}");
    Ok(())
}

pub fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Colorpatch { file } => colorpatch(file)?,
    }

    Ok(())
}
