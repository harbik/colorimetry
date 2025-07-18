// cargo run --bin color --features="csv" bin/test.csv1

use colorimetry::lab::CieLab;
use colorimetry::{
    colorant::Colorant, illuminant::D65, observer::Observer, rgb::RgbSpace,
    spectrum::Spectrum,
};

use std::env;
use std::fs::File;
use std::error::Error;
use csv::ReaderBuilder;

fn main() -> Result<(), Box<dyn Error>> {
    let path = env::args().nth(1).expect("Usage: color <file.csv>");
    let file = File::open(path)?;
    let mut rdr = ReaderBuilder::new()
        .comment(Some(b'#'))  // 👈 ignore lines starting with '#'
        .from_reader(file);

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
    
    println!("{xyz:?} {lab:?} {rgb:?}");
    Ok(())
}