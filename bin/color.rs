use colorimetry::{
    colorant::Colorant, illuminant::D65, lab::CieLCh, observer::Observer, rgb::RgbSpace,
    spectrum::Spectrum,
};

use std::env;
use std::fs::File;
use std::error::Error;
use csv::Reader;

    /*
pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    let wavelengths = [
        379.6304167,
        385.5067911,
        385.5067911,
        389.7377807,
        390.6780006,
        396.2801441,
        399.6100896,
        406.0741015,
        413.2432782,
        420.7650374,
        428.1524795,
        452.4974591,
        478.3535064,
        504.2095537,
        530.065601,
        555.9216483,
    ];
    let reflectivities = [
        0.5317849104,
        0.578869226,
        0.6245618571,
        0.667396702,
        0.7116690332,
        0.7618164653,
        0.8117336725,
        0.8532283814,
        0.8941943853,
        0.9399369657,
        0.9887871218,
        0.9997185357,
        0.9997889018,
        0.9998592679,
        0.9999296339,
        1.0,
    ];
    let spectrum = Spectrum::linear_interpolate(&wavelengths, &reflectivities)?;
    let colorant = Colorant::new(spectrum)?;
    let lab = Observer::Cie1931.lab(&D65, &colorant);
    let white = Observer::Cie1931.lab(&D65, &Colorant::white());
    let de = lab.ciede2000(&white)?;
    let xyz = Observer::Cie1931.rel_xyz(&D65, &colorant);
    let rgb: [u8; 3] = xyz.xyz().rgb(RgbSpace::SRGB).compress().into();
    println!("{rgb:?} {de}");
    Ok(())
}
     */


fn main() -> Result<(), Box<dyn Error>> {
    let path = env::args().nth(1).expect("Usage: color <file.csv>");
    let file = File::open(path)?;
    let mut rdr = Reader::from_reader(file);

    let mut wavelengths = Vec::new();
    let mut values = Vec::new();

    for result in rdr.records() {
        let record = result?;
        let wl: f64 = record[0].trim().parse()?;
        let val: f64 = record[1].trim().parse()?;
        wavelengths.push(wl);
        values.push(val);
    }

    let spectrum = Spectrum::linear_interpolate(&wavelengths, &reflectivities)?;
    let colorant = Colorant::new(spectrum)?;
    let xyz = Observer::Cie1931.rel_xyz(&D65, &colorant);
    let lab = Lab::from_xyz(&xyz);
    let rgb: [u8; 3] = xyz.xyz().rgb(RgbSpace::SRGB).compress().into();
    
    println!("{rgb:?} {de}");
    Ok(())
}