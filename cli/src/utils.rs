use colored::Colorize;
use colorimetry::illuminant::CieIlluminant;
use colorimetry::lab::CieLab;
use colorimetry::{colorant::Colorant, observer::Observer, rgb::RgbSpace, spectrum::Spectrum};
use csv::ReaderBuilder;
use std::error::Error;
use std::fs::File;
use strum::IntoEnumIterator;

/// The input file should be a CSV file with two columns:
///
/// - Wavelength (in nm)
/// - Reflectivity (as a fraction between 0 and 1)
///
/// The output will be the XYZ, CIE Lab, and sRGB values computed from the reflectivity data.
/// The CSV file can have comments starting with `#`.
pub fn colorvalue(file: String, illuminant: String, format: String) -> Result<(), Box<dyn Error>> {
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
    let illuminant = parse_illuminant(&illuminant)?;

    match format.to_lowercase().as_str() {
        "xyz" => print_xyz(&illuminant, &colorant),
        "yxy" => print_yxy(&illuminant, &colorant),
        "lab" | "cielab" => print_lab(&illuminant, &colorant),
        "srgb" => print_rgb(&illuminant, &colorant),
        _ => return Err(format!("Unsupported format: {format}").into()),
    }
    Ok(())
}

fn print_yxy(illuminant: &CieIlluminant, colorant: &Colorant) {
    let xyz = Observer::Cie1931.rel_xyz(illuminant, colorant);
    let (x, y) = xyz.xyz().chromaticity().to_tuple();
    let [_x, yy, _z] = xyz.xyz().values();
    println!("{}", "Yxy (Y,x,y) chromaticity values".bold());
    println!("{yy:.2}\t{x:.4}\t{y:.4}");
}

fn print_xyz(illuminant: &CieIlluminant, colorant: &Colorant) {
    let xyz = Observer::Cie1931.rel_xyz(illuminant, colorant);
    let [x, y, z] = xyz.xyz().values();
    println!("{}", "YYZ (X,Y,Z) Tristimulus Values".bold());
    println!("{x:.2}\t{y:.2}\t{z:.2}");
}

fn print_lab(illuminant: &CieIlluminant, colorant: &Colorant) {
    let xyz = Observer::Cie1931.rel_xyz(illuminant, colorant);
    let lab = CieLab::from_rxyz(xyz);
    println!("{}", "CIELAB (L*, a*, b*) values".bold());
    println!("{:.2}\t{:.2}\t{:.2}", lab.l(), lab.a(), lab.b());
}

fn print_rgb(illuminant: &CieIlluminant, colorant: &Colorant) {
    let xyz = Observer::Cie1931.rel_xyz(illuminant, colorant);
    let rgb: [u8; 3] = xyz.xyz().rgb(RgbSpace::SRGB).compress().into();
    println!("{}", "sRGB (R,G,B) values".bold());
    println!("{:3}\t{:3}\t{:3}", rgb[0], rgb[1], rgb[2]);
}

fn parse_illuminant(illuminant: &str) -> Result<CieIlluminant, Box<dyn Error>> {
    for ill in CieIlluminant::iter() {
        if format!("{ill}").to_lowercase() == illuminant.to_lowercase() {
            return Ok(ill);
        }
    }
    eprintln!("No exact match found for illuminant: {illuminant}");
    eprintln!("Supported illuminants are:");
    for ill in CieIlluminant::iter() {
        eprintln!("- {ill}");
    }
    Err(format!("Unsupported illuminant: {illuminant}").into())
}
