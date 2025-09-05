// to run this example use:
//  `cargo run --features cct --example cct -- -x 0.333 -y 0.333`

use clap::Parser;
use colored::Colorize;

/// Calculate Correlated Color Temperature, Planckian Distance, and Tint, for a given set of (x,y)
/// chromaticity values.
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// CIE 1931 x Chromaticity coordinate
    #[arg(short, long)]
    x: f64,

    /// CIE 1931 y Chromaticity coordinate
    #[arg(short, long)]
    y: f64,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let chromaticity = colorimetry::xyz::Chromaticity::new(args.x, args.y);
    let xyz: colorimetry::xyz::XYZ =
        colorimetry::xyz::XYZ::from_chromaticity(chromaticity, None, None)?;
    let cct: colorimetry::illuminant::CCT = xyz.try_into()?;
    let [t, d] = cct.into();
    let tint = d * 1000.0;
    println!("\n{}\n", "Result".bold().underline());
    println!("CCT  = {t:.1} Kelvin");
    println!("Duv  = {d:.5}");
    println!("Tint = {tint:.1}");
    println!();

    Ok(())
}
