// to run this example use:
//  `cargo run --features cie-illuminants,cri --example cri`
use colored::Colorize;
use colorimetry::{illuminant::CCT, prelude::*};
use strum::IntoEnumIterator;

/// Prints the standard illuminants in the library, with their elated color temperature, with
/// parameters distance to the Planckian, the general Color Rendering Index Ra, and the spectal
/// color rendering index R9.
fn main() -> Result<(), Box<dyn std::error::Error>> {
    for spc in CieIlluminant::iter() {
        // Calculate CRI parameters
        let ill: &Illuminant = spc.as_ref();
        let cri: colorimetry::illuminant::CRI = ill.try_into()?;
        let ra = cri.ra();
        let r9 = cri[9];

        // Calculate Correlated Color Temperature
        let xyz = Cie1931.xyz(&spc, None);
        let cct: CCT = xyz.try_into()?;
        let [t, d] = cct.into();
        let tint = d * 1000.0;

        // Output Results
        let s = format!("{}", spc);
        println!("\n{}", s.bold().underline());
        println!("CCT  = {t:.0} Kelvin");
        println!("Duv  = {d:.5}");
        println!("Tint = {tint:.1}");
        println!("CRI  = {ra:.0}");
        println!("R9   = {r9:.0}");
        println!();
    }

    Ok(())
}
