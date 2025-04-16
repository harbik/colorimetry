// to run this example use: 
//  `cargo run --example spectral_locus`
use colorimetry::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>>{
    let nm_min = CIE1931.spectral_locus_nm_min();
    let nm_max = CIE1931.spectral_locus_nm_max();
    println!("Spectral locus");
    println!("nm\tx\t y");
    for nm in nm_min..=nm_max {
        // unwrap OK because nm is in range
        let xyz = CIE1931.spectral_locus_by_nm(nm).unwrap();
        let [x, y] = xyz.chromaticity();
        println!("{nm}\t{x:.4}\t{y:.4}");
    }
    Ok(())

}