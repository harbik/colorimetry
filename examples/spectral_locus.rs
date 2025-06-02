// to run this example use:
//  `cargo run --example spectral_locus`
use colorimetry::observer::Observer;
use strum::IntoEnumIterator;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for observer in Observer::iter() {
        let wavelength_range = observer.data().spectral_locus_wavelength_range();
        println!(
            "Spectral locus for {observer} goes from {} to {}",
            wavelength_range.start(),
            wavelength_range.end()
        );
        println!("nm\tx\t y");
        for wavelength in wavelength_range {
            // unwrap OK because nm is in range
            let xyz = observer.data().xyz_at_wavelength(wavelength).unwrap();
            let chromaticity = xyz.chromaticity();
            println!(
                "{wavelength}\t{:.4}\t{:.4}",
                chromaticity.x(),
                chromaticity.y()
            );
        }
        println!("==========================");
    }

    Ok(())
}
