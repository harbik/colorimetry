// to run this example use:
//  `cargo run --example spectral_locus`
use colorimetry::observer::Observer;
use strum::IntoEnumIterator;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for observer in Observer::iter() {
        let nm_min = observer.data().spectral_locus_nm_min();
        let nm_max = observer.data().spectral_locus_nm_max();
        println!("Spectral locus for {observer:?} goes from {nm_min} to {nm_max}:");
        println!("nm\tx\t y");
        for nm in nm_min..=nm_max {
            // unwrap OK because nm is in range
            let xyz = observer.data().spectral_locus_by_nm(nm).unwrap();
            let [x, y] = xyz.chromaticity();
            println!("{nm}\t{x:.4}\t{y:.4}");
        }
        println!("==========================");
    }

    Ok(())
}
