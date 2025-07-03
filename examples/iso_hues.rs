use colorimetry::{
    lab::CieLChGamut,
    observer::{
        Observer::{self, Cie1931},
        SpectralLocus,
    },
    rgb::RgbSpace,
};

use plotters::{coord::types::RangedCoordf64, prelude::*};

const OUT_FILE_NAME: &str = "tmp/iso_hues.svg";

#[derive()]
struct Plot<'a> {
    chart: ChartContext<'a, SVGBackend<'a>, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
}

impl<'a> Plot<'a> {
    fn new(caption: &str) -> Result<Self, Box<dyn std::error::Error>> {
        std::fs::create_dir_all("plots")?;
        let root = SVGBackend::new(OUT_FILE_NAME, (800, 800)).into_drawing_area();
        let mut chart = ChartBuilder::on(&root)
            .caption(caption, ("sans-serif", 20))
            .margin(10)
            .x_label_area_size(30)
            .y_label_area_size(30)
            .build_cartesian_2d(0f64..0.8f64, 0f64..0.9f64)?;

        chart.configure_mesh().draw()?;

        Ok(Plot { chart })
    }

    fn add_spectral_locus(&mut self, observer: Observer) -> Result<(), Box<dyn std::error::Error>> {
        let spectral_locus = SpectralLocus::new(observer);
        self.chart
            .draw_series(LineSeries::new(spectral_locus.into_iter(), &BLACK))?;

        Ok(())
    }

    fn plot_iso_hue_lines_adobe(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        let gamut = CieLChGamut::new(Cie1931, RgbSpace::Adobe);
        for h in (0..360).step_by(5) {
            //  let h = 85;
            let mut points = Vec::with_capacity(100);
            for l in (0..1000).rev() {
                //  let lch = gamut.full(l as f64, h as f64);
                if let Some(lch) = gamut.rgb_gamut(l as f64 / 10.0, h as f64) {
                    let [x, y] = lch.rxyz().xyz().chromaticity().to_array();
                    points.push((x, y));
                    let [rr, gg, bb] = lch.rxyz().xyz().rgb(RgbSpace::Adobe).values();
                    println!(
                        " {:3.0}, {:4.2}, {:5.2}, {x:.5}, {y:.5}, {rr:.3}, {gg:.3}, {bb:.3}",
                        lch.h(),
                        lch.l(),
                        lch.c()
                    );
                } else {
                    break;
                }
            }
            self.chart.draw_series(LineSeries::new(points, &RED))?;
        }
        Ok(())
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    if !cfg!(feature = "gamut-tables") {
        eprintln!("This example requires the 'gamut-tables' feature to be enabled.");
        return Ok(());
    }
    // Create a new plot
    let mut plot = Plot::new("Iso-hue lines in CIE 1931 chromaticity diagram")?;

    // Plot the spectral locus
    plot.add_spectral_locus(Cie1931).unwrap();

    plot.plot_iso_hue_lines_adobe()?;
    Ok(())
}
