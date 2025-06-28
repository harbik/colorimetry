use colorimetry::{
    illuminant::{CieIlluminant, CieIlluminant::D65},
    lab::CieLChGamut,
    observer::{Observer, Observer::Cie1931},
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
            .caption(&caption, ("sans-serif", 20))
            .margin(10)
            .x_label_area_size(30)
            .y_label_area_size(30)
            .build_cartesian_2d(0f64..0.8f64, 0f64..0.9f64)?;

        chart.configure_mesh().draw()?;

        Ok(Plot { chart })
    }

    fn add_spectral_locus(&mut self, observer: Observer) -> Result<(), Box<dyn std::error::Error>> {
        let spectral_locus = observer.monochromes(D65);
        let mut points: Vec<(f64, f64)> = spectral_locus
            .iter()
            .map(|&(_wl, rxyz)| {
                let chromaticity = rxyz.xyz().chromaticity();
                (chromaticity.x(), chromaticity.y())
            })
            .collect();
        // add purple line
        let first_point = spectral_locus[0].1.xyz().chromaticity();
        points.push((first_point.x(), first_point.y()));
        self.chart.draw_series(LineSeries::new(points, &BLACK))?;

        Ok(())
    }

    fn plot_iso_hue_lines(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        let gamut = CieLChGamut::new(Observer::Cie1931, CieIlluminant::D65);

        gamut.colors().iter().for_each(|([l_bin, h_bin], c_bin)| {
            let lab = gamut.bins_to_cielab(*l_bin, *c_bin, *h_bin);
            if !lab.is_valid() {
                eprintln!("Invalid color at L*={}, C*={}, h*={}", l_bin, c_bin, h_bin);
            };
        });

        for h in 0..72 {
            let points: Vec<(f64, f64)> = (1..=100)
                .filter_map(|l| {
                    gamut.max_chroma(l, h).map(|c| {
                        let lab = gamut.bins_to_cielab(l, c, h);
                        let [x, y] = lab.xyz().xyz().chromaticity().to_array();
                        (x, y)
                    })
                })
                .collect();
            self.chart.draw_series(LineSeries::new(points, &BLACK))?;
        }
        Ok(())
    }

    // fn save_svg(&mut self) -> Result<(), Box<dyn std::error::Error>> {
    //     self.root.present()?;
    //     Ok(())
    // }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create a new plot
    let mut plot = Plot::new("Iso-hue lines in CIE 1931 chromaticity diagram")?;

    // Plot the spectral locus
    plot.add_spectral_locus(Cie1931).unwrap();

    plot.plot_iso_hue_lines()
}
