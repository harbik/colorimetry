use std::fs;

use plotters::style::ShapeStyle;

use colorimetry::{
    lab::CieLChGamut,
    observer::{
        Observer::{self, Cie1931},
        SpectralLocus,
    },
    rgb::RgbSpace,
};

use plotters::{coord::types::RangedCoordf64, prelude::*};

const SVG_FILE_NAME: &str = "tmp/iso_hues.svg";
const DATA_FILE_NAME: &str = "tmp/iso_hues.csv";

#[derive()]
struct Plot<'a> {
    chart: ChartContext<'a, SVGBackend<'a>, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    data: Vec<String>,
}

impl<'a> Plot<'a> {
    fn new(caption: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let root = SVGBackend::new(SVG_FILE_NAME, (820, 920)).into_drawing_area();
        let mut chart = ChartBuilder::on(&root)
            .caption(caption, ("sans-serif", 15))
            .margin(20)
            .x_label_area_size(30)
            .y_label_area_size(30)
            .build_cartesian_2d(0f64..0.8f64, 0f64..0.9f64)?;

        chart.configure_mesh().draw()?;

        Ok(Plot {
            chart,
            data: Vec::new(),
        })
    }

    fn add_spectral_locus(&mut self, observer: Observer) -> Result<(), Box<dyn std::error::Error>> {
        let spectral_locus = SpectralLocus::new(observer);
        self.chart.draw_series(LineSeries::new(
            spectral_locus.into_iter(),
            BLACK.stroke_width(2),
        ))?;

        Ok(())
    }

    fn plot_iso_hue_lines(
        &mut self,
        l_step: f64, // Lightness range is 0.0 to 100.0
        h_step: f64, // Hue range is from 0.0 to 360.0
        rgb_space: RgbSpace,
        style: &ShapeStyle,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let gamut = CieLChGamut::new(Cie1931, rgb_space);
        self.data
            .push("Hue, Lightness, Chroma, x, y, R, G, B".to_string());

        let mut h = 0.0;
        while h < 360.0 {
            let mut points = Vec::new();

            let mut l = 0.0;
            while l < 100.0 {
                if let Some(lch) = gamut.max_chroma(l, h) {
                    // ^^ replace max_chroma_in_gamut with max_chroma if you want to plot the full chromatricity diagram
                    let [x, y] = lch.rxyz().xyz().chromaticity().to_array();
                    points.push((x, y));
                    self.data.push(format!(
                        "{:3.0}, {:4.2}, {:5.2}, {x:.5}, {y:.5}",
                        h,
                        l,
                        lch.c()
                    ));
                }
                l += l_step;
            }
            self.chart.draw_series(LineSeries::new(points, *style))?;
            h += h_step;
        }
        Ok(())
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create a new plot
    let mut plot = Plot::new("Iso-hue lines in CIE 1931 chromaticity diagram")?;

    // Plot the spectral locus
    plot.add_spectral_locus(Cie1931).unwrap();

    plot.plot_iso_hue_lines(0.1, 5.0, RgbSpace::Adobe, &BLACK.stroke_width(1))?;

    let data = plot.data.join("\n");
    fs::write(DATA_FILE_NAME, data)?;
    Ok(())
}
