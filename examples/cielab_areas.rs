use plotters::style::ShapeStyle;

use colorimetry::{
    lab::{CieLCh, CieLChGamut},
    observer::{
        Observer::{self, Cie1931},
        SpectralLocus,
    },
    rgb::RgbSpace,
};

use plotters::{coord::types::RangedCoordf64, prelude::*};

const SVG_FILE_NAME: &str = "tmp/cielab_areas.svg";

#[derive()]
struct Plot<'a> {
    chart: ChartContext<'a, SVGBackend<'a>, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
}

impl<'a> Plot<'a> {
    fn new(caption: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let root = SVGBackend::new(SVG_FILE_NAME, (800, 900)).into_drawing_area();
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

    fn hue_lines(
        &mut self,
        l_step: f64, // Lightness range is 0.0 to 100.0
        h_step: f64, // Hue range is from 0.0 to 360.0
        rgb_space: RgbSpace,
        style: &ShapeStyle,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let gamut = CieLChGamut::new(Cie1931, rgb_space);
        let mut h = 0.0;
        while h < 360.0 {
            let mut points = Vec::new();

            let mut l = 0.0;
            while l < 100.0 {
                if let Some(lch) = gamut.max_chroma_in_gamut(l, h) {
                    // ^^ replace max_chroma_in_gamut with max_chroma if you want to plot the full chromatricity diagram
                    let [x, y] = lch.rxyz().xyz().chromaticity().to_array();
                    points.push((x, y));
                }
                l += l_step;
            }
            self.chart.draw_series(LineSeries::new(points, *style))?;
            h += h_step;
        }
        Ok(())
    }

    fn cielab_area(
        &mut self,
        l_step: f64,
        h_step: f64,
        rgb_space: RgbSpace,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let to_point = |lch: CieLCh| {
            let [x, y] = lch.rxyz().xyz().chromaticity().to_array();
            (x, y)
        };

        let mut l = 0.0;
        let gamut = CieLChGamut::new(Cie1931, rgb_space);
        while l + l_step <= 100.0 {
            let mut h = 0.0;
            while h + h_step <= 360.0 {
                let mut points = Vec::new();
                if let Some(lch) = gamut.max_chroma_in_gamut(l + l_step, h + h_step) {
                    points.push(to_point(lch));
                    let lch_close = lch;
                    if let Some(lch) = gamut.max_chroma_in_gamut(l, h + h_step) {
                        points.push(to_point(lch));
                        if let Some(lch) = gamut.max_chroma_in_gamut(l, h) {
                            points.push(to_point(lch));
                            if let Some(lch) = gamut.max_chroma_in_gamut(l + l_step, h) {
                                points.push(to_point(lch));
                                points.push(to_point(lch_close));
                                let [r, g, b] = gamut
                                    .max_chroma_in_gamut(l, h)
                                    .unwrap()
                                    .rgb(rgb_space)
                                    .compress()
                                    .values();
                                let fill_color = RGBColor(
                                    (r * 255.0).clamp(0.0, 255.0) as u8,
                                    (g * 255.0).clamp(0.0, 255.0) as u8,
                                    (b * 255.0).clamp(0.0, 255.0) as u8,
                                );
                                self.chart.draw_series(AreaSeries::new(
                                    points,
                                    0.0,
                                    ShapeStyle {
                                        color: fill_color.to_rgba(),
                                        filled: true,
                                        stroke_width: 0,
                                    },
                                ))?;
                            }
                        }
                    }
                }
                h += h_step;
            }
            l += l_step;
        }
        Ok(())
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create a new plot
    let mut plot = Plot::new("Iso-hue lines in CIE 1931 chromaticity diagram")?;

    // Plot the spectral locus
    plot.add_spectral_locus(Cie1931).unwrap();

    let l_step = 0.2;
    let h_step = 5.0;
    let colorspace = RgbSpace::SRGB;

    plot.cielab_area(l_step, h_step, colorspace)?;
    plot.hue_lines(
        l_step,
        h_step,
        colorspace,
        &ShapeStyle {
            color: RGBAColor(255, 255, 255, 0.5),
            filled: false,
            stroke_width: 1,
        },
    )?;
    // plot.plot_iso_lightness_lines(RgbSpace::Adobe, &BLACK.stroke_width(1))?;

    Ok(())
}
