use colorimetry::{observer::Observer, rgb::RgbSpace::SRGB};
use colorimetry_plot::{chart::XYChromaticity, style_attr, svgdoc::SvgDocument};

const STYLE: &str = "
    :root {
        --plot-background-color: #888;
        --spectral-locus-color: #DDD;
        --axis-color: #DDD;
        --grid-color: #AAA;
    }
    .fine-grid {
        stroke: var(--grid-color);
        stroke-width: 0.5;
    }
    .grid {
        stroke: var(--grid-color); 
        stroke-width: 1.0;
    }
    .chart {
        fill: var(--plot-background-color);
        stroke: none;
        stroke-width: 0;
    }
    .spectral-locus {
        fill: var(--spectral-locus-color);
        stroke: none;
        stroke-width: 0;
        stroke-linecap: round;,
    }
    .spectral-locus-ticks {
        fill:none;
        stroke:var(--plot-background-color);
        stroke-width: 1;
        stroke-linecap: round;
    }
    .planckian {
        stroke: var(--plot-background-color);
        fill: none;
        stroke-width: 2;
        stroke-linecap: round;
    }
    text {
        fill: black;
        stroke: none;
        font-size: 12pt;
        font-family: sans-serif;
    }
    text.spectral-locus-labels {
        fill: var(--spectral-locus-color);
        stroke: none;
        stroke-width: 0;
    }
    .white-point {
        stroke: black;
        stroke-width: 1;
    }
    text.white-point {
        fill:black;
        stroke:none;
        stroke-width:0;
    }
";

pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    let observer = Observer::default();

    let xy_chromaticity = XYChromaticity::new(
        observer,
        (775, 875),
        (-0.025..=0.75, 0.0..=0.875),
        style_attr!(class: "chart"),
    )
    .ticks(0.01, 0.01, 5, style_attr!(class: "fine-grid"))
    .ticks(0.1, 0.1, 10, style_attr!(class: "grid"))
    .x_labels(0.1, 10)
    .x_axis_description("CIE 1931 x Chromaticity")
    .y_labels(0.1, 10)
    .y_axis_description("CIE 1931 y Chromaticity")
    .plot_spectral_locus(style_attr!(class: "spectral-locus"))
    .draw_spectral_locus_ticks(440..=650, 10, 15, style_attr!(class: "spectral-locus-ticks"))
    .draw_spectral_locus_ticks(460..=630, 1, 7, style_attr!(class:"spectral-locus-ticks"))
    .draw_spectral_locus_labels(460..=620, 10, 2, style_attr!(class:"spectral-locus-labels"))
    .plot_rgb_gamut(SRGB, style_attr!())
    .plot_planckian_locus(style_attr!(class:"planckian"))
    .plot_grid(0.01, 0.01, style_attr!(class: "fine-grid"))
    .plot_grid(0.1, 0.1, style_attr!(class: "grid"));

    SvgDocument::new(1000, 1100, STYLE)
        .add_svg(Box::new(xy_chromaticity))
        .save("tmp/srgb_gamut.svg")
}
