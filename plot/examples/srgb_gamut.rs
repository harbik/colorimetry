use colorimetry::{observer::Observer, rgb::RgbSpace::SRGB};
use colorimetry_plot::{chromaticity::XYChromaticity, svgdoc::SvgDocument};

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
    .chart-area {
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
        "cie1931_chromaticity_diagram",
        observer,
        (775, 875),
        (-0.025..=0.75, 0.0..=0.875),
        (Some("chart-area"), None),
    )
    .add_ticks(0.01, 0.01, 5, Some("fine-grid"), None)
    .add_ticks(0.1, 0.1, 10, Some("grid"), None)
    /*
    .add_axis(
        Some("CIE 1931 x Chromaticity"),
        AxisSide::Bottom,
        0.1,
        6,
        true,
        Some("grid"),
        None
    )
     */
    .add_x_labels(0.1, 10)
    .x_axis_description("CIE 1931 x Chromaticity")
    .add_y_labels(0.1, 10)
    .y_axis_description("CIE 1931 y Chromaticity")
    /*
    .add_axis(
        Some("CIE 1931 y Chromaticity"),
        AxisSide::Left,
        0.1,
        6,
        true,
        Some("grid"),
        None
    )
     */
    .draw_spectral_locus(Some("spectral-locus"), None)
    .draw_spectral_locus_ticks(440..651, 10, 15, Some("spectral-locus-ticks"), None)
    .draw_spectral_locus_ticks(460..631, 1, 7, Some("spectral-locus-ticks"), None)
    .draw_spectral_locus_labels(460..=620, 10, 2, Some("spectral-locus-labels"), None)
    .draw_rgb_gamut(SRGB, None, None)
    .draw_planckian_locus(Some("planckian"), None)
    .draw_grid(0.01, 0.01, Some("fine-grid"), None)
    .draw_grid(0.1, 0.1, Some("grid"), None);

    Ok(SvgDocument::new(1000, 1100, STYLE)
        .add_svg(Box::new(xy_chromaticity))
        .save("tmp/srgb_gamut.svg")?)
}
