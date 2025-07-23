use colorimetry::observer::Observer;
use colorimetry_plot::{
    chart::XYChart,
    svgdoc::{SvgDocument, NORTH_WEST, SOUTH_EAST},
};

const STYLE: &str = "
    :root {
        --chart-background-color: #888;
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
        fill: var(--chart-background-color);
        stroke: black;
        stroke-width: 5;
    }
    .spectral-locus {
        fill: white;
        stroke: black;
        stroke-width: 2;
        stroke-linecap: round;
    }
    .planckian-locus {
        fill: none;
        stroke: gray;
        stroke-width: 2;
        stroke-linecap: round;
    }
    text {
        fill: black;
        stroke: none;
        font-size: 12pt;
        font-family:sans-serif;
    }
    .white-point {
        stroke: black;
        stroke-width: 1;
    }
    text.white-point {
        fill: black;
        stroke: none;
        stroke-width: 0;
    }
";

pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    let observer = Observer::default();
    let d65 = observer.xyz_d65().chromaticity().to_tuple();
    let d50 = observer.xyz_d50().chromaticity().to_tuple();

    let chart = XYChart::new(
        "cie1931_chromaticity_diagram",
        (500, 500),
        (0.25..=0.45, 0.25..=0.45),
        (Some("chart"), None),
    )
    .add_ticks(0.01, 0.01, 4, Some("fine-grid"), None)
    .add_ticks(0.1, 0.1, 6, Some("fine-grid"), None)
    .add_x_labels(0.1, 3)
    .add_y_labels(0.1, 3)
    .x_axis_description("CIE 1931 x Chromaticity")
    .y_axis_description("CIE 1931 y Chromaticity")
    .draw_shape(
        observer.spectral_locus().into_iter().take(330),
        Some("spectral-locus"),
        None,
    )
    .draw_grid(0.01, 0.01, Some("fine-grid"), None)
    .draw_grid(0.1, 0.1, Some("grid"), None)
    .draw_line(observer.planckian_locus(), Some("planckian-locus"), None)
    .annotate(
        (1. / 3., 1. / 3.),
        3.0,
        (SOUTH_EAST, 20),
        "E",
        Some("white-point"),
        None,
    )
    .annotate(d65, 3.0, (NORTH_WEST, 20), "D65", Some("white-point"), None)
    .annotate(d50, 3.0, (NORTH_WEST, 20), "D50", Some("white-point"), None);

    SvgDocument::new(800, 800, STYLE)
        .add_svg(Box::new(chart.clone()))
        .save("tmp/blackbody_locus.svg")
}
