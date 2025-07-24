use colorimetry::observer::Observer;
use colorimetry_plot::{
    chart::XYChart, style_attr, svgdoc::{SvgDocument, NORTH_WEST, SOUTH_EAST}
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
        (500, 500),
        (0.25..=0.45, 0.25..=0.45),
        style_attr!(class: "chart"),
    )
    .ticks(0.01, 0.01, 4, style_attr!(class:"fine-grid"))
    .ticks(0.1, 0.05, 8, style_attr!(class:"grid"))
    .ticks(0.1, 0.1, 12, style_attr!(class:"grid"))
    .x_labels(0.1, 3)
    .y_labels(0.1, 3)
    .x_axis_description("CIE 1931 x Chromaticity")
    .y_axis_description("CIE 1931 y Chromaticity")
    .plot_shape(
        observer.spectral_locus().into_iter().take(330),
        style_attr!(class:"spectral-locus"),
    )
    .plot_grid(0.01, 0.01, style_attr!(class: "fine-grid"))
    .plot_grid(0.1, 0.1, style_attr!(class: "grid"))
    .plot_poly_line(observer.planckian_locus(), style_attr!(class:"planckian-locus"))
    .label_pin(
        (1. / 3., 1. / 3.),
        3.0,
        (SOUTH_EAST, 20),
        "E",
        style_attr!(class: "white-point"),
    )
    .label_pin(d65, 3.0, (NORTH_WEST, 20), "D65", style_attr!(class: "white-point"))
    .label_pin(d50, 3.0, (NORTH_WEST, 20), "D50", style_attr!(class: "white-point"));

    SvgDocument::new(800, 800, STYLE)
        .add_svg(Box::new(chart.clone()))
        .save("tmp/blackbody_locus.svg")
}
