use colorimetry::observer::Observer;
use colorimetry_plot::{
    axis::AxisSide,
    chart::XYChart,
    svgdoc::{SvgDocument, NORTH_WEST, SOUTH_EAST},
};

pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    let observer = Observer::default();
    let d65 = observer.xyz_d65().chromaticity().to_tuple();
    let d50 = observer.xyz_d50().chromaticity().to_tuple();

    let mut svgdoc = SvgDocument::new(800, 800)
        .add_css_rule(".fine-grid", "stroke: #88888888; stroke-width: 0.5;")
        .add_css_rule(".grid", "stroke: #88888888; stroke-width: 1.0;")
        .add_css_rule(
            ".chart-area",
            "fill: #DDDDDD; stroke: none; stroke-width: 0;",
        )
        .add_css_rule(
            ".spectral-locus",
            "fill: white; stroke: black; stroke-width: 2; stroke-linecap: round;",
        )
        .add_css_rule(
            ".planckian-locus",
            "fill: none; stroke: gray; stroke-width: 2; stroke-linecap: round;",
        )
        .add_css_rule(
            "text",
            "fill: black; stroke: none; font-size: 12pt; font-family: sans-serif;",
        )
        .add_css_rule(".white-point", "stroke: black; stroke-width: 1;")
        .add_css_rule(
            "text.white-point",
            "fill: black; stroke: None; stroke-width: 0;",
        );

    let chart = XYChart::new(
        "cie1931_chromaticity_diagram",
        (500, 500),
        (0.25..=0.45, 0.25..=0.45),
        (Some("chart-area"), None),
    )
        .add_axis(
            Some("CIE 1931 x Chromaticity"),
            AxisSide::Bottom,
            0.1,
            6,
            true,
            Some("grid"),
        )
        .add_axis(None, AxisSide::Bottom, 0.01, 4, false, Some("fine-grid"))
        .add_axis(
            Some("y Chromaticity"),
            AxisSide::Left,
            0.1,
            6,
            true,
            Some("grid"),
        )
        .add_axis(None, AxisSide::Left, 0.01, 4, false, Some("fine-grid"))
        .add_axis(
            Some("CIE 1931 x Chromaticity"),
            AxisSide::Top,
            0.1,
            6,
            true,
            Some("grid"),
        )
        .add_axis(None, AxisSide::Top, 0.01, 4, false, Some("fine-grid"))
        .add_axis(
            Some("y Chromaticity"),
            AxisSide::Right,
            0.1,
            6,
            true,
            Some("grid"),
        )
        .add_axis(None, AxisSide::Right, 0.01, 4, false, Some("fine-grid"))
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

    svgdoc.add_svg(Box::new(chart.clone()));
    //    svgdoc.place_center(chart);

    svgdoc.save("tmp/blackbody_locus.svg").unwrap();
    Ok(())
}
