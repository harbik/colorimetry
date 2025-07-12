use colorimetry::observer::Observer;
use colorimetry_plot::{axis::AxisSide, chart::Chart, svgdoc::SvgDocument};

pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    let observer = Observer::default();

    let svgdoc = SvgDocument::new(900, 1000)
        .add_css_rule(".fine-grid", "stroke: #88888888; stroke-width: 0.5;")
        .add_css_rule(".grid", "stroke: #88888888; stroke-width: 1.0;")
        .add_css_rule(
            ".chart-area",
            "fill: #DDDDDD; stroke: none; stroke-width: 0;",
        )
        .add_css_rule(".spectral-locus", "fill: white; stroke: black; stroke-width: 2; stroke-linecap: round;")
        .add_css_rule(".planckian-locus", "fill: none; stroke: gray; stroke-width: 4; stroke-linecap: round;")
        .add_css_rule("text", "fill: black; stroke: none; font-size: 12pt; font-family: sans-serif;");
    

    Chart::new(&svgdoc,[50, 50, 750, 850],
            [[0.0, 0.75], [0.0, 0.85]],
            Some("chart-area"),
            None,
        )
        .add_axis(Some("CIE 1931 x Chromaticity"), AxisSide::Bottom, 0.1, 6, true, Some("grid"))
        .add_axis(None, AxisSide::Bottom, 0.01,4, false, Some("fine-grid"))
        .add_axis(Some("y Chromaticity"), AxisSide::Left, 0.1, 6,  true, Some("grid"))
        .add_axis(None, AxisSide::Left, 0.01, 4, false, Some("fine-grid"))
        .draw_shape(observer.spectral_locus().into_iter().take(330), Some("spectral-locus"), None)
        .draw_grid(0.01, 0.01, Some("fine-grid"), None)
        .draw_grid(0.1, 0.1, Some("grid"), None)
        .draw_line(observer.planckian_locus(), Some("planckian-locus"), None)
        .draw_dot(1./3., 1./3., 5.0, None, None)
        .render();

    svgdoc.save("tmp/blackbody_locus.svg").unwrap();
    Ok(())
}