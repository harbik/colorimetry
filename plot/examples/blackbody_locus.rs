use colorimetry::observer::Observer;
use colorimetry_plot::{axis::AxisSide, canvas::Canvas};

pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    let observer = Observer::default();

    let canvas = Canvas::new(900, 1000)
        .add_style(".fine-grid", "stroke: #88888888; stroke-width: 0.5;")
        .add_style(".grid", "stroke: #88888888; stroke-width: 1.0;")
        .add_style(
            ".chart-area",
            "fill: #DDDDDD; stroke: none; stroke-width: 0;",
        )
        .add_style(".spectral-locus", "fill: white; stroke: black; stroke-width: 2; stroke-linecap: round;")
        .add_style(".planckian-locus", "fill: none; stroke: gray; stroke-width: 4; stroke-linecap: round;")
        .add_style("text", "fill: black; stroke: none; font-size: 12pt; font-family: sans-serif;");
    

    canvas
        .add_chart(
            [50, 50, 750, 850],
            [[0.0, 0.75], [0.0, 0.85]],
            Some("chart-area"),
            None,
        )
//        .add_x_axis(0.1, 10, None, None)
        .add_axis(Some("CIE 1931 x Chromaticity"), AxisSide::Bottom, 0.1, 6, true, Some("grid"))
        .add_axis(None, AxisSide::Bottom, 0.01,4, false, Some("fine-grid"))
        .add_axis(Some("CIE 1931 y Chromaticity"), AxisSide::Left, 0.1, 6,  true, Some("grid"))
        .add_axis(None, AxisSide::Left, 0.01, 4, false, Some("fine-grid"))
        .draw_area(observer.spectral_locus().into_iter().take(330), Some("spectral-locus"), None)
        .draw_grid(0.01, 0.01, Some("fine-grid"), None)
        .draw_grid(0.1, 0.1, Some("grid"), None)
        .draw_line(observer.planckian_locus(), Some("planckian-locus"), None)
        .draw_dot(1./3., 1./3., 5.0, None, None)
        .render();

    canvas.save("tmp/blackbody_locus.svg").unwrap();
    Ok(())
}