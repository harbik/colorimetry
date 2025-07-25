use colorimetry_plot::{chart::XYChart, style_attr, svgdoc::SvgDocument};

const STYLE: &str = include_str!("sine_curve.scss");

pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    let pi = std::f64::consts::PI;
    let pi2 = 2.0 * pi;

    const N: i32 = 100;
    let sin: Vec<(f64, f64)> = (0..N)
        .map(|i| {
            let x = i as f64 * pi2 / (N as f64 - 1.0);
            (x, x.sin())
        })
        .collect();

    let chart = XYChart::new((600, 300), (..pi2, -1.1..1.1))
        .x_labels(1.0, 10)
        .y_labels(0.5, 10)
        .x_axis_description("x")
        .y_axis_description("sin(x)")
        .plot_grid(0.2, 0.1, style_attr!(class: "fine-grid"))
        .plot_grid(1.0, 0.5, style_attr!(class: "grid"))
        .plot_poly_line(sin, style_attr!(class:"curve"))
        .plot_poly_line(vec![(0.0, 0.0), (pi2, 0.0)], style_attr!(class:"base-line"));

    SvgDocument::new()
        .append_scss(STYLE)
        .add_svg(Box::new(chart))
        .save("tmp/sine_curve.svg")
}
