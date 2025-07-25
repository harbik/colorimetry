use colorimetry_plot::{style_attr::class, chart::XYChart, svgdoc::SvgDocument};

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
        .x_labels(1.0, 10, None)
        .y_labels(0.5, 10, None)
        .x_axis_description("x", None)
        .y_axis_description("sin(x)", None)
        .plot_grid(0.2, 0.1, class("fine-grid"))
        .plot_grid(1.0, 0.5, class("grid"))
        .plot_poly_line(sin, class("curve"))
        .plot_poly_line(vec![(0.0, 0.0), (pi2, 0.0)], class("base-line"));

    SvgDocument::new()
        .append_scss(STYLE)
        .add_svg(Box::new(chart))
        .save("docs/img/sine_curve.svg")
}
