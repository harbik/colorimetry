use colorimetry::{observer::Observer, rgb::RgbSpace::SRGB};
use colorimetry_plot::{axis::AxisSide, chromaticity::XYChromaticity, svgdoc::SvgDocument};

pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    let observer = Observer::default();

    let mut svgdoc = SvgDocument::new(1000, 1100)
        .add_css_rule(".fine-grid", "stroke: #AAA; stroke-width: 0.5;")
        .add_css_rule(".grid", "stroke: #AAA; stroke-width: 1.0;")
        .add_css_rule(".chart-area", "fill: #888; stroke: none; stroke-width: 0;")
        .add_css_rule(
            ".spectral-locus",
            "fill: #DDD; stroke: none; stroke-width: 0; stroke-linecap: round;",
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

    let mut xy_chromaticity = XYChromaticity::new(
        "cie1931_chromaticity_diagram",
        observer,
        (750, 850),
        (0.0..=0.75, 0.0..=0.85),
        (Some("chart-area"), None),
    );

    xy_chromaticity
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
        .draw_spectral_locus(Some("spectral-locus"), None)
        .draw_rgb_gamut(SRGB, None, None)
        .draw_grid(0.01, 0.01, Some("fine-grid"), None)
        .draw_grid(0.1, 0.1, Some("grid"), None);

    svgdoc.add_svg(Box::new(xy_chromaticity));

    svgdoc.save("tmp/srgb_gamut.svg").unwrap();
    Ok(())
}
