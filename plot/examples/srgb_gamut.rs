use colorimetry::{observer::Observer, rgb::RgbSpace::SRGB};
use colorimetry_plot::{chart::XYChromaticity, rendable::Rendable, style_attr, svgdoc::SvgDocument};

/// Includes the style for the SVG document from an external SCSS file.
/// This is a SCSS stylesheet that styles the sRGB gamut plot and is embedded into the SVG output.
const STYLE: &str = include_str!("srgb_gamut.scss");

pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    let observer = Observer::default();

    // Create an XYChromaticity chart with the specified observer and ranges
    let xy_chromaticity = XYChromaticity::new(
        observer,
        (775, 875),
        (-0.025..=0.75, 0.0..=0.875),
    )
    .ticks(0.01, 0.01, 5, style_attr!(class: "fine-grid"))
    .ticks(0.1, 0.1, 10, style_attr!(class: "grid"))
    .x_labels(0.1, 10)
    .x_axis_description("CIE 1931 x Chromaticity")
    .y_labels(0.1, 10)
    .y_axis_description("CIE 1931 y Chromaticity")
    .plot_spectral_locus(style_attr!(class: "spectral-locus"))
    .plot_spectral_locus_ticks(
        440..=650,
        10,
        15,
        style_attr!(class: "spectral-locus-ticks"),
    )
    .plot_spectral_locus_ticks(460..=630, 1, 7, style_attr!(class:"spectral-locus-ticks"))
    .plot_spectral_locus_labels(460..=620, 10, 2, style_attr!(class:"spectral-locus-labels"))
    .plot_rgb_gamut(SRGB, style_attr!())
    .plot_planckian_locus(style_attr!(class:"planckian"))
    .plot_grid(0.01, 0.01, style_attr!(class: "fine-grid"))
    .plot_grid(0.1, 0.1, style_attr!(class: "grid"));

    let margin = 50;
    let width = xy_chromaticity.width() + margin;
    let height = xy_chromaticity.height() + margin;

    // crate the plot
    SvgDocument::new(width, height, STYLE)
        .add_svg(Box::new(xy_chromaticity))
        .save("tmp/srgb_gamut.svg")
}
