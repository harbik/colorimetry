use colorimetry::rgb::RgbSpace::SRGB;
use colorimetry_plot::{style_attr::class, chart::XYChromaticity, svgdoc::SvgDocument,
};

/// Includes the style for the SVG document from an external SCSS file.
/// This is a SCSS stylesheet that styles the sRGB gamut plot and is embedded into the SVG output.
const STYLE: &str = include_str!("srgb_gamut.scss");

pub fn main() -> Result<(), Box<dyn std::error::Error>> {

    // Create an XYChromaticity chart with the specified observer and ranges
    let xy_chromaticity = XYChromaticity::new((775, 875), (-0.025..=0.75, 0.0..=0.875))
        .ticks(0.01, 0.01, 5, class("fine-grid"))
        .ticks(0.1, 0.1, 10, class("grid"))
        .x_labels(0.1, 10, None)
        .x_axis_description("CIE 1931 x Chromaticity",None)
        .y_labels(0.1, 10, class("y-labels"))
        .y_axis_description("CIE 1931 y Chromaticity",None) 
        .plot_spectral_locus(class("spectral-locus"))
        .plot_spectral_locus_ticks(
            440..=650,
            10,
            15,
            class("spectral-locus-ticks"),
        )
        .plot_spectral_locus_ticks(460..=630, 1, 7, class("spectral-locus-ticks"))
        .plot_spectral_locus_labels(460..=620, 10, 2, class("spectral-locus-labels"))
        .plot_rgb_gamut(SRGB, None)
        .plot_planckian_locus(class("planckian"))
        .plot_grid(0.01, 0.01, class("fine-grid"))
        .plot_grid(0.1, 0.1, class("grid"));

    // create the plot
    SvgDocument::new()
        .append_scss(STYLE)
        .add_svg(Box::new(xy_chromaticity))
        .save("docs/img/srgb_gamut.svg")
}
