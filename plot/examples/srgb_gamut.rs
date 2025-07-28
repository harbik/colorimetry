use colorimetry::{observer::Observer, rgb::RgbSpace};
use colorimetry_plot::{chart::XYChromaticity, style_attr::class, svgdoc::SvgDocument};

/// Includes the style for the SVG document from an external SCSS file.
/// This is a SCSS stylesheet that styles the sRGB gamut plot and is embedded into the SVG output.
const STYLE: &str = include_str!("srgb_gamut.scss");
const PLANCKIAN_LABELS_AT: &[u32] = &[
    2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6500, 7500, 9300,
];

pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    let observer = Observer::Cie2015_10; // Use the CIE 1931 observer
    let space = RgbSpace::DisplayP3; // Use the sRGB color space

    // Create an XYChromaticity chart with the specified observer and ranges
    let xy_chromaticity = XYChromaticity::new((775, 875), (-0.025..=0.75, 0.0..=0.875))
        .ticks(0.01, 0.01, 5, class("fine-grid"))
        .ticks(0.1, 0.1, 10, class("grid"))
        .x_labels(0.1, 10, None)
        .x_axis_description(&format!("{} x", observer.name()), None)
        .y_labels(0.1, 10, class("y-labels"))
        .y_axis_description(&format!("{} y", observer.name()), None)
        .plot_spectral_locus(class("spectral-locus"))
        .plot_spectral_locus_ticks(440..=650, 10, 15, class("spectral-locus-ticks"))
        .plot_spectral_locus_ticks(460..=630, 1, 7, class("spectral-locus-ticks"))
        .plot_spectral_locus_labels(460..=620, 10, 3, class("spectral-locus-labels"))
        .plot_rgb_gamut(space, None)
        .plot_planckian_locus(class("planckian"))
        .plot_planckian_locus_ticks((1500..=7500).step_by(100), 7, class("planckian-ticks-fine"))
        .plot_planckian_locus_ticks(PLANCKIAN_LABELS_AT.to_vec(), 15, class("planckian-ticks"))
        .plot_planckian_locus_labels(PLANCKIAN_LABELS_AT.to_vec(), 18, class("planckian-labels"))
        .plot_grid(0.01, 0.01, class("fine-grid"))
        .plot_grid(0.1, 0.1, class("grid"))
    ;

    // create the plot
    SvgDocument::new()
        .append_scss(STYLE)
        .add_svg(Box::new(xy_chromaticity))
        .save("docs/img/srgb_gamut.svg")
}
