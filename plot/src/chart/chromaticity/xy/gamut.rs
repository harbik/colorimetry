use std::fmt::Display;
use std::fs::File;

use base64::engine::general_purpose;
use base64::engine::Engine;
use colorimetry::rgb::Rgb;
use colorimetry::rgb::WideRgb;
use colorimetry::xyz::XYZ;
use colorimetry::{math::Triangle, observer::Observer, rgb::RgbSpace};
use image::ImageEncoder;
use image::{codecs::png::PngEncoder, Rgba, RgbaImage};

use crate::chart::CoordinateTransform;
use svg::node::element::Image;

pub struct PngImageData {
    png: String,
    position: (i32, i32),
    dimensions: (u32, u32),
}

impl PngImageData {
    pub fn new(png: String, position: (i32, i32), dimensions: (u32, u32)) -> Self {
        Self {
            png,
            position,
            dimensions,
        }
    }

    /// Creates a new PngImageData from an RGB space and observer.
    /// # Arguments
    /// * `observer` - The observer for which the RGB space is defined.
    /// * `space` - The RGB space to generate the PNG for.
    /// * `on_canvas` - A function that transforms chromaticity coordinates to pixel coordinates on the canvas.
    /// * `position` - The position of the image on the canvas.
    /// * `dimensions` - The dimensions of the image.
    pub fn from_rgb_space(
        observer: Observer,
        space: RgbSpace,
        to_plot: CoordinateTransform,
        to_world: CoordinateTransform,
    ) -> Self {
        png_from_rgb_space(observer, space, to_plot, to_world)
    }

    pub fn png(&self) -> &str {
        &self.png
    }

    pub fn position(&self) -> (i32, i32) {
        self.position
    }

    pub fn dimensions(&self) -> (u32, u32) {
        self.dimensions
    }
}

/// Converts PngImageData to an SVG Image element.
impl From<PngImageData> for Image {
    fn from(data: PngImageData) -> Self {
        Image::new()
            .set("x", data.position.0)
            .set("y", data.position.1)
            .set("width", data.dimensions.0)
            .set("height", data.dimensions.1)
            .set("href", format!("data:image/png;base64,{}", data.png))
    }
}

/// Generates a PNG image of the RGB gamut for a given observer and RGB space,
fn png_from_rgb_space(
    observer: Observer,
    space: RgbSpace,
    to_plot: CoordinateTransform,
    to_world: CoordinateTransform,
) -> PngImageData {
    // get chromaticity coordinate of the primaries
    let chromaticities = space.chromaticities(observer);

    // calculate the pixel coordinates of the triangle vertices, using the CoordinateTransform function
    let gamut_plot_coordinates = chromaticities
        .iter()
        .map(|c| c.to_tuple())
        .map(|xy| to_plot(xy))
        .map(|(x, y)| [x as u32, y as u32])
        .collect::<Vec<_>>();

    let mut v_min = u32::MAX;
    let mut h_min = u32::MAX;
    let mut v_max = 0;
    let mut h_max = 0;

    // Find the minimum and maximum coordinates to determine the image size
    for &[h, v] in gamut_plot_coordinates.iter() {
        if h < h_min {
            h_min = h;
        }
        if h > h_max {
            h_max = h;
        }
        if v < v_min {
            v_min = v;
        }
        if v > v_max {
            v_max = v;
        }
    }
    let width = h_max - h_min;
    let height = v_max - v_min;
    let mut image = RgbaImage::new(width, height);

    // Create triangle from the plot coordinates, to only calculate the pixels inside the triangle
    // The coordinates are in the form [x, y] where x and y are the pixel coordinates
    // Triange needs the coordinates as f64
    let gamut_plot_triangle = Triangle::new(
        gamut_plot_coordinates[0].map(|c| c as f64),
        gamut_plot_coordinates[1].map(|c| c as f64),
        gamut_plot_coordinates[2].map(|c| c as f64),
    )
    .unwrap();

    // Fill the image with the triangle gradient
    for v in 0..height { // vertical pixel index in the png image
        for h in 0..width { // horizontal pixel index in the png image
            if gamut_plot_triangle.contains((h + h_min) as f64, (v + v_min) as f64) {
                let (x, y) = to_world((h as f64 + h_min as f64, v as f64 + v_min as f64));
                let xyz = XYZ::new([x, y, 1.0 - x - y], observer).set_illuminance(100.0);
                let wrgb = xyz.rgb(space);
                let [r, g, b] = wrgb.compress().into();
                image.put_pixel(h, v, Rgba([r, g, b, 255]));
            } else {
                image.put_pixel(h, v, Rgba([0, 0, 0, 0])); // fully transparent pixel
            }
        }
    }

    // Save PNG to memory
    let mut png_data = Vec::new();
    PngEncoder::new(&mut png_data)
        .write_image(
            image.as_raw(),
            image.width(),
            image.height(),
            image::ExtendedColorType::Rgba8,
        )
        .unwrap();

    // Encode PNG as base64
    let base64_png = general_purpose::STANDARD.encode(png_data);

    // Create PngImageData with the base64 string and dimensions
    PngImageData::new(base64_png, (h_min as i32, v_min as i32), (width, height))
}
