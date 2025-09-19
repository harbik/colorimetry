// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2025, Harbers Bik LLC

//! # Chromaticity XY Plot Module
//!
//! This module provides the [`XYChromaticity`] struct and related functionality for plotting chromaticity diagrams
//! in the CIE xy color space. It enables visualization of color gamuts, spectral loci, Planckian loci, and
//! white points, as well as the embedding of RGB gamut images within chromaticity plots.
//!
//! ## Features
//! - Plot CIE xy chromaticity diagrams for any observer
//! - Visualize spectral locus, Planckian locus, and standard illuminants
//! - Annotate white points and color gamuts
//! - Embed RGB gamut images as raster overlays
//! - Full access to all [`XYChart`] methods via delegation
//!
//! This module is intended for scientific visualization and color science applications.

mod gamut;
use gamut::PngImageData;

use std::ops::RangeBounds;

use crate::{
    chart::{ScaleRangeWithStep, XYChart},
    delegate_xy_chart_methods,
    rendable::Rendable,
    StyleAttr,
};
use colorimetry::{
    illuminant::{CieIlluminant, CCT},
    observer::Observer,
    rgb::RgbSpace,
    xyz::XYZ,
};
use svg::{
    node::element::{path::Data, Group, Path, Text, SVG},
    Node,
};

#[derive(Clone)]
pub struct XYChromaticity {
    pub(crate) observer: Observer,
    pub(crate) xy_chart: XYChart,
}

delegate_xy_chart_methods!(XYChromaticity, xy_chart);

impl XYChromaticity {
    pub const ANNOTATE_SEP: u32 = 2;

    pub fn new(
        width_and_height: (u32, u32),
        ranges: (impl RangeBounds<f64>, impl RangeBounds<f64>),
    ) -> XYChromaticity {
        let xy_chart = XYChart::new(width_and_height, ranges);
        let observer = Observer::default();
        XYChromaticity { observer, xy_chart }
    }
    pub fn set_observer(mut self, observer: Observer) -> Self {
        self.observer = observer;
        self
    }

    pub fn plot_spectral_locus(self, style_attr: Option<StyleAttr>) -> Self {
        let obs = self.observer;
        let locus = obs.spectral_locus();
        self.plot_shape(locus, style_attr)
    }

    /// Plots spectral locus ticks perpendicular to the spectral locus for
    /// a range for wavelengths, specficied by a start wavelength and an end wavelength,
    /// in `usize` units of nanometer, and with a `step` in nanometers.
    ///
    /// The `length` parameter specifies the length of the tick lines in pixels.
    /// A positive `length` will draw the ticks pointing inwards, towards the white point,
    /// and negative `length` will draw the ticks pointing outwards, away from the white point.
    ///
    /// Typically you will use this method multiple times to draw ticks for different ranges,
    /// and different lengths, with the finest ticks plotted first, and the coarsest ticks plotted last,
    /// plotting over the fine lines.
    ///
    /// The `style_attr` parameter allows you to specify the style of the ticks,
    /// such as stroke color, width, and other SVG style attributes.
    pub fn plot_spectral_locus_ticks(
        self,
        range: impl RangeBounds<usize>,
        step: usize,
        length: usize,
        style_attr: Option<StyleAttr>,
    ) -> Self {
        let length = length as f64;
        let this = self;
        let locus = this.observer.spectral_locus();
        let to_plot = this.xy_chart.to_plot.clone();
        let mut data = Data::new();
        for (xy, angle) in locus.iter_range_with_slope(range, step) {
            let pxy1 = to_plot(xy);
            data = data.move_to(pxy1);
            let pxy2 = (pxy1.0 + length * angle.sin(), pxy1.1 + length * angle.cos());
            data = data.line_to(pxy2);
        }
        this.draw_data("plot", data, style_attr)
    }

    /// Plots spectral locus labels for the specified range of wavelengths,
    /// with a step size in nanometers.
    ///
    /// The labels are rotated to align with the spectral locus slope at each point.
    /// The `distance` parameter specifies the distance from the spectral locus to the label in pixels.
    ///
    /// The `style_attr` parameter allows you to specify the style of the labels,
    pub fn plot_spectral_locus_labels(
        self,
        range: impl RangeBounds<usize> + Clone,
        step: usize,
        distance: usize,
        style_attr: Option<StyleAttr>,
    ) -> Self {
        let d = distance as f64;
        let mut self_as_mut = self;
        let locus = self_as_mut.observer.spectral_locus();
        let range_f64 = (
            range.start_bound().map(|&x| x as f64),
            range.end_bound().map(|&x| x as f64),
        );
        let values: ScaleRangeWithStep = (range_f64, step as f64).into();
        let to_plot = self_as_mut.xy_chart.to_plot.clone();
        let mut group = Group::new();
        for ((xy, angle), v) in locus.iter_range_with_slope(range, step).zip(values.iter()) {
            let pxy1 = to_plot(xy);
            let pxy2 = (pxy1.0 - d * angle.sin(), pxy1.1 - d * angle.cos());
            let label = format!("{v:.0}");
            let rotation_angle = -angle.to_degrees();
            let text = Text::new(label)
                .set("x", pxy2.0)
                .set("y", pxy2.1)
                .set("text-anchor", "middle")
                .set("dominant-baseline", "before-edge")
                .set(
                    "transform",
                    format!("rotate({:.3} {:.3} {:.3})", rotation_angle, pxy2.0, pxy2.1),
                );
            group.append(text);
        }
        style_attr.unwrap_or_default().assign(&mut group);
        self_as_mut
            .xy_chart
            .layers
            .get_mut("plot")
            .unwrap()
            .append(group);
        self_as_mut
    }

    pub fn plot_planckian_locus(self, style_attr: Option<StyleAttr>) -> Self {
        let locus = self.observer.planckian_locus();
        self.plot_poly_line(locus, style_attr)
    }

    /// Calculates the slope angle of the Planckian locus at a given CCT.
    /// Returns a tuple containing the chromaticity coordinates (x, y) and the slope angle in degrees.
    /// The slope angle is the angle of the tangent to the Planckian locus in the xy chromatity chart at the specified CCT,
    /// and points left with increasing CCT.
    /// Angles in SVG are measured clockwise from the positive x-axis.
    ///
    /// # Arguments
    /// * `cct` - The correlated color temperature in Kelvin.
    /// # Returns
    /// A tuple containing the chromaticity coordinates (x, y) and the slope angle in degrees.
    ///
    /// # Example
    /// ```
    /// use colorimetry_plot::chart::XYChromaticity;
    /// let xy_chromaticity = XYChromaticity::new((800, 600), (0.0..=0.75, 0.0..=0.875));
    /// let (xy, angle) = xy_chromaticity.planckian_xy_slope_angle(2300.0);
    /// approx::assert_abs_diff_eq!(angle.to_degrees(), -178.0, epsilon = 0.5);
    /// let (xy, angle) = xy_chromaticity.planckian_xy_slope_angle(6500.0);
    /// approx::assert_abs_diff_eq!(angle.to_degrees(), -136.0, epsilon = 0.5); // Check if the angle is a finite number
    /// ```
    pub fn planckian_xy_slope_angle(&self, cct: f64) -> ((f64, f64), f64) {
        // Get the XYZ coordinates and their derivatives at the given CCT
        let xyz = self.observer.xyz_planckian_locus(cct);
        let [x, y, z] = xyz.to_array();
        let nom = x + y + z;
        let [dxdt, dydt, dzdt] = self.observer.xyz_planckian_locus_slope(cct).to_array();
        let dnom_dt = dxdt + dydt + dzdt;

        // convert XYZ derivatives to xy derivatives using the quotient rule
        // omit the division by nom^2, since we are only interested in the angle
        let dx_chromaticity = dxdt * nom - x * dnom_dt; //  /(nom * nom);
        let dy_chromaticity = dydt * nom - y * dnom_dt; // / (nom * nom);
        let angle = dy_chromaticity.atan2(dx_chromaticity);
        (xyz.chromaticity().to_tuple(), angle)
    }

    /// The iso-temperature lines are defined as the lines that are perpendicular to the slope of
    /// the Planckian locus at a given CCT in the CIE 1960 (u,v) chromaticity diagram, using
    /// the CIE 1931 standard observer.
    /// # References
    /// - CIE 015:2018, "Colorimetry, 4th Edition", Section 9.4
    pub fn planckian_uvp_normal_angle(&self, cct: f64) -> ((f64, f64), f64) {
        // Get the XYZ coordinates and their derivatives at the given CCT
        let cct_observer = colorimetry::observer::Observer::Cie1931;
        let xyz = cct_observer.xyz_planckian_locus(cct);
        let [x, y, z] = xyz.to_array();
        let sigma_t = x + 15.0 * y + 3.0 * z;
        let [dxdt, dydt, dzdt] = cct_observer.xyz_planckian_locus_slope(cct).to_array();
        let dsigma_dt = dxdt + 15.0 * dydt + 3.0 * dzdt;

        // convert XYZ derivatives to xy derivatives using the quotient rule
        // omit the division by sigma_t^2, since we are only interested in the angle
        let du_dt = 4.0 * dxdt * sigma_t - 4.0 * x * dsigma_dt; // / (sigma_t * sigma_t)
        let dv_dt = 9.0 * dydt * sigma_t - 9.0 * y * dsigma_dt; // / (sigma_t * sigma_t
        let angle = dv_dt.atan2(du_dt);
        (
            xyz.chromaticity().to_tuple(),
            angle + std::f64::consts::FRAC_PI_2,
        )
    }

    /// Transform a (CCT, duv) pair to a plot point using the CIE 1931 observer.
    ///
    /// # Arguments
    /// * `cct_duv` - A tuple containing the correlated color temperature (CCT) in Kelvin and the chromaticity deviation from the Planckian locus (duv).
    ///
    /// # Returns
    /// A Result containing a tuple of the transformed chromaticity coordinates (x, y) in the plot space.
    /// # Errors
    /// Returns an error if the CCT is invalid or cannot be converted to XYZ coordinates.
    pub fn cct_transform(&self, cct_duv: (f64, f64)) -> Result<(f64, f64), colorimetry::Error> {
        let (cct, duv) = cct_duv;
        let xyz: XYZ = CCT::new(cct, duv)?.try_into()?;
        let xy = xyz.chromaticity().to_tuple();
        Ok((self.xy_chart.to_plot)(xy))
    }

    /// Plots the ANSI C78.377-2017 step 7 Quadrangles in the plot,
    /// filled with the color corresponding to the CCT and duv target.
    pub fn plot_ansi_step7(mut self, rgb_space: RgbSpace, style_attr: Option<StyleAttr>) -> Self {
        // ANSI C78.377-2017, Table 1, Basic Nominal CCT Specification.
        const DATA: &[(&str, (i32, i32), f64)] = &[
            ("2200", (2238, 102), 0.0000),
            ("2500", (2460, 120), 0.0000),
            ("2700", (2725, 145), 0.0000),
            ("3000", (3045, 175), 0.0001),
            ("3500", (3465, 245), 0.0005),
            ("4000", (3985, 275), 0.0010),
            ("4500", (4503, 243), 0.0015),
            ("5000", (5029, 283), 0.0020),
            ("5700", (5667, 355), 0.0025),
            ("6500", (6532, 510), 0.0031),
        ];

        let duv = |cct: i32| {
            if cct < 2780 {
                0.0
            } else {
                let r = 1.0 / cct as f64;
                57_700.0 * r * r - 44.6 * r + 0.00854
            }
        };
        const DUV_TOLERANCE: f64 = 0.006; // tolerance for duv

        let mut ansi = Group::new();
        style_attr.unwrap_or_default().assign(&mut ansi);
        for &(_, (cct, tol), duv_target) in DATA {
            let mut data = Data::new();

            let (px0, py0) = self
                .cct_transform(((cct - tol) as f64, duv(cct - tol) - DUV_TOLERANCE))
                .unwrap();
            let (px1, py1) = self
                .cct_transform(((cct - tol) as f64, duv(cct - tol) + DUV_TOLERANCE))
                .unwrap();
            let (px2, py2) = self
                .cct_transform(((cct + tol) as f64, duv(cct + tol) + DUV_TOLERANCE))
                .unwrap();
            let (px3, py3) = self
                .cct_transform(((cct + tol) as f64, duv(cct + tol) - DUV_TOLERANCE))
                .unwrap();

            data = data
                .move_to((px0, py0))
                .line_to((px1, py1))
                .line_to((px2, py2))
                .line_to((px3, py3))
                .close();
            let xyz: XYZ = CCT::new(cct as f64, duv_target)
                .unwrap()
                .try_into()
                .unwrap();
            let [r, g, b]: [u8; 3] = xyz.rgb(rgb_space).compress().into();
            let path = Path::new()
                .set("d", data.clone())
                .set("style", format!("fill: rgb({r:.0}, {g:.0}, {b:.0})"));
            ansi.append(path);
        }
        self.xy_chart.layers.get_mut("plot").unwrap().append(ansi);
        self
    }

    /// Get the normal to the Planckian locus, which is perpendicular to the slope at the given CCT.
    /// Returns a tuple containing the chromaticity coordinates (x, y) and the normal angle in radians.
    pub fn planckian_xy_normal_angle(&self, cct: f64) -> ((f64, f64), f64) {
        let (xy, slope_angle) = self.planckian_xy_slope_angle(cct);
        (xy, slope_angle + std::f64::consts::FRAC_PI_2)
    }

    /// Plots the Planckian locus ticks at specified CCT values.
    /// Typically, you will use this method several times to create a nice looking scale,
    /// with different ranges and lengths.
    /// The scale if very non-linear, and requires small steps for the lower CCT values,
    /// and larger steps for the higher CCT values.
    pub fn plot_planckian_locus_ticks(
        self,
        values: impl IntoIterator<Item = u32>,
        length: usize,
        style_attr: Option<StyleAttr>,
    ) -> Self {
        let mut data = Data::new();
        let to_plot = self.xy_chart.to_plot.clone();
        for cct in values {
            let (xy, angle) = self.planckian_xy_normal_angle(cct as f64);
            let (px, py) = to_plot(xy);
            let pdx = length as f64 * angle.cos();
            let pdy = length as f64 * angle.sin();
            data = data
                .move_to((px - pdx, py + pdy)) // top in plot
                .line_to((px + pdx, py - pdy)); // bottom
        }
        self.draw_data("plot", data, style_attr)
    }

    /// Plots labels for the Planckian locus at specified CCT values, divided by 100.
    /// The labels are rotated to align with the normal to the the Planckian locus at each point,
    /// centered away from the white point.
    pub fn plot_planckian_locus_labels(
        mut self,
        values: impl IntoIterator<Item = u32>,
        distance: usize,
        style_attr: Option<StyleAttr>,
    ) -> Self {
        let mut planckian_labels = Group::new();
        let to_plot = self.xy_chart.to_plot.clone();
        for cct in values {
            let (xy, angle) = self.planckian_xy_normal_angle(cct as f64);
            let (px, py) = to_plot(xy);
            let pdx = distance as f64 * angle.cos();
            let pdy = distance as f64 * angle.sin();

            let px2 = px - pdx;
            let py2 = py + pdy;
            let text = Text::new(format!("{}", cct / 100))
                .set("x", px2)
                .set("y", py2)
                .set("text-anchor", "end")
                .set("dominant-baseline", "middle")
                .set(
                    "transform",
                    format!("rotate({:.3} {px2:.3} {py2:.3}) ", -angle.to_degrees()),
                );
            planckian_labels.append(text);
        }
        style_attr.unwrap_or_default().assign(&mut planckian_labels);
        self.xy_chart
            .layers
            .get_mut("plot")
            .unwrap()
            .append(planckian_labels);
        self
    }

    /// Plots the sRGB gamut as an embeded image overlay in the plot layer.
    /// At this point, only correct colors are shown in case the RGB space is sRGB.
    /// For other color spaces, the colors are clipped to the sRGB gamut.
    pub fn plot_rgb_gamut(self, rgb_space: RgbSpace, style_attr: Option<StyleAttr>) -> Self {
        // **TODO**
        // Include the color profiles in the embeded image data, for requested colors space.  This will allow
        // the colors to be displayed correctly for displays that have a wide enough gamut, and
        // browser which support SVG embeded PNG images with color profiles (such as recent versions
        // of Chrome, Firefox, and Safari).
        let gamut_fill = PngImageData::from_rgb_space(
            self.observer,
            rgb_space,
            self.xy_chart.to_plot.clone(),
            self.xy_chart.to_world.clone(),
        );
        self.plot_image(gamut_fill, style_attr)
    }

    /// Draw white points on the chromaticity diagram as an iterator of CieIlluminant, and i32 angle and length pairs.
    #[allow(unused)]
    pub fn annotate_white_points(
        &mut self,
        point: impl IntoIterator<Item = (CieIlluminant, (i32, i32))>,
    ) -> &mut Self {
        todo!()
    }
}

/// Implements the XYChromaticity as a Rendable object, allowing it to be rendered as an SVG.
impl Rendable for XYChromaticity {
    fn render(&self) -> SVG {
        self.xy_chart.render()
    }

    fn view_parameters(&self) -> crate::view::ViewParameters {
        self.xy_chart.view_parameters()
    }

    fn set_view_parameters(&mut self, view_box: crate::view::ViewParameters) {
        self.xy_chart.set_view_parameters(view_box);
    }
}
