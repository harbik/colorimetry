#![allow(unused)]
mod gamut;

use std::ops::{Deref, DerefMut, RangeBounds};

use crate::{
    axis::AxisSide,
    chart::XYChart,
    chromaticity::xy::{self, gamut::PngImageData},
    delegate_xy_chart_methods,
    rendable::Rendable,
    svgdoc::SvgDocument,
};
use colorimetry::{
    observer::{self, Observer},
    prelude::CieIlluminant,
    rgb::RgbSpace,
};
use svg::node::element::{Image, SVG};

#[derive(Clone)]
pub struct XYChromaticity {
    pub(crate) observer: Observer,
    pub(crate) xy_chart: XYChart,
}

delegate_xy_chart_methods!(XYChromaticity, xy_chart);

impl XYChromaticity {
    pub const ANNOTATE_SEP: u32 = 2;

    pub fn new(
        id: impl AsRef<str>,
        observer: Observer,
        width_and_height: (u32, u32),
        ranges: (impl RangeBounds<f64>, impl RangeBounds<f64>),
        class_and_style: (Option<&str>, Option<&str>),
    ) -> XYChromaticity {
        let (class, style) = class_and_style;
        let xy_chart = XYChart::new(id.as_ref(), width_and_height, ranges, class_and_style);

        XYChromaticity { observer, xy_chart }
    }

    pub fn draw_spectral_locus(mut self, class: Option<&str>, style: Option<&str>) -> Self {
        let locus = self.observer.spectral_locus();
        self.xy_chart = self.xy_chart.draw_shape(locus, class, style);
        self
    }
    
    pub fn draw_planckian_locus(mut self, class: Option<&str>, style: Option<&str>) -> Self {
        let locus = self.observer.planckian_locus();
        self.xy_chart = self.xy_chart.draw_line(locus, class, style);
        self
    }

    /// Draw white points on the chromaticity diagram as an iterator of CieIlluminant, and i32 angle and length pairs.
    pub fn annotate_white_points(
        &mut self,
        point: impl IntoIterator<Item = (CieIlluminant, (i32, i32))>,
    ) -> &mut Self {
        todo!()
    }

    pub fn draw_rgb_gamut(
        mut self,
        rgb_space: RgbSpace,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let gamut_fill = PngImageData::from_rgb_space(
            self.observer,
            rgb_space,
            self.xy_chart.to_plot.clone(),
            self.xy_chart.to_world.clone(),
        );
        self = self.draw_image(gamut_fill, class, style);
        self
    }
}

/*
/// Implements Deref and DerefMut traits for XYChromaticity, allowing it to use the methods from XYChart.
impl Deref for XYChromaticity {
    type Target = XYChart;

    fn deref(&self) -> &Self::Target {
        &self.xy_chart
    }
}

impl DerefMut for XYChromaticity {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.xy_chart
    }
}
 */

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
