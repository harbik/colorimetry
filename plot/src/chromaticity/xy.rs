#![allow(unused)]
use std::ops::RangeBounds;

use crate::{axis::AxisSide, chart::XYChart};
use colorimetry::{
    observer::{self, Observer},
    prelude::CieIlluminant,
    rgb::RgbSpace,
};

#[derive(Clone)]
pub struct XYChromaticity {
    pub(crate) observer: Observer,
    pub(crate) chart: XYChart,
}

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
        XYChromaticity {
            observer,
            chart: xy_chart
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
                .draw_grid(0.01, 0.01, Some("fine-grid"), None)
                .draw_grid(0.1, 0.1, Some("grid"), None),
        }
    }

    pub fn draw_spectral_locus(mut self, class: Option<&str>, style: Option<&str>) -> Self {
        let locus = self.observer.spectral_locus();
        self.chart = self.chart.draw_shape(locus, class, style);
        self
    }
    pub fn draw_planckian_locus(mut self, class: Option<&str>, style: Option<&str>) -> Self {
        let locus = self.observer.planckian_locus();
        self.chart = self.chart.draw_line(locus, class, style);
        self
    }

    /// Draw white points on the chromaticity diagram as an iterator of CieIlluminant, and i32 angle and length pairs.
    pub fn annotate_white_points(
        mut self,
        point: impl IntoIterator<Item = (CieIlluminant, (i32, i32))>,
    ) -> Self {
        todo!()
    }

    pub fn draw_rgb_gamut(
        mut self,
        rgb_space: RgbSpace,
        rgb_fill: bool,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        todo!()
    }
}
