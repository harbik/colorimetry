#![allow(unused)]
mod gamut;
use gamut::PngImageData;

use std::{
    ops::Range,
    ops::{Deref, DerefMut, RangeBounds},
};

use crate::{
    chart::{ScaleRange, ScaleRangeWithStep, XYChart},
    delegate_xy_chart_methods,
    rendable::Rendable,
    svgdoc::SvgDocument,
    StyleAttr,
};
use colorimetry::{
    observer::{self, Observer},
    prelude::CieIlluminant,
    rgb::RgbSpace,
};
use svg::{
    node::element::{path::Data, Group, Image, Line, Text, SVG},
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
    pub fn set_observer(mut self, observer: Observer) ->  Self {
        self.observer = observer;
        self
    }

    pub fn plot_spectral_locus(self, style_attr: Option<StyleAttr>) -> Self {
        let obs = self.observer;
        let locus = obs.spectral_locus();
        self.plot_shape(locus, style_attr)
    }

    pub fn plot_spectral_locus_ticks(
        self,
        range: impl RangeBounds<usize>,
        step: usize,
        length: usize,
        style_attr: Option<StyleAttr>,
    ) -> Self {
        let length = length as f64;
        let mut self_as_mut = self;
        let locus = self_as_mut.observer.spectral_locus();
        let to_plot = self_as_mut.xy_chart.to_plot.clone();
        let mut data = Data::new();
        for (xy, angle) in locus.iter_range_with_slope(range, step) {
            let pxy1 = to_plot(xy);
            data = data.move_to(pxy1);
            let pxy2 = (pxy1.0 + length * angle.sin(), pxy1.1 + length * angle.cos());
            data = data.line_to(pxy2);
        }
        self_as_mut.draw_data("plot", data, style_attr)
    }

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
            let mut text = Text::new(label)
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

    pub fn plot_rgb_gamut(self, rgb_space: RgbSpace, style_attr: Option<StyleAttr>) -> Self {
        let gamut_fill = PngImageData::from_rgb_space(
            self.observer,
            rgb_space,
            self.xy_chart.to_plot.clone(),
            self.xy_chart.to_world.clone(),
        );
        self.plot_image(gamut_fill, style_attr)
    }

    /// Draw white points on the chromaticity diagram as an iterator of CieIlluminant, and i32 angle and length pairs.
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
