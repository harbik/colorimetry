use svg::{
    node::element::{path::Data, Group, Path, Text},
    Node,
};

pub mod tick;
pub use tick::Tick;
mod range;
pub use range::{ScaleRange, ScaleRangeWithStep};

#[derive(Debug, Clone, Copy, strum::Display, PartialEq, Eq)]
pub enum AxisSide {
    Left,
    Right,
    Top,
    Bottom,
}

pub struct Axis {
    pub(super) target: (i32, i32, u32, u32), // left, top, width, height
    pub(super) side: AxisSide,
    pub(super) range_with_step: ScaleRangeWithStep, // Range of the axis
    pub(super) show_labels: bool,                   // Whether to show labels or not
    pub(super) class: Option<String>,
    pub(super) tick_length: i32,            // Length of the ticks
    pub(super) description: Option<String>, // Description of the axis
}

impl Axis {
    pub fn new(
        description: Option<&str>,
        target: (i32, i32, u32, u32), // left, top, width, height
        range_with_step: ScaleRangeWithStep,
        side: AxisSide,
        tick_length: i32,
        show_labels: bool,
        class: Option<&str>, // invisible if not set
    ) -> Self {
        Axis {
            description: description.map(|d| d.to_string()),
            target,
            side,
            range_with_step,
            show_labels,
            class: class.map(|c| c.to_string()),
            tick_length,
        }
    }

    pub fn to_group(&self) -> Group {
        let (left, top, width, height) = self.target;
        let range = self.range_with_step.range;
        let step = self.range_with_step.step;
        // This closure converts the axis value to canvas coordinates based on the axis side
        let to_canvas: Box<dyn Fn(f64) -> f64> = {
            let length = match self.side {
                AxisSide::Left | AxisSide::Right => height as f64, // height
                AxisSide::Top | AxisSide::Bottom => width as f64,  // width
            };
            match self.side {
                AxisSide::Left | AxisSide::Right => {
                    Box::new(move |x: f64| (range.scale_descent(x) * length) + top as f64)
                }
                AxisSide::Top | AxisSide::Bottom => {
                    Box::new(move |y: f64| (range.scale(y) * length) + left as f64)
                }
            }
        };

        let mut ticks = Group::new().set("class", "ticks");
        let mut labels = Group::new().set("class", "labels");

        let mut data = Data::new();
        for x in self.range_with_step.iter() {
            let tick = Tick(x, step);
            let (start, end) = tick.tick(self.tick_length, self.target, to_canvas(x), self.side);
            data = data.move_to(start).line_to(end);
            if self.show_labels {
                labels.append(tick.label(
                    self.tick_length,
                    self.target,
                    to_canvas(x),
                    round_to_precision(x, 4),
                    self.side,
                ));
            }
        }
        let path = Path::new().set("d", data.clone());
        ticks.append(path);

        // Create the axis group with the ticks and optional labels
        let mut group = Group::new()
            .set(
                "id",
                format!("axis-{}-{}", self.side.to_string().to_lowercase(), new_id()),
            )
            .set(
                "class",
                format!("axis-{}", self.side.to_string().to_lowercase()),
            )
            .add(ticks);

        if let Some(class) = self.class.as_ref() {
            group = group.set("class", class.clone());
        }

        if self.show_labels {
            group.append(labels);
        }

        if let Some(description) = self.description.as_ref() {
            let x = match self.side {
                AxisSide::Left => left as f64 - XYChart::DESCRIPTION_OFFSET as f64, // Position left of the axis
                AxisSide::Right => left as f64 + width as f64 + XYChart::DESCRIPTION_OFFSET as f64, // Position right of the axis
                AxisSide::Top | AxisSide::Bottom => left as f64 + width as f64 / 2.0, // Centered above the axis
            };
            let y = match self.side {
                AxisSide::Left | AxisSide::Right => top as f64 + height as f64 / 2.0, // Centered vertically
                AxisSide::Top => top as f64 - XYChart::DESCRIPTION_OFFSET as f64, // Position above the axis
                AxisSide::Bottom => top as f64 + height as f64 + XYChart::DESCRIPTION_OFFSET as f64, // Position below the axis
            };
            let mut text = Text::new(description.clone())
                .set("x", x)
                .set("y", y) // Position above the axis
                .set("text-anchor", "middle")
                .set("class", "axis-title");
            match self.side {
                AxisSide::Left => {
                    text = text
                        .set("dominant-baseline", "text-after-edge")
                        .set("transform", format!("rotate(-90, {x}, {y})"));
                }
                AxisSide::Right => {
                    text = text
                        .set("dominant-baseline", "text-after-edge")
                        .set("transform", format!("rotate(90, {x}, {y})"));
                }
                AxisSide::Top => {
                    text = text.set("dominant-baseline", "text-after-edge");
                }
                AxisSide::Bottom => {
                    text = text.set("dominant-baseline", "text-before-edge");
                }
            }

            group.append(text);
        }
        group
    }
}

impl From<Axis> for Group {
    fn from(axis: Axis) -> Self {
        axis.to_group()
    }
}

impl Display for Axis {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "{}", Axis::to_group(self))
    }
}

use std::fmt::{Display, Formatter, Result as FmtResult};

use crate::{chart::XYChart, new_id, round_to_precision};
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_axis_tick() {
        fn format_tick(value: f64, step: f64) -> String {
            format!("{}", Tick(value, step))
        }
        assert_eq!(format_tick(1.2345, 0.1), "1.2");
        assert_eq!(format_tick(1.2345, 0.01), "1.23");
        assert_eq!(format_tick(1.2345, 0.001), "1.235");
        assert_eq!(format_tick(1.2345, 1.0), "1");
        assert_eq!(format_tick(1.0, 0.001), "1.000");
    }
}
