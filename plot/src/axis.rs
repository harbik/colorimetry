use svg::{node::element::{Group, Line, Text}, Node};

#[derive(Debug, Clone, Copy, strum::Display, PartialEq, Eq)]
pub enum AxisSide {
    Left,
    Right,
    Top,
    Bottom,
}

impl AxisSide {
    pub fn geo(&self, target: [u32; 4]) -> (String, &str, &str, &str) {
        let [left, top, width, height] = target;
        let offset = Axis::MARGIN;
        match self {
            AxisSide::Bottom => (
                format!("translate({} {})", left, top + height + offset),
                "horizontal",
                "middle",
                "hanging",
            ),
            AxisSide::Top => (
                format!("translate({} {}) scale(1 -1)", left, top - offset),
                "horizontal",
                "middle",
                "hanging",
            ),
            AxisSide::Left => (
                format!("translate({} {}) rotate(-90)", left - offset, top),
                "vertical",
                "end",
                "middle",
            ),
            AxisSide::Right => (
                format!("translate({} {}) rotate(90)", left + width + offset, top),
                "vertical",
                "start",
                "middle",
            ),
        }
    }
}

pub struct Axis {
    pub(super) target: [u32; 4], // left, top, width, height
    pub(super) side: AxisSide,
    pub(super) scale: [f64; 2], // min, max, step
    pub(super) step: f64,
    pub(super) show_labels: bool, // Whether to show labels or not
    pub(super) class: Option<String>,
}

impl Axis {
    const MARGIN: u32 = 10; // Margin for axis labels
   // const HEIGHT: u32 = 30; // Default height for axis labels
    const TICK_LENGTH: u32 = 5; // Length of the tick marks
    
    pub fn new(
        target: [u32; 4], // left, top, width, height
        min_max: [f64; 2], // [min, max]
        step: f64,
        side: AxisSide,
        show_labels: bool,
        class: Option<&str> // invisible if niot 
    ) -> Self {

        Axis {
            target,
            side,
            scale: min_max,
            step,
            show_labels,
            class: class.map(|c| c.to_string()),
        }
    }

    pub fn to_group(&self) -> Group {
        let [x_min, x_max] = self.scale;

        let to_canvas = {
            let length = match self.side {
                AxisSide::Left | AxisSide::Right => self.target[3] as f64, // height
                AxisSide::Top | AxisSide::Bottom => self.target[2] as f64, // width
            };
            let left = self.target[0] as f64;
            match self.side {
                AxisSide::Left => move |x: f64| (left + (x - x_min) / (x_max - x_min) * length) as u32,
                AxisSide::Right => move |x: f64| (left + (x - x_min) / (x_max - x_min) * length) as u32,
                AxisSide::Top => move |x: f64| (left + (x - x_min) / (x_max - x_min) * length) as u32,
                AxisSide::Bottom => move |x : f64| (left as f64 + (x - x_min) / (x_max - x_min) * length) as u32
        };

        let round = |x: f64| {
            (x * 1000.0).round() / 1000.0 // Round to 3 decimal places
        };

        let mut ticks = Group::new()
            .set("class", "ticks");
        let mut labels = Group::new()
            .set("class", "labels");

        let mut x = (x_min/ self.step).ceil() * self.step; // Start from the first tick
        while x < x_max {
            let tick = Tick(x, self.step);
            ticks.append(tick.tick(self.target, to_canvas(x), self.side));
            if self.show_labels {
                labels.append(tick.label(self.target, to_canvas(x), round(x), self.side));
            }
            x += self.step;
        }

        let mut group = Group::new()
            .set("id", format!("axis-{}-{}", self.side.to_string().to_lowercase(), new_id()))
            .set("class", format!("axis-{}", self.side.to_string().to_lowercase()))
            .add(ticks);

        if let Some(class) = self.class.as_ref() {
            group = group.set("class", class.clone());

        }

        if self.show_labels {
            group.append(labels);

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


use std::fmt::{self, Display, Formatter, Result as FmtResult};

use crate::new_id;
struct Tick(f64, f64); // (value, step)

impl Tick {

    pub fn tick(&self, target: [u32;4], pos: u32, side: AxisSide) -> Line {
        let [left, top, width, height] = target;

        let (x_start, y_start) = match side {
            AxisSide::Bottom => (pos, top + height),
            AxisSide::Top => (pos, top),
            AxisSide::Left => (left, pos),
            AxisSide::Right => (left + width, pos),
        };  

        let (x_end, y_end) = match side {
            AxisSide::Bottom => (pos, top + height + Axis::TICK_LENGTH),
            AxisSide::Top => (pos, top - Axis::TICK_LENGTH),
            AxisSide::Left => (left - Axis::TICK_LENGTH, pos),
            AxisSide::Right => (left + width + Axis::TICK_LENGTH, pos),
        };  
        
        Line::new()
            .set("x1", x_start)
            .set("y1", y_start)
            .set("x2", x_end)
            .set("y2", y_end)
    }


    pub fn label(&self, target: [u32;4], pos:u32, value: f64, side: AxisSide) -> Text {
        let [left, top, width, height] = target;

        let (x_pos, y_pos) = match side {
            AxisSide::Bottom => (pos, top + height + Axis::TICK_LENGTH),
            AxisSide::Top => (pos, top - Axis::TICK_LENGTH),
            AxisSide::Left => (left - Axis::TICK_LENGTH, pos),
            AxisSide::Right => (left + width + Axis::TICK_LENGTH, pos),
        };  
        
        let mut text = Text::new(format!("{}", value))
            .set("x", x_pos)
            .set("y", y_pos)
            .set("text-anchor", "middle")
            .set("dominant-baseline", "text-after-edge");

        match side {
            AxisSide::Bottom => {
                text = text
                    .set("dominant-baseline", "text-before-edge")
            },
            AxisSide::Top => {
                text = text
                    .set("dominant-baseline", "text-after-edge")
                    .set("text-anchor", "middle")
                    .set("transform", format!("rotate(180, {x_pos}, {y_pos})"));
            },
            AxisSide::Left => {
                text = text
                    .set("dominant-baseline", "text-before-edge")
                    .set("transform", format!("rotate(-90, {x_pos}, {y_pos}')"));
            },
            AxisSide::Right => {
                text = text
                    .set("dominant-baseline", "text-before-edge")
                    .set("transform", format!("rotate(90, {x_pos}, {y_pos})"));
            },
        }
        text
    }
}

impl fmt::Display for Tick {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        format_axis_tick(self.0, self.1, f)
    }
}

/// Round value to the nearest tick step and format using the formatter's precision.
pub fn format_axis_tick(value: f64, step: f64, f: &mut Formatter<'_>) -> FmtResult {
    let rounded = (value / step).round() * step;

    // Determine precision from tick step if not explicitly set
    let precision = match f.precision() {
        Some(p) => p,
        None => {
            if step >= 1.0 {
                0
            } else {
                // Count decimal digits in step: e.g. 0.01 → 2
                let mut digits = 0;
                let mut s = step;
                while s < 1.0 {
                    s *= 10.0;
                    digits += 1;
                    if digits > 10 { break; } // prevent infinite loops
                }
                digits
            }
        }
    };

    write!(f, "{:.*}", precision, rounded)
}

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
