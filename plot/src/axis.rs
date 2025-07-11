use svg::{node::element::{Group, Line, Text}, Node};

#[derive(Debug, Clone, Copy, strum::Display, PartialEq, Eq)]
pub enum AxisSide {
    Left,
    Right,
    Top,
    Bottom,
}

pub struct Axis {
    pub(super) target: [u32; 4], // left, top, width, height
    pub(super) side: AxisSide,
    pub(super) scale: [f64; 2], // min, max, step
    pub(super) step: f64,
    pub(super) show_labels: bool, // Whether to show labels or not
    pub(super) class: Option<String>,
    pub(super) tick_length: u32, // Length of the ticks
    pub(super) description: Option<String>, // Description of the axis
}

impl Axis {
    
    pub fn new(
        description: Option<&str>,
        target: [u32; 4], // left, top, width, height
        min_max: [f64; 2], // [min, max]
        step: f64,
        side: AxisSide,
        tick_length: u32,
        show_labels: bool,
        class: Option<&str> // invisible if niot 
    ) -> Self {

        Axis {
            description: description.map(|d| d.to_string()),
            target,
            side,
            scale: min_max,
            step,
            show_labels,
            class: class.map(|c| c.to_string()),
            tick_length
        }
    }

    pub fn to_group(&self) -> Group {
        let [x_min, x_max] = self.scale;
        let [left, top, width, height] = self.target;
        let to_canvas: Box<dyn Fn(f64) -> f64> = {
            let length = match self.side {
                AxisSide::Left | AxisSide::Right => height as f64, // height
                AxisSide::Top | AxisSide::Bottom => width as f64, // width
            };
            match self.side {
                AxisSide::Left | AxisSide::Right => Box::new(move |x: f64| ((x_max - x ) / (x_max - x_min) * length) + top as f64),
                AxisSide::Top | AxisSide::Bottom => Box::new(move |x: f64| ((x - x_min) / (x_max - x_min) * length) + left as f64),
            }
        };

        let round = |x: f64| {
            (x * 10000.0).round() / 10000.0 // Round to 4 decimal places
        };

        let mut ticks = Group::new()
            .set("class", "ticks");
        let mut labels = Group::new()
            .set("class", "labels");

        let mut x = (x_min/ self.step).ceil() * self.step; // Start from the first tick
        while x < x_max {
            let tick = Tick(x, self.step);
            ticks.append(tick.tick(self.tick_length, self.target, to_canvas(x), self.side));
            if self.show_labels {
                labels.append(tick.label(self.tick_length, self.target, to_canvas(x), round(x), self.side));
            }
            x += self.step;
        }

        // Create the axis group with the ticks and optional labels
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

        const OFFSET: f64 = 35.0;
        if let Some(description) = self.description.as_ref() {
            let x = match self.side {
                AxisSide::Left => left as f64 - OFFSET, // Position left of the axis
                AxisSide::Right => left as f64 + width as f64 + OFFSET, // Position right of the axis
                AxisSide::Top | AxisSide::Bottom => left as f64 + width as f64 / 2.0, // Centered above the axis
            };
            let y = match self.side {
                AxisSide::Left | AxisSide::Right => top as f64 + height as f64 / 2.0, // Centered vertically
                AxisSide::Top => top as f64 - OFFSET, // Position above the axis
                AxisSide::Bottom => top as f64 + height as f64 + OFFSET, // Position below the axis
            };
            let mut text = Text::new(description.clone())
                .set("x", x)
                .set("y", y) // Position above the axis
                .set("text-anchor", "middle")
                .set("class", "axis-title");
            match self.side {
                AxisSide::Left => {
                    text = text.set("dominant-baseline", "text-after-edge")
                        .set("transform", format!("rotate(-90, {x}, {y})"));
                },
                AxisSide::Right => {
                    text = text.set("dominant-baseline", "text-after-edge")
                        .set("transform", format!("rotate(90, {x}, {y})"));
                },
                AxisSide::Top => {
                    text = text.set("dominant-baseline", "text-after-edge");
                },
                AxisSide::Bottom => {
                    text = text.set("dominant-baseline", "hanging");
                },
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


use std::fmt::{self, Display, Formatter, Result as FmtResult};

use crate::new_id;
struct Tick(f64, f64); // (value, step)

impl Tick {

    pub fn tick(&self, tick_length: u32, target: [u32;4], pos: f64, side: AxisSide) -> Line {
        let [left, top, width, height] = target.map(|x| x as f64);
        let tick_length = tick_length as f64;
        let (x_start, y_start) = match side {
            AxisSide::Bottom => (pos, top + height),
            AxisSide::Top => (pos, top),
            AxisSide::Left => (left, pos),
            AxisSide::Right => (left + width, pos),
        };  

        let (x_end, y_end) = match side {
            AxisSide::Bottom => (pos, top + height + tick_length),
            AxisSide::Top => (pos, top - tick_length),
            AxisSide::Left => (left - tick_length, pos),
            AxisSide::Right => (left + width + tick_length, pos),
        };  
        
        Line::new()
            .set("x1", x_start)
            .set("y1", y_start)
            .set("x2", x_end)
            .set("y2", y_end)
    }


    pub fn label(&self, tick_length: u32, target: [u32;4], pos:f64, value: f64, side: AxisSide) -> Text {
        let [left, top, width, height] = target.map(|x| x as f64);
        let tick_length = tick_length as f64;
        let (x_pos, y_pos) = match side {
            AxisSide::Bottom => (pos, top + height + tick_length),
            AxisSide::Top => (pos, top - tick_length),
            AxisSide::Left => (left - tick_length, pos),
            AxisSide::Right => (left + width + tick_length, pos),
        };  
        
        let mut text = Text::new(format!("{}", value))
            .set("x", x_pos)
            .set("y", y_pos)
            .set("text-anchor", "middle");

        match side {
            AxisSide::Bottom => {
                text = text
                    .set("dominant-baseline", "text-before-edge");
            },
            AxisSide::Top => {
                text = text
                    .set("dominant-baseline", "text-after-edge");
            },
            AxisSide::Left => {
                text = text
                    .set("dominant-baseline", "text-after-edge")
                    .set("transform", format!("rotate(-90, {x_pos}, {y_pos})"));
            },
            AxisSide::Right => {
                text = text
                    .set("dominant-baseline", "text-after-edge")
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
