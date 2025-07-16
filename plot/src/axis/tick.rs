use std::fmt::{self, Formatter, Result as FmtResult};

use svg::node::element::Text;

use crate::axis::AxisSide;


pub struct Tick(pub(crate) f64, pub(crate) f64); // (value, step)

impl Tick {

    pub fn tick(&self, tick_length: i32, target: (i32, i32, u32, u32), pos: f64, side: AxisSide) -> ((f64,f64), (f64,f64)) {
        let left = target.0 as f64;
        let top = target.1 as f64;
        let width = target.2 as f64;
        let height = target.3 as f64;
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
        
        ((x_start, y_start), (x_end, y_end))
    }


    pub fn label(&self, tick_length: i32, target: (i32, i32, u32, u32), pos:f64, value: f64, side: AxisSide) -> Text {
        let left = target.0 as f64;
        let top = target.1 as f64;
        let width = target.2 as f64;
        let height = target.3 as f64;
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
