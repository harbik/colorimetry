// Copyright (c) 2025, Harbers Bik LLC
// SPDX-License-Identifier: Apache-2.0 OR MIT
//
use std::fmt::{self, Formatter, Result as FmtResult};
use std::ops::{Bound, Range, RangeBounds};

pub struct ScaleValue(pub f64, pub f64);

impl fmt::Display for ScaleValue {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        format_axis_tick(self.0, self.1, f)
    }
}

/// Round value to the nearest tick step and format using the formatter's precision.
pub fn format_axis_tick(value: f64, step: f64, f: &mut Formatter<'_>) -> FmtResult {
    let rounded = (value / step).round() * step;
    if approx::abs_diff_eq!(rounded, 0.0, epsilon = 1e-10) {
        // Handle zero case explicitly to avoid formatting issues
        return write!(f, "0");
    }
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
                    if digits > 10 {
                        break;
                    } // prevent infinite loops
                }
                digits
            }
        }
    };

    write!(f, "{rounded:.precision$}")
}

#[derive(Debug, Clone, Copy)]
pub struct ScaleRange {
    pub start: f64,
    pub end: f64,
    pub end_included: bool,
}

impl ScaleRange {
    pub fn new<R: RangeBounds<f64>>(range: R) -> Self {
        let start = match range.start_bound() {
            Bound::Included(&s) => s,
            Bound::Excluded(&s) => s,
            Bound::Unbounded => 0.0,
        };
        let (end, end_included) = match range.end_bound() {
            Bound::Included(&e) => (e, true),
            Bound::Excluded(&e) => (e, false),
            Bound::Unbounded => panic!("Unbounded end not allowed"),
        };

        Self {
            start,
            end,
            end_included,
        }
    }

    pub fn as_range(&self) -> Range<f64> {
        self.start..self.end
    }

    pub fn span(&self) -> f64 {
        self.end - self.start
    }

    pub fn scale(&self, value: f64) -> f64 {
        (value - self.start) / self.span()
    }

    pub fn unscale(&self, value: f64) -> f64 {
        self.start + value * self.span()
    }

    pub fn scale_descent(&self, value: f64) -> f64 {
        (self.end - value) / self.span()
    }

    pub fn unscale_descent(&self, value: f64) -> f64 {
        self.end - value * self.span()
    }

    /// Creates an iterator over the range with a specified step size.
    /// Produces values aligned with the step size, for use to draw grid lines and ticks.
    pub fn iter_with_step(&self, step: f64) -> ScaleRangeIterator {
        let start = (self.start / step).ceil() * step; // Start from the first tick

        ScaleRangeIterator {
            start,
            end: self.end,
            step,
            end_included: self.end_included,
        }
    }
}

pub struct ScaleRangeIterator {
    start: f64,
    end: f64,
    end_included: bool,
    step: f64,
}

impl Iterator for ScaleRangeIterator {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start < self.end
            || (self.end_included && approx::abs_diff_eq!(self.start, self.end, epsilon = 1e-10))
        {
            let current = self.start;
            self.start += self.step;
            Some(current)
        } else {
            None
        }
    }
}

pub struct ScaleRangeWithStep {
    pub range: ScaleRange,
    pub step: f64,
}

impl ScaleRangeWithStep {
    pub fn new(range: ScaleRange, step: f64) -> Self {
        Self { range, step }
    }

    pub fn iter(&self) -> ScaleRangeIterator {
        self.range.iter_with_step(self.step)
    }
}

impl<R: RangeBounds<f64>> From<(R, f64)> for ScaleRangeWithStep {
    fn from((range, step): (R, f64)) -> Self {
        Self {
            range: ScaleRange::new(range),
            step,
        }
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn test_chart_range() {
        let range = ScaleRange::new(0.0..1.0);
        assert_abs_diff_eq!(range.end, 1.0);
        assert_abs_diff_eq!(range.span(), 1.0);
        assert_abs_diff_eq!(range.scale(0.5), 0.5);
        assert_abs_diff_eq!(range.scale_descent(0.5), 0.5);
    }

    #[test]
    fn test_iter_with_step() {
        let range = ScaleRange::new(0.0..=1.0);
        let mut iter = range.iter_with_step(0.2);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.0);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.2);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.4);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.6);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.8);
        assert_abs_diff_eq!(iter.next().unwrap(), 1.0);
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_iter_with_step_excluding_end() {
        let range = ScaleRange::new(0.0..1.0);
        let mut iter = range.iter_with_step(0.2);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.0);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.2);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.4);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.6);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.8);
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_iter_with_open_start() {
        let range = ScaleRange::new(..1.0);
        let mut iter = range.iter_with_step(0.2);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.0);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.2);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.4);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.6);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.8);
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_range_with_step() {
        let rws: ScaleRangeWithStep = (0.0..=1.0, 0.2).into();
        let mut iter = rws.iter();
        assert_abs_diff_eq!(iter.next().unwrap(), 0.0);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.2);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.4);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.6);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.8);
        assert_abs_diff_eq!(iter.next().unwrap(), 1.0);
        assert_eq!(iter.next(), None);
    }
}
