use std::ops::{Bound, Range, RangeBounds};

#[derive(Debug, Clone, Copy)]
pub struct ChartRange {
    pub start: f64,
    pub end: f64,
}

impl ChartRange {
    pub fn new<R: RangeBounds<f64>>(range: R) -> Self {
        let start = match range.start_bound() {
            Bound::Included(&s) | Bound::Excluded(&s) => s,
            Bound::Unbounded => panic!("Unbounded start not allowed"),
        };
        let end = match range.end_bound() {
            Bound::Included(&e) | Bound::Excluded(&e) => e,
            Bound::Unbounded => panic!("Unbounded end not allowed"),
        };

        Self { start, end }
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

    pub fn scale_descent(&self, value: f64) -> f64 {
        (self.end - value) / self.span()
    }

    /// Creates an iterator over the range with a specified step size.
    pub fn iter_with_step(&self, step: f64) -> ChartRangeIterator {
        let start = (self.start/ step).ceil() * step; // Start from the first tick
        ChartRangeIterator {
            start,
            end: self.end,
            step,
        }
    }
}

pub struct ChartRangeIterator {
    start: f64,
    end: f64,
    step: f64,
}

impl Iterator for ChartRangeIterator {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start < self.end {
            let current = self.start;
            self.start += self.step;
            Some(current)
        } else {
            None
        }
    }
}
