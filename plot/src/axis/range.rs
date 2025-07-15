use std::ops::{Bound, Range, RangeBounds};

#[derive(Debug, Clone, Copy)]
pub struct ChartRange {
    pub start: f64,
    pub end: f64,
    pub end_included: bool,
}

impl ChartRange {
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

        Self { start, end, end_included }
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
    /// Produces values aligned with the step size, for use to draw grid lines and ticks.
    pub fn iter_with_step(&self, step: f64) -> ChartRangeIterator {
        let start = (self.start/ step).ceil() * step; // Start from the first tick

        ChartRangeIterator {
            start,
            end: self.end,
            step,
            end_included: self.end_included,
        }
    }
}

pub struct ChartRangeIterator {
    start: f64,
    end: f64,
    end_included: bool,
    step: f64,
}

impl Iterator for ChartRangeIterator {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start < self.end || ( self.end_included && approx::abs_diff_eq!(self.start, self.end, epsilon = 1e-10)) {
            let current = self.start;
            self.start += self.step;
            Some(current)
        } else {
            None
        }
    }
}


#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn test_chart_range() {
        let range = ChartRange::new(0.0..1.0);
        assert_abs_diff_eq!(range.end, 1.0);
        assert_abs_diff_eq!(range.span(), 1.0);
        assert_abs_diff_eq!(range.scale(0.5), 0.5);
        assert_abs_diff_eq!(range.scale_descent(0.5), 0.5);
    }

    #[test]
    fn test_iter_with_step() {
        let range = ChartRange::new(0.0..=1.0);
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
        let range = ChartRange::new(0.0..1.0);
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
        let range = ChartRange::new(..1.0);
        let mut iter = range.iter_with_step(0.2);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.0);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.2);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.4);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.6);
        assert_abs_diff_eq!(iter.next().unwrap(), 0.8);
        assert_eq!(iter.next(), None);
    }
}