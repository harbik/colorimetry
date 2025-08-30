// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2024-2025, Harbers Bik LLC

//! # Utility Mathematical functions
//!
//! ## Usage
//! Use these primitives to generate accurate spectral power distributions for illuminants,
//! stimuli, colorants, and to perform photometric and radiometric calculations.
//!

use approx::abs_diff_eq;

use crate::Error;
use core::f64;
use std::f64::consts::PI;

/// Constant for converting a Gaussian distribution’s standard deviation (σ)
/// to its full width at half maximum (FWHM).
const STDEV2FWHM: f64 = 2.3548200450309493; // (8.0 * 2f64.ln()).sqrt())

/// A Gaussian distribution, defined by its mean (μ) and standard deviation (σ).
/// This struct provides methods to create a Gaussian from its parameters,
/// calculate its full width at half maximum (FWHM), evaluate the Gaussian at a point,
/// and compute the peak value of the Gaussian at a given point.
pub struct Gaussian {
    mu: f64,
    sigma: f64,
}

impl Gaussian {
    /// Creates a new Gaussian distribution with the specified mean (μ) and standard deviation (σ).
    pub fn new(mu: f64, sigma: f64) -> Self {
        Self { mu, sigma }
    }

    /// Creates a new Gaussian distribution from its mean (μ) and full width at half maximum (FWHM).
    pub fn from_fwhm(mu: f64, fwhm: f64) -> Self {
        let sigma = fwhm / STDEV2FWHM;
        Self { mu, sigma }
    }

    /// Returns the mean (μ) of the Gaussian distribution.
    pub fn mu(&self) -> f64 {
        self.mu
    }

    /// Returns the standard deviation (σ) of the Gaussian distribution.
    pub fn sigma(&self) -> f64 {
        self.sigma
    }

    /// Calculates the full width at half maximum (FWHM) of the Gaussian distribution.
    pub fn fwhm(&self) -> f64 {
        let sigma = self.sigma;
        sigma * STDEV2FWHM
    }

    /// Computes the Gaussian function at `x`, normalized so that its peak (at the mean) has a value
    /// of 1.0.
    pub fn peak_one(&self, x: f64) -> f64 {
        let mu = self.mu;
        let sigma = self.sigma;
        let exponent = -((x - mu).powi(2)) / (2.0 * sigma.powi(2));
        exponent.exp()
    }

    /// Evaluates the normalized Gaussian function at `x`, such that the total area under the curve
    /// is 1.0.  In spectral terms, this represents a distribution whose integrated power is 1 W.
    pub fn normalized(&self, x: f64) -> f64 {
        let mu = self.mu;
        let sigma = self.sigma;
        let exponent = -((x - mu).powi(2)) / (2.0 * sigma.powi(2));
        (1.0 / (sigma * (2.0 * PI).sqrt())) * exponent.exp()
    }
}

#[test]
fn gaussian_peak_one_test() {
    use approx::assert_ulps_eq;
    let sigma = 10E-9;
    let mu = 500E-9;
    let gauss = Gaussian::new(mu, sigma);
    let x = mu - sigma;
    let v = gauss.peak_one(x);
    assert_ulps_eq!(v, 0.60653065971, epsilon = 1E-10);
}

#[derive(Clone, Copy)]
pub struct LineAB {
    xa: f64,
    ya: f64,
    xb: f64,
    yb: f64,
    l: f64,
    angle: f64,
}

/// Orientation of a point relative to an infinite line through points A and B,
/// when moving in the direction from point A to point B.
#[derive(Debug, PartialEq)]
pub enum Orientation {
    Left,
    Right,
    Colinear,
}

impl LineAB {
    pub fn new(a: [f64; 2], b: [f64; 2]) -> Result<Self, Error> {
        let [[xa, ya], [xb, yb]] = [a, b];
        let l = (xb - xa).hypot(yb - ya);
        let angle = (yb - ya).atan2(xb - xa);
        if l > 1E-10 {
            //if l>f64::EPSILON {
            Ok(Self {
                xa,
                ya,
                xb,
                yb,
                l,
                angle,
            })
        } else {
            Err(Error::RequiresDistinctPoints)
        }
    }

    pub fn len(&self) -> f64 {
        self.l
    }

    pub fn angle(&self) -> f64 {
        self.angle
    }

    pub fn angle_diff(&self, other: LineAB) -> f64 {
        other.angle - self.angle
    }

    pub fn angle_deg(&self) -> f64 {
        self.angle.to_degrees()
    }

    fn nom(&self, x: f64, y: f64) -> f64 {
        //(self.ya - self.yb) * x + (self.xb - self.xa) * y + (self.xa * self.yb - self.xb * self.ya)
        // see [wikipedia](https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points)
        -((self.yb - self.ya) * x - (self.xb - self.xa) * y + self.xb * self.ya - self.yb * self.xa)
    }

    pub fn orientation(&self, x: f64, y: f64) -> Orientation {
        match self.nom(x, y) {
            d if d > f64::EPSILON => Orientation::Left,
            d if d < -f64::EPSILON => Orientation::Right,
            _ => Orientation::Colinear,
        }
    }

    /// Get distance from a point to a line with a negative value for a point on
    /// the left, and positive value for a point on the right of the line,
    /// looking at the direction from A to B.
    pub fn distance_with_sign(&self, x: f64, y: f64) -> f64 {
        -self.nom(x, y) / self.l
    }

    pub fn distance(&self, x: f64, y: f64) -> f64 {
        self.distance_with_sign(x, y).abs()
    }

    /// Intersection of two lines, if not parallel or colinear.
    ///
    /// If there is an intersection, the intersection point is given, and two
    /// parameters t and u.
    /// The parameters indicated the location of the intersection on the line
    /// segments, and have a value between 0 and 1 if the intersection is
    /// between the two points used to define the lineAB.
    /// See [Wikipedia](https://en.wikipedia.org/wiki/Line–line_intersection#Given_two_points_on_each_line_segment) for the algorithm used.
    pub fn intersect(&self, line: &LineAB) -> Result<([f64; 2], f64, f64), Error> {
        if (self.angle() - line.angle()).abs() < 2.0 * f64::EPSILON {
            Err(Error::NoIntersection)
        } else {
            let [x1, y1, x2, y2] = [self.xa, self.ya, self.xb, self.yb];
            let [x3, y3, x4, y4] = [line.xa, line.ya, line.xb, line.yb];
            let den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
            let t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den;
            let u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / den;
            Ok(([x1 + t * (x2 - x1), (y1 + t * (y2 - y1))], t, u))
        }
    }
}

/// Distance of point (x,y) to line going through points (x0,y0) having slope m
///
/// See Robertson68 formula (4), but included the multiplication with the sign of m to extend the Robertson method to
/// lower color temperatures to deal with the change of the sign value of the Blackbody locus normal: it
/// is negative for temperatures higher than 1667K, and positive otherwise.
pub fn distance_to_line(x: f64, y: f64, x0: f64, y0: f64, m: f64) -> f64 {
    m.signum() * ((y - y0) - m * (x - x0)) / m.hypot(1.0)
}

#[test]
fn lineab() {
    use approx::assert_ulps_eq;
    // line pointing North
    let abup = LineAB::new([0.0, 0.0], [0.0, 1.0]).unwrap();
    assert_eq!(abup.len(), 1.0);
    let orientation = abup.orientation(-0.5, 10.0);
    assert_eq!(orientation, Orientation::Left);
    let orientation = abup.orientation(0.5, 0.5);
    assert_eq!(orientation, Orientation::Right);
    let orientation = abup.orientation(0.0, 0.5);
    assert_eq!(orientation, Orientation::Colinear);
    assert_ulps_eq!(abup.distance(0.0, 10.0), 0.0);
    assert_ulps_eq!(abup.distance_with_sign(10.0, 0.0), 10.0); // right is positive
    assert_ulps_eq!(abup.distance_with_sign(-10.0, 0.0), -10.0); // left is negative
    assert_ulps_eq!(abup.angle().to_degrees(), 90.0);

    //line point South
    let abdown = LineAB::new([0.0, 1.0], [0.0, 0.0]).unwrap();
    let orientation = abdown.orientation(-0.5, 0.5);
    assert_eq!(orientation, Orientation::Right);
    let orientation = abdown.orientation(0.5, 0.5);
    assert_eq!(orientation, Orientation::Left);
    let orientation = abdown.orientation(0.0, 0.5);
    assert_eq!(orientation, Orientation::Colinear);
    assert_ulps_eq!(abdown.angle().to_degrees(), -90.0);

    //line point East
    let abeast = LineAB::new([0.0, 0.0], [1.0, 0.0]).unwrap();
    let orientation = abeast.orientation(0.0, 0.5);
    assert_eq!(orientation, Orientation::Left);
    let orientation = abeast.orientation(0.0, -0.5);
    assert_eq!(orientation, Orientation::Right);
    let orientation = abeast.orientation(0.0, 0.0);
    assert_eq!(orientation, Orientation::Colinear);
    assert_ulps_eq!(abeast.angle().to_degrees(), 0.0);

    //line point North West
    let abnw = LineAB::new([0.0, 0.0], [-1.0, 1.0]).unwrap();
    assert_eq!(abnw.len(), 2f64.sqrt());
    let orientation = abnw.orientation(-0.5, 0.0);
    assert_eq!(orientation, Orientation::Left);
    let orientation = abnw.orientation(0.0, 1.0);
    assert_eq!(orientation, Orientation::Right);
    let orientation = abnw.orientation(-2.0, 2.0);
    assert_eq!(orientation, Orientation::Colinear);
    assert_ulps_eq!(abnw.angle().to_degrees(), 135.0);

    //line point South West
    let absw = LineAB::new([0.0, 0.0], [-1.0, -1.0]).unwrap();
    let orientation = absw.orientation(-0.5, -1.0);
    assert_eq!(orientation, Orientation::Left);
    let orientation = absw.orientation(0.0, 1.0);
    assert_eq!(orientation, Orientation::Right);
    let orientation = absw.orientation(-2.0, -2.0);
    assert_eq!(orientation, Orientation::Colinear);
    assert_ulps_eq!(absw.angle().to_degrees(), -135.0);
}

#[test]
fn lineab_intersect_test() {
    use approx::assert_ulps_eq;
    // line pointing North
    let v = LineAB::new([0.0, 0.0], [0.0, 1.0]).unwrap();
    let h = LineAB::new([-1.0, 0.0], [1.0, 0.0]).unwrap();
    let ([x, y], t, u) = v.intersect(&h).unwrap();
    assert_ulps_eq!(x, 0.0);
    assert_ulps_eq!(y, 0.0);
    assert_ulps_eq!(t, 0.0);
    assert_ulps_eq!(u, 0.5);

    let v2 = LineAB::new([0.5, 0.5], [1.0, 1.0]).unwrap();
    let ([x, y], t, u) = v2.intersect(&h).unwrap();
    assert_ulps_eq!(x, 0.0);
    assert_ulps_eq!(y, 0.0);
    assert_ulps_eq!(t, -1.0);
    assert_ulps_eq!(u, 0.5);
}

pub struct Triangle {
    xa: f64,
    ya: f64,
    xb: f64,
    yb: f64,
    xc: f64,
    yc: f64,
    area: f64,
    nom: f64,
}
impl Triangle {
    pub fn new(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Result<Self, Error> {
        let [[xa, ya], [xb, yb], [xc, yc]] = [a, b, c];
        let la = (xc - xb).hypot(yc - yb);
        let lb = (xc - xa).hypot(yc - ya);
        let lc = (xb - xa).hypot(yb - ya);
        let s = (la + lb + lc) / 2.0;
        let area2 = s * (s - la) * (s - lb) * (s - lc);
        let nom = (yb - yc) * (xa - xc) + (xc - xb) * (ya - yc);
        if s > f64::EPSILON {
            Ok(Self {
                xa,
                ya,
                xb,
                yb,
                xc,
                yc,
                area: area2.sqrt(),
                nom,
            })
        } else {
            Err(Error::RequiresDistinctPoints)
        }
    }

    pub fn area(&self) -> f64 {
        self.area
    }

    pub fn barycentric_coordinates(&self, x: f64, y: f64) -> [f64; 3] {
        let a =
            ((self.yb - self.yc) * (x - self.xc) + (self.xc - self.xb) * (y - self.yc)) / self.nom;
        let b =
            ((self.yc - self.ya) * (x - self.xc) + (self.xa - self.xc) * (y - self.yc)) / self.nom;
        let c = 1.0 - a - b;
        [a, b, c]
    }

    pub fn contains(&self, x: f64, y: f64) -> bool {
        let abc = self.barycentric_coordinates(x, y);
        !abc.iter().any(|v| !(0.0..=1.0).contains(v))
    }
}

/// A quadratic curve passing through three points.
///
/// This implementation uses a quadratic Bézier curve representation, which can handle cases
/// where the points are not functions of x (e.g., vertical tangents). The three points
/// provided at creation, `p0`, `p1`, and `p2`, are assumed to correspond to the curve at
/// parameter `t` values of 0, 0.5, and 1, respectively.
pub struct QuadraticCurve {
    // Control points for the Bézier curve
    c0: (f64, f64),
    c1: (f64, f64),
    c2: (f64, f64),
}

impl QuadraticCurve {
    /// Creates a new quadratic curve that passes through three given points.
    ///
    /// The points `p0`, `p1`, and `p2` are mapped to the curve at `t=0`, `t=0.5`, and `t=1`.
    ///
    /// # Errors
    ///
    /// Returns `Error::RequiresDistinctPoints` if the points are not distinct, which would
    /// result in a degenerate curve.
    pub fn new(p0: (f64, f64), p1: (f64, f64), p2: (f64, f64)) -> Result<Self, Error> {
        if (abs_diff_eq!(p0.0, p1.0, epsilon = f64::EPSILON)
            && abs_diff_eq!(p0.1, p1.1, epsilon = f64::EPSILON))
            || (abs_diff_eq!(p1.0, p2.0, epsilon = f64::EPSILON)
                && abs_diff_eq!(p1.1, p2.1, epsilon = f64::EPSILON))
            || (abs_diff_eq!(p0.0, p2.0, epsilon = f64::EPSILON)
                && abs_diff_eq!(p0.1, p2.1, epsilon = f64::EPSILON))
        {
            return Err(Error::RequiresDistinctPoints);
        }
        // Determine control points for a quadratic Bézier curve passing through p0, p1, p2
        // at t=0, t=0.5, and t=1.
        // p0 -> c0
        // p2 -> c2
        // p1 = (1-0.5)^2*c0 + 2*(1-0.5)*0.5*c1 + 0.5^2*c2
        // p1 = 0.25*c0 + 0.5*c1 + 0.25*c2
        // 0.5*c1 = p1 - 0.25*c0 - 0.25*c2
        // c1 = 2*p1 - 0.5*c0 - 0.5*c2
        let c0 = p0;
        let c2 = p2;
        let c1 = (
            2.0 * p1.0 - 0.5 * (p0.0 + p2.0),
            2.0 * p1.1 - 0.5 * (p0.1 + p2.1),
        );

        Ok(Self { c0, c1, c2 })
    }

    /// Evaluates the curve at a parameter `t`.
    ///
    /// The parameter `t` typically ranges from 0 to 1.
    /// - `t=0` returns the start point `p0`.
    /// - `t=0.5` returns the mid-point `p1`.
    /// - `t=1` returns the end point `p2`.
    pub fn value(&self, t: f64) -> (f64, f64) {
        let omt = 1.0 - t;
        let x = omt.powi(2) * self.c0.0 + 2.0 * omt * t * self.c1.0 + t.powi(2) * self.c2.0;
        let y = omt.powi(2) * self.c0.1 + 2.0 * omt * t * self.c1.1 + t.powi(2) * self.c2.1;
        (x, y)
    }

    /// Calculates the derivative of the curve with respect to `t`.
    ///
    /// This gives the tangent vector `(dx/dt, dy/dt)` at parameter `t`.
    pub fn derivative(&self, t: f64) -> (f64, f64) {
        let dx_dt = 2.0 * (1.0 - t) * (self.c1.0 - self.c0.0) + 2.0 * t * (self.c2.0 - self.c1.0);
        let dy_dt = 2.0 * (1.0 - t) * (self.c1.1 - self.c0.1) + 2.0 * t * (self.c2.1 - self.c1.1);
        (dx_dt, dy_dt)
    }

    /// Calculates the slope angle of the tangent to the curve at parameter `t`.
    pub fn slope_angle(&self, t: f64) -> f64 {
        let (dx_dt, dy_dt) = self.derivative(t);
        dy_dt.atan2(dx_dt)
    }
}

#[test]
fn quadratic_curve_test() {
    use approx::assert_ulps_eq;

    let p0 = (0.0, 0.0);
    let p1 = (1.0, 2.0);
    let p2 = (2.0, 0.0);
    let curve = QuadraticCurve::new(p0, p1, p2).unwrap();

    // Evaluate at t=0, t=0.5, t=1
    assert_eq!(curve.value(0.0), p0);
    let mid_point = curve.value(0.5);
    assert_ulps_eq!(mid_point.0, 1.0);
    assert_ulps_eq!(mid_point.1, 2.0);
    assert_eq!(curve.value(1.0), p2);

    // Derivative at t=0 and t=1 should be horizontal
    let deriv_0 = curve.derivative(0.0);
    assert_ulps_eq!(deriv_0.0, 2.0);
    assert_ulps_eq!(deriv_0.1, 8.0);
    let deriv_1 = curve.derivative(1.0);
    assert_ulps_eq!(deriv_1.0, 2.0);
    assert_ulps_eq!(deriv_1.1, -8.0);

    // Slope angle at t=0 and t=1
    assert_ulps_eq!(curve.slope_angle(0.0), 4f64.atan());
    assert_ulps_eq!(curve.slope_angle(1.0), -4f64.atan());
}

#[test]
fn triangle_test() {
    use approx::assert_ulps_eq;

    let t1 = Triangle::new([0.0, 0.0], [0.0, 1.0], [1.0, 0.0]).unwrap();
    assert_ulps_eq!(t1.area(), 0.5);

    let [a, b, c] = t1.barycentric_coordinates(0.0, 0.0);
    assert_ulps_eq!(a, 1.0);
    assert_ulps_eq!(b, 0.0);
    assert_ulps_eq!(c, 0.0);

    let [a, b, c] = t1.barycentric_coordinates(0.0, 1.0);
    assert_ulps_eq!(a, 0.0);
    assert_ulps_eq!(b, 1.0);
    assert_ulps_eq!(c, 0.0);

    let [a, b, c] = t1.barycentric_coordinates(1.0, 0.0);
    assert_ulps_eq!(a, 0.0);
    assert_ulps_eq!(b, 0.0);
    assert_ulps_eq!(c, 1.0);

    let h = 3f64.sqrt() / 2.0; // height equilateral rectangle with side length of 1
    let t2 = Triangle::new([0.0, h], [-0.5, 0.0], [0.5, 0.0]).unwrap();
    assert_ulps_eq!(t2.area(), h / 2.0); // https://en.wikipedia.org/wiki/Equilateral_triangle
    let [a, b, c] = t2.barycentric_coordinates(0.0, h / 3.0); // Apothem is h/3
    println!("{a:.4} {b:.4} {c:.4}");
    assert_ulps_eq!(a, 1.0 / 3.0);
    assert_ulps_eq!(b, 1.0 / 3.0);
    assert_ulps_eq!(c, 1.0 / 3.0);
}

#[inline]
/// Calculates the Euclidian Distance between two points,
/// in arbitrary space dimension. Typically used in two or three
/// dimensional spaces.
pub fn distance(p: &[f64], q: &[f64]) -> f64 {
    let s2 = p
        .iter()
        .zip(q)
        .fold(0.0, |s, (&pi, &qi)| s + (pi - qi).powi(2));
    s2.sqrt()
}
#[test]
fn test_distance() {
    approx::assert_abs_diff_eq!(distance(&[0.0], &[1.0]), 1.0);
    approx::assert_abs_diff_eq!(distance(&[0.0], &[0.0]), 0.0);
    approx::assert_abs_diff_eq!(distance(&[0.0], &[2.0]), 2.0);
    approx::assert_abs_diff_eq!(distance(&[0.0, 0.0], &[3.0, 4.0]), 5.0);
    approx::assert_abs_diff_eq!(distance(&[0.0, 0.0, 0.0], &[1.0, 1.0, 1.0]), 3.0f64.sqrt());
}
