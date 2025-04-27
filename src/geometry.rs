use core::f64;

use nalgebra::ComplexField;
use wasm_bindgen::prelude::wasm_bindgen;

use crate::error::CmtError;

/// Distance of point (x,y) to line going through points (x0,y0) having slope m
///
/// See Robertson68 formula (4), but included the multiplication with the sign of m to extend the Robertson method to
/// lower color temperatures to deal with the change of the sign value of the Blackbody locus normal: it
/// is negative for temperatures higher than 1667K, and positive otherwise.
pub fn distance_to_line(x: f64, y: f64, x0: f64, y0: f64, m: f64) -> f64 {
    m.signum() * ((y - y0) - m * (x - x0)) / m.hypot(1.0)
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
    pub fn try_new(a: [f64; 2], b: [f64; 2]) -> Result<Self, CmtError> {
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
            Err(CmtError::RequiresDistinctPoints)
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
        let n = self.nom(x, y);
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
    pub fn intersect(&self, line: &LineAB) -> Result<([f64; 2], f64, f64), CmtError> {
        if (self.angle() - line.angle()).abs() < 2.0 * f64::EPSILON {
            Err(CmtError::NoIntersection)
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

#[test]
fn lineab() {
    use approx::assert_ulps_eq;
    // line pointing North
    let abup = LineAB::try_new([0.0, 0.0], [0.0, 1.0]).unwrap();
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
    let abdown = LineAB::try_new([0.0, 1.0], [0.0, 0.0]).unwrap();
    let orientation = abdown.orientation(-0.5, 0.5);
    assert_eq!(orientation, Orientation::Right);
    let orientation = abdown.orientation(0.5, 0.5);
    assert_eq!(orientation, Orientation::Left);
    let orientation = abdown.orientation(0.0, 0.5);
    assert_eq!(orientation, Orientation::Colinear);
    assert_ulps_eq!(abdown.angle().to_degrees(), -90.0);

    //line point East
    let abeast = LineAB::try_new([0.0, 0.0], [1.0, 0.0]).unwrap();
    let orientation = abeast.orientation(0.0, 0.5);
    assert_eq!(orientation, Orientation::Left);
    let orientation = abeast.orientation(0.0, -0.5);
    assert_eq!(orientation, Orientation::Right);
    let orientation = abeast.orientation(0.0, 0.0);
    assert_eq!(orientation, Orientation::Colinear);
    assert_ulps_eq!(abeast.angle().to_degrees(), 0.0);

    //line point North West
    let abnw = LineAB::try_new([0.0, 0.0], [-1.0, 1.0]).unwrap();
    assert_eq!(abnw.len(), 2f64.sqrt());
    let orientation = abnw.orientation(-0.5, 0.0);
    assert_eq!(orientation, Orientation::Left);
    let orientation = abnw.orientation(0.0, 1.0);
    assert_eq!(orientation, Orientation::Right);
    let orientation = abnw.orientation(-2.0, 2.0);
    assert_eq!(orientation, Orientation::Colinear);
    assert_ulps_eq!(abnw.angle().to_degrees(), 135.0);

    //line point South West
    let absw = LineAB::try_new([0.0, 0.0], [-1.0, -1.0]).unwrap();
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
    let v = LineAB::try_new([0.0, 0.0], [0.0, 1.0]).unwrap();
    let h = LineAB::try_new([-1.0, 0.0], [1.0, 0.0]).unwrap();
    let ([x, y], t, u) = v.intersect(&h).unwrap();
    assert_ulps_eq!(x, 0.0);
    assert_ulps_eq!(y, 0.0);
    assert_ulps_eq!(t, 0.0);
    assert_ulps_eq!(u, 0.5);

    let v2 = LineAB::try_new([0.5, 0.5], [1.0, 1.0]).unwrap();
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
    pub fn try_new(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Result<Self, CmtError> {
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
            Err(CmtError::RequiresDistinctPoints)
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

    pub fn within(&self, x: f64, y: f64) -> bool {
        let abc = self.barycentric_coordinates(x, y);
        !abc.iter().any(|v| !(0.0..=1.0).contains(v))
    }
}

#[test]
fn triangle_test() {
    use approx::assert_ulps_eq;

    let t1 = Triangle::try_new([0.0, 0.0], [0.0, 1.0], [1.0, 0.0]).unwrap();
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
    let t2 = Triangle::try_new([0.0, h], [-0.5, 0.0], [0.5, 0.0]).unwrap();
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
