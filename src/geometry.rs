
use wasm_bindgen::prelude::wasm_bindgen;


pub struct LineAB {
    xa: f64,
    ya: f64,
    xb: f64,
    yb: f64
}

/// Orientation of a point relative to an infinite line through points A and B,
/// when moving in the direction from point A to point B.
#[derive(Debug, PartialEq)]
pub enum Orientation {
    Left,
    Right,
    Colinear
}

impl LineAB {
    pub fn try_new(a: [f64;2], b: [f64;2]) -> Result<Self, CmError> {
        if a!=b { Err(CmError::RequiresDistinctPoints)}
        else {
        let [[xa, ya], [xb, yb]] = [a, b];
        Ok( Self { xa, ya, xb, yb }) 
    }
    }

    pub fn orientation(&self, x:f64 ,y:f64) -> Orientation {
        let d = (self.ya - self.yb) * x + (self.xb - self.xa) * y + (self.xa * self.yb - self.xb * self.ya);
        match d {
            _ if d>0.0 => Orientation::Left,
            _ if d<0.0 => Orientation::Right,
            _ => Orientation::Colinear
        }
    }
}

#[test]
fn lineab() {
    // line pointing North
    let abup = LineAB::try_new([0.0, 0.0], [0.0, 1.0]).unwrap();
    let orientation = abup.orientation(-0.5, 10.0);
    assert_eq!(orientation, Orientation::Left);
    let orientation = abup.orientation(0.5, 0.5);
    assert_eq!(orientation, Orientation::Right);
    let orientation = abup.orientation(0.0, 0.5);
    assert_eq!(orientation, Orientation::Colinear);
    
    //line point South
    let abdown = LineAB::try_new([0.0, 1.0], [0.0, 0.0]).unwrap();
    let orientation = abdown.orientation(-0.5, 0.5);
    assert_eq!(orientation, Orientation::Right);
    let orientation = abdown.orientation(0.5, 0.5);
    assert_eq!(orientation, Orientation::Left);
    let orientation = abdown.orientation(0.0, 0.5);
    assert_eq!(orientation, Orientation::Colinear);

    //line point East
    let abeast = LineAB::try_new([0.0, 0.0], [1.0, 0.0]).unwrap();
    let orientation = abeast.orientation(0.0, 0.5);
    assert_eq!(orientation, Orientation::Left);
    let orientation = abeast.orientation(0.0, -0.5);
    assert_eq!(orientation, Orientation::Right);
    let orientation = abeast.orientation(0.0, 0.0);
    assert_eq!(orientation, Orientation::Colinear);
}