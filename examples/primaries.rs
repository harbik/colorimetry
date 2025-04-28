use argmin::{
    core::{CostFunction, Executor, State, TerminationReason},
    solver::neldermead::NelderMead,
};
use colored::Colorize;
use colorimetry::prelude::*;
use strum::IntoEnumIterator as _;

#[allow(dead_code)]
struct Gauss {
    x: f64,
    y: f64,
    d: StdIlluminant,
}

impl Gauss {
    fn new(xyz: XYZ, d: StdIlluminant) -> Self {
        let [x, y] = xyz.chromaticity();
        Self { x, y, d }
    }
}

impl CostFunction for Gauss {
    type Param = Vec<f64>;

    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        let [l, w] = param.clone().try_into().unwrap();
        let [xt, yt] = CIE1931
            .xyz(&StdIlluminant::D65, Some(&Colorant::gaussian(l, w)))
            .chromaticity();
        //  println!("({l},{w}) cost: {xt:.4}, {yt:.4}");
        Ok((xt - self.x).hypot(yt - self.y))
    }
}
#[allow(dead_code)]
struct GaussWithAnchor {
    x: f64,
    y: f64,
    d: StdIlluminant,
    anchor: XYZ,
}

impl GaussWithAnchor {
    fn new(xyz: XYZ, anchor: XYZ, d: StdIlluminant) -> Self {
        let [x, y] = xyz.chromaticity();
        Self {
            x,
            y,
            anchor: anchor.set_illuminance(100.0),
            d,
        }
    }
}

impl CostFunction for GaussWithAnchor {
    type Param = Vec<f64>;

    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        let [l, w, c]: [f64; 3] = param.clone().try_into().unwrap();
        let r = CIE1931
            .xyz(&StdIlluminant::D65, Some(&Colorant::gaussian(l, w)))
            .set_illuminance(100.0);
        let t = c * self.anchor + (1.0 - c) * r;
        let [xt, yt] = t.chromaticity();
        Ok((xt - self.x).hypot(yt - self.y))
    }
}

fn gauss(space: RgbSpace, i: usize) -> Result<Vec<f64>, String> {
    let [x, y] = space.primaries_chromaticity()[i];
    let xyz = XYZ::try_from_chromaticity(x, y, None, None).unwrap();
    let d = space.white();
    let problem = Gauss::new(xyz, d);

    let l = xyz.dominant_wavelength(CIE1931.xyz_d65()).unwrap();
    let w = 40.0;

    let solver = NelderMead::new(vec![vec![l, w], vec![l + 5.0, w], vec![l, w + 5.0]]);
    let res = Executor::new(problem, solver)
        .configure(|state| state.max_iters(1000).target_cost(1E-5))
        .run()
        .unwrap();
    let tr = res.state.get_termination_reason().unwrap();
    if tr == &TerminationReason::TargetCostReached {
        Ok(res.state.best_param.unwrap())
    } else {
        Err(res
            .state
            .get_termination_reason()
            .unwrap()
            .text()
            .to_string())
    }
}

fn gauss_with_anchor(space: RgbSpace, i: usize, j: usize) -> Result<Vec<f64>, String> {
    let [x, y] = space.primaries_chromaticity()[i];
    let xyz = XYZ::try_from_chromaticity(x, y, None, None).unwrap();
    let [xb, yb] = space.primaries_chromaticity()[j];
    let xyzb = XYZ::try_from_chromaticity(xb, yb, None, None).unwrap();
    let d = space.white();
    let problem = GaussWithAnchor::new(xyz, xyzb, d);

    // start parameters
    let l = xyz.dominant_wavelength(CIE1931.xyz_d65()).unwrap();
    let w = 40.0;
    let c = 0.1;

    let solver = NelderMead::new(vec![
        vec![l, w, c],
        vec![l + 5.0, w, c],
        vec![l, w + 5.0, c],
        vec![l, w, c + 0.05],
    ]);
    let res = Executor::new(problem, solver)
        .configure(|state| state.max_iters(1000).target_cost(1E-5))
        .run()
        .unwrap();
    let tr = res.state.get_termination_reason().unwrap();
    if tr == &TerminationReason::TargetCostReached {
        Ok(res.state.best_param.unwrap())
    } else {
        Err(res
            .state
            .get_termination_reason()
            .unwrap()
            .text()
            .to_string())
    }
}

fn report(desc: &str, r: [f64; 3], g: [f64; 2], b: [f64; 2]) {
    println!("\n{}", desc.bold().underline());
    println!("{}", format!("{:?}", r).red());
    println!("{}", format!("{:?}", g).green());
    println!("{}", format!("{:?}", b).blue());
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for rgbspace in RgbSpace::iter() {
        report(
            rgbspace.name(),
            gauss_with_anchor(rgbspace, 0, 2)?.try_into().unwrap(),
            gauss(rgbspace, 1)?.try_into().unwrap(),
            gauss(rgbspace, 2)?.try_into().unwrap(),
        );
    }

    Ok(())
}
