use argmin::{
    core::{CostFunction, Executor, State, TerminationReason},
    solver::neldermead::NelderMead,
};
use colored::Colorize;
use colorimetry::{
    colorant::Colorant,
    illuminant::CieIlluminant,
    observer::Observer::Cie1931,
    rgb::RgbSpace,
    xyz::{Chromaticity, XYZ},
};
use strum::IntoEnumIterator;

#[allow(dead_code)]
struct Gauss {
    chromaticity: Chromaticity,
    d: CieIlluminant,
}

impl Gauss {
    fn new(xyz: XYZ, d: CieIlluminant) -> Self {
        let chromaticity = xyz.chromaticity();
        Self { chromaticity, d }
    }
}

impl CostFunction for Gauss {
    type Param = Vec<f64>;

    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        let [l, w] = param.clone().try_into().unwrap();
        let param_chromaticity = colorimetry::observer::Observer::Cie1931
            .xyz(&CieIlluminant::D65, Some(&Colorant::gaussian(l, w)))
            .chromaticity();
        //  println!("({l},{w}) cost: {xt:.4}, {yt:.4}");
        Ok((param_chromaticity.x() - self.chromaticity.x())
            .hypot(param_chromaticity.y() - self.chromaticity.y()))
    }
}

#[allow(dead_code)]
struct GaussWithAnchor {
    chromaticity: Chromaticity,
    d: CieIlluminant,
    anchor: XYZ,
}

impl GaussWithAnchor {
    fn new(xyz: XYZ, anchor: XYZ, d: CieIlluminant) -> Self {
        let chromaticity = xyz.chromaticity();
        Self {
            chromaticity,
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
        let r = colorimetry::observer::Observer::Cie1931
            .xyz(&CieIlluminant::D65, Some(&Colorant::gaussian(l, w)))
            .set_illuminance(100.0);
        let t = c * self.anchor + (1.0 - c) * r;
        let t_chromaticity = t.chromaticity();
        Ok((t_chromaticity.x() - self.chromaticity.x())
            .hypot(t_chromaticity.y() - self.chromaticity.y()))
    }
}

fn gauss(space: RgbSpace, rgb_channel_i: usize) -> Result<Vec<f64>, String> {
    let chromaticity = space.target_chromaticities_cie1931()[rgb_channel_i];
    let xyz = XYZ::from_chromaticity(chromaticity, None, None).unwrap();
    let d = space.white();
    let problem = Gauss::new(xyz, d);

    let l = xyz.dominant_wavelength(Cie1931.xyz_d65()).unwrap();
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
    let chromaticity_i = space.target_chromaticities_cie1931()[i];
    let xyz = XYZ::from_chromaticity(chromaticity_i, None, None).unwrap();
    let chromaticity_j = space.target_chromaticities_cie1931()[j];
    let xyzb = XYZ::from_chromaticity(chromaticity_j, None, None).unwrap();
    let d = space.white();
    let problem = GaussWithAnchor::new(xyz, xyzb, d);

    // start parameters
    let l = xyz.dominant_wavelength(Cie1931.xyz_d65()).unwrap();
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

/// Report results
/// * `desc` - Description of the RGB space
/// * `r` - Result for red channel, printed in red in the terminal
/// * `g` - Result for green channel, printed in green in the terminal
/// * `b` - Result for blue channel, printed in blue in the terminal
fn report(desc: &str, r: [f64; 3], g: [f64; 2], b: [f64; 2]) {
    println!("\n{}", desc.bold().underline());
    println!("{}", format!("{r:?}").red());
    println!("{}", format!("{g:?}").green());
    println!("{}", format!("{b:?}").blue());
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
