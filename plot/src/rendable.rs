
use svg::node::element::SVG;

use crate::view::ViewParameters;

pub trait Rendable {
    // Required
    fn view_parameters(&self) -> ViewParameters;

    fn set_view_parameters(&mut self, view_box: ViewParameters); 

    fn render(&self) -> SVG;

    // Optional

    fn set_x(&mut self, x: i32) {
        let mut vb = self.view_parameters();
        vb.set_x(x);
        self.set_view_parameters(vb);
    }
    
    fn set_y(&mut self, y: i32) {
        let mut vb = self.view_parameters();
        vb.set_y(y);
        self.set_view_parameters(vb);
    }

    fn width(&self) -> u32 {
        self.view_parameters().width()
    }

    fn set_width(&mut self, width: u32) {
        let mut vb = self.view_parameters();
        vb.set_width(width);
        self.set_view_parameters(vb);
    }

    fn height(&self) -> u32 {
        self.view_parameters().height()
    }

    fn set_height(&mut self, height: u32) {
        let mut vb = self.view_parameters();
        vb.set_height(height);
        self.set_view_parameters(vb);
    }

    fn view_box(&mut self) -> (i32, i32, u32, u32) {
        let vp = self.view_parameters();
        vp.view_box()
    }

    fn set_view_box(&mut self, vx:i32, vy: i32, vw: u32, vh:u32) {
        let mut vp = self.view_parameters();
        vp.set_view_box(vx, vy, vw, vh);
        self.set_view_parameters(vp);
    }
}
