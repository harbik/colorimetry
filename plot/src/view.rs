use std::fmt::Display;

#[derive(Clone, Debug, Default)]
pub struct ViewParameters {
    x: Option<i32>,
    y: Option<i32>,
    width: u32,
    height: u32,
    vx: i32,
    vy: i32,
    vw: u32,
    vh: u32,
}

impl ViewParameters {
    pub fn new(vx: i32, vy: i32, vw: u32, vh: u32, width: u32, height: u32) -> Self {
        ViewParameters {
            vx,
            vy,
            vw,
            vh,
            x: None,
            y: None,
            width,
            height,
        }
    }

    pub fn view_box(&self) -> (i32, i32, u32, u32) {
        (self.vx, self.vy, self.vw, self.vh)
    }

    pub fn view_box_str(&self) -> String {
        format!("{} {} {} {}", self.vx, self.vy, self.vw, self.vh)
    }

    pub fn set_view_box(&mut self, vx: i32, vy: i32, vw: u32, vh: u32) {
        self.vx = vx;
        self.vy = vy;
        self.vw = vw;
        self.vh = vh;

    }

    // Optional SVG x value
    pub fn x(&self) ->  Option<i32> {
        self.x
    }

    pub fn set_x(&mut self, x: i32) {
        self.x = Some(x);
    }

    // Optinal SVG y value
    pub fn y(&self) ->  Option<i32> {
        self.y
    }

    pub fn set_y(&mut self, y: i32) {
        self.y = Some(y);
    }


    pub fn width(&self) -> u32 {
        self.width
    }

    pub fn set_width(&mut self, width: u32) {
        self.width = width;
    }

    pub fn height(&self) -> u32 {
        self.height
    }

    pub fn set_height(&mut self, height: u32) {
        self.height = height;
    }

    
    pub fn extend(&mut self, other: &ViewParameters) {
        self.vx = self.vx.min(other.vx);
        self.vy = self.vy.min(other.vy);
        self.vw = self.vw.max(other.vw);
        self.vh = self.vh.max(other.vh);
    }

    pub fn extend_with_pos(&mut self, other: &ViewParameters, x: i32, y: i32) {
        self.vx = self.vx.min(other.vx);
        self.vy = self.vy.min(other.vy);
        self.vw = self.vw.max((other.vw as i32 + x) as u32);
        self.vh = self.vh.max((other.vh as i32 + y) as u32);
    }

    pub fn extend_with_margin(&mut self, other: &ViewParameters, margin: i32) {
        self.vx = self.vx.min(other.vx);
        self.vy = self.vy.min(other.vy);
        self.vw = self.vw.max((other.vw as i32 + 2 * margin) as u32);
        self.vh = self.vh.max((other.vh as i32 + 2 * margin) as u32);
    }
}

impl Display for ViewParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {} {} {}", self.vx, self.vy, self.vw, self.vh)
    }
}

impl From<ViewParameters> for (i32, i32, u32, u32) {
    fn from(value: ViewParameters) -> Self {
        (value.vx, value.vy, value.vw, value.vh)
    }
}
