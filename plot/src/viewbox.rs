use std::fmt::Display;

pub trait GetViewBox {
    fn view_box(&self) -> ViewBox;
}
#[derive(Clone, Debug, Default)]
pub struct ViewBox(i32, i32, u32, u32);

impl ViewBox {
    pub fn new(x: i32, y: i32, width: u32, height: u32) -> Self {
        ViewBox(x, y, width, height)
    }
    
    pub fn x(&self) -> i32 {
        self.0
    }

    pub fn y(&self) -> i32 {
        self.1
    }

    pub fn width(&self) -> u32 {
        self.2
    }
    
    pub fn height(&self) -> u32 {
        self.3
    }
    
    pub fn extend(&mut self, other: &ViewBox) {
        self.0 = self.0.min(other.0);
        self.1 = self.1.min(other.1);
        self.2 = self.2.max(other.2);
        self.3 = self.3.max(other.3);
    }

    pub fn extend_with_pos(&mut self, other: &ViewBox, x: i32, y: i32) {
        self.0 = self.0.min(other.0 + x);
        self.1 = self.1.min(other.1 + y);
        self.2 = self.2.max((other.2 as i32 + x) as u32);
        self.3 = self.3.max((other.3 as i32 + y) as u32);
    }
}

impl Display for ViewBox {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {} {} {}", self.0, self.1, self.2, self.3)
    }
}

impl From<ViewBox> for (i32, i32, u32, u32) {
    fn from(value: ViewBox) -> Self {
        (value.0, value.1, value.2, value.3)
    }
}
