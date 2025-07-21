use std::ops::Deref;

use svg::{node::element::Group, Node};

#[derive(Debug, Clone, Default)]
pub struct Layer(pub(crate) Group);

impl Layer {
    pub fn new() -> Self {
        Layer(Group::new())
    }

    pub fn add_layer(self, node: impl Into<Box<dyn Node>>) -> Self {
        Layer(self.0.add(node))
    }

    pub fn set(self, key: &str, value: impl Into<String>) -> Self {
        Layer(self.0.set(key, value.into()))
    }

    pub fn set_class_and_style(self, class: Option<&str>, style: Option<&str>) -> Self {
        let mut layer = self;
        if let Some(c) = class {
            layer = layer.set("class", c);
        }
        if let Some(s) = style {
            layer = layer.set("style", s);
        }
        layer
    }
}


impl From<Layer> for Box<dyn Node> {
    fn from(layer: Layer) -> Self {
        Box::new(layer.0)
    }
}

impl std::fmt::Display for Layer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Deref for Layer {
    type Target = Group;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for Layer {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
