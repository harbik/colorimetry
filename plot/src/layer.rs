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

    /*
    pub fn to_group(self) -> Group {
        self.0
    }
     */
}

/*
impl From<Layer> for Group {
    fn from(layer: Layer) -> Self {
        layer.to_group()
    }
}

impl From<Group> for Layer {
    fn from(group: Group) -> Self {
        Layer(group)
    }
}
*/

impl From<Layer> for Box<dyn Node> {
    fn from(layer: Layer) -> Self {
        Box::new(layer.0)
    }
}

/*
impl From<Layer> for Element {
    fn from(layer: Layer) -> Self {
        layer.to_group().into()
    }
}
 */

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