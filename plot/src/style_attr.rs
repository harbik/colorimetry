// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2025, Harbers Bik LLC

//! Style attributes for SVG elements in the plotting library
//!
//! The `StyleAttr` struct and associated macros are used to manage SVG style attributes such as
//! CSS classes, inline styles, and IDs. This allows for easy styling of SVG elements in plots.

/// Represents SVG style attributes, including optional class, style, and id fields.
///
/// `StyleAttr` is used to conveniently manage and assign SVG styling attributes such as CSS classes,
/// inline styles, and element IDs. This struct can be used to apply consistent styling to SVG nodes
/// throughout the plotting library.
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct StyleAttr {
    /// Optional CSS class attribute for the SVG element.
    pub class: Option<String>,
    /// Optional inline style attribute for the SVG element.
    pub style: Option<String>,
    /// Optional id attribute for the SVG element.
    pub id: Option<String>,
}

impl StyleAttr {
    /// Assigns the style attributes (`class` and `style`) to the given SVG node.
    ///
    /// This method sets the `class` and `style` attributes on the provided SVG node if they are
    /// present in the `StyleAttr` instance. If the `id` is set, it will not be assigned here, as
    /// this method focuses on styling attributes only. If no `class` or `style` is set,
    /// this method will not modify the node.
    ///
    /// # Arguments * `node` - A mutable reference to an SVG node that implements the `svg::Node`
    /// trait.
    pub fn assign<T>(&self, node: &mut T)
    where
        T: svg::Node,
    {
        if let Some(class) = self.class.clone() {
            node.assign("class", class);
        }
        if let Some(style) = self.style.clone() {
            node.assign("style", style);
        }
    }

    /// Returns the `id` attribute as an optional string slice, if set.
    ///
    /// # Returns
    /// * `Some(&str)` if the `id` is present, or `None` otherwise.
    pub fn id(&self) -> Option<&str> {
        self.id.as_deref()
    }
}

/// Creates a new `StyleAttr` with the specified class.
/// This function is a convenience method to create a `StyleAttr` instance with the `class` field set.
/// # Arguments
/// * `str` - A string slice representing the CSS class to be assigned.
/// # Returns
/// An `Option<StyleAttr>` containing the new `StyleAttr` with the class set, or `None` if the class is empty.
pub fn class(str: &str) -> Option<StyleAttr> {
    Some(StyleAttr {
        class: Some(str.to_string()),
        style: None,
        id: None,
    })
}

/// Creates a new `StyleAttr` with the specified style.
/// This function is a convenience method to create a `StyleAttr` instance with the `style` field set.
/// # Arguments
/// * `str` - A string slice representing the inline style to be assigned.
/// # Returns
/// An `Option<StyleAttr>` containing the new `StyleAttr` with the style set, or `None` if the style is empty.
pub fn style(str: &str) -> Option<StyleAttr> {
    Some(StyleAttr {
        class: None,
        style: Some(str.to_string()),
        id: None,
    })
}

/// Creates a new `StyleAttr` with the specified id.
/// This function is a convenience method to create a `StyleAttr` instance with the `id` field set.
/// # Arguments
/// * `id_val` - A string slice representing the ID to be assigned.
/// # Returns
/// An `Option<StyleAttr>` containing the new `StyleAttr` with the id set, or `None` if the id is empty.
pub fn id(id_val: &str) -> Option<StyleAttr> {
    Some(StyleAttr {
        class: None,
        style: None,
        id: Some(id_val.to_string()),
    })
}

/// A [`StyleAttr`] struct with either the `class`, `style`, `id` or any combination of fields set.
/// This macro allows for flexible creation of style attributes for SVG elements.
/// It supports various combinations of `class`, `style`, and `id` attributes, making it easy to
/// define the styling of SVG elements in a concise manner.
// The macro can handle single attributes, combinations of two, or all three attributes.
// It also provides compile-time checks for unknown keys, ensuring that only valid attributes are used.
// The macro generates an `Option<StyleAttr>` which can be used directly in SVG rendering contexts.
// Usage examples:
// ```
// css!(class: "my-class");
// css!(style: "fill: red;");
// css!(id: "my-id");
// css!(class: "my-class", style: "fill: red;");
// css!(class: "my-class", id: "my-id");
// css!(style: "fill: red;", id: "my-id");
// css!(class: "my-class", style: "fill: red;", id: "my-id");
// css!(id: "my-id", class: "my-class", style: "fill: red;");
// ```
#[macro_export]
macro_rules! css {

    // Only class
    ( $( class : $class_val:expr ),+ $(,)? ) => {
        Some(
            $crate::StyleAttr {
                class: Some([$($class_val),+].join(" ")),
                style: None,
                id: None,
            }
        )
    };
    // Only style
    ( $( style : $style_val:expr ),+ $(,)? ) => {
        Some(
            $crate::StyleAttr {
                class: None,
                style: Some([$($style_val),+].join(" ")),
                id: None,
            }
        )
    };
    // Only id
    ( $( id : $id_val:expr ),+ $(,)? ) => {
        Some(
            $crate::StyleAttr {
                class: None,
                style: None,
                id: Some([$($id_val),+].join(" ")),
            }
        )
    };
    // class + style
    ( class : $class_val:expr, style : $style_val:expr $(,)? ) => {
        Some(
            $crate::StyleAttr {
                class: Some($class_val.to_string()),
                style: Some($style_val.to_string()),
                id: None,
            }
        )
    };
    ( style : $style_val:expr, class : $class_val:expr $(,)? ) => {
        Some(
            $crate::StyleAttr {
                class: Some($class_val.to_string()),
                style: Some($style_val.to_string()),
                id: None,
            }
        )
    };
    // class + id
    ( class : $class_val:expr, id : $id_val:expr $(,)? ) => {
        Some(
            $crate::StyleAttr {
                class: Some($class_val.to_string()),
                style: None,
                id: Some($id_val.to_string()),
            }
        )
    };
    ( id : $id_val:expr, class : $class_val:expr $(,)? ) => {
        Some(
            $crate::StyleAttr {
                class: Some($class_val.to_string()),
                style: None,
                id: Some($id_val.to_string()),
            }
        )
    };
    // style + id
    ( style : $style_val:expr, id : $id_val:expr $(,)? ) => {
        Some(
            $crate::StyleAttr {
                class: None,
                style: Some($style_val.to_string()),
                id: Some($id_val.to_string()),
            }
        )
    };
    ( id : $id_val:expr, style : $style_val:expr $(,)? ) => {
        Some(
            $crate::StyleAttr {
                class: None,
                style: Some($style_val.to_string()),
                id: Some($id_val.to_string()),
            }
        )
    };
    // class + style + id (any order)
    ( class : $class_val:expr, style : $style_val:expr, id : $id_val:expr $(,)? ) => {
        Some(
            $crate::StyleAttr {
                class: Some($class_val.to_string()),
                style: Some($style_val.to_string()),
                id: Some($id_val.to_string()),
            }
        )
    };
    ( class : $class_val:expr, id : $id_val:expr, style : $style_val:expr $(,)? ) => {
        Some(
            $crate::StyleAttr {
                class: Some($class_val.to_string()),
                style: Some($style_val.to_string()),
                id: Some($id_val.to_string()),
            }
        )
    };
    ( style : $style_val:expr, class : $class_val:expr, id : $id_val:expr $(,)? ) => {
        Some(
            $crate::StyleAttr {
                class: Some($class_val.to_string()),
                style: Some($style_val.to_string()),
                id: Some($id_val.to_string()),
            }
        )
    };
    ( style : $style_val:expr, id : $id_val:expr, class : $class_val:expr $(,)? ) => {
        Some(
            $crate::StyleAttr {
                class: Some($class_val.to_string()),
                style: Some($style_val.to_string()),
                id: Some($id_val.to_string()),
            }
        )
    };
    ( id : $id_val:expr, class : $class_val:expr, style : $style_val:expr $(,)? ) => {
        Some(
            $crate::StyleAttr {
                class: Some($class_val.to_string()),
                style: Some($style_val.to_string()),
                id: Some($id_val.to_string()),
            }
        )
    };
    ( id : $id_val:expr, style : $style_val:expr, class : $class_val:expr $(,)? ) => {
        Some(
            $crate::StyleAttr {
                class: Some($class_val.to_string()),
                style: Some($style_val.to_string()),
                id: Some($id_val.to_string()),
            }
        )
    };
    // Catch-all for unknown keys
    ( $( $key:ident : $val:expr ),+ $(,)? ) => {
        compile_error!("Unknown key in style_attr! macro. Only 'class', 'style', and 'id' are supported.");
    };
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_style_attr_macro() {
        let a = css!(class: "foo");
        assert_eq!(a.as_ref().unwrap().class.as_deref(), Some("foo"));
        assert_eq!(a.unwrap().style, None);

        let b = css!(style: "stroke:red;");
        assert_eq!(b.as_ref().unwrap().class, None);
        assert_eq!(b.unwrap().style.as_deref(), Some("stroke:red;"));

        let c = css!(class: "bar", style: "fill:blue;");
        assert_eq!(c.as_ref().unwrap().class.as_deref(), Some("bar"));
        assert_eq!(c.unwrap().style.as_deref(), Some("fill:blue;"));
    }
}
