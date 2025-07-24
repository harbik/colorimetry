
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StyleAttr {
    pub class: Option<String>,
    pub style: Option<String>,
    pub id: Option<String>,
}

impl StyleAttr {
    /// Converts to a single string like: `class="foo" style="stroke:red;"`
    pub fn assign<T>(&self, node: &mut T)
    where T: svg::Node,
    {
        if let Some(class) = self.class.clone() {
            node.assign("class", class);
        }
        if let Some(style) = self.style.clone() {
            node.assign("style", style);
        }
        // If neither class nor style is set, assign a default class
        // This is useful for elements that need a default style as else they are hidden on the plot
        // In this libarary we use "default" as a fallback class, showing the otherwise hidden elements
        // with a bright green color ('chartreuse') in the plot.
        if self.class.is_none() && self.style.is_none() {
            node.assign("class", "default");
        }
    }

    pub fn id(&self) -> Option<&str> {
        self.id.as_deref()
    }
}



/// Creates a [`StyleAttr`] struct for use in HTML or SVG element attributes.
///
/// This macro supports multiple invocation patterns:
///
/// ```rust
/// use colorimetry_plot::{style_attr, StyleAttr};
/// let style = style_attr!(class: "button", style: "fill:red;");
/// assert_eq!(style.class.as_deref(), Some("button"));
/// assert_eq!(style.style.as_deref(), Some("fill:red;"));
/// ```
///
/// This is useful for freeform style attributes that don’t follow key-value
/// syntax or when you want to include semicolons without quoting.
///
/// # Returns
///
/// A [`Style`] struct with either the `class`, `style`, or both fields set.
///
/// # Notes
///
/// - If both `class:` and `style:` keys are used, both fields are set.
/// - If the macro receives a single identifier (e.g. `highlight`), it sets only the `class`.
/// - Any other token sequence (e.g. semicolon-delimited) is treated as a full `style` string.
///
/// # Panics
///
/// If an unknown key is used in key-value form (e.g. `color: "red"`), this will cause a compile-time error.
///
/// # Examples
///
/// ```rust
/// use colorimetry_plot::{style_attr, StyleAttr};
/// let a = style_attr!(); // No style specfied
/// let b = style_attr!(class: "foo");
/// let c = style_attr!(style: "stroke:red;");
/// let d = style_attr!(class: "bar", style: "fill:blue;");
/// let e = style_attr!(id: "line3");
/// ```
#[macro_export]
macro_rules! style_attr {
    // No arguments, defaults to empty fields
    () => {
        $crate::StyleAttr {
            class: None,
            style: None,
            id: None,
        }
    };

    // Only class
    ( $( class : $class_val:expr ),+ $(,)? ) => {
        $crate::StyleAttr {
            class: Some([$($class_val),+].join(" ")),
            style: None,
            id: None,
        }
    };
    // Only style
    ( $( style : $style_val:expr ),+ $(,)? ) => {
        $crate::StyleAttr {
            class: None,
            style: Some([$($style_val),+].join(" ")),
            id: None,
        }
    };
    // Only id
    ( $( id : $id_val:expr ),+ $(,)? ) => {
        $crate::StyleAttr {
            class: None,
            style: None,
            id: Some([$($id_val),+].join(" ")),
        }
    };
    // class + style
    ( class : $class_val:expr, style : $style_val:expr $(,)? ) => {
        $crate::StyleAttr {
            class: Some($class_val.to_string()),
            style: Some($style_val.to_string()),
            id: None,
        }
    };
    ( style : $style_val:expr, class : $class_val:expr $(,)? ) => {
        $crate::StyleAttr {
            class: Some($class_val.to_string()),
            style: Some($style_val.to_string()),
            id: None,
        }
    };
    // class + id
    ( class : $class_val:expr, id : $id_val:expr $(,)? ) => {
        $crate::StyleAttr {
            class: Some($class_val.to_string()),
            style: None,
            id: Some($id_val.to_string()),
        }
    };
    ( id : $id_val:expr, class : $class_val:expr $(,)? ) => {
        $crate::StyleAttr {
            class: Some($class_val.to_string()),
            style: None,
            id: Some($id_val.to_string()),
        }
    };
    // style + id
    ( style : $style_val:expr, id : $id_val:expr $(,)? ) => {
        $crate::StyleAttr {
            class: None,
            style: Some($style_val.to_string()),
            id: Some($id_val.to_string()),
        }
    };
    ( id : $id_val:expr, style : $style_val:expr $(,)? ) => {
        $crate::StyleAttr {
            class: None,
            style: Some($style_val.to_string()),
            id: Some($id_val.to_string()),
        }
    };
    // class + style + id (any order)
    ( class : $class_val:expr, style : $style_val:expr, id : $id_val:expr $(,)? ) => {
        $crate::StyleAttr {
            class: Some($class_val.to_string()),
            style: Some($style_val.to_string()),
            id: Some($id_val.to_string()),
        }
    };
    ( class : $class_val:expr, id : $id_val:expr, style : $style_val:expr $(,)? ) => {
        $crate::StyleAttr {
            class: Some($class_val.to_string()),
            style: Some($style_val.to_string()),
            id: Some($id_val.to_string()),
        }
    };
    ( style : $style_val:expr, class : $class_val:expr, id : $id_val:expr $(,)? ) => {
        $crate::StyleAttr {
            class: Some($class_val.to_string()),
            style: Some($style_val.to_string()),
            id: Some($id_val.to_string()),
        }
    };
    ( style : $style_val:expr, id : $id_val:expr, class : $class_val:expr $(,)? ) => {
        $crate::StyleAttr {
            class: Some($class_val.to_string()),
            style: Some($style_val.to_string()),
            id: Some($id_val.to_string()),
        }
    };
    ( id : $id_val:expr, class : $class_val:expr, style : $style_val:expr $(,)? ) => {
        $crate::StyleAttr {
            class: Some($class_val.to_string()),
            style: Some($style_val.to_string()),
            id: Some($id_val.to_string()),
        }
    };
    ( id : $id_val:expr, style : $style_val:expr, class : $class_val:expr $(,)? ) => {
        $crate::StyleAttr {
            class: Some($class_val.to_string()),
            style: Some($style_val.to_string()),
            id: Some($id_val.to_string()),
        }
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
        let a = style_attr!(class: "foo");
        assert_eq!(a.class.as_deref(), Some("foo"));
        assert_eq!(a.style, None);

        let b = style_attr!(style: "stroke:red;");
        assert_eq!(b.class, None);
        assert_eq!(b.style.as_deref(), Some("stroke:red;"));

        let c = style_attr!(class: "bar", style: "fill:blue;");
        assert_eq!(c.class.as_deref(), Some("bar"));
        assert_eq!(c.style.as_deref(), Some("fill:blue;"));

    }
}