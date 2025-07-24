
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StyleAttr {
    pub class: Option<String>,
    pub style: Option<String>,
}

impl StyleAttr {
    /// Converts to a single string like: `class="foo" style="stroke:red;"`
    pub fn node_assign<T>(&self, mut node: T)
    where T: svg::Node,
    {
        if let Some(class) = self.class.clone() {
            node.assign("class", class);
        }
        if let Some(style) = self.style.clone() {
            node.assign("style", style);
        }
    }
}



/// Creates a [`StyleAttr`] struct for use in HTML or SVG element attributes.
///
/// This macro supports multiple invocation patterns:
///
/// ### 1. Shorthand identifier for `class`
///
/// ```rust
/// use colorimetry-plot::style_attr;
/// let style = style_attr!(highlight);
/// assert_eq!(style.class.as_deref(), Some("highlight"));
/// assert_eq!(style.style, None);
/// ```
///
/// ### 2. Key-value pairs for `class` and/or `style`
///
/// ```rust
/// let style = style_attr!(class: "button", style: "fill:red;");
/// assert_eq!(style.class.as_deref(), Some("button"));
/// assert_eq!(style.style.as_deref(), Some("fill:red;"));
/// ```
///
/// ### 3. Raw style string as token tree
///
/// ```rust
/// let style = style_attr!(stroke:blue; fill:none;);
/// assert_eq!(style.class, None);
/// assert_eq!(style.style.as_deref(), Some("stroke:blue; fill:none;"));
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
/// let a = style_attr!(highlight);
/// let b = style_attr!(class: "foo");
/// let c = style_attr!(style: "stroke:red;");
/// let d = style_attr!(class: "bar", style: "fill:blue;");
/// let e = style_attr!(stroke:black; fill:none;);
/// ```
#[macro_export]
macro_rules! style_attr {
    // Case 1: shorthand identifier -> class
    ($ident:ident) => {
        $crate::StyleAttr {
            class: Some(stringify!($ident).to_string()),
            style: None,
        }
    };

    // Case 2: named key-value pairs for class/style only
    ( $( class : $class_val:expr ),+ $(,)? ) => {
        $crate::StyleAttr {
            class: Some([$($class_val),+].join(" ")),
            style: None,
        }
    };
    ( $( style : $style_val:expr ),+ $(,)? ) => {
        $crate::StyleAttr {
            class: None,
            style: Some([$($style_val),+].join(" ")),
        }
    };
    ( class : $class_val:expr, style : $style_val:expr $(,)? ) => {
        $crate::StyleAttr {
            class: Some($class_val.to_string()),
            style: Some($style_val.to_string()),
        }
    };
    ( style : $style_val:expr, class : $class_val:expr $(,)? ) => {
        $crate::StyleAttr {
            class: Some($class_val.to_string()),
            style: Some($style_val.to_string()),
        }
    };
    // Catch-all for unknown keys
    ( $( $key:ident : $val:expr ),+ $(,)? ) => {
        compile_error!("Unknown key in style_attr! macro. Only 'class' and 'style' are supported.");
    };

    // Case 3: fallback token tree (e.g. style_attr!(color:red; stroke:blue;))
    ($($style:tt)+) => {
        $crate::StyleAttr {
            class: None,
            style: Some(format!($($style)+)),
        }
    };
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_style_attr_macro() {
        let a = style_attr!(highlight);
        assert_eq!(a.class.as_deref(), Some("highlight"));
        assert_eq!(a.style, None);

        let b = style_attr!(class: "foo");
        assert_eq!(b.class.as_deref(), Some("foo"));
        assert_eq!(b.style, None);

        let c = style_attr!(style: "stroke:red;");
        assert_eq!(c.class, None);
        assert_eq!(c.style.as_deref(), Some("stroke:red;"));

        let d = style_attr!(class: "bar", style: "fill:blue;");
        assert_eq!(d.class.as_deref(), Some("bar"));
        assert_eq!(d.style.as_deref(), Some("fill:blue;"));

        let e = style_attr!("stroke:black; fill:none;");
        assert_eq!(e.class, None);
        assert_eq!(e.style.as_deref(), Some("stroke:black; fill:none;"));
    }
}