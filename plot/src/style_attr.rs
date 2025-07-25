#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct StyleAttr {
    pub class: Option<String>,
    pub style: Option<String>,
    pub id: Option<String>,
}

impl StyleAttr {
    /// Converts to a single string like: `class="foo" style="stroke:red;"`
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

    pub fn id(&self) -> Option<&str> {
        self.id.as_deref()
    }
}

pub fn class(str: &str) -> Option<StyleAttr> {
    Some(StyleAttr {
        class: Some(str.to_string()),
        style: None,
        id: None,
    })
}

pub fn style(str: &str) -> Option<StyleAttr> {
    Some(StyleAttr {
        class: None,
        style: Some(str.to_string()),
        id: None,
    })
}

pub fn id(id_val: &str) -> Option<StyleAttr> {
    Some(StyleAttr {
        class: None,
        style: None,
        id: Some(id_val.to_string()),
    })
}

/// A [`StyleAttr`] struct with either the `class`, `style`, `id` or any combination of fields set.
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
