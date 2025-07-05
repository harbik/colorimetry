use chrono::Utc;
use colorimetry::observer::Observer;
use colorimetry::rgb::RgbSpace;
use handlebars::Handlebars;
use serde::Serialize;
use strum::IntoEnumIterator;

#[derive(Serialize)]
struct Entry {
    observer: String,
    rgbspace: String,
    identifier: String,
    matrix: String,
    is_inv: bool,
    supplemental_observers: bool,
}

#[derive(Serialize)]
struct TemplateContext {
    date: String,
    entries: Vec<Entry>,
}

use std::fs;
use std::io;
use std::path::Path;

pub fn matrices() -> io::Result<()> {
    let contents = render_rgb_template();
    let dir = Path::new("src/observer");
    fs::create_dir_all(dir)?;

    let file_path = dir.join("rgbxyz.rs");
    fs::write(file_path, contents)?;

    Ok(())
}
pub fn render_rgb_template() -> String {
    let mut entries = Vec::new();

    for observer in Observer::iter() {
        let supplemental_observers = observer != Observer::Cie1931;
        for rgbspace in RgbSpace::iter() {
            let m_rgb2xyz = observer.calc_rgb2xyz_matrix(rgbspace);
            let m_xyz2rgb = observer.calc_xyz2rgb_matrix(rgbspace);
            let observer = format!("{:?}", observer);
            let rgbspace = format!("{:?}", rgbspace);
            let matrix = format!("{:?}", m_rgb2xyz);
            let matrix_inv = format!("{:?}", m_xyz2rgb);
            let identifier_inv = format!(
                "XYZ2RGB_{}_{}",
                observer.to_uppercase(),
                rgbspace.to_uppercase()
            );
            let identifier = format!(
                "RGB2XYZ_{}_{}",
                observer.to_uppercase(),
                rgbspace.to_uppercase()
            );

            entries.push(Entry {
                matrix,
                identifier,
                observer: observer.clone(),
                rgbspace: rgbspace.clone(),
                is_inv: false,
                supplemental_observers,
            });

            entries.push(Entry {
                matrix: matrix_inv,
                identifier: identifier_inv,
                observer,
                rgbspace,
                is_inv: true,
                supplemental_observers,
            });
        }
    }

    let context = TemplateContext {
        date: Utc::now().format("%Y-%m-%d").to_string(),
        entries,
    };

    let mut handlebars = Handlebars::new();
    // 👇 Reads your rgb.hbs file at runtime
    handlebars
        .register_template_file("rgbxyz", "xtask/templates/rgbxyz.rs.hbs")
        .expect("Failed to read rgb.hbs");

    handlebars
        .render("rgbxyz", &context)
        .expect("Failed to render template")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_print_rgb_observers() {
        render_rgb_template();
    }
}
