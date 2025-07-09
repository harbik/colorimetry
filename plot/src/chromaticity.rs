
#[cfg(test)]
mod tests {
    use crate::{canvas::Canvas, chart::Chart};

    use colorimetry::observer::Observer;
    

    #[test]
    fn test_draw_spectral_locus() {
        let observer = Observer::default();

        let canvas = Canvas::new(850, 950);

        let mut xy_chart = Chart::new(
            &canvas,
            "xy_chart",
            [50, 50, 750, 850],
            [[0.0, 0.75], [0.0, 0.85]],
            true,
            Some("fill: #CCCCCC"),
        );

        // Draw the spectral locus
        xy_chart.path("locus", observer.spectral_locus(), false, Some("stroke: black; stroke-width: 1; fill: #DDDDDD;"));
        xy_chart.draw_grid(0.01, 0.01, "stroke: #AAAAAA; stroke-width: 0.5;");
        xy_chart.draw_grid(0.1, 0.1, "stroke: #AAAAAA; stroke-width: 1.0;");


        xy_chart.render();
        /*
        chart.draw_labels();
        chart.draw_grid();
         */
        canvas.save("tmp/test_chromaticity_diagram.svg").unwrap();
    }
}
