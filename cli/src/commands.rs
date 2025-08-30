// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2025, Harbers Bik LLC

use crate::{args::Commands, utils::colorvalue};

pub fn dispatch(cmd: Commands) {
    match cmd {
        Commands::Sample {
            file,
            model,
            illuminant,
        } => {
            println!(
                "Converting Color Spectrum in '{file}' to '{model}'-values using '{illuminant}' as illumant ...");
            if let Err(e) = colorvalue(file, illuminant, model) {
                eprintln!("Error processing file: {e}");
            }
        }
        Commands::Illuminants => {
            println!("Available illuminants: D65, A, C...");
            println!("This feature is not yet implemented. Please check back later.");
        }

        Commands::Smooth {
            file,
            output,
            width,
            plot,
        } => {
            println!(
                "{file} {output} {width} {plot} Smoothing spectral data...",
                output = output.as_deref().unwrap_or("stdout"),
                plot = plot.as_deref().unwrap_or("no plot")
            );
            // Placeholder for smoothing logic
            println!("This feature is not yet implemented. Please check back later.");
        }

        Commands::XYPlot {
            observer,
            file,
            rgbspace,
            planck,
        } => {
            println!("Generating XY plot with observer: {observer}, output: {output}, RGB space: {space}, Planckian locus: {planck}",
                output = file.as_deref().unwrap_or("default"), 
                space = rgbspace.as_deref().unwrap_or("default"), 
            );
            println!("This feature is not yet implemented. Please check back later.");
        }
    }
}
