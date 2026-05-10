//! Convert a SpectraShop `.txt` export to the `colorimetry-io` JSON format.
//!
//! # Usage
//!
//! ```text
//! spectrashop_to_json [-c <copyright>] <input.txt> [output.json]
//! ```
//!
//! If no output path is given the JSON is written to stdout.

use colorimetry_io::SpectrumFile;
use std::{env, fs, path::PathBuf, process};

fn usage(prog: &str) -> ! {
    eprintln!("Usage: {prog} [-c <copyright>] <input.txt> [output.json]");
    eprintln!();
    eprintln!("Options:");
    eprintln!("  -c <string>   Copyright notice embedded in every spectrum record.");
    eprintln!();
    eprintln!("If no output file is given the JSON is written to stdout.");
    process::exit(1);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let prog = &args[0];

    let mut copyright: Option<String> = None;
    let mut positional: Vec<String> = Vec::new();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-c" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("error: -c requires an argument");
                    usage(prog);
                }
                copyright = Some(args[i].clone());
            }
            arg if arg.starts_with('-') => {
                eprintln!("error: unknown option '{arg}'");
                usage(prog);
            }
            _ => positional.push(args[i].clone()),
        }
        i += 1;
    }

    if positional.is_empty() {
        usage(prog);
    }

    let input = PathBuf::from(&positional[0]);
    let output: Option<PathBuf> = positional.get(1).map(PathBuf::from);

    let mut file = SpectrumFile::from_spectrashop_path(&input).unwrap_or_else(|e| {
        eprintln!("error: {e}");
        process::exit(1);
    });

    if let Some(ref cr) = copyright {
        match &mut file {
            SpectrumFile::Single { spectrum, .. } => {
                spectrum.metadata.copyright = Some(cr.clone());
            }
            SpectrumFile::Batch { spectra, .. } => {
                for sp in spectra.iter_mut() {
                    sp.metadata.copyright = Some(cr.clone());
                }
            }
        }
    }

    let json = serde_json::to_string_pretty(&file).unwrap_or_else(|e| {
        eprintln!("error serialising to JSON: {e}");
        process::exit(1);
    });

    match output {
        Some(ref path) => {
            fs::write(path, &json).unwrap_or_else(|e| {
                eprintln!("error writing to {}: {e}", path.display());
                process::exit(1);
            });
            eprintln!("Written to {}", path.display());
        }
        None => print!("{json}"),
    }
}
