//! This is the main entry point for the colorimetry CLI tool.
//! It parses command-line arguments and dispatches to the appropriate command handler.
use colorimetry_cli::args::parse_args;
use colorimetry_cli::commands;

fn main() {
    // cargo run -- arguments data/test.csv --format lab
    let cli = parse_args();
    commands::dispatch(cli.command);
}
