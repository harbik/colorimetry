use clap::{Parser, Subcommand};
use std::process::Command;

#[derive(Parser)]
#[command(name = "xtask", about = "Custom build tasks")]
struct XtaskArgs {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Check formatting, clippy, and build
    Check,
    /// Run tests
    Test,
    /// Builds the documentation, and opens it
    Doc,
}

fn main() {
    let args = XtaskArgs::parse();

    match args.command {
        Commands::Check => {
            run("cargo", &["fmt", "--", "--check"]);
            run(
                "cargo",
                &[
                    "clippy",
                    "--all-targets",
                    "--all-features",
                    "--",
                    "-D",
                    "warnings",
                ],
            );
            run("cargo", &["check", "--all-targets"]);
            run("cargo", &["rdme", "--check"]);
        }
        Commands::Test => {
            run("cargo", &["test", "--all-features"]);
            run("cargo", &["test", "--no-default-features"]);
        }
        Commands::Doc => {
            run("cargo", &["doc", "--all-features", "--no-deps", "--open"]);
        }
    }
}

fn run(cmd: &str, args: &[&str]) {
    println!("$ {} {}", cmd, args.join(" "));
    let status = Command::new(cmd)
        .args(args)
        .status()
        .expect("failed to run command");
    if !status.success() {
        std::process::exit(status.code().unwrap_or(1));
    }
}
