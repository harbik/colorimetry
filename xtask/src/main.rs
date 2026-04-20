use clap::{Parser, Subcommand};
use std::process::Command;

mod gen_observer;
mod gen_rgbspace;
mod gen_rgbxyz;

/// Represents the command-line arguments for the xtask utility.
#[derive(Parser)]
#[command(
    name = "xtask",
    about = "Custom build tasks for the colorimetry library"
)]
struct XtaskArgs {
    #[command(subcommand)]
    command: Commands,
}

// Update Commands enum to include Gen with subcommand
#[derive(Subcommand)]
enum Commands {
    /// Check formatting, clippy, and build
    Check,
    /// Run tests
    Test,
    /// Builds the documentation, and opens it
    Doc,
    /// Builds the wasm-bindgen output using wasm-pack
    Wasm,
    /// Builds and publishes the WASM package to npm
    ///
    /// Requires `wasm-pack`, `wasm-opt`, and an active `npm login` session.
    /// The package version in pkg/package.json is set by wasm-pack from Cargo.toml.
    PublishWasm,
    /// Generates data structures which are expensive to generate,
    /// and which are used in the library.
    Gen {
        #[command(subcommand)]
        subcommand: GenCommands,
    },
}

#[derive(Subcommand)]
enum GenCommands {
    /// Generate RGB-XYZ matrices
    RgbTransforms,
    /// Generate Observer Data
    Observers,
    /// Generate RGB Space Data
    RgbSpaces,
}

impl Commands {
    fn handle(&self) {
        match self {
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
                check_or_force_rdme();
                run_env(
                    "cargo",
                    &["doc", "--all-features", "--no-deps"],
                    &[("RUSTDOCFLAGS", "--deny warnings")],
                )
                .expect("failed to run cargo doc");
                println!("✅ Documentation build complete");
            }
            Commands::Wasm => {
                build_wasm();
            }
            Commands::PublishWasm => {
                build_wasm();
                publish_wasm();
            }
            Commands::Gen { subcommand } => match subcommand {
                GenCommands::RgbTransforms => {
                    gen_rgbxyz::main().unwrap();
                    run("rustfmt", &["src/observer/rgbxyz.rs"]);
                }
                GenCommands::Observers => {
                    gen_observer::main();
                }
                GenCommands::RgbSpaces => {
                    gen_rgbspace::main().expect("Failed to generate rgbspace data");
                    // Use glob to format all .rs files in src/rgb/rgbspace
                    for entry in
                        glob::glob("src/rgb/rgbspace/*.rs").expect("Failed to read glob pattern")
                    {
                        match entry {
                            Ok(path) => {
                                run("rustfmt", &[path.to_str().unwrap()]);
                            }
                            Err(e) => eprintln!("Glob error: {}", e),
                        }
                    }
                }
            },
        }
    }
}

fn main() {
    let args = XtaskArgs::parse();
    args.command.handle();
}

fn build_wasm() {
    let status = Command::new("wasm-pack")
        .args(["build", "--target", "web", "--release", "--out-dir", "pkg"])
        .status()
        .expect("failed to run wasm-pack");

    if !status.success() {
        std::process::exit(status.code().unwrap_or(1));
    }

    println!("✅ wasm-bindgen build complete");
}

fn publish_wasm() {
    // Read the version from pkg/package.json and print it so the user can
    // confirm before the publish reaches npm.
    if let Ok(json) = std::fs::read_to_string("pkg/package.json") {
        if let Some(version) = json
            .lines()
            .find(|l| l.contains("\"version\""))
            .and_then(|l| {
                let (_, after) = l.split_once("\"version\"")?;
                let (_, after) = after.split_once('"')?;
                let (version, _) = after.split_once('"')?;
                Some(version.to_string())
            })
        {
            println!("📦 Publishing colorimetry@{version} to npm");
        }
    }

    let status = Command::new("wasm-pack")
        .args(["publish", "--pkg-dir", "pkg"])
        .status()
        .expect("failed to run wasm-pack publish");

    if !status.success() {
        std::process::exit(status.code().unwrap_or(1));
    }

    println!("✅ Published to npm");
}

fn check_or_force_rdme() {
    let status = Command::new("cargo")
        .args(["rdme", "--check"])
        .status()
        .expect("faid to run cargo rdme --check");

    if status.success() {
        println!("✅ README is up to date")
    } else {
        println!("⚠️ README is out of date, regenerating...");
        let force_status = Command::new("cargo")
            .args(["rdme", "--force"])
            .status()
            .expect("failed to run cargo rdme --force");

        if force_status.success() {
            println!("✅ README regenerated successfully");
        } else {
            eprintln!("❌ Failed to regenerate README");
            std::process::exit(force_status.code().unwrap_or(1));
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

fn run_env(cmd: &str, args: &[&str], env: &[(&str, &str)]) -> Result<(), std::io::Error> {
    let mut c = Command::new(cmd);
    c.args(args);
    for (k, v) in env {
        c.env(k, v);
    }
    let status = c.status()?;
    if !status.success() {
        std::process::exit(status.code().unwrap_or(1));
    } else {
        Ok(())
    }
}
