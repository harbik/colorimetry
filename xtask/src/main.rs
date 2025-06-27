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
    /// Builds the wasm-bindgen output using wasm-pack
    Wasm,
    /// Generates data structures which are expensive to generate,
    /// and which are used in the library.
    Gen,
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
        Commands::Gen => {
            todo!()
        }
    }
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

fn check_or_force_rdme() {
    let status = Command::new("cargo")
        .args(["rdme", "--check"])
        .status()
        .expect("failed to run cargo rdme --check");

    if status.success() {
        println!("✅ README is up to date");
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
    }
    Ok(())
}
