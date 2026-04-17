# Claude Code instructions for colorimetry

## Pre-release / CI checks

Whenever code is modified or tests are run, always execute the full xtask pipeline in this order and report any errors before marking work complete:

```
cargo xtask check   # fmt, clippy, build, rdme
cargo xtask test    # test --all-features and --no-default-features
cargo xtask doc     # rustdoc with --deny warnings
```

These three commands together are the definition of "the build is green" for this project.
