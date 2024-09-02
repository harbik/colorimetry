
# Standard Illuminants

Many of the CIE Standard Illuminants are made available in this module through
`StdIlluminant` for data-defined illuminants and the
`Spectrum::cie_d_illuminant(cct: f64)` function for generic D-illuminant.

The `StdIlluminant` object gives access to all the CIE illuminants defined in this library, obtained from the datasets published in the
CIE15::2018 standard, downloaded from the [CIE Website](https://web.archive.org/web/20240314231650/https://cie.co.at/data-tables) August 2024.

As the data is compiled into the library, you can choose only to include the two basic illuminants `D65` and `D50` to limit the size of your
recent `F3_X` series included here,
executable by using the `--no-default-features` in the compiler, or the `default-features = false` in its dependence declaration of this crate in
the `cargo.toml` file of your application. The library uses the "cie-illuminants" feature flag to select the inclusion of these illuminants.

If you use this library in a JavaScript application, the default `colorimetry` package excludes all the features that allow fast load times for
lightweight applications.  Use the `colorimetry-all` package to use the library with all its features enabled.

For more detailed information on the CIE Standard Illuminant Datasets, see
[Standard illuminant](https://en.wikipedia.org/wiki/Standard_illuminant#White_points_of_standard_illuminants)
on Wikipedia.  Instead of a dash, use the `_` character to access these
illuminants by their name here, so use `StdIlluminant::LED_BH1` to use the
phosphor-converted Blue LED and Red LED standard illuminant.
The Fluorescent `F3_X` series is included here, with X ranging from 1 to 15.
