use clap::{Parser, Subcommand};

/// Top-level CLI struct.
#[derive(Parser, Debug)]
#[command(name = "color")]
#[command(about = "A commandline tool for colorimetric calculations", long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Convert spectral data of measured color samples to various color models.
    Sample {
        #[arg()] // required
        file: String,

        #[arg(short, long, default_value = "xyz")]
        model: String,

        // Optional illuminant to use for conversion
        #[arg(short, long, default_value = "D65")]
        illuminant: String,
    },

    /// Smooth spectral data using a Gaussian filter.
    /// The input file should contain spectral data in CSV format with two columns:

    /// - Wavelength (in nm), no restrictions on the number of spectral values.
    /// - Spectral Values (no restriction in range)
    ///
    /// The output will be the smoothed spectral data, over a range from 380 to 780 nm.
    Smooth {
        #[arg()] // required
        /// Input file containing spectral data in CSV format
        /// The file should have two columns: Wavelength (in nm) and Spectral Values.
        /// There are no restrictions on the number of spectral values, or their wavelength range.
        file: String,

        #[arg()]
        /// Optional output file for the smoothed data.
        /// If not specified, it will write the data to standard output.
        output: Option<String>,

        #[arg(short, long, default_value = "5.0")]
        /// FWHM maximum of the Gaussian smoothing filter.
        width: f64,

        #[arg(short, long)]
        /// Create a plot of the smoothed data in form of an SVG, if specified.
        /// If not specified, no plot is created.
        plot: Option<String>,
    },

    /// Show available standard illuminants and their key characteristics.
    Illuminants,

    /// Generate a general XY chromaticity diagram with a selection of graphical elements,
    /// such as the Planckian locus, RGB color space gamut, and more.
    /// The output is an SVG file that can be viewed in a web browser or edited with an SVG editor such as Inkscape.
    XYPlot {
        //
        #[arg(short, long, default_value = "cie1931")]
        observer: String,

        /// Output file for the XY plot, writes to standard output if not specified
        #[arg(short, long)]
        file: Option<String>,
        /// Input file containing spectral data in CSV format

        /// Shows the sRGB gamut in the XY chromaticity diagram, if specified.
        #[arg(short, long)]
        rgbspace: Option<String>,

        /// Shows the the Planckian Locus in the XY chromaticity diagram, if specified.
        #[arg(short, long)]
        planck: bool,
    },
}

/// Public function to parse CLI arguments.
pub fn parse_args() -> Cli {
    Cli::parse()
}
