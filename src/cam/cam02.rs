/// CIECAM02 Color Appearance Model
///
/// Implements the CIECAM02 appearance model, which predicts how we perceive color under different
/// viewing conditions.  It computes the JCh correlates:
/// - **J**: Lightness 
/// - **C**: Chroma (colorfulness)
/// - **h**: Hue angle  
/// 
/// It also performs chromatic adaptation, allowing colors to be accurately transformed between
/// different illuminants and viewing environments.
///
#[derive(Debug)]
pub struct CieCam02(CamJCh);
