
use std::sync::OnceLock;

use wasm_bindgen::prelude::wasm_bindgen;


use crate::Spectrum;

static CTS: OnceLock<[Spectrum;15]> = OnceLock::new();
#[wasm_bindgen]
pub struct CRI([f64;15]);

impl CRI {
    pub async fn init() {
        // use reqwest library
        todo!()
    }
}

// JS-WASM Interface code

#[cfg(target_arch="wasm32")] 
#[wasm_bindgen]
impl CRI {

    pub async fn init_js() {
        // use fetch WEB API
        if CTS.get().is_none() {
            
            // get cts's from harbik.github.io
        };
    }

}
