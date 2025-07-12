pub mod svgdoc;
pub mod chart;
pub mod spectrum;
pub mod transforms;
pub mod axis;
pub mod layer;

use std::sync::atomic::{AtomicUsize, Ordering};

static COUNTER: AtomicUsize = AtomicUsize::new(0);

pub fn new_id() -> String {
    format!("id{}", COUNTER.fetch_add(1, Ordering::Relaxed))
}

pub fn last_id() -> String {
    COUNTER.load(Ordering::Relaxed).to_string()
}
