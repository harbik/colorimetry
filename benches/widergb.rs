use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

use colorimetry::widergb::WideRgb;

fn criterion_benchmark(c: &mut Criterion) {
    let wide_rgb = WideRgb::new(0.1, 0.3, 1.0, None, None);
    c.bench_function("WideRgb::compress in-gamut", |b| {
        b.iter(|| black_box(wide_rgb).compress())
    });

    let wide_rgb = WideRgb::new(-0.5, 0.1, 1.2, None, None);
    c.bench_function("WideRgb::compress out-of-gamut", |b| {
        b.iter(|| black_box(wide_rgb).compress())
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
