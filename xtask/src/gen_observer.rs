use colorimetry::observer::Observer;
use strum::IntoEnumIterator;

pub fn main() {
    // for each observer
    for observer in Observer::iter() {
        println!("{}", observer);
    }
}
