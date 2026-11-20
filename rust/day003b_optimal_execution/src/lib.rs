#![forbid(unsafe_code)]
//! Eksekusi optimal dalam Rust murni: Almgren-Chriss, dampak transien Obizhaeva-Wang, dan DP likuiditas stokastik.

pub mod ac;
pub mod dp_liquidity;
pub mod json;
pub mod linalg;
pub mod parallel;
pub mod rng;
pub mod stats;
pub mod svg;
pub mod transient;
