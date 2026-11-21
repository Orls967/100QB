//! Lightweight 17-digit precision JSON writer with zero external crates.

use std::fmt::Write;
use std::fs;
use std::path::Path;

#[derive(Clone, Debug)]
pub enum JsonVal {
    Null,
    Bool(bool),
    Int(i64),
    UInt(usize),
    Float(f64),
    Str(String),
    Arr(Vec<JsonVal>),
    Obj(Vec<(&'static str, JsonVal)>),
}

impl From<bool> for JsonVal {
    fn from(v: bool) -> Self {
        JsonVal::Bool(v)
    }
}
impl From<usize> for JsonVal {
    fn from(v: usize) -> Self {
        JsonVal::UInt(v)
    }
}
impl From<i64> for JsonVal {
    fn from(v: i64) -> Self {
        JsonVal::Int(v)
    }
}
impl From<f64> for JsonVal {
    fn from(v: f64) -> Self {
        JsonVal::Float(v)
    }
}
impl From<&str> for JsonVal {
    fn from(v: &str) -> Self {
        JsonVal::Str(v.to_string())
    }
}
impl From<String> for JsonVal {
    fn from(v: String) -> Self {
        JsonVal::Str(v)
    }
}
impl<T: Into<JsonVal>> From<Vec<T>> for JsonVal {
    fn from(v: Vec<T>) -> Self {
        JsonVal::Arr(v.into_iter().map(Into::into).collect())
    }
}

fn escape_str(s: &str, out: &mut String) {
    out.push('"');
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

fn write_val(val: &JsonVal, indent: usize, out: &mut String) {
    match val {
        JsonVal::Null => out.push_str("null"),
        JsonVal::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
        JsonVal::Int(i) => {
            let _ = write!(out, "{i}");
        }
        JsonVal::UInt(u) => {
            let _ = write!(out, "{u}");
        }
        JsonVal::Float(f) => {
            if f.is_nan() || f.is_infinite() {
                out.push_str("null");
            } else {
                let _ = write!(out, "{f:.17}");
            }
        }
        JsonVal::Str(s) => escape_str(s, out),
        JsonVal::Arr(arr) => {
            if arr.is_empty() {
                out.push_str("[]");
                return;
            }
            let primitive = arr
                .iter()
                .all(|v| !matches!(v, JsonVal::Obj(_) | JsonVal::Arr(_)));
            if primitive && arr.len() <= 16 {
                out.push('[');
                for (i, item) in arr.iter().enumerate() {
                    if i > 0 {
                        out.push_str(", ");
                    }
                    write_val(item, 0, out);
                }
                out.push(']');
                return;
            }
            out.push_str("[\n");
            for (i, item) in arr.iter().enumerate() {
                for _ in 0..(indent + 2) {
                    out.push(' ');
                }
                write_val(item, indent + 2, out);
                if i + 1 < arr.len() {
                    out.push(',');
                }
                out.push('\n');
            }
            for _ in 0..indent {
                out.push(' ');
            }
            out.push(']');
        }
        JsonVal::Obj(fields) => {
            if fields.is_empty() {
                out.push_str("{}");
                return;
            }
            out.push_str("{\n");
            for (i, (k, v)) in fields.iter().enumerate() {
                for _ in 0..(indent + 2) {
                    out.push(' ');
                }
                escape_str(k, out);
                out.push_str(": ");
                write_val(v, indent + 2, out);
                if i + 1 < fields.len() {
                    out.push(',');
                }
                out.push('\n');
            }
            for _ in 0..indent {
                out.push(' ');
            }
            out.push('}');
        }
    }
}

pub fn write_json_file(path: &Path, val: &JsonVal) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut out = String::with_capacity(65536);
    write_val(val, 0, &mut out);
    out.push('\n');
    fs::write(path, out)
}
