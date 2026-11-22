use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Clone, Debug)]
pub enum JsonVal {
    Null,
    Bool(bool),
    Int(i64),
    UInt(usize),
    Num(f64),
    Str(String),
    Arr(Vec<JsonVal>),
    Obj(Vec<(String, JsonVal)>),
}

impl From<bool> for JsonVal {
    fn from(v: bool) -> Self {
        JsonVal::Bool(v)
    }
}
impl From<i64> for JsonVal {
    fn from(v: i64) -> Self {
        JsonVal::Int(v)
    }
}
impl From<usize> for JsonVal {
    fn from(v: usize) -> Self {
        JsonVal::UInt(v)
    }
}
impl From<f64> for JsonVal {
    fn from(v: f64) -> Self {
        JsonVal::Num(v)
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

pub fn obj<const N: usize>(pairs: [(&str, JsonVal); N]) -> JsonVal {
    JsonVal::Obj(pairs.into_iter().map(|(k, v)| (k.to_string(), v)).collect())
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
                out.push_str(&format!("\\u{:04x}", c as u32));
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
        JsonVal::Int(i) => out.push_str(&i.to_string()),
        JsonVal::UInt(u) => out.push_str(&u.to_string()),
        JsonVal::Num(x) => {
            if x.is_finite() {
                out.push_str(&format!("{x:.17}"));
            } else {
                out.push_str("null");
            }
        }
        JsonVal::Str(s) => escape_str(s, out),
        JsonVal::Arr(arr) => {
            if arr.is_empty() {
                out.push_str("[]");
                return;
            }
            out.push_str("[\n");
            let pad = " ".repeat(indent + 2);
            for (idx, item) in arr.iter().enumerate() {
                out.push_str(&pad);
                write_val(item, indent + 2, out);
                if idx + 1 < arr.len() {
                    out.push(',');
                }
                out.push('\n');
            }
            out.push_str(&" ".repeat(indent));
            out.push(']');
        }
        JsonVal::Obj(pairs) => {
            if pairs.is_empty() {
                out.push_str("{}");
                return;
            }
            out.push_str("{\n");
            let pad = " ".repeat(indent + 2);
            for (idx, (k, v)) in pairs.iter().enumerate() {
                out.push_str(&pad);
                escape_str(k, out);
                out.push_str(": ");
                write_val(v, indent + 2, out);
                if idx + 1 < pairs.len() {
                    out.push(',');
                }
                out.push('\n');
            }
            out.push_str(&" ".repeat(indent));
            out.push('}');
        }
    }
}

pub fn write_json_file(path: &Path, val: &JsonVal) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut s = String::with_capacity(32_768);
    write_val(val, 0, &mut s);
    s.push('\n');
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    w.write_all(s.as_bytes())?;
    w.flush()
}
