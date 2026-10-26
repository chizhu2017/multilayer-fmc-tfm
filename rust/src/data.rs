use std::{
    fs,
    path::{Path, PathBuf},
};

pub fn read_matrix(path: &Path) -> Result<Vec<Vec<f32>>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut rows = Vec::new();
    let mut width = 0;
    for line in text.trim_start_matches('\u{feff}').lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('%') || line.starts_with('#') {
            continue;
        }
        let row: Vec<f32> = line
            .replace(',', " ")
            .split_whitespace()
            .map(|s| {
                // Parse through f64, matching MATLAB readmatrix -> single.
                let v = s
                    .parse::<f64>()
                    .map_err(|_| format!("{}: invalid number {s}", path.display()))?;
                let f = v as f32;
                if f.is_finite() {
                    Ok(f)
                } else {
                    Err(format!("{}: nonfinite float32", path.display()))
                }
            })
            .collect::<Result<_, String>>()?;
        if width == 0 {
            width = row.len();
        }
        if width != row.len() {
            return Err("inconsistent row widths".into());
        }
        rows.push(row);
    }
    if rows.len() < 2 || width < 2 {
        return Err("need >=2 samples and >=1 receiver".into());
    }
    Ok(rows)
}

pub struct Dataset {
    pub n: usize,
    pub nt: usize,
    pub dt: f32,
    time: Vec<f32>,
    paths: Vec<PathBuf>,
}

impl Dataset {
    pub fn open(directory: &Path) -> Result<Self, String> {
        let mut found: Vec<PathBuf> = fs::read_dir(directory)
            .map_err(|e| e.to_string())?
            .map(|e| e.map(|v| v.path()).map_err(|e| e.to_string()))
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .filter(|p| p.extension().is_some_and(|e| e.eq_ignore_ascii_case("txt")))
            .collect();
        let n = found.len();
        if n == 0 {
            return Err("no TXT files found".into());
        }
        let paths: Vec<_> = (1..=n)
            .map(|i| directory.join(format!("{i}.txt")))
            .collect();
        found.sort();
        let mut expected = paths.clone();
        expected.sort();
        if found != expected {
            return Err("input must contain exactly 1.txt through N.txt".into());
        }
        let first = read_matrix(&paths[0])?;
        if first[0].len() != n + 1 {
            return Err("file count must match receiver count".into());
        }
        let time: Vec<_> = first.iter().map(|r| r[0]).collect();
        let delta: Vec<_> = time.windows(2).map(|w| w[1] - w[0]).collect();
        if delta.iter().any(|d| *d <= 0.) {
            return Err("time must strictly increase".into());
        }
        // f64 reduction avoids sequential f32 sum error; final mean stays f32.
        let mean =
            (delta.iter().map(|d| (1f32 / d) as f64).sum::<f64>() / delta.len() as f64) as f32;
        let fs = mean.round();
        let dt = 1f32 / fs;
        if !dt.is_finite() || dt <= 0. {
            return Err("invalid sampling interval".into());
        }
        if delta.iter().any(|d| (*d - dt).abs() > 0.01f32 * dt) {
            return Err("time step must be uniform within 1% of dt".into());
        }
        Ok(Self {
            n,
            nt: first.len(),
            dt,
            time,
            paths,
        })
    }
    pub fn read(&self, k: usize) -> Result<Vec<Vec<f32>>, String> {
        let matrix = read_matrix(&self.paths[k])?;
        if matrix.len() != self.nt || matrix[0].len() != self.n + 1 {
            return Err(format!(
                "{}: inconsistent dimensions",
                self.paths[k].display()
            ));
        }
        if matrix
            .iter()
            .zip(&self.time)
            .any(|(r, t)| (r[0] - t).abs() > 0.01f32 * self.dt)
        {
            return Err("inconsistent time axis".into());
        }
        // Channel-major traces, no full FMC cube.
        Ok((0..self.n)
            .map(|j| matrix.iter().map(|r| r[j + 1]).collect())
            .collect())
    }
    pub fn peaks(&self) -> Result<Vec<f32>, String> {
        let mut peaks = Vec::with_capacity(self.n * self.n);
        for k in 0..self.n {
            for trace in self.read(k)? {
                let mut best = 0;
                for i in 1..trace.len() {
                    if trace[i] > trace[best] {
                        best = i;
                    }
                }
                peaks.push((best + 1) as f32);
            }
        }
        Ok(peaks)
    }
}
