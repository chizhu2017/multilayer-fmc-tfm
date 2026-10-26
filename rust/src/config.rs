use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct Config {
    pub nx: usize,
    pub h1: usize,
    pub h2: usize,
    pub h3: usize,
    pub v1: f32,
    pub v2: f32,
    pub v3: f32,
    pub x_loc: usize,
    pub gap: usize,
    pub cut_point: f64,
    pub t0_fine: f32,
    pub color_ratio: f32,
    pub show_interface: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            nx: 600,
            h1: 350,
            h2: 200,
            h3: 50,
            v1: 4200.,
            v2: 2800.,
            v3: 4600.,
            x_loc: 50,
            gap: 50,
            cut_point: 20.,
            t0_fine: 0.,
            color_ratio: 0.7,
            show_interface: false,
        }
    }
}

impl Config {
    pub fn nz(&self) -> usize {
        self.h1 + self.h2 + self.h3
    }
    pub fn validate(&self) -> Result<(), String> {
        if self.nx == 0
            || self.h1 < 2
            || self.h2 == 0
            || self.h3 == 0
            || self.x_loc == 0
            || self.gap == 0
        {
            return Err("positive integer dimensions/probe coordinates required; h1 >= 2".into());
        }
        if [self.v1, self.v2, self.v3, self.color_ratio]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.)
            || !self.cut_point.is_finite()
            || !self.t0_fine.is_finite()
        {
            return Err("invalid speed, display ratio or time correction".into());
        }
        let nz = self
            .h1
            .checked_add(self.h2)
            .and_then(|n| n.checked_add(self.h3))
            .ok_or("depth overflow")?;
        nz.checked_mul(self.nx).ok_or("image size overflow")?;
        Ok(())
    }
    pub fn probes(&self, n: usize) -> Result<Vec<f32>, String> {
        (0..n)
            .map(|i| {
                let x = i
                    .checked_mul(self.gap)
                    .and_then(|d| d.checked_add(self.x_loc))
                    .ok_or("probe overflow")?;
                if x > self.nx {
                    Err("probe coordinates must lie in 1..nx mm".into())
                } else {
                    Ok(x as f32)
                }
            })
            .collect()
    }
}
