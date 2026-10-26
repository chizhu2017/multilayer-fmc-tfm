use crate::config::Config;

pub struct Cache {
    pub dist: Vec<f32>,
    pub t2: Vec<f32>,
    pub t3: Vec<f32>,
}

pub fn tables(c: &Config, max_dx: usize) -> (Vec<f32>, Vec<f32>) {
    let width = max_dx + 1;
    let mut t12 = vec![0f32; width];
    let mut lut2 = vec![0f32; c.h2 * width];
    let mut lut3 = vec![0f32; (c.h3 + 1) * width];
    for d in 0..width {
        let mut best = f32::INFINITY;
        for u in 0..=d {
            let a = u as f32;
            let b = (d - u) as f32;
            let t = (a * a + (c.h1 as f32).powi(2)).sqrt() / c.v1
                + (b * b + (c.h2 as f32).powi(2)).sqrt() / c.v2;
            best = best.min(t);
        }
        t12[d] = best;
        for z in 0..c.h2 {
            let mut best = f32::INFINITY;
            for u in 0..=d {
                let a = u as f32;
                let b = (d - u) as f32;
                let t = (a * a + (c.h1 as f32).powi(2)).sqrt() / c.v1
                    + (b * b + (z as f32).powi(2)).sqrt() / c.v2;
                best = best.min(t);
            }
            lut2[z * width + d] = best;
        }
        for z in 0..=c.h3 {
            let mut best = f32::INFINITY;
            for (u, &pre) in t12.iter().enumerate().take(d + 1) {
                let b = (d - u) as f32;
                best = best.min(pre + (b * b + (z as f32).powi(2)).sqrt() / c.v3);
            }
            lut3[z * width + d] = best;
        }
    }
    (lut2, lut3)
}

pub fn cache(c: &Config, probes: &[f32]) -> Cache {
    let max_dx = probes
        .iter()
        .map(|p| (p - 1.).abs().max((p - c.nx as f32).abs()) as usize)
        .max()
        .unwrap();
    let (lut2, lut3) = tables(c, max_dx);
    let width = max_dx + 1;
    let mut dist = Vec::with_capacity(probes.len() * (c.h1 - 1) * c.nx);
    let mut t2 = Vec::with_capacity(probes.len() * c.h2 * c.nx);
    let mut t3 = Vec::with_capacity(probes.len() * (c.h3 + 1) * c.nx);
    for &p in probes {
        for z in 1..c.h1 {
            for x in 1..=c.nx {
                let dx = x as f32 - p;
                dist.push((dx * dx + (z as f32).powi(2)).sqrt());
            }
        }
        for z in 0..c.h2 {
            for x in 1..=c.nx {
                t2.push(lut2[z * width + (x as f32 - p).abs() as usize]);
            }
        }
        for z in 0..=c.h3 {
            for x in 1..=c.nx {
                t3.push(lut3[z * width + (x as f32 - p).abs() as usize]);
            }
        }
    }
    Cache { dist, t2, t3 }
}
