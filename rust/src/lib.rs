pub mod config;
pub mod data;
pub mod travel;
use config::Config;
use data::Dataset;
use rustfft::{num_complex::Complex32, FftPlanner};
use std::path::Path;

pub fn interpolate(signal: &[f32], index: f32) -> f32 {
    let base = index.floor();
    if !base.is_finite() || base < 1. || base >= signal.len() as f32 {
        return 0.;
    }
    let i = base as usize - 1;
    let alpha = index - base;
    (1. - alpha) * signal[i] + alpha * signal[i + 1]
}

pub fn time_zero(peaks: &[f32], probes: &[f32], dt: f32, c: &Config) -> Result<f32, String> {
    let n = probes.len();
    let mut values = Vec::new();
    for k in 0..n {
        for j in 0..n {
            let direct = (probes[k] - probes[j]).abs() * 1e-3f32 / c.v1;
            let value = peaks[k * n + j] - direct / dt;
            if value.is_finite() && value > 0. {
                values.push(value);
            }
        }
    }
    if values.is_empty() {
        return Err("no valid time-zero candidates".into());
    }
    values.sort_by(|a, b| a.total_cmp(b));
    let m = values.len() / 2;
    let median = if values.len() % 2 == 1 {
        values[m]
    } else {
        (values[m - 1] + values[m]) / 2.
    };
    let t0 = median + c.t0_fine;
    if !t0.is_finite() {
        return Err("invalid time zero".into());
    }
    Ok(t0)
}

pub fn envelope(image: &[f32], nz: usize, nx: usize) -> Vec<f32> {
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(nz);
    let ifft = planner.plan_fft_inverse(nz);
    let mut column = vec![Complex32::new(0., 0.); nz];
    let mut out = vec![0.; image.len()];
    for x in 0..nx {
        for z in 0..nz {
            column[z] = Complex32::new(image[z * nx + x], 0.);
        }
        fft.process(&mut column);
        for (k, v) in column.iter_mut().enumerate() {
            let multiplier = if k == 0 || (nz.is_multiple_of(2) && k == nz / 2) {
                1.
            } else if k <= (nz - 1) / 2 {
                2.
            } else {
                0.
            };
            *v *= multiplier;
        }
        ifft.process(&mut column);
        for z in 0..nz {
            out[z * nx + x] = column[z].norm() / nz as f32;
        }
    }
    out
}

pub struct ResultImage {
    pub image: Vec<f32>,
    pub envelope: Vec<f32>,
    pub t0_index: f32,
    pub dt: f32,
    pub probes: Vec<f32>,
}

pub fn reconstruct(directory: &Path, c: &Config) -> Result<ResultImage, String> {
    c.validate()?;
    let ds = Dataset::open(directory)?;
    let probes = c.probes(ds.n)?;
    let peaks = ds.peaks()?;
    let t0 = time_zero(&peaks, &probes, ds.dt, c)?;
    let cache = travel::cache(c, &probes);
    let s1 = (c.h1 - 1) * c.nx;
    let s2 = c.h2 * c.nx;
    let s3 = (c.h3 + 1) * c.nx;
    let mut image = vec![0f32; c.nz() * c.nx];
    for k in 0..ds.n {
        let signals = ds.read(k)?;
        for j in 0..ds.n {
            let mut signal = signals[j].clone();
            let cut = (peaks[k * ds.n + j] as f64 + c.cut_point)
                .floor()
                .clamp(0., ds.nt as f64) as usize;
            signal[..cut].fill(0.);
            for (i, pixel) in image[..s1].iter_mut().enumerate() {
                let a = cache.dist[j * s1 + i];
                let b = cache.dist[k * s1 + i];
                let index = (a + b) * 1e-3f32 / c.v1 / ds.dt + t0;
                let z = (i / c.nx + 1) as f32;
                *pixel += interpolate(&signal, index) * ((z / a) * (z / b));
            }
            for (i, pixel) in image[s1..s1 + s2].iter_mut().enumerate() {
                let index = (cache.t2[j * s2 + i] + cache.t2[k * s2 + i]) * 1e-3f32 / ds.dt + t0;
                *pixel += interpolate(&signal, index);
            }
            for (i, pixel) in image[s1 + s2..].iter_mut().enumerate() {
                let index = (cache.t3[j * s3 + i] + cache.t3[k * s3 + i]) * 1e-3f32 / ds.dt + t0;
                *pixel += interpolate(&signal, index);
            }
        }
        eprintln!("FMC {}/{}", k + 1, ds.n);
    }
    let env = envelope(&image, c.nz(), c.nx);
    if env.iter().any(|v| !v.is_finite()) || !env.iter().any(|v| *v > 0.) {
        return Err("no finite nonzero image signal".into());
    }
    Ok(ResultImage {
        image,
        envelope: env,
        t0_index: t0,
        dt: ds.dt,
        probes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn interpolation_preserves_matlab_edges() {
        let s = [10., 20., 40.];
        assert_eq!(interpolate(&s, 1.), 10.);
        assert_eq!(interpolate(&s, 1.5), 15.);
        assert_eq!(interpolate(&s, 2.5), 30.);
        assert_eq!(interpolate(&s, 3.), 0.);
        assert_eq!(interpolate(&s, 0.99), 0.);
        assert_eq!(interpolate(&s, f32::NAN), 0.);
    }
    #[test]
    fn hilbert_cosine_odd_even_and_dc() {
        for n in [7, 8] {
            let input: Vec<_> = (0..n)
                .flat_map(|z| [(2. * std::f32::consts::PI * z as f32 / n as f32).cos(), 3.])
                .collect();
            let env = envelope(&input, n, 2);
            for z in 0..n {
                assert!((env[2 * z] - 1.).abs() < 2e-6);
                assert!((env[2 * z + 1] - 3.).abs() < 2e-6);
            }
        }
    }
    #[test]
    fn zero_offset_has_vertical_time() {
        let c = Config {
            h1: 12,
            h2: 10,
            h3: 9,
            ..Config::default()
        };
        let (a, b) = travel::tables(&c, 0);
        assert!((a[3] - (12. / c.v1 + 3. / c.v2)).abs() < 1e-8);
        assert!((b[5] - (12. / c.v1 + 10. / c.v2 + 5. / c.v3)).abs() < 1e-8);
    }
}
