use multilayer_fmc_tfm::{config::Config, reconstruct};
use std::{
    env, fs,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};

fn csv(path: &Path, data: &[f32], nx: usize) -> Result<(), String> {
    let mut out = BufWriter::new(fs::File::create(path).map_err(|e| e.to_string())?);
    for row in data.chunks(nx) {
        for (i, v) in row.iter().enumerate() {
            if i > 0 {
                write!(out, ",").map_err(|e| e.to_string())?;
            }
            write!(out, "{v:.9e}").map_err(|e| e.to_string())?;
        }
        writeln!(out).map_err(|e| e.to_string())?;
    }
    out.flush().map_err(|e| e.to_string())
}
fn run() -> Result<(), String> {
    let mut args = env::args().skip(1);
    let input = args
        .next()
        .ok_or("usage: tfm-rs INPUT [--config CONFIG.json] [--output DIRECTORY]")?;
    if input == "--help" || input == "-h" {
        println!("tfm-rs INPUT [--config CONFIG.json] [--output DIRECTORY]");
        return Ok(());
    }
    let mut c = Config::default();
    let mut output = PathBuf::from("results/rust");
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--config" => {
                let p = args.next().ok_or("missing config path")?;
                c = serde_json::from_str(&fs::read_to_string(p).map_err(|e| e.to_string())?)
                    .map_err(|e| e.to_string())?;
            }
            "--output" => output = PathBuf::from(args.next().ok_or("missing output path")?),
            _ => return Err(format!("unknown argument: {arg}")),
        }
    }
    let start = Instant::now();
    let result = reconstruct(Path::new(&input), &c)?;
    fs::create_dir_all(&output).map_err(|e| e.to_string())?;
    csv(&output.join("image.csv"), &result.image, c.nx)?;
    csv(&output.join("envelope.csv"), &result.envelope, c.nx)?;
    let metadata = serde_json::json!({"config":c,"shape":[c.nz(),c.nx],"t0_index":result.t0_index,
        "dt":result.dt,"probe_x_mm":result.probes,"elapsed_seconds":start.elapsed().as_secs_f64(),
        "units":{"coordinates":"mm","dt":"s","t0_index":"samples (MATLAB 1-based convention)"}});
    fs::write(
        output.join("metadata.json"),
        serde_json::to_string_pretty(&metadata).unwrap(),
    )
    .map_err(|e| e.to_string())?;
    println!("{metadata}");
    Ok(())
}
fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}
