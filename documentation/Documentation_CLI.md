# FragHub - Command Line Interface (CLI) Documentation

FragHub now includes a fully featured **Headless CLI mode**. This mode allows you to run the highly optimized Rust parsing and filtering engine directly from your terminal, without launching the Vue.js/Electron graphical user interface.

This is particularly useful for:
- **Server-side automation** (CI/CD, cron jobs, etc.)
- **Batch processing** of massive datasets without GUI overhead
- **Integration** into broader data-science pipelines

---

## 🚀 How to Run

To trigger the CLI mode, you **must** append the `--cli` flag when executing FragHub.

### Running from source (Python)
```bash
uv run python scripts/FragHub.py --cli --input_directory "data/input" --output_directory "data/output"
```
> [!IMPORTANT]
> Running from source requires **RDKit C++** (installed with mamba/micromamba in `~/rdkit_cpp`, or in the folder given by `RDKIT_PREFIX`) and the Rust module compiled with `python build_backend.py --rust-only`. See section 7 of the developer documentation for your platform. The compiled executable needs nothing: RDKit is bundled inside it.

### Running the compiled binary
If FragHub is packaged as an executable (e.g., via PyInstaller/Electron), you can pass the arguments directly to the executable:
```bash
./fraghub.exe --cli --input_directory "data/input" --output_directory "data/output"
```

> [!NOTE]
> When running in CLI mode, FragHub outputs real-time progress bars and cleaning reports directly to your terminal (`stdout`), instead of writing to the debug log file.

---

## 📋 Required Arguments

These arguments must be provided for FragHub to run.

| Argument | Type | Description |
| :--- | :--- | :--- |
| `--input_directory` | `List[str]` | One or multiple paths to input directories or files (MGF, MSP, JSON, CSV). |
| `--output_directory` | `str` | The directory where the cleaned and converted datasets will be saved. |

---

## 📁 Input Resolution

The `--input_directory` argument supports recursive directory parsing and accepts multiple inputs:
- **Files**: Absolute or relative paths to individual files (`--input_directory file1.msp file2.mgf`).
- **Directories**: Absolute or relative paths to directories (`--input_directory /my/data/folder`). FragHub will recursively scan all subdirectories to find valid MS files (`.msp`, `.mgf`, `.csv`, `.json`).
- **Mixed**: Files and directories can be provided simultaneously.

---

## ⏱️ Execution Metrics

By default, the CLI outputs real-time execution metrics to `stdout`:
- **Elapsed Time**: Time spent on the current process.
- **ETA**: Estimated time of arrival for the current process (`< 00:15`).
- **Speed**: Processing rate in items per second (`it/s`).

*Example output:*
`▶ derivation and calculation (RDKit via Rust):`
`  [█████████████░░░░░░░░░] 50.0% | 500/1000 [00:15 < 00:15, 33.3 it/s]`

## 🔇 Quiet Mode

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--quiet` | `False` | Silence all progress bars, step indicators, and deletion reports. The program will only output a single line upon successful completion with the total elapsed time. |

---

## 🎛️ Filter Arguments (Optional)

FragHub uses standard `yes` or `no` arguments for enabling/disabling filters and outputs. All filters are **enabled by default (`yes`)** unless specified otherwise, matching the GUI's default behavior.

### General Intensity
| Argument | Default | Description |
| :--- | :--- | :--- |
| `--normalize_intensity` | `yes` | Normalize the intensity of peaks to 100%. |
| `--remove_peak_above_precursormz` | `yes` | Remove all peaks with an m/z greater than the precursor m/z. |

### Minimum Peaks Rule
| Argument | Default | Description |
| :--- | :--- | :--- |
| `--check_minimum_peak_requiered` | `yes` | Enable the minimum peaks threshold filter. |
| `--check_minimum_peak_requiered_n_peaks` | `3.0` | Minimum number of peaks a spectrum must have to be kept. |

### Peak List Reduction
| Argument | Default | Description |
| :--- | :--- | :--- |
| `--reduce_peak_list` | `yes` | Enable reducing the peak list to the top most intense peaks. |
| `--reduce_peak_list_max_peaks` | `500.0` | The maximum number of top intense peaks to keep per spectrum. |

### Entropy Score Filtering
| Argument | Default | Description |
| :--- | :--- | :--- |
| `--remove_spectrum_under_entropy_score` | `yes` | Enable spectrum deletion based on its calculated Shannon entropy. |
| `--remove_spectrum_under_entropy_score_value` | `0.5` | Threshold below which the spectrum will be deleted. |

### M/Z Range Filtering
| Argument | Default | Description |
| :--- | :--- | :--- |
| `--keep_mz_in_range` | `yes` | Enable deleting peaks outside a specific m/z range. |
| `--keep_mz_in_range_from_mz` | `50.0` | The minimum m/z allowed. |
| `--keep_mz_in_range_to_mz` | `2000.0` | The maximum m/z allowed. |

### High Peaks Rule
| Argument | Default | Description |
| :--- | :--- | :--- |
| `--check_minimum_of_high_peaks_requiered` | `yes` | Enable minimum high peaks rule. |
| `--check_minimum_of_high_peaks_requiered_intensity_percent` | `5.0` | What percentage of the base peak defines a "high peak". |
| `--check_minimum_of_high_peaks_requiered_no_peaks` | `2.0` | Minimum number of "high peaks" a spectrum must contain. |

---

## 🔬 De Novo Calculations

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--calculate_de_novo` | `no` | Enable De Novo peptide fragment annotation (Requires heavy computation). |
| `--de_novo_ppm_tolerance` | `10.0` | PPM tolerance used during De Novo calculations. |

---

## 💾 Output Formats & Project Management

Specify which formats you want FragHub to generate at the end of the processing pipeline.

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--csv` | `yes` | Output data as tabular `.csv` files. |
| `--msp` | `yes` | Output data as `.msp` spectral files. |
| `--json` | `yes` | Output data as highly compacted `.json` files. |
| `--reset_updates` | `no` | If `yes`, forces FragHub to ignore previously cleaned caches and run everything from scratch in the output directory. |

---

## 💡 Examples

### 1. Simple Execution (Default Settings)
Runs FragHub with all the default filters enabled (mirroring the GUI defaults).
```bash
uv run python scripts/FragHub.py \
  --cli \
  --input_directory "/path/to/my_library.msp" \
  --output_directory "/path/to/output_folder"
```

### 2. Disabling Specific Outputs
Runs FragHub but only outputs JSON files, disabling MSP and CSV generation.
```bash
uv run python scripts/FragHub.py \
  --cli \
  --input_directory "/path/to/my_library.mgf" \
  --output_directory "/path/to/output_folder" \
  --csv no \
  --msp no \
  --json yes
```

### 3. Strict Custom Filtering
Runs FragHub with very strict entropy and high-peak rules, and enforces a reset of previous project caches.
```bash
uv run python scripts/FragHub.py \
  --cli \
  --input_directory "/data/raw/batch1/" "/data/raw/batch2/" \
  --output_directory "/data/processed/" \
  --reset_updates yes \
  --remove_spectrum_under_entropy_score_value 1.2 \
  --check_minimum_of_high_peaks_requiered_intensity_percent 10.0 \
  --check_minimum_of_high_peaks_requiered_no_peaks 5.0
```
