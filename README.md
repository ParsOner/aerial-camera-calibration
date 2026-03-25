# 📷 Aerial Camera Calibration — MATLAB Toolset

Camera calibration pipeline for aerial/UAV imaging systems.  
Implements checkerboard-based intrinsic calibration with radial, tangential, and skew distortion models using MATLAB's Computer Vision Toolbox.

---

## 📁 Repository Structure

```
aerial-camera-calibration/
│
├── README.md
├── .gitignore
│
├── calibration/
│   ├── run_calibration_model.m       # Main calibration script (run this first)
│   ├── N_Calibration.m               # Image count sweep & robustness analysis
│   └── calibration_analyze.m         # Post-calibration visualization & analysis
│
├── results/
│   ├── M1_4th_noTan_noSkew_calibration.mat   # Model 1 — calibration output
│   ├── M2_4th_Tan_noSkew_calibration.mat     # Model 2 — calibration output
│   ├── M3_6th_Tan_noSkew_calibration.mat     # Model 3 — calibration output
│   ├── M4_6th_Tan_Skew_calibration.mat       # Model 4 — calibration output (recommended)
│   └── calibrationSession.mat                # Full session backup
│
└── docs/
    └── Project_Presentation_Aerial_Part_3.pptx   # Project presentation
```

---

## 🔧 Requirements

| Requirement | Details |
|---|---|
| MATLAB | R2022b or later recommended |
| Computer Vision Toolbox | `estimateCameraParameters`, `detectCheckerboardPoints` |
| Checkerboard target | Physical calibration pattern |
| Calibration images | Not included — see note below |

> 📸 **Calibration images are not included** in this repository.  
> To re-run calibration from scratch, collect checkerboard images and update `imagesDir` in `run_calibration_model.m`.  
> Minimum recommended: **15–20 images** from varied angles and distances.  
> To use pre-computed results directly, load any `.mat` file from the `results/` folder.

---

## 🚀 Quick Start

### Option A — Use pre-computed results (no images needed)

```matlab
% Load a result and run the analysis/visualization script
load("results/M4_6th_Tan_Skew_calibration.mat");   % loads 'calib' struct
run('calibration/calibration_analyze.m')
```

> ⚠️ Before running `calibration_analyze.m`, update the `load()` path at the top of the script to point to your local `results/` folder.

### Option B — Re-run calibration from scratch

**1. Configure paths and model**

Open `calibration/run_calibration_model.m` and set:

```matlab
imagesDir  = "path/to/your/calibration/images";
resultsDir = "path/to/output/results";
squareSize = 25;   % checkerboard square size in mm
```

**2. Choose a model configuration**

```matlab
% Model 1: 4th-order radial, no tangential, no skew
cfg = makeModelCfg("M1_4th_noTan_noSkew", 2, false, false);

% Model 2: 4th-order radial, with tangential, no skew
cfg = makeModelCfg("M2_4th_Tan_noSkew", 2, true, false);

% Model 3: 6th-order radial, with tangential, no skew
cfg = makeModelCfg("M3_6th_Tan_noSkew", 3, true, false);

% Model 4: 6th-order radial, with tangential, with skew  <- recommended
cfg = makeModelCfg("M4_6th_Tan_Skew", 3, true, true);
```

**3. Run**

```matlab
run('calibration/run_calibration_model.m')
```

Output workspace variables: `calib`, `params`, `errors`, `report`, `analysis`, `score`

### Option C — Analyze optimal image count

```matlab
run('calibration/N_Calibration.m')
```

Sweeps N = {5, 8, 10, 15, 20} images with K=10 random subsets each.  
Outputs score, reprojection error, and RMS components vs. N.

---

## 📊 Calibration Models

| Model | Radial Order | Tangential | Skew | Parameters |
|---|---|---|---|---|
| M1 | 4th (k1, k2) | No | No | 6 |
| M2 | 4th (k1, k2) | Yes (p1, p2) | No | 8 |
| M3 | 6th (k1, k2, k3) | Yes (p1, p2) | No | 9 |
| M4 | 6th (k1, k2, k3) | Yes (p1, p2) | Yes | 10 |

### Distortion Model

**Radial:**

```
x' = x * L(r),   y' = y * L(r)
L(r) = 1 + k1*r^2 + k2*r^4 + k3*r^6
```

**Tangential:**

```
dx_t = 2*p1*x*y + p2*(r^2 + 2*x^2)
dy_t = p1*(r^2 + 2*y^2) + 2*p2*x*y
```

---

## 📈 Scoring Metric

A robust assignment-style score (lower is better) is computed after each calibration run.  
It combines normalized parameter uncertainties with the RMS reprojection error:

```
Score = |σ_fx|/|fx| + |σ_fy|/|fy| + |σ_cx|/|cx| + |σ_cy|/|cy|
      + Σ |σ_ki| / max(|ki|, floor_i)
      + |σ_alpha| / max(|alpha|, 1.0)
      + sqrt(err_x^2 + err_y^2)
```

Floors prevent division blow-ups when parameters (e.g. skew, tangential) are near zero.

---

## 💾 Results — `.mat` File Contents

Each `*_calibration.mat` contains a single struct `calib` with the following fields:

| Field | Description |
|---|---|
| `calib.params` | `cameraParameters` object (focal length, principal point, distortion) |
| `calib.errors` | `cameraCalibrationErrors` object (std errors for all parameters) |
| `calib.report` | Struct with `intrinsicsTable`, `distortionTable`, `extrinsicsTable` |
| `calib.analysis` | Reprojection error statistics (per-image mean, RMS, overall) |
| `calib.score` | Scalar quality score (lower is better) |
| `calib.cfg` | Model configuration used |

---

## 📬 Contact

m.oner@studenti.unina.it

---

## 📄 License

This project is released under the [MIT License](LICENSE).
