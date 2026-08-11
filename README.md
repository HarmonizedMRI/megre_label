# 3D Multi-Echo Gradient Echo (ME-GRE) Pulseq Sequence

This repository contains a MATLAB script (`script_writeGradientEcho3D_label_spoil_github_v0.m`) used to generate a 3D Multi-Echo Gradient Echo (ME-GRE) MRI sequence using the open-source **Pulseq** framework. 

This sequence is optimized for phase imaging, Quantitative Susceptibility Mapping (QSM), and parallel imaging reconstruction, featuring steady-state gradient/RF spoiling and built-in mapVBVD labels for GRAPPA reconstruction.

## 🚀 Key Features

* **Multi-Echo Acquisition:** Acquires 5 echoes per TR using a unipolar flyback readout gradient design.
* **GRAPPA Acceleration:** Implements 1D phase-encode (Y-axis) undersampling with integrated Auto-Calibration Signal (ACS) lines.
* The sequence also allows for R=1 fully sampled acquisition, enabled by setting Ry = 1;
* **Siemens mapVBVD Labels:** Injects inline Pulseq labels (`LIN`, `PAR`, `ECO`, `REF`, `IMA`, `NOISE`) directly into the ADC blocks. This allows seamless data sorting and standard pipeline reconstruction using tools like `mapVBVD`.
* **Spoiling & Steady State:** * Continuous quadratic RF phase spoiling globally across all TRs (using an 84° increment for optimal steady-state stabilization).
  * Z-axis gradient phase-encode blips are perfectly rewound and combined with the spoiler (`gzRewindAndSpoil`) to maintain a constant net gradient moment per TR.
* **Dummy Scans:** Plays out 50 dummy TRs prior to acquisition with fully matched gradient structures to gracefully drive the magnetization into the steady state.
* Gx, Gy, and Gz gradients are flipped so that the online reconstructed image matches those reconstructed using the native Siemens GRE sequence.

## ⚙️ Sequence Parameters

By default, the script generates a protocol with the following parameters (which can be easily modified in the script):

* **Matrix Size:** 264 (RO) × 184 (PE) × 144 (PAR)
* **FOV:** 264 × 184 × 144 mm³
* **Resolution:** 1 mm isotropic
* **TEs (Echo Times):** 5.0, 11.0, 17.0, 23.0, 29.0 ms
* **TR (Repetition Time):** 35.0 ms
* **Flip Angle:** 15°
* **Acceleration (Ry):** 2
* **ACS Lines:** 32

## 🛠 Prerequisites

To run this script, you will need:
1. **MATLAB** (R2019a or newer recommended)
2. **Pulseq MATLAB Toolbox (version newer than 2026-05-19):** You must have the official [Pulseq repository](https://github.com/pulseq/pulseq) downloaded and added to your MATLAB path.
3. (This is only needed when running the sequence on a GE scanner) Download and add pulceq toolbox to matlab path [PulCeq v2.5.2.0](https://github.com/HarmonizedMRI/PulCeq/releases/tag/v2.5.2.0).
4. **On the scanner console:** from the resolution tab, set GRAPPA as the acceleration method using R=2 and 32 acs lines (integrated). From the sequence>special tab, select data handling -> ICE_STD to get the online recon to work.
5. Set Orientation to "Transversal" and phase encoding dir. to "R>>L". The slab position can also rotated/translated.
6. **if fully sampled acquisition is desired** set acceleration to "None", still using ICE_STD for data handling.

## 🏃 Getting Started

1. Clone this repository and its submodules using
   `git clone --recurse-submodules`.
2. Open `script_writeGradientEcho3D_label_spoil_github.m` in MATLAB.
3. Update the Pulseq directory path at the top of the script to point to your local installation:
   ```matlab
   addpath('pulseq-version-newer-than-20260519/matlab/')
   addpath('path/to/pulceq/v2.5.2.0/matlab') % only needed for GE scanner
   ```

## Raw-data reconstruction for Siemens scanners

The [`reconstruction`](reconstruction/) package reconstructs the labeled
Siemens raw data produced by this sequence. It includes the first-party
GRAPPA, ESPIRiT, centered-FFT, and coil-compression utilities and links
mapVBVD as a Git submodule.

The pipeline:

1. Reads the `.dat` image and reference acquisitions with mapVBVD.
2. Inserts the integrated ACS lines and compresses the receive coils.
3. Performs slice-parallel 2-D GRAPPA reconstruction.
4. Estimates ESPIRiT receive maps and produces complex multi-echo images.
5. Writes separate 4-D magnitude and phase NIfTI files.

To run it, edit the raw-data path in
[`reconstruction/example_reconstruct.m`](reconstruction/example_reconstruct.m),
then execute that script in MATLAB. The default reconstruction settings match
the sequence defaults: `Ry=2`, `Rz=1`, 32 ACS lines, and 32 virtual coils.
See the [reconstruction README](reconstruction/README.md) for configuration,
requirements, outputs, and assumptions.

Generated raw data and reconstruction outputs are intentionally excluded from
version control.

## Raw-data reconstruction for GE scanners
