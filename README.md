# 3D Multi-Echo Gradient Echo (ME-GRE) Pulseq Sequence

This repository contains a MATLAB script (`script_writeGradientEcho3D_label_spoil_v1.m`) used to generate a 3D Multi-Echo Gradient Echo (ME-GRE) MRI sequence using the open-source **Pulseq** framework. 

This sequence is optimized for phase imaging, Quantitative Susceptibility Mapping (QSM), and parallel imaging reconstruction, featuring steady-state gradient/RF spoiling and built-in mapVBVD labels for GRAPPA reconstruction.

## 🚀 Key Features

* **Multi-Echo Acquisition:** Acquires 5 echoes per TR using a unipolar flyback readout gradient design.
* **GRAPPA Acceleration:** Implements 1D phase-encode (Y-axis) undersampling with integrated Auto-Calibration Signal (ACS) lines.
* **Siemens mapVBVD Labels:** Injects inline Pulseq labels (`LIN`, `PAR`, `ECO`, `REF`, `IMA`, `NOISE`) directly into the ADC blocks. This allows seamless data sorting and standard pipeline reconstruction using tools like `mapVBVD`.
* **Spoiling & Steady State:** * Continuous quadratic RF phase spoiling globally across all TRs (using an 84° increment for optimal steady-state stabilization).
  * Z-axis gradient phase-encode blips are perfectly rewound and combined with the spoiler (`gzRewindAndSpoil`) to maintain a constant net gradient moment per TR.
* **Dummy Scans:** Plays out 50 dummy TRs prior to acquisition with fully matched gradient structures to gracefully drive the magnetization into the steady state.

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
2. **Pulseq MATLAB Toolbox (v1.5.1):** You must have the official [Pulseq repository](https://github.com/pulseq/pulseq) downloaded and added to your MATLAB path.
3. On the scanner console, from the resolution tab, set GRAPPA as the acceleration method using R=2 and 32 acs lines (integrated). From the sequence>special tab, select data handling -> ICE_STD to get the online recon to work.

## 🏃 Getting Started

1. Clone this repository to your local machine.
2. Open `script_writeGradientEcho3D_label_spoil_v1.m` in MATLAB.
3. Update the Pulseq directory path at the top of the script to point to your local installation:
   ```matlab
   addpath('path/to/your/pulseq-1.5.1/matlab/')
