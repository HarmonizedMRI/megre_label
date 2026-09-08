# QSM 3-D GRE GRAPPA reconstruction toolbox

This folder packages the raw-data reconstruction beginning at line 450 of
`script_writeGradientEcho3D_label_spoil_github_recon_bay4_v0.m`. It reads a
Siemens `.dat` file, restores ACS lines from the reference scan, compresses
the receive coils, reconstructs missing phase/partition samples with 2-D
GRAPPA, estimates ESPIRiT receive maps, and returns magnitude and complex
multi-echo volumes.

## Requirements

- MATLAB R2019b or newer.
- Parallel Computing Toolbox only when `opts.numWorkers > 0`.
- Enough memory for the raw, compressed, RSS, and complex volumes.
- A Siemens raw file containing both `image` and `refscan` acquisitions.

Clone this repository with submodules so that the upstream
[mapVBVD](https://github.com/pehses/mapVBVD) dependency is available:

```bash
git clone --recurse-submodules https://github.com/HarmonizedMRI/megre_label.git
```

If the repository was cloned without submodules, run
`git submodule update --init --recursive`. The first-party ESPIRiT,
centered-FFT, coil-compression, and 2-D GRAPPA utilities are included under
`external/recon_utils` and distributed under the repository's MIT license.
The reconstruction function adds both dependency directories to the MATLAB
path automatically.

## Example Siemens raw data

The Siemens raw data are available in this
[Dropbox folder](https://www.dropbox.com/scl/fo/xk1q7wmbca5mirlaau4ab/AHAkgtLHZW7qOyfU7SVh5DQ?rlkey=xqrpv1w7cgtbqjxd0rvlb9pp0&st=qt6xvqr9&dl=0).
Download the following file for reconstruction:

```text
meas_MID00048_FID27124_pulseq151fix_gre3d_label_spoil_xyzflip_2.dat
```

`example_reconstruct.m` is configured to use this filename. Set `dataPath`
in that script to the local folder containing the downloaded file.

## Quick start

Set `dataPath` in `example_reconstruct.m`, then run the script. To reconstruct
a different file, also update `filename` without the `.dat` extension.
Alternatively:

```matlab
addpath('/path/to/qsm_grappa_recon_toolbox')
opts = struct('acceleration',[2 1], 'acsLines',32, ...
              'phaseEncodingLines',184, ...
              'numVirtualCoils',32, 'numWorkers',16);
recon = reconstruct_qsm_grappa('/path/to/meas_file.dat', opts);
```

The supplied example writes
`*_combined_magnitude_oriented.nii.gz` and
`*_combined_phase_oriented.nii.gz` under `output/`. It swaps readout/phase
and flips the partition dimension to match the original script's NIfTI
orientation. Every echo is retained in the fourth dimension; phase is in
radians over `[-pi, pi]`.

The main outputs are:

- `recon.rssImage`: real `[Nx Ny Nz Necho]` GRAPPA magnitude volume.
- `recon.combinedImage`: complex `[Nx Ny Nz Necho]` ESPIRiT-combined volume.
- `recon.kSpace`: ACS-filled, even-padded multi-coil k-space.
- `recon.kSpaceCompressed`: virtual-coil k-space used by GRAPPA.
- `recon.compressionMatrix`, `recon.options`, and `recon.sourceFile`.

## Parameters and assumptions

The defaults reproduce the source script: acceleration `[Ry Rz]=[2 1]`, 32
ACS lines, 32 virtual coils, a `[3 3]` GRAPPA kernel, and the first echo for
coil-map calibration. Set `phaseEncodingLines` to the prescribed nominal
matrix size (184 for the supplied sequence). A value of zero infers the next
even size, which handles the 183-line accelerated raw extent. `acsLines` must equal the phase dimension of the
reference scan. Odd spatial dimensions are zero-padded at the high-index
edge. A reference scan with fewer trailing partitions is zero-padded, as in
the original Bay 4 handling; a reference with more partitions is rejected.

Set `numWorkers=0` for serial execution (the `parfor` loop runs as a normal
loop). If the function starts a pool, it also closes that pool. It does not
close a pool that was already running.

## Notes

The source block allocated GRAPPA coil arrays using the original coil count
after SVD compression. This packaged version uses the actual virtual-coil
count, preventing a dimension mismatch when compression reduces the number
of channels. Hard-coded data paths, slice selection, figures, and scanner
matrix size were removed; dimensions are inferred from the raw data.

`mapVBVD` remains an independently maintained upstream project and is linked
as a Git submodule rather than relicensed by this repository. Full end-to-end
reconstruction requires a compatible, de-identified Siemens raw-data file and
is not exercised by the repository's dependency smoke test.
