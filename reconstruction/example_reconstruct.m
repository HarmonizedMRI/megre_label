%% Example: multi-echo 3-D GRE raw-data reconstruction
toolboxRoot = fileparts(mfilename('fullpath'));
addpath(toolboxRoot);

dataPath = '/path/to/raw/data';
filename = 'meas_MID00048_FID27124_pulseq151fix_gre3d_label_spoil_xyzflip_2';
rawFile = fullfile(dataPath, [filename '.dat']);
opts = struct;
opts.acceleration = [2 1];
opts.acsLines = 32;
opts.phaseEncodingLines = 184;
opts.numVirtualCoils = 32;
opts.numWorkers = 16;       % use 0 to run without opening a pool
opts.verbose = true;

recon = reconstruct_qsm_grappa(rawFile, opts);

% Preserve all echoes as the fourth NIfTI dimension. niftiwrite does not
% support complex arrays, so magnitude and phase are stored separately.
outputDir = fullfile(toolboxRoot, 'output');
if ~isfolder(outputDir), mkdir(outputDir); end
% Match the orientation used at the end of the original reconstruction:
% [readout phase partition echo] -> [phase readout partition echo], then
% reverse the partition direction.
combinedOriented = flip(permute(recon.combinedImage, [2 1 3 4]), 3);
magnitudeFile = fullfile(outputDir, [filename '_combined_magnitude_oriented.nii']);
phaseFile = fullfile(outputDir, [filename '_combined_phase_oriented.nii']);
niftiwrite(single(abs(combinedOriented)), magnitudeFile, 'Compressed', true);
niftiwrite(single(angle(combinedOriented)), phaseFile, 'Compressed', true);
fprintf('Wrote magnitude: %s.gz\n', magnitudeFile);
fprintf('Wrote phase:     %s.gz\n', phaseFile);

% Magnitude and phase for the first echo:
figure; imagesc(squeeze(recon.rssImage(:,:,round(end/2),1))); axis image off;
colormap gray; title('GRAPPA root-sum-of-squares');
figure; imagesc(angle(squeeze(recon.combinedImage(:,:,round(end/2),1))));
axis image off; colorbar; clim([-pi pi]); title('ESPIRiT-combined phase');

% save('reconstruction.mat', 'recon', '-v7.3');
