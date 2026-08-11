function result = reconstruct_qsm_grappa(datFile, opts)
%RECONSTRUCT_QSM_GRAPPA Reconstruct multi-echo 3-D GRE Siemens raw data.
%   RESULT = RECONSTRUCT_QSM_GRAPPA(DATFILE, OPTS) reads DATFILE with
%   mapVBVD2, inserts the reference/ACS data, compresses coils, performs
%   slice-wise 2-D GRAPPA in the phase/partition plane, and combines coils
%   with ESPIRiT receive maps.

if nargin < 2, opts = struct; end
assert(ischar(datFile) || (isstring(datFile) && isscalar(datFile)), ...
    'datFile must be a character vector or string scalar.');
datFile = char(datFile);
defaults = struct('acceleration',[2 1], 'acsLines',32, ...
    'phaseEncodingLines',0, ...
    'numVirtualCoils',32, 'coilCalibrationEcho',1, ...
    'grappaKernel',[3 3], 'lambdaPercent',1e-3, ...
    'espiritKernel',[6 6], 'espiritEigenThreshold',0.5, ...
    'numWorkers',0, 'verbose',true);
opts = merge_options(defaults, opts);
validateattributes(opts.acceleration, {'double'}, ...
    {'row','numel',2,'integer','positive'});
validateattributes(opts.acsLines, {'double'}, {'scalar','integer','positive'});
validateattributes(opts.phaseEncodingLines, {'double'}, ...
    {'scalar','integer','nonnegative'});
validateattributes(opts.numVirtualCoils, {'double'}, {'scalar','integer','positive'});
validateattributes(opts.coilCalibrationEcho, {'double'}, {'scalar','integer','positive'});
validateattributes(opts.grappaKernel, {'double'}, ...
    {'row','numel',2,'integer','positive'});
validateattributes(opts.lambdaPercent, {'double'}, {'scalar','nonnegative'});
validateattributes(opts.espiritKernel, {'double'}, ...
    {'row','numel',2,'integer','positive'});
validateattributes(opts.espiritEigenThreshold, {'double'}, {'scalar','>',0,'<=',1});
validateattributes(opts.numWorkers, {'double'}, {'scalar','integer','nonnegative'});

requiredFunctions = {'mapVBVD', 'ifft3call', 'svd_compress3d', ...
    'svd_apply3d', 'ifftc', 'ecalib_soft', 'dot_mult', ...
    'grappa_gfactor_2d_jvc3', 'fft2call', 'ifft2call', ...
    'coil_combine', 'rsos'};
missingFunctions = requiredFunctions(cellfun( ...
    @(name) isempty(which(name)), requiredFunctions));
assert(isempty(missingFunctions), ...
    ['Missing reconstruction dependencies on the MATLAB path: %s. ', ...
     'See reconstruction/README.md for installation instructions.'], ...
    strjoin(missingFunctions, ', '));

assert(isfile(datFile), 'Raw data file not found: %s', datFile);
if opts.verbose, fprintf('Reading %s\n', datFile); end
twix = mapVBVD(datFile);
if iscell(twix), twix = twix{end}; end
assert(isfield(twix, 'image') || isprop(twix, 'image'), ...
    'The selected RAID has no image data.');

k = squeeze(twix.image());
kref = squeeze(twix.refscan());
assert(ndims(k) >= 4, 'Unexpected image data dimensions.');
assert(ndims(kref) >= 4, 'No usable reference scan was found.');

% mapVBVD order [readout coil phase partition echo] -> [x y z coil echo].
k = ensure_five_dims(k);
kref = ensure_five_dims(kref);
k = permute(k, [1 3 4 2 5]);
kref = permute(kref, [1 3 4 2 5]);

if size(kref,3) < size(k,3)
    kref(:,:,end+1:size(k,3),:,:) = 0;
elseif size(kref,3) > size(k,3)
    error('Reference scan has more partitions than imaging data.');
end

% Accelerated mapVBVD data can stop at line Ny-1 when the final nominal
% phase-encode line is unsampled. Establish the nominal matrix before ACS
% insertion; otherwise the reference block is shifted by one k-space line.
nyRaw = size(k,2);
if opts.phaseEncodingLines == 0
    nominalNy = nyRaw + mod(nyRaw, 2);
else
    nominalNy = opts.phaseEncodingLines;
end
assert(nominalNy >= nyRaw, ...
    'phaseEncodingLines (%d) is smaller than raw line extent (%d).', ...
    nominalNy, nyRaw);
if nominalNy > nyRaw
    k(:,nominalNy,:,:,:) = 0;
end
assert(size(kref,2) == opts.acsLines, ...
    'Reference scan has %d phase lines; opts.acsLines is %d.', ...
    size(kref,2), opts.acsLines);
acsY = centered_indices(nominalNy, opts.acsLines);
k(:,acsY,:,:,:) = kref;

% Make spatial dimensions even, matching the source reconstruction.
spatialSize = size(k); spatialSize = spatialSize(1:3);
for d = 1:3
    if mod(spatialSize(d),2)
        idx = repmat({':'}, 1, 5); idx{d} = spatialSize(d)+1;
        k(idx{:}) = 0;
    end
end
[nx, ny, nz, nc, ne] = size(k);
assert(opts.coilCalibrationEcho <= ne, 'Calibration echo exceeds echo count.');
nvc = min(opts.numVirtualCoils, nc);

refImage = ifft3call(kref);
[~, compression] = svd_compress3d( ...
    refImage(:,:,:,:,opts.coilCalibrationEcho), nvc, 1);
kCompressed = complex(zeros(nx,ny,nz,nvc,ne, 'like', k));
for echo = 1:ne
    kCompressed(:,:,:,:,echo) = svd_apply3d(k(:,:,:,:,echo), compression);
end

% Hybrid space: image along readout, k-space along phase and partition.
kHybrid = ifftc(kCompressed, 1);
rssImage = zeros(nx,ny,nz,ne, 'like', real(k));
combinedImage = complex(zeros(nx,ny,nz,ne, 'like', k));

numAcs = floor([size(kref,2), nz] / 2) * 2 - 2;
assert(all(numAcs >= opts.grappaKernel .* opts.acceleration), ...
    'ACS region is too small for the selected GRAPPA kernel.');

poolStartedHere = false;
if opts.numWorkers > 0
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', opts.numWorkers); poolStartedHere = true;
    elseif pool.NumWorkers ~= opts.numWorkers
        error(['An existing pool has %d workers. Close it or set ', ...
            'opts.numWorkers=%d.'], pool.NumWorkers, pool.NumWorkers);
    end
end

parfor x = 1:nx
    slab = squeeze(kHybrid(x,:,:,:,:));
    acs = complex(zeros(size(slab), 'like', slab));
    iy = centered_indices(ny, numAcs(1));
    iz = centered_indices(nz, numAcs(2));
    acs(iy,iz,:,:) = slab(iy,iz,:,:);

    [maps, weights] = ecalib_soft(acs(:,:,:,opts.coilCalibrationEcho), ...
        opts.acsLines, opts.espiritKernel, opts.espiritEigenThreshold);
    receive = dot_mult(maps, weights >= opts.espiritEigenThreshold);
    receive = abs(receive) .* exp(1i * angle(receive .* ...
        repmat(conj(receive(:,:,1)), [1 1 nvc])));

    coilImages = complex(zeros(ny,nz,nvc,ne, 'like', slab));
    combo = complex(zeros(ny,nz,ne, 'like', slab));
    for echo = 1:ne
        coilImages(:,:,:,echo) = grappa_gfactor_2d_jvc3( ...
            slab(:,:,:,echo), acs(:,:,:,echo), ...
            opts.acceleration(1), opts.acceleration(2), numAcs, ...
            opts.grappaKernel, opts.lambdaPercent, 1);
    end
    kRecon = fft2call(coilImages);
    kRecon(acsY,:,:,:) = slab(acsY,:,:,:);
    coilImages = ifft2call(kRecon);
    for echo = 1:ne
        combo(:,:,echo) = coil_combine(coilImages(:,:,:,echo), receive, 3);
    end
    rssImage(x,:,:,:) = rsos(coilImages, 3);
    combinedImage(x,:,:,:) = combo;
end

if poolStartedHere, delete(gcp('nocreate')); end
result = struct('kSpace', k, 'kSpaceCompressed', kCompressed, ...
    'rssImage', rssImage, 'combinedImage', combinedImage, ...
    'compressionMatrix', compression, 'options', opts, ...
    'sourceFile', datFile);
end

function x = ensure_five_dims(x)
sz = size(x); sz(end+1:5) = 1; x = reshape(x, sz(1:5));
end

function idx = centered_indices(n, width)
assert(width <= n && mod(width,2)==0, ...
    'Centered width must be even and no larger than the dimension.');
idx = (floor(n/2)-width/2+1):(floor(n/2)+width/2);
end

function out = merge_options(defaults, supplied)
assert(isstruct(supplied) && isscalar(supplied), 'opts must be a scalar struct.');
out = defaults;
names = fieldnames(supplied);
unknown = setdiff(names, fieldnames(defaults));
assert(isempty(unknown), 'Unknown option: %s', strjoin(unknown, ', '));
for n = 1:numel(names), out.(names{n}) = supplied.(names{n}); end
end
