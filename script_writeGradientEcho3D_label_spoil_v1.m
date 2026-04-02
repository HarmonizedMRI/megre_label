% -------------------------------------------------------------------------
%%
% -------------------------------------------------------------------------

clear
restoredefaultpath

addpath /cluster/berkin/berkin/Matlab_Code_New/PULSEQ/QSM_harmonization/pulseq-1.5.1/matlab/
% addpath /cluster/berkin/berkin/Matlab_Code_New/PULSEQ/QSM_harmonization/pulseq-1.4.2/matlab/
addpath /autofs/cluster/berkin/xingwang/share/seqeyes/

% Define FOV and resolution
fov = [264e-3 184e-3 144e-3];    

Nx = 264;
Ny = 184/1; 
Nz = 144/1;            

adc_dwell = 15e-6;

Tread = adc_dwell * Nx;

disp(['Tread ', num2str(Tread*1e3), ' ms'])
% Tread = 4e-3;     % adc duration

Ndummy = 50;

% define system properties
Gmax = 20;
Smax = 100;

sys = mr.opts('MaxGrad', Gmax, 'GradUnit', 'mT/m', ...
    'MaxSlew', Smax, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

seq = mr.Sequence(sys);           % Create a new sequence object

alpha = 15;     % flip angle
rf_duration = 1.e-3;
rf_tbp = 20;    % RF time bandwidth
thickness = 128e-3;   % slab thickness to excite

% Create alpha-degree slice selection pulse and gradient
[rf, gz, gz_rf_reph] = mr.makeSincPulse(alpha*pi/180, 'Duration', rf_duration,...
    'SliceThickness',thickness,'apodization',0.42,'timeBwProduct',rf_tbp,'use','excitation','system',sys);


% -------------------------------------------------------------------------
% Define other gradients and ADC events
% -------------------------------------------------------------------------

deltak = 1./fov;
spoil_cycles = 2;

Tspoil = 3.2e-3;    % allow longer time for spoiler
Tpre = 1.6e-3;    % prewinder duration


gx = mr.makeTrapezoid('x',sys,'FlatArea',Nx*deltak(1),'FlatTime',Tread);

% adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime);
% sys has to be provided otherwise adc does not have dead time by default,
% which causes an error during noise adc -> other adc events are saved by
% the gx ramp up/down times
adc = mr.makeAdc(Nx,sys,'Duration',gx.flatTime,'Delay',gx.riseTime);

gxPre = mr.makeTrapezoid('x',sys,'Area',-gx.area/2,'Duration',Tpre);

gxSpoil = mr.makeTrapezoid('x',sys,'Area',spoil_cycles * Nx * deltak(1),'Duration',Tspoil);
gzSpoil = mr.makeTrapezoid('z',sys,'Area',spoil_cycles * Nz * deltak(3),'Duration',Tspoil);

gxFlyback = mr.makeTrapezoid('x',sys,'Area',-gx.area,'Duration',Tpre);

areaY = ((0:Ny-1)-Ny/2)*deltak(2);
areaZ = ((0:Nz-1)-Nz/2)*deltak(3);

% Calculate timing -> assume fixed delta_TE
TE = [5, 11, 17, 23, 29] * 1e-3;
delta_TE = TE(2) - TE(1);    % for delay calculation

num_TE = length(TE);
TR = 35e-3;

Ry = 2;
num_acs = 32;

ky_indices_acs = 1+Ny/2-num_acs/2:Ny/2+num_acs/2;   % acs indices
ky_indices_accl = 1:Ry:Ny;
 
% all of the sampled ky lines:
ky_indices = union(ky_indices_acs, ky_indices_accl);


% assume that gxPre, gyPre, gzPre can be played within Tpre simultaneously
delayTE1 = ceil((TE(1) - mr.calcDuration(gz) + mr.calcRfCenter(rf) + rf.delay - mr.calcDuration(gxPre)  ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;

assert(delayTE1>0, 'delayTE1 needs to be non-negative')

% delay for the later TEs
delayTE_next = ceil( (delta_TE - Tpre - mr.calcDuration(gx)) / seq.gradRasterTime ) * seq.gradRasterTime;

assert(delayTE_next>0, 'delayTE_next needs to be non-negative')


% also assume gzPre (including Gz blip and rf gradient rephaser) can be played within the same duration as gxPre

% assume that both prewinders and flyback gradients are played within Tpre,
% and spoiler requires Tspoil
delayTR = ceil((TR - mr.calcDuration(gz) - Tpre * num_TE ...
    - mr.calcDuration(gx) * num_TE - mr.calcDuration(gxSpoil) - delayTE1 - delayTE_next * (num_TE-1) )/seq.gradRasterTime)*seq.gradRasterTime;

assert(delayTR>0, 'delayTR needs to be non-negative')

dTE1 = mr.makeDelay(delayTE1);
dTEn = mr.makeDelay(delayTE_next);

dTR = mr.makeDelay(delayTR);


% -------------------------------------------------------------------------
% PE Lines logic for online recon
% -------------------------------------------------------------------------

accelFactorPE = Ry;
ACSnum = num_acs;
centerLineIdx = floor(Ny/2) + 1 ; % index of the center k-space line, starting from 1.

count = 1 ;
clear PEsamp_u ;

for i = 1:Ny
    if ( mod(i-centerLineIdx, accelFactorPE)==0 )
        PEsamp_u(count) = i ;
        count = count + 1 ;
    end
end

minPATRefLineIdx = centerLineIdx - ACSnum/2 ; % mininum PAT line starting from 1
maxPATRefLineIdx = centerLineIdx + floor(ACSnum-1)/2 ; % maximum PAT line starting from 1
PEsamp_ACS = minPATRefLineIdx : maxPATRefLineIdx ; % GRAPPA autocalibration lines

PEsamp = union(PEsamp_u, PEsamp_ACS) ; % actually sampled lines
nPEsamp = length(PEsamp) ; % number of actually sampled
PEsamp_INC = diff([PEsamp, PEsamp(end)]) ;


% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------


% Make trapezoids for inner loop to save computation
clear gyPre gyReph;

for iY = ky_indices
    gyPre(iY) = mr.makeTrapezoid('y','Area',areaY(iY),'Duration',Tpre);
    gyReph(iY) = mr.makeTrapezoid('y','Area',-areaY(iY),'Duration',Tpre);
end


% preregister constant objects to accelerate computations
% this is not necessary, but accelerates the sequence creation by up to a factor of 2
% there is one more place in the second loop
gxPre.id = seq.registerGradEvent(gxPre);
gx.id = seq.registerGradEvent(gx);

gxFlyback.id = seq.registerGradEvent(gxFlyback);

gxSpoil.id = seq.registerGradEvent(gxSpoil);
gzSpoil.id = seq.registerGradEvent(gzSpoil);


[~, rf.shapeIDs] = seq.registerRfEvent(rf); % the phase of the RF object will change, therefore we only per-register the shapes 

for iY = ky_indices
    gyPre(iY).id = seq.registerGradEvent(gyPre(iY));
    gyReph(iY).id = seq.registerGradEvent(gyReph(iY));
end

% Create a Z-gradient that ONLY does the RF rephasing for the dummy loop
gzReph_dummy = mr.makeTrapezoid('z', sys, 'Area', gz_rf_reph.area, 'Duration', Tpre);


% rf spoiling parameters
rfSpoilingInc = 84;

rf_phase = 0;
rf_inc = 0;


% Drive magnetization to the steady state
for iY = 1:Ndummy
    % RF
    % RF spoiling phase increment = 84° for smoother transient decay, https://doi.org/10.1002/mrm.1910350216, 169° for diffusion independent rf spoiling in steady-state https://doi.org/10.1371/journal.pone.0324455

    rf.phaseOffset = rf_phase/180*pi;
    adc.phaseOffset = rf_phase/180*pi;

    rf_inc = mod(rf_inc + rfSpoilingInc, 360.0);
    rf_phase = mod(rf_phase + rf_inc, 360.0);
    
    seq.addBlock(rf,gz);

    % Gradients
    seq.addBlock(gxPre,gyPre(centerLineIdx),gzReph_dummy);
    seq.addBlock(dTE1);
    seq.addBlock(gx);
    
    for t = 2:num_TE
        seq.addBlock(gxFlyback);
        seq.addBlock(dTEn);
        seq.addBlock(gx);
    end

    seq.addBlock(gyReph(centerLineIdx),gxSpoil,gzSpoil);
    seq.addBlock(dTR);
end

% define labels
lblSetRefScan = mr.makeLabel('SET','REF', true) ;
lblSetRefAndImaScan = mr.makeLabel('SET','IMA', true) ;
lblResetRefScan = mr.makeLabel('SET','REF', false) ;
lblResetRefAndImaScan = mr.makeLabel('SET','IMA', false) ;

% register labels
lblSetRefScan.id=seq.registerLabelEvent(lblSetRefScan);
lblSetRefAndImaScan.id=seq.registerLabelEvent(lblSetRefAndImaScan);
lblResetRefScan.id=seq.registerLabelEvent(lblResetRefScan);
lblResetRefAndImaScan.id=seq.registerLabelEvent(lblResetRefAndImaScan);

% Add noise scans.
seq.addBlock(mr.makeLabel('SET', 'LIN', 0),mr.makeLabel('SET','PAR', 0)) ;
seq.addBlock(adc, mr.makeLabel('SET', 'NOISE', true),lblResetRefScan,lblResetRefAndImaScan) ;
seq.addBlock(mr.makeLabel('SET', 'NOISE', false)) ;


% pre-make labels
% lbl_ref_true = mr.makeLabel('SET', 'REF', 1);
% lbl_ref_false = mr.makeLabel('SET', 'REF', 0);

lbl_eco = [];
for t = 1:num_TE
    lbl_eco{t} = mr.makeLabel('SET', 'ECO', t-1);
end


lbl_lin = [];
for iY = ky_indices
    lbl_lin{iY} = mr.makeLabel('SET', 'LIN', iY - 1); 
end


cnt_adc = 0;

% Loop over phase encodes and define sequence blocks
tic
for iZ = 1:Nz
    disp(['iZ: ', num2str(iZ)])

    gzPre = mr.makeTrapezoid('z','Area',areaZ(iZ) + gz_rf_reph.area, 'Duration', Tpre);  % combine gz phase encode blip and rf rephaser

    % NEW: Combine the Z rewinder (-areaZ(iZ)) and Z spoiler (spoil_cycles * Nz * deltak(3))
    gzRewindAndSpoil = mr.makeTrapezoid('z', sys, 'Area', -areaZ(iZ) + gzSpoil.area, 'Duration', Tspoil);

    % optional pre-registration for acceleration
    gzPre.id = seq.registerGradEvent(gzPre);
    gzRewindAndSpoil.id = seq.registerGradEvent(gzRewindAndSpoil);

    lbl_par = mr.makeLabel('SET', 'PAR', iZ - 1);

    for iY = ky_indices
        % RF spoiling
        % RF spoiling phase increment = 84° for smoother transient decay, https://doi.org/10.1002/mrm.1910350216, 169° for diffusion independent rf spoiling in steady-state https://doi.org/10.1371/journal.pone.0324455
        
        rf.phaseOffset = rf_phase/180*pi;
        adc.phaseOffset = rf_phase/180*pi;

        rf_inc = mod(rf_inc + rfSpoilingInc, 360.0);
        rf_phase = mod(rf_phase + rf_inc, 360.0);

        % Excitation
        seq.addBlock(rf,gz);
        
        % Encoding
        seq.addBlock(gxPre,gyPre(iY),gzPre);
        seq.addBlock(dTE1);


        if ismember(iY,PEsamp_ACS)
            if ismember(iY,PEsamp_u)
                seq.addBlock(lblSetRefAndImaScan, lblSetRefScan) ;
            else
                seq.addBlock(lblResetRefAndImaScan, lblSetRefScan) ;
            end
        else
            seq.addBlock(lblResetRefAndImaScan, lblResetRefScan) ;
        end


        % seq.addBlock(gx,adc);
        seq.addBlock(gx, adc, lbl_lin{iY}, lbl_par, lbl_eco{1});

        cnt_adc = cnt_adc+1;
        
        for t = 2:num_TE
            seq.addBlock(gxFlyback);
            seq.addBlock(dTEn);

            % seq.addBlock(gx,adc);
            seq.addBlock(gx, adc, lbl_lin{iY}, lbl_par, lbl_eco{t});

            cnt_adc = cnt_adc+1;
        end

        seq.addBlock(gyReph(iY),gzRewindAndSpoil,gxSpoil);  % gxSpoil duration dominates
        seq.addBlock(dTR)
    end
end
toc

disp(['num adc:', num2str(cnt_adc)])

fprintf('Sequence ready\n');


% -------------------------------------------------------------------------
%% check whether the timing of the sequence is correct
% -------------------------------------------------------------------------

use_v141 = 1;

[ok, error_report] = seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

phaseResolution = fov(1)/Nx / (fov(2)/Ny) ;

seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'gre3d');
seq.setDefinition('AccelerationFactor', Ry);
seq.setDefinition('AccelerationFactorPE', Ry);

seq.setDefinition('kSpaceCenterLine', centerLineIdx-1) ;
seq.setDefinition('PhaseResolution', phaseResolution) ;

if use_v141
    seq.write_v141('gre3d_test_label_spoil_v0_v141.seq');
else
    seq.write('gre3d_test_lable.seq');
end


% -------------------------------------------------------------------------
%% plot sequence
% -------------------------------------------------------------------------

seqeyes(seq)

% === Compute and Display Total Scan Time ===
total_time_sec = seq.duration();
minutes = floor(total_time_sec / 60);
seconds = mod(total_time_sec, 60);

% Calculate the theoretical full box TRs for comparison
total_box_TRs = Nz * length(ky_indices);

disp('-----------------------------------------');
disp('Sequence successfully generated!');
fprintf('Total Scan Time: %d min %.1f sec\n', minutes, seconds);
fprintf('Total Dummies: %d\n', Ndummy);
disp('-----------------------------------------');




% -------------------------------------------------------------------------
%% Reconstruct 3D Multi-Echo GRE Pulseq Data with GRAPPA
% Requires mapVBVD to be in your MATLAB path.
% -------------------------------------------------------------------------

addpath(genpath('/cluster/berkin/berkin/Matlab_Code_New/LIBRARY/'))

% 1. Load the raw data using mapVBVD
data_path = '/autofs/space/marduk_001/users/berkin/2026_03_31_bay5_phantom_qsm_harmonization/';

filename = 'meas_MID00234_FID12436_pulseq151fix_gre3d_test_lable'; % Replace with your .dat file

disp('Loading data with mapVBVD...');
twix = mapVBVD2([data_path, filename, '.dat']);

% Handle multi-raid files (often the actual data is the last cell)
if iscell(twix)
    twix = twix{end};
end

% 2. Extract k-space data
% mapVBVD default squeezed dimensions are usually: [Col, Cha, Lin, Par, Eco]
% Col = Nx (Readout)
% Cha = Nc (Coils)
% Lin = Ny (Phase Encoding Y)
% Par = Nz (Partition/Phase Encoding Z)
% Eco = Ne (Echoes)
disp('Extracting k-space...');

tic
    kSpace = twix.image(); 
    kSpace_ref = twix.refscan(); 
toc

kSpace = squeeze(kSpace); 
kSpace_ref = squeeze(kSpace_ref); 

% -------------------------------------------------------------------------
%%
% -------------------------------------------------------------------------

kSpace = permute(kSpace, [1, 3, 4, 2, 5]);
[~, ~, ~, Nc, Ne] = size(kSpace);

kSpace_ref = permute(kSpace_ref, [1, 3, 4, 2, 5]);
img_ref = ifft3call(kSpace_ref);    


% -------------------------------------------------------------------------
%% Reorder dimensions to [Nx, Ny, Nz, Nc, Ne] for easier mapping
% From [Col, Cha, Lin, Par, Eco] -> [1, 3, 4, 2, 5]
% -------------------------------------------------------------------------

% Sequence Parameters (Must match your Pulseq script)
R = 2; 
acs_lines = 32; 

Ny = 184;
ky_indices_acs = 1+Ny/2-acs_lines/2:Ny/2+acs_lines/2;

kSpace_all = kSpace;

kSpace_all(:,ky_indices_acs,:,:,:) = kSpace_ref; %

kSpace_all = padarray(kSpace_all, rem(size(kSpace_all(:,:,:,1,1)),2), 'post');

[Nx,Ny,Nz] = size(kSpace_all(:,:,:,1,1));  % update k-space size after zero padding to even size

img_zf = ifft3call(kSpace_all);

for t = 1:Ne
    imagesc3d2( rsos(kSpace_all(:,:,:,:,t),4), s(img_zf)/2, t, [90,90,90], [0.,1e-3]), setGcf(.6)
    imagesc3d2( rsos(img_zf(:,:,:,:,t),4), s(img_zf)/2, 10+t, [90,90,0], [0.,5e-4]), setGcf(.6)
end


%--------------------------------------------------------------------------
%% SVD coil compression: use 1st chan as body coil reference
%--------------------------------------------------------------------------

num_chan = 32;
select_acs_echo = 1;

[img_ref_svd, cmp_mtx] = svd_compress3d(img_ref(:,:,:,:,select_acs_echo), num_chan, 1);

rmse(rsos(img_ref_svd,4),rsos(img_ref(:,:,:,:,select_acs_echo),4))

% apply svd
kSpace_all_svd = zeross([Nx,Ny,Nz,num_chan,Ne]);

for t = 1:Ne
    tmp = svd_apply3d(kSpace_all(:,:,:,:,t),cmp_mtx);
    kSpace_all_svd(:,:,:,:,t) = tmp; 
end

rmse(rsos(kSpace_all_svd,4), rsos(kSpace_all,4))


%--------------------------------------------------------------------------
%% grappa: recon single slice
%--------------------------------------------------------------------------

addpath '/autofs/cluster/berkin/berkin/Matlab_Code_New/TOOLBOXES/Espirit_matlab_only_toolbox/utils'

substitute_acs = 1;

Ry = 2;
Rz = 1;


lambda_percent = 1e-3;      % percentage of sigma_min to use as regularizer

num_acs = floor([size(kSpace_ref,2), Nz] / 2) * 2 - 2;
kernel_size = [3,3];        % odd kernel size

% espirit parameters
num_acs_espirit = acs_lines;
kernel_size_espirit = [6,6];
eigen_thresh = 0.5;

k_hybrid = ifftc(kSpace_all_svd,1);

Img_Grappa_rsos = zeross([Nx,Ny,Nz,Ne]);
Img_Grappa_combo = zeross([Nx,Ny,Nz,Ne]);

tic
for slc_select = 132%1:Nx
    disp(slc_select)

    k_hybrid_slc = sq(k_hybrid(slc_select,:,:,:,:));

    k_hybrid_acs = zeross(size(k_hybrid_slc));
    k_hybrid_acs(1+end/2-num_acs(1)/2:end/2+num_acs(1)/2, 1+end/2-num_acs(2)/2:end/2+num_acs(2)/2,:,:) = ...
            k_hybrid_slc(1+end/2-num_acs(1)/2:end/2+num_acs(1)/2, 1+end/2-num_acs(2)/2:end/2+num_acs(2)/2,:,:);

    % espirit
    [maps, weights] = ecalib_soft( k_hybrid_acs(:,:,:,select_acs_echo), num_acs_espirit, kernel_size_espirit, eigen_thresh );

    receive_slc = dot_mult(maps, weights >= eigen_thresh);


    Im_Grappa = zeros([Ny, Nz, Nc, Ne]);    
    Im_Grappa_combo = zeros([Ny, Nz, Ne]);    
    
    for t = 1:Ne
        Im_Grappa(:,:,:,t) = grappa_gfactor_2d_jvc3( k_hybrid_slc(:,:,:,t), k_hybrid_acs(:,:,:,t), Ry, Rz, num_acs, kernel_size, lambda_percent, substitute_acs );
    end
    
    % substitute acs
    k_Grappa = fft2call(Im_Grappa);
    k_Grappa(ky_indices_acs,:,:,:) = k_hybrid_slc(ky_indices_acs,:,:,:);
    Im_Grappa = ifft2call(k_Grappa);

    for t = 1:Ne
        Im_Grappa_combo(:,:,t) = coil_combine(Im_Grappa(:,:,:,t),receive_slc,3);
    end

    Img_Grappa_rsos(slc_select,:,:,:) = rsos(Im_Grappa,3);
    Img_Grappa_combo(slc_select,:,:,:) = Im_Grappa_combo;
end
toc


close all

mosaic( rsos(ifft2call(k_hybrid_acs),3), 1, 5, 20, '', [0,3e-4], 90 ), setGcf(.7)
mosaic( rsos(ifft2call(k_hybrid_slc),3), 1, 5, 21, '', [0,3e-4], 90 ), setGcf(.7)
mosaic( Img_Grappa_rsos(slc_select,:,:,:), 1, 5, 22, '', [0,3e-4], 90 ), setGcf(.7)

mosaic( Img_Grappa_combo(slc_select,:,:,:), 1, 5, 23, '', [0,3e-4], 90 ), setGcf(.7)
mosaic( angle(Img_Grappa_combo(slc_select,:,:,:)), 1, 5, 24, '', [-pi,pi], 90 ), setGcf(.7)

mosaic( rsos(k_hybrid_acs,3), 1, 5, 30, '', [0,1e-4] )
mosaic( rsos(k_hybrid_slc,3), 1, 5, 31, '', [0,1e-4] )



%--------------------------------------------------------------------------
%% grappa: recon single slice with parfor
%--------------------------------------------------------------------------

% --- 1. Parallel Pool Setup ---
num_cores = 16; % Specify the number of cores you want to use

% Check if a parallel pool already exists
poolobj = gcp('nocreate'); 

% If no pool exists, or if the current pool has the wrong number of workers, start a new one
if isempty(poolobj)
    parpool('local', num_cores);
elseif poolobj.NumWorkers ~= num_cores
    disp('Restarting parallel pool with the requested number of cores...');
    delete(poolobj);
    parpool('local', num_cores);
end
% ------------------------------

substitute_acs = 1;

Ry = 2;
Rz = 1;

lambda_percent = 1e-3;      % percentage of sigma_min to use as regularizer

num_acs = floor([size(kSpace_ref,2), Nz] / 2) * 2 - 2;
kernel_size = [3,3];        % odd kernel size

% espirit parameters
num_acs_espirit = acs_lines;
kernel_size_espirit = [6,6];
eigen_thresh = 0.5;

k_hybrid = ifftc(kSpace_all,1);

Img_Grappa_rsos = zeros([Nx, Ny, Nz, Ne]);
Img_Grappa_combo = zeros([Nx, Ny, Nz, Ne]);

tic
parfor slc_select = 1:Nx
    % disp() inside parfor prints out of order, which is normal
    disp(['Reconstructing slice: ', num2str(slc_select)]);

    % Extract slice using squeeze to ensure proper slicing and save RAM
    k_hybrid_slc = sq(k_hybrid(slc_select,:,:,:,:));

    k_hybrid_acs = zeross(size(k_hybrid_slc));
    k_hybrid_acs(1+end/2-num_acs(1)/2:end/2+num_acs(1)/2, 1+end/2-num_acs(2)/2:end/2+num_acs(2)/2,:,:) = ...
            k_hybrid_slc(1+end/2-num_acs(1)/2:end/2+num_acs(1)/2, 1+end/2-num_acs(2)/2:end/2+num_acs(2)/2,:,:);
        
    % espirit
    [maps, weights] = ecalib_soft( k_hybrid_acs(:,:,:,select_acs_echo), num_acs_espirit, kernel_size_espirit, eigen_thresh );

    receive_slc = dot_mult(maps, weights >= eigen_thresh);

    receive_slc = abs(receive_slc) .* exp(1i * angle( receive_slc .* repmat(conj(receive_slc(:,:,1)), [1,1,num_chan]) ));

    % Temporary variable inside the parfor loop for the current slice
    Im_Grappa = zeros([Ny, Nz, Nc, Ne]);    
    Im_Grappa_combo = zeros([Ny, Nz, Ne]);    
    
    for t = 1:Ne
        % Perform GRAPPA on the current echo/time-point
        Im_Grappa(:,:,:,t) = grappa_gfactor_2d_jvc3(...
            k_hybrid_slc(:,:,:,t), ...
            k_hybrid_acs(:,:,:,t), ...
            Ry, Rz, num_acs, kernel_size, lambda_percent, substitute_acs);
    end

    % substitute acs
    k_Grappa = fft2call(Im_Grappa);
    k_Grappa(ky_indices_acs,:,:,:) = k_hybrid_slc(ky_indices_acs,:,:,:);
    Im_Grappa = ifft2call(k_Grappa);

    for t = 1:Ne
        Im_Grappa_combo(:,:,t) = coil_combine(Im_Grappa(:,:,:,t),receive_slc,3);
    end


    % Root sum of squares over the coil dimension (dim 3)
    Img_Grappa_rsos(slc_select, :, :, :) = rsos(Im_Grappa, 3);

    % espirit coil combination
    Img_Grappa_combo(slc_select,:,:,:) = Im_Grappa_combo;
end
toc

% --- 2. Release Cores ---
disp('Shutting down parallel pool and releasing cores...');
delete(gcp('nocreate'));
% ------------------------

for t = 1:Ne
    imagesc3d2( Img_Grappa_rsos(:,:,:,t), s(Img_Grappa_rsos)/2, t, [90,90,90], [0.,5e-4]), setGcf(.6)
    imagesc3d2( Img_Grappa_combo(:,:,:,t), s(Img_Grappa_rsos)/2, 10+t, [90,90,90], [0.,5e-4]), setGcf(.6)
    imagesc3d2( angle(Img_Grappa_combo(:,:,:,t)), s(Img_Grappa_rsos)/2, 20+t, [90,90,90], [-pi,pi]), setGcf(.6)
end

tic
    save([data_path, filename, 'Img_Grappa_rsos.mat'], 'Img_Grappa_rsos', '-v7.3')
    save([data_path, filename, 'Img_Grappa_combo.mat'], 'Img_Grappa_combo', '-v7.3')
toc

genNiigz(abs(Img_Grappa_combo(:,:,:,end)),[1,1,1],[data_path, filename, 'Img_Grappa_combo.nii'])