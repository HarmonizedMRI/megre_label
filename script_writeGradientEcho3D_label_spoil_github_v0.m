% -------------------------------------------------------------------------
%%
% -------------------------------------------------------------------------

clear
restoredefaultpath

% here add pulseq matlab directory to path e.g.
addpath pulseq-version-newer-than-20260519/matlab/
addpath path/to/pulceq/v2.5.2.0/matlab % only needed for GE scanner




% Define FOV and resolution
fov = [264e-3 184e-3 144e-3];    

Nx = 264;
Ny = 184/1; 
Nz = 144/1;            

adc_dwell = 14e-6;

Tread_adc = adc_dwell * Nx;
Tread = ceil(Tread_adc / 20e-6) * 20e-6;

disp(['Tread ', num2str(Tread*1e3), ' ms'])
disp(['ADC duration ', num2str(Tread_adc*1e3), ' ms'])

Ndummy = 50;

% define system properties
Gmax = 20;
Smax = 100;

sys = mr.opts('MaxGrad', Gmax, 'GradUnit', 'mT/m', ...
    'MaxSlew', Smax, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, ...
    'adcRasterTime', 2e-6, 'rfRasterTime', 2e-6, ...
    'gradRasterTime', 20e-6, 'blockDurationRaster', 20e-6);

seq = mr.Sequence(sys);           % Create a new sequence object

assert(ismethod(seq, 'addTRID'), 'seq must have an addTRID method. Please download a new version of matlab pulseq toolbox');

alpha = 15;     % flip angle
rf_duration = 1.e-3;
rf_tbp = 20;    % RF time bandwidth
thickness = 128e-3;   % slab thickness to excite

% Create alpha-degree slice selection pulse and gradient
[rf, gz, gz_rf_reph] = mr.makeSincPulse(alpha*pi/180, 'Duration', rf_duration,...
    'SliceThickness',thickness,'apodization',0.42,'timeBwProduct',rf_tbp,'use','excitation','system',sys);


% -------------------------------------------------------------------------
%% Define other gradients and ADC events
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
adc = mr.makeAdc(Nx,sys,'Duration',Tread_adc,'Delay',gx.riseTime);

gxPre = mr.makeTrapezoid('x',sys,'Area',-gx.area/2,'Duration',Tpre);

gxSpoil = mr.makeTrapezoid('x',sys,'Area',spoil_cycles * Nx * deltak(1),'Duration',Tspoil);
gzSpoil = mr.makeTrapezoid('z',sys,'Area',spoil_cycles * Nz * deltak(3),'Duration',Tspoil);

gxFlyback = mr.makeTrapezoid('x',sys,'Area',-gx.area,'Duration',Tpre);

areaZ = ((0:Nz-1)-Nz/2)*deltak(3);

% NEW: Flip Y-axis (phase encode)
areaY = -((0:Ny-1)-Ny/2)*deltak(2);


% Calculate timing -> assume fixed delta_TE
TE = [5, 11, 17, 23, 29] * 1e-3;
delta_TE = TE(2) - TE(1);    % for delay calculation

num_TE = length(TE);
TR = 35e-3;

Ry = 2;
% set to Ry = 1 if fully sampled acquisition is desired
% Ry = 1;

if Ry > 1
    num_acs = 32;

    ky_indices_acs = 1+Ny/2-num_acs/2:Ny/2+num_acs/2;   % acs indices
    ky_indices_accl = 1:Ry:Ny;
     
    % all of the sampled ky lines:
    ky_indices = union(ky_indices_acs, ky_indices_accl);
else
    ky_indices = 1:Ny;    
end



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
%% PE Lines logic for online recon
% -------------------------------------------------------------------------

accelFactorPE = Ry;
centerLineIdx = floor(Ny/2) + 1 ; % index of the center k-space line, starting from 1.

count = 1 ;
clear PEsamp_u ;

for i = 1:Ny
    if ( mod(i-centerLineIdx, accelFactorPE)==0 )
        PEsamp_u(count) = i ;
        count = count + 1 ;
    end
end

if Ry > 1
    ACSnum = num_acs;

    minPATRefLineIdx = centerLineIdx - ACSnum/2 ; % mininum PAT line starting from 1
    maxPATRefLineIdx = centerLineIdx + floor(ACSnum-1)/2 ; % maximum PAT line starting from 1
    PEsamp_ACS = minPATRefLineIdx : maxPATRefLineIdx ; % GRAPPA autocalibration lines
    
    PEsamp = union(PEsamp_u, PEsamp_ACS) ; % actually sampled lines
else
    PEsamp = PEsamp_u;
end

nPEsamp = length(PEsamp) ; % number of actually sampled
PEsamp_INC = diff([PEsamp, PEsamp(end)]) ;


% -------------------------------------------------------------------------
%%
% -------------------------------------------------------------------------


% GE conversion expects repeated waveform shapes. Define one base PE
% trapezoid and scale it inside the scan loop.
gyPreBase = mr.makeTrapezoid('y', sys, 'Area', max(abs(areaY(ky_indices))), 'Duration', Tpre);
gyRephBase = mr.makeTrapezoid('y', sys, 'Area', max(abs(areaY(ky_indices))), 'Duration', Tpre);
gyPreScales = areaY / gyPreBase.area;
gyRephScales = -areaY / gyRephBase.area;




[~, rf.shapeIDs] = seq.registerRfEvent(rf); % the phase of the RF object will change, therefore we only per-register the shapes 

% Create a Z-gradient that ONLY does the RF rephasing for the dummy loop
gzReph_dummy = mr.makeTrapezoid('z', sys, 'Area', gz_rf_reph.area, 'Duration', Tpre);


% -------------------------------------------------------------------------
%% change orientation to match the siemens product sequence for ICE recon
% -------------------------------------------------------------------------


gz = mr.scaleGrad(gz, -1) ;
gz_rf_reph = mr.scaleGrad(gz_rf_reph, -1) ;
gzReph_dummy = mr.scaleGrad(gzReph_dummy, -1) ;
gzSpoil = mr.scaleGrad(gzSpoil, -1) ;


% NEW: Flip X-axis (Readout)
gx = mr.scaleGrad(gx, -1);
gxPre = mr.scaleGrad(gxPre, -1);
gxFlyback = mr.scaleGrad(gxFlyback, -1);
gxSpoil = mr.scaleGrad(gxSpoil, -1);


% correction: pre-registration of gradient events needs to be done after the scaling, 
% otherwise the pre-registered events will not reflect the flipped gradient waveforms and cause errors during sequence creation
% preregister constant objects to accelerate computations
% this is not necessary, but accelerates the sequence creation by up to a factor of 2
% there is one more place in the second loop
gxPre.id = seq.registerGradEvent(gxPre);
gx.id = seq.registerGradEvent(gx);

gxFlyback.id = seq.registerGradEvent(gxFlyback);

gxSpoil.id = seq.registerGradEvent(gxSpoil);
gzSpoil.id = seq.registerGradEvent(gzSpoil);



% -------------------------------------------------------------------------
%% build sequence
% -------------------------------------------------------------------------


% rf spoiling parameters
rfSpoilingInc = 84;

rf_phase = 0;
rf_inc = 0;

pislquant_adc_count = 0;
siemensOnlineReconOff = mr.makeLabel('SET', 'OFF', true);
siemensOnlineRecon = mr.makeLabel('SET', 'OFF', false);

% Drive magnetization to the steady state
for iY = 1:Ndummy
    % RF
    % RF spoiling phase increment = 84° for smoother transient decay, https://doi.org/10.1002/mrm.1910350216, 169° for diffusion independent rf spoiling in steady-state https://doi.org/10.1371/journal.pone.0324455

    rf.phaseOffset = rf_phase/180*pi;
    adc.phaseOffset = rf_phase/180*pi;

    rf_inc = mod(rf_inc + rfSpoilingInc, 360.0);
    rf_phase = mod(rf_phase + rf_inc, 360.0);
    
    seq.addTRID('receive_gain_calib');
    seq.addBlock(rf,gz);

    % Gradients
    centerScalePre = gyPreScales(centerLineIdx);
    centerScalePre = centerScalePre + (centerScalePre == 0) * eps;
    centerScaleReph = gyRephScales(centerLineIdx);
    centerScaleReph = centerScaleReph + (centerScaleReph == 0) * eps;

    seq.addBlock(gxPre,mr.scaleGrad(gyPreBase, centerScalePre),gzReph_dummy);
    seq.addBlock(dTE1);
    seq.addBlock(gx, adc, siemensOnlineReconOff);
    pislquant_adc_count = pislquant_adc_count + 1;
    
    for t = 2:num_TE
        seq.addBlock(gxFlyback);
        seq.addBlock(dTEn);
        seq.addBlock(gx, adc);
        pislquant_adc_count = pislquant_adc_count + 1;
    end

    seq.addBlock(mr.scaleGrad(gyRephBase, centerScaleReph),gxSpoil,gzSpoil);
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
adc_dur_round_up = ceil(mr.calcDuration(adc)/seq.blockDurationRaster)*seq.blockDurationRaster;
seq.addTRID('noise_scan');
seq.addBlock(mr.makeLabel('SET', 'LIN', 0),mr.makeLabel('SET','PAR', 0)) ;
% seq.addBlock(adc, mr.makeLabel('SET', 'NOISE', true),lblResetRefScan,lblResetRefAndImaScan) ;
seq.addBlock(adc_dur_round_up,siemensOnlineRecon,adc, mr.makeLabel('SET', 'NOISE', true),lblResetRefScan,lblResetRefAndImaScan) ;
seq.addBlock(mr.makeLabel('SET', 'NOISE', false)) ;


% pre-make labels

lbl_eco = [];
for t = 1:num_TE
    lbl_eco{t} = mr.makeLabel('SET', 'ECO', t-1);
end


lbl_lin = [];
for iY = ky_indices
    lbl_lin{iY} = mr.makeLabel('SET', 'LIN', iY - 1); 
end


cnt_adc = 0;

gzPreAreas = -areaZ + gz_rf_reph.area;
gzRewindAndSpoilAreas = areaZ + gzSpoil.area;
gzPreBase = mr.makeTrapezoid('z', sys, 'Area', max(abs(gzPreAreas)), 'Duration', Tpre);
gzRewindAndSpoilBase = mr.makeTrapezoid('z', sys, 'Area', max(abs(gzRewindAndSpoilAreas)), 'Duration', Tspoil);

% Loop over phase encodes and define sequence blocks
tic
for iZ = 1:Nz
    disp(['iZ: ', num2str(iZ)])

    gzPreScale = gzPreAreas(iZ) / gzPreBase.area;
    gzPreScale = gzPreScale + (gzPreScale == 0) * eps;
    gzRewindAndSpoilScale = gzRewindAndSpoilAreas(iZ) / gzRewindAndSpoilBase.area;
    gzRewindAndSpoilScale = gzRewindAndSpoilScale + (gzRewindAndSpoilScale == 0) * eps;

    lbl_par = mr.makeLabel('SET', 'PAR', iZ - 1);

    for iY = ky_indices
        % RF spoiling
        % RF spoiling phase increment = 84° for smoother transient decay, https://doi.org/10.1002/mrm.1910350216, 169° for diffusion independent rf spoiling in steady-state https://doi.org/10.1371/journal.pone.0324455
        
        rf.phaseOffset = rf_phase/180*pi;
        adc.phaseOffset = rf_phase/180*pi;

        rf_inc = mod(rf_inc + rfSpoilingInc, 360.0);
        rf_phase = mod(rf_phase + rf_inc, 360.0);

        % Excitation
        seq.addTRID('imaging');
        seq.addBlock(rf,gz);
        
        % Encoding
        gyPreScale = gyPreScales(iY);
        gyPreScale = gyPreScale + (gyPreScale == 0) * eps;
        gyRephScale = gyRephScales(iY);
        gyRephScale = gyRephScale + (gyRephScale == 0) * eps;

        seq.addBlock(gxPre,mr.scaleGrad(gyPreBase, gyPreScale),mr.scaleGrad(gzPreBase, gzPreScale));
        seq.addBlock(dTE1);


        if Ry > 1
            if ismember(iY,PEsamp_ACS)
                if ismember(iY,PEsamp_u)
                    seq.addBlock(lblSetRefAndImaScan, lblSetRefScan) ;
                else
                    seq.addBlock(lblResetRefAndImaScan, lblSetRefScan) ;
                end
            else
                seq.addBlock(lblResetRefAndImaScan, lblResetRefScan) ;
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

        seq.addBlock(mr.scaleGrad(gyRephBase, gyRephScale),mr.scaleGrad(gzRewindAndSpoilBase, gzRewindAndSpoilScale),gxSpoil);  % gxSpoil duration dominates
        seq.addBlock(dTR)
    end
end
toc

disp(['num adc:', num2str(cnt_adc)])

% Add noise scans after pislquant and imaging ADCs.
seq.addTRID('noise');
seq.addBlock(mr.makeLabel('SET', 'LIN', 0),mr.makeLabel('SET','PAR', 0)) ;
seq.addBlock(adc, mr.makeDelay(3.8e-3), mr.makeLabel('SET', 'NOISE', true),lblResetRefScan,lblResetRefAndImaScan) ;
seq.addBlock(mr.makeLabel('SET', 'NOISE', false)) ;

fprintf('Sequence ready\n');


% -------------------------------------------------------------------------
%% check whether the timing of the sequence is correct
% -------------------------------------------------------------------------

use_v141 = 0;  % set to 1 to use old interpreter

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
seq.setDefinition('num_adcs_pislquant', pislquant_adc_count);
seq.setDefinition('pislquant', pislquant_adc_count);

seq.setDefinition('kSpaceCenterLine', centerLineIdx-1) ;
seq.setDefinition('PhaseResolution', phaseResolution) ;

seq.setDefinition('TridIdName', strjoin(seq.tridId2Name, ',')); % for seqeyes



% Generate a date string in the format YYYYMMDD (e.g., 20260324)
dateString = char(datetime('today', 'Format', 'yyyyMMdd'));

file_path = [pwd, '/'];  % use current folder

if ~isfolder([file_path, dateString])
    mkdir([file_path, dateString])
end

% Create the full filename
fileName = [file_path, dateString, '/gre3d_test_label_spoil_v2_xyzflip_', num2str(Ry)];


% Save the sequence
if use_v141
    seq.write_v141([fileName, '_v141.seq']);
else
    seq.write([fileName, '.seq']);
end
writeceq(seq2ceq(seq), [fileName, '.pge'], 'pislquant', seq.getDefinition('pislquant')); % for GE scanner

% -------------------------------------------------------------------------
%%  sequence duration
% -------------------------------------------------------------------------


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
