%% === PARAMETERS ===
clear;
clc;

% Radar parameters
fc = 3e9;                 % 3 GHz
c = physconst('LightSpeed');
lambda = c/fc;
prf = 1000;               % PRF
nPulses = 16;             % Pulses per CPI
nElements = 8;            % ULA elements
fs = 1e6;                 % Sampling frequency
rng(2025);                % For reproducibility

% Target parameters (multiple)
tgtAngles = [20, -10];    % degrees
tgtVels   = [50, -30];    % m/s
tgtRCS    = [1, 1];       % m^2
tgtRanges = [500, 700];   % meters

% Clutter parameters
clutterCNR = 40;          % dB

% Jammer parameters
jamAngle = -40;           % degrees
jamJNR   = 30;            % dB

% Noise
noisePower = -10;         % dB

%% === ARRAY & WAVEFORM ===
array = phased.ULA('NumElements',nElements,'ElementSpacing',lambda/2);
waveform = phased.RectangularWaveform('SampleRate',fs,'PulseWidth',2e-6,'PRF',prf);
sig = waveform();
nSamples = length(sig);

% Steering helper
pos = getElementPosition(array)/lambda;

%% === SIMULATION ===
nRanges = 128; % range bins for simulation
rxCube = zeros(nRanges,nElements,nPulses);

for m = 1:nPulses
    pulseData = zeros(nRanges, nElements);
    
    % ---- Targets ----
    for k = 1:length(tgtAngles)
        % Steering in space & Doppler phase shift
        sv_space = steervec(pos, tgtAngles(k));
        dopplerPhase = exp(1j*2*pi*m*tgtVels(k)/(lambda*prf));
        
        rangeBin = round(tgtRanges(k) / (c/(2*fs))) + 1;
        pulseData(rangeBin,:) = pulseData(rangeBin,:) + ...
            sqrt(db2pow(tgtRCS(k))) * dopplerPhase * sv_space.';
    end
    
    % ---- Clutter (simplified random phase model) ----
    clutterAngle = linspace(-90,90,10);
    for cl = 1:length(clutterAngle)
        sv_space = steervec(pos, clutterAngle(cl));
        rangeBin = randi([10, nRanges-10]);
        pulseData(rangeBin,:) = pulseData(rangeBin,:) + ...
            sqrt(db2pow(clutterCNR)) * (randn(1)+1j*randn(1)) * sv_space.';
    end
    
    % ---- Jammer ----
    sv_jam = steervec(pos, jamAngle);
    pulseData(:,:) = pulseData(:,:) + ...
        sqrt(db2pow(jamJNR)) * (randn(nRanges,1)+1j*randn(nRanges,1)) * sv_jam.';
    
    % ---- Noise ----
    noise = sqrt(db2pow(noisePower)/2) * ...
        (randn(nRanges,nElements) + 1j*randn(nRanges,nElements));
    
    rxCube(:,:,m) = pulseData + noise;
end

%% === RANGE-DOPPLER BEFORE STAP ===
RD_before = zeros(nRanges, nPulses);
for r = 1:nRanges
    % Doppler FFT across pulses
    dopplerSig = fftshift(fft(squeeze(rxCube(r,:,:)), nPulses, 2), 2);
    RD_before(r,:) = sum(abs(dopplerSig),1); % non-coherent sum over elements
end

%% === STAP PROCESSING ===
RD_after = zeros(nRanges, nPulses);
cutBin = 64; % middle range bin for covariance estimate

% Training bins (exclude CUT ± guard)
guard = 2;
trainBins = setdiff(1:nRanges, cutBin-guard:cutBin+guard);

% Build training covariance
Xtrain = [];
for r = trainBins
    snapshot = squeeze(rxCube(r,:,:)); % elements x pulses
    Xtrain = [Xtrain, snapshot(:)];
end
Rhat = (Xtrain*Xtrain')/size(Xtrain,2);

% Apply STAP to each range bin
for r = 1:nRanges
    for d = 1:nPulses
        % Space-time snapshot for this bin & Doppler
        snapshot = squeeze(rxCube(r,:,:));
        x = snapshot(:);
        
        % Assume target steering vector for demonstration (20 deg, 50 m/s)
        sv_space = steervec(pos, 20);
        sv_dopp = exp(1j*2*pi*(0:nPulses-1)'*50/(lambda*prf));
        s = kron(sv_dopp, sv_space);
        
        % SMI weights
        w = Rhat \ s / (s'*(Rhat\s));
        
        % Apply
        y = w' * x;
        RD_after(r,d) = abs(y);
    end
end

%% === PLOTS ===
dopplerAxis = (-nPulses/2:nPulses/2-1)*(prf/nPulses);
rangeAxis = (0:nRanges-1)*c/(2*fs);

figure;
subplot(1,2,1);
imagesc(dopplerAxis, rangeAxis, 20*log10(abs(RD_before)));
title('Range-Doppler BEFORE STAP');
xlabel('Doppler (Hz)'); ylabel('Range (m)'); colorbar; caxis([0 60]);

subplot(1,2,2);
imagesc(dopplerAxis, rangeAxis, 20*log10(abs(RD_after)));
title('Range-Doppler AFTER STAP');
xlabel('Doppler (Hz)'); ylabel('Range (m)'); colorbar; caxis([0 60]);
