clc; clear; close all;

%% Radar and Array Parameters
numPulses = 64;     % Number of pulses in the data cube
numChannels= 8;  
C=physconst('LightSpeed');
fc=1e9;
antenna=phased.IsotropicAntennaElement('FrequencyRange',[800e6 1.2e9]);
Array=phased.ULA('NumElements',numChannels,'ElementSpacing',0.5*(C/fc),'Element',antenna);

fs= 20e6;      % Sampling rate
PRF=10e3;      % Pulse repetition frequency
pw=33e-6;      % Pulse width

%% Waveform
waveform=phased.RectangularWaveform('SampleRate',fs,'PRF',PRF,'PulseWidth',pw);
refPulse = waveform();  % Reference pulse for matched filter

% Matched Filter
matchedFilter = phased.MatchedFilter('Coefficients',conj(flipud(refPulse)));

% Plot transmitted pulse
t   = (0:length(refPulse)-1)/fs;
figure;
plot(t*1e6, abs(refPulse));
xlabel('Time (\mus)'); ylabel('Amplitude');
title('Transmitted Rectangular Pulse');

%% Transmitter, Target, Channel, Collector, Receiver
tx_gain = 20; % dB
transmitter = phased.Transmitter('Gain', tx_gain, 'PeakPower', 5000);

target_range = 10000; % 10 km
target_rcs = 5; % m²
target = phased.RadarTarget('MeanRCS', target_rcs, 'OperatingFrequency', fc);

channel = phased.FreeSpace('SampleRate', fs,'TwoWayPropagation',false,'OperatingFrequency', fc);

collector = phased.Collector('Sensor',Array,'OperatingFrequency',fc);

receiver = phased.ReceiverPreamp('Gain', 0, 'NoiseFigure', 10, 'SampleRate', fs);

%% Initialize Data Cube
samplesPerPulse = round(fs/PRF); % Range bins per PRI
dataCube = complex(zeros(samplesPerPulse, numPulses, numChannels));
mfCube   = complex(zeros(samplesPerPulse+length(refPulse)-1, numPulses, numChannels));

%% Platform Setup
radarPos = [0; 0; 0]; radarVel = [0; 0; 0];
targetPos = [target_range; 0; 0]; targetVel = [0; 0; 0];

%% Loop over pulses
for pulseIdx = 1:numPulses
    % Transmit pulse
    sig = waveform();
    txsig = transmitter(sig);

    % Propagation -> Target
    prop_signal = channel(txsig, radarPos, targetPos, radarVel, targetVel);

    % Reflection
    refl_signal = target(prop_signal);

    % Back to radar
    rx_signal = channel(refl_signal, targetPos, radarPos, targetVel, radarVel);

    % Collect across array
    rx_array = collector(rx_signal, [0; 0]); % broadside

    % Receiver noise
    rx_array = receiver(rx_array);

    % Store raw received data
    dataCube(:,pulseIdx,:) = rx_array(1:samplesPerPulse,:);

    % Apply matched filter (pulse compression)
    mfSig = matchedFilter(rx_array);
    mfCube(:,pulseIdx,:) = mfSig;
end

%% Plots
figure;
imagesc(abs(dataCube(:,:,1)));
xlabel('Pulse Index (Slow Time)');
ylabel('Range Bin (Fast Time)');
title('Raw Radar Data Cube (Channel 1)');
colorbar;

figure;
imagesc(abs(mfCube(:,:,1)));
xlabel('Pulse Index (Slow Time)');
ylabel('Range Bin (Fast Time)');
title('Matched Filtered Data Cube (Channel 1)');
colorbar;

%% Range Profile after compression (one pulse, one channel)
figure;
plot(abs(mfCube(:,10,1)));
xlabel('Range Bin');
ylabel('Amplitude');
title('Range Profile after Matched Filtering');
grid on;
