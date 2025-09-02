% === Radar Parameters ===
fc = 1e9;          % Operating frequency = 1 GHz
c = physconst('LightSpeed');  % Speed of light

PRF = 10e3;         % Pulse Repetition Frequency (10 kHz) 
pw = 2e-6;         % Pulse width = 2 µs
fs = 20e6;         % Sampling rate = 20 MHz

% === Additional Parameters for Data Cube ===
numPulses = 64;     % Number of pulses in the data cube
numChannels = 4;    % Number of receive channels (for array processing)
target_range = 10000; % 10 km
target_rcs = 1;     % 1 m² RCS

% === Define Rectangular Waveform ===
waveform = phased.RectangularWaveform(...
    'SampleRate', fs, ...
    'PulseWidth', pw, ...
    'PRF', PRF);

% === Transmitter Setup ===
tx_gain = 20; % Transmitter gain in dB
transmitter = phased.Transmitter('Gain', tx_gain, 'PeakPower', 5000);

% === Target ===
target = phased.RadarTarget('MeanRCS', target_rcs, 'OperatingFrequency', fc);

% === Free Space Channel ===
channel = phased.FreeSpace(...
    'SampleRate', fs, ...
    'TwoWayPropagation', false, ...
    'OperatingFrequency', fc);

% === Receive Array (for multiple channels) ===
array = phased.ULA('NumElements', numChannels, 'ElementSpacing', 0.5*(c/fc));
collector = phased.Collector('Sensor', array, 'OperatingFrequency', fc);
receiver = phased.ReceiverPreamp('Gain', 0, 'NoiseFigure', 10, 'SampleRate', fs);

% === Initialize Data Cube ===
samplesPerPulse = round(pw * fs); % Samples per pulse
dataCube = zeros(samplesPerPulse, numPulses, numChannels);

% === Platform Positions ===
radarPos = [0; 0; 0];
radarVel = [0; 0; 0];
targetPos = [target_range; 0; 0];
targetVel = [0; 0; 0];

% === Simulate multiple pulses ===
for pulseIdx = 1:numPulses   % 1:2:3:4:......,63:64
    % Generate the transmitted pulse
    sig = waveform();
    
    % Transmit the pulse
    transmitted_signal = transmitter(sig);
    
    % Calculate the time for this pulse
    t = (0:length(sig)-1)/fs + (pulseIdx-1)*(1/PRF);
    
    % Propagate to target
    prop_signal = channel(transmitted_signal, radarPos, targetPos, radarVel, targetVel);
    
    % Reflect off target
    reflected_signal = target(prop_signal);
    
    % Propagate back to radar
    received_signal = channel(reflected_signal, targetPos, radarPos, targetVel, radarVel);
    
    % Collect the signal with the array (simulate multiple channels)
    % CORRECTED: Using collector instead of direct array call
    received_signal = collector(received_signal, [0; 0]); % Assuming target at broadside
    
    % Add receiver noise
    received_signal = receiver(received_signal);
    
    % Store in data cube (only keep the pulse duration)
    dataCube(:, pulseIdx, :) = received_signal(1:samplesPerPulse, :);
end

% === Display Data Cube Information ===
fprintf('Data Cube Dimensions: %d (Range) x %d (Pulses) x %d (Channels)\n', ...
    size(dataCube, 1), size(dataCube, 2), size(dataCube, 3));
fprintf('Range Resolution: %.2f m\n', c/(2*fs));
fprintf('Maximum Unambiguous Range: %.2f km\n', c/(2*PRF)/1000);

% === Plot the data cube slices ===
figure;

% Range-Time plot for first channel
subplot(2,2,1);
range_axis = (0:samplesPerPulse-1)*c/(2*fs);
imagesc(1:numPulses, range_axis, 20*log10(abs(dataCube(:,:,1))));
xlabel('Pulse Number');
ylabel('Range (m)');
title('Range-Time Plot (Channel 1)');
colorbar;
grid on;

% Range-Channel plot for first pulse
subplot(2,2,2);
% FIXED: Use squeeze to convert 3D to 2D for plotting
range_channel_data = squeeze(dataCube(:,1,:));
imagesc(1:numChannels, range_axis, 20*log10(abs(range_channel_data)));
xlabel('Channel Number');
ylabel('Range (m)');
title('Range-Channel Plot (Pulse 1)');
colorbar;
grid on;

% Pulse-Channel plot for target range bin
targetBin = round(2*target_range/c * fs); % Range bin where target is located
targetBin = min(max(targetBin, 1), samplesPerPulse); % Ensure within bounds
subplot(2,2,3);
% FIXED: Use squeeze and transpose for proper 2D plotting
pulse_channel_data = squeeze(dataCube(targetBin,:,:));
imagesc(1:numChannels, 1:numPulses, 20*log10(abs(pulse_channel_data')));
xlabel('Channel Number');
ylabel('Pulse Number');
title('Pulse-Channel Plot (Target Range Bin)');
colorbar;
grid on;

% 3D visualization of the data cube
subplot(2,2,4);
slice = squeeze(20*log10(abs(dataCube(:,:,1))));
surf(1:numPulses, range_axis, slice);
xlabel('Pulse Number');
ylabel('Range (m)');
zlabel('Amplitude (dB)');
title('3D View of Data Cube (Channel 1)');
grid on;
shading interp;

% === Doppler Processing (Example) ===
% Extract data from target range bin across all pulses and channels
targetData = squeeze(dataCube(targetBin, :, :));

% Apply Doppler processing (FFT across pulses)
dopplerFFT = fftshift(fft(targetData, [], 1), 1);
dopplerAxis = (-numPulses/2:numPulses/2-1) * PRF / numPulses;

figure;
plot(dopplerAxis, 20*log10(abs(dopplerFFT(:,1))));
xlabel('Doppler Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Doppler Spectrum (Channel 1)');
grid on;