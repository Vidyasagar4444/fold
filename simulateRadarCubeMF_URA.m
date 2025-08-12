function y = simulateRadarCubeMF_URA()
%#codegen
% Returns matched-filtered radar data cube [samples x elements x pulses] for URA

% === Parameters ===
nRows = 8;              % URA rows (elevation)
nCols = 8;              % URA cols (azimuth)
nElements = nRows * nCols;
nPulses = 64;
fs = 1e6;
prf = 1000;
pw = 1e-4;
fc = 300e6;
lambda = physconst('LightSpeed') / fc;

% === Waveform ===
waveform = phased.LinearFMWaveform( ...
    'SampleRate', fs, ...
    'PRF', prf, ...
    'PulseWidth', pw);
nSamples = fs / prf;
refPulse = getMatchedFilter(waveform);

% === URA Antenna ===
dy = lambda / 2;
dx = lambda / 2;
antenna = phased.URA('Size', [nRows, nCols], ...
                     'ElementSpacing', [dy, dx]);

% === Radar and Target Components ===
TX = phased.Transmitter('Gain', 20);
TgtModel = phased.RadarTarget;

tgtPos = [10e3 * sqrt(3); 10e3; 0];
tgtVel = [75 * sqrt(3); 75; 0];
PlatformModel = phased.Platform('InitialPosition', tgtPos, ...
                                'Velocity', tgtVel);
ChannelModel = phased.FreeSpace('TwoWayPropagation', true, ...
                                'SampleRate', fs, ...
                                'OperatingFrequency', fc);

txArray = phased.Radiator('Sensor', antenna, ...
                          'OperatingFrequency', fc);
rxArray = phased.Collector('Sensor', antenna, ...
                           'OperatingFrequency', fc);
rxPreamp = phased.ReceiverPreamp('Gain', 10, ...
                                 'NoiseFigure', 5, ...
                                 'SampleRate', fs);

radarPos = [0;0;0];
radarVel = [0;0;0];

% === Preallocate ===
dataCube = complex(zeros(nSamples, nElements, nPulses));

% === Loop over pulses ===
for ii = 1:nPulses
    wf = waveform();
    [tgtPos, tgtVel] = PlatformModel(1/prf);
    [~, tgtAng] = rangeangle(tgtPos, radarPos);

    s0 = TX(wf);
    s1 = txArray(s0, tgtAng);
    s2 = ChannelModel(s1, radarPos, tgtPos, radarVel, tgtVel);
    s3 = TgtModel(s2);
    s4 = rxArray(s3, tgtAng);
    s5 = rxPreamp(s4);

    % === Matched Filter per element ===
    for el = 1:nElements
        dataCube(:, el, ii) = conv(s5(:,el), refPulse, 'same');
    end
end

y = dataCube;  % Return 3D matched-filtered cube
end
