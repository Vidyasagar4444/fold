% Parameters
N = 50;                  % Total range bins
CUT = 25;                % Cell Under Test index
numGuard = 2;            % Guard cells on each side
numTraining = 8;         % Training cells on each side

% Create synthetic power profile
signal = 10*rand(1,N);   % noise
signal(CUT) = 50;        % target at CUT

% Identify guard and training cell indices
guardIdx = (CUT-numGuard):(CUT+numGuard);      % 23:27
trainIdxLeft = (CUT-numGuard-numTraining):(CUT-numGuard-1); % 15:22
trainIdxRight = (CUT+numGuard+1):(CUT+numGuard+numTraining); % 28:35

% Plot
figure;
stem(1:N, signal, 'k', 'LineWidth', 1.5); hold on;
stem(CUT, signal(CUT), 'r', 'LineWidth', 2); % CUT in red
stem(guardIdx, signal(guardIdx), 'b', 'LineWidth', 2); % Guard cells in blue
stem(trainIdxLeft, signal(trainIdxLeft), 'g', 'LineWidth', 2); % Training cells in green
stem(trainIdxRight, signal(trainIdxRight), 'g', 'LineWidth', 2); % Training cells in green

% Labels
xlabel('Range Bin');
ylabel('Amplitude');
title('Guard Cells and Training Cells around CUT');
legend('Noise/Clutter', 'CUT', 'Guard Cells', 'Training Cells');
grid on;
