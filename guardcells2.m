% Parameters
N = 50;                   % Total range bins
numGuard = 2;             % Guard cells each side
numTraining = 8;          % Training cells each side

% Create synthetic power profile (random noise + few targets)
signal = 10*rand(1,N);
targetPositions = [15, 30, 45]; % target locations
signal(targetPositions) = [40, 55, 50]; % target strengths

figure;
for CUT = (numTraining+numGuard+2):(N-numTraining-numGuard-1) % 12:39
    
    % Identify guard and training cell indices
    guardIdx = (CUT-numGuard):(CUT+numGuard); 
    trainIdxLeft = (CUT-numGuard-numTraining):(CUT-numGuard-1);
    trainIdxRight = (CUT+numGuard+1):(CUT+numGuard+numTraining);
    
    % Plot all bins
    stem(1:N, signal, 'k', 'LineWidth', 1.5); hold on;
    
    % Highlight CUT
    stem(CUT, signal(CUT), 'r', 'LineWidth', 2);
    
    % Highlight Guard Cells
    stem(guardIdx, signal(guardIdx), 'b', 'LineWidth', 2);
    
    % Highlight Training Cells
    stem(trainIdxLeft, signal(trainIdxLeft), 'g', 'LineWidth', 2);
    stem(trainIdxRight, signal(trainIdxRight), 'g', 'LineWidth', 2);
    
    % Labels
    xlabel('Range Bin');
    ylabel('Amplitude');
    title(sprintf('Guard Cells & Training Cells | CUT = %d', CUT));
    legend('Noise/Clutter', 'CUT', 'Guard Cells', 'Training Cells');
    grid on;
    hold off;
    
    pause(1.0); % animation speed
end
