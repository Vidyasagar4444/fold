

Rubber = simulateRadarCubeMF_URA();


figure;
elementIdx = 1;
rdMap = squeeze(abs(Rubber(:, elementIdx, :))); % size: [range × pulses]
imagesc(rdMap);
xlabel('Pulse Index'); ylabel('Range Bin');
title('Range-Doppler Map for 1 Element');
colorbar;

volshow(abs(Rubber));