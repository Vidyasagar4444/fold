%% STAP Angle-Doppler Visualization: Before & After Processing
clear; clc; close all;

% --- Grid parameters ---
az = linspace(-1,1,200);           % sin(azimuth)
doppler = linspace(-0.5,0.5,200);  % normalized Doppler
[AZ, DOP] = meshgrid(az, doppler);

%% --- Generate interference environment ---
% Clutter ridge (centered at Doppler=0)
clutter = 30 * exp(-((AZ-0).^2 / 0.1) - ((DOP).^2 / 0.02));

% Jammer (wall in azimuth)
jammer = zeros(size(AZ));
jammer(abs(AZ - 0.8) < 0.02) = 40;  % strong, wide in Doppler

% Target (small peak away from clutter/jammer)
target = 15 * exp(-((AZ - (-0.2)).^2 / 0.005) - ((DOP - 0.25).^2 / 0.005));

% Combined SNR surface (Before STAP)
SNR_before = clutter + jammer + target;

%% --- Simulated STAP filtering ---
% Create suppression masks (notch filter in clutter & jammer zones)
clutter_mask = 1 - exp(-((AZ-0).^2 / 0.08) - ((DOP).^2 / 0.015));   % notch near clutter ridge
jammer_mask  = 1 - exp(-((AZ-0.8).^2 / 0.001));                     % notch at jammer azimuth

% Total STAP suppression filter
stap_filter = clutter_mask .* jammer_mask;

% Apply filter to suppress clutter & jammer
SNR_after = SNR_before .* stap_filter;

%% --- Plot BEFORE STAP ---
figure
subplot(1,2,1)
surf(AZ, DOP, SNR_before, 'EdgeColor','none');
xlabel('SIN(Azimuth)');
ylabel('Normalized Doppler');
zlabel('SNR (dB)');
title('Before STAP Processing');
colormap jet
colorbar
view([-45 30])
grid on
hold on
text(-0.05, 0, max(clutter(:))+5, 'CLUTTER', 'FontWeight','bold', 'Color','k');
text(0.82, 0, max(jammer(:))+5, 'JAMMING', 'FontWeight','bold', 'Color','k');
text(-0.2, 0.25, max(target(:))+5, 'TARGET', 'FontWeight','bold', 'Color','k');

%% --- Plot AFTER STAP ---
subplot(1,2,2)
surf(AZ, DOP, SNR_after, 'EdgeColor','none');
xlabel('SIN(Azimuth)');
ylabel('Normalized Doppler');
zlabel('SNR (dB)');
title('After STAP Processing');
colormap jet
colorbar
view([-45 30])
grid on
hold on
text(-0.2, 0.25, max(target(:))+5, 'TARGET', 'FontWeight','bold', 'Color','k');
