% Clear workspace
clear;
close;

%% Filter structure design
    % 6th order filter composed of 2 sequential ANF's (Adaptive Notch Filters)
    % and a 1st order LPF.
    % The 1st Notch filter is set to a constant frequency os ~120Hz. Other
    % ANF follows a linear trace based on thrust (see RollTestAnalysis.m for
    % that).
    
% Create Notch Filters
fs = 1000;
res = 200;
thrustArray = linspace(1, 100, res);
NotchC_Array = zeros(1, res) * tf('s');
NotchD_Array = NotchC_Array;

for thrust = 1:1:res
    [NotchC_Array(1, thrust), NotchD_Array(1, thrust)] = FilterTF_Gen(thrustArray(thrust));
end

% Frequency range in Hz
f = linspace(0, 500, 2000);    % Linear spacing 0–500 Hz
w = 2 * pi * f;                % Convert to rad/s for 'bode' or 'freqresp'

% Get frequency response at each thrust level
magC = zeros(res, size(f, 2));
phaseC = magC;
for thrust = 1:1:res
    [mag, phase] = bode(NotchC_Array(1, thrust), w);
    mag = squeeze(mag);
    phase = squeeze(phase);
    phase = wrapTo180(phase);
    magC(thrust, :) = mag; 
    phaseC(thrust, :) = phase;
end

% Plots
% --- 1. Fix the Phase Wrap Artifact ---

% 'unwrap' works in radians, so we convert, unwrap, and convert back
% We unwrap along dimension 2 (the rows), which is your frequency axis
magC_dB = 20 * log10(magC);
phaseC_rad = deg2rad(phaseC); 
phaseC_unwrapped_rad = unwrap(phaseC_rad, [], 2);
phaseC_unwrapped_deg = rad2deg(phaseC_unwrapped_rad);

% --- 3. Create the 3D Surface Plots (Convention: X=Thrust, Y=Freq) ---
figure;

% --- Top plot for Magnitude ---
subplot(2,1,1);
colormap turbo
% Swap axes (f, thrustArray -> thrustArray, f) and transpose data (magC_dB -> magC_dB')
surf(thrustArray, f, magC_dB', 'EdgeColor', 'none');
shading interp; 
xlabel('Thrust [%]');           % <-- Swapped
ylabel('Frequency [Hz]');      % <-- Swapped
zlabel('Magnitude [dB]');
title('ANF Magnitude (3D Surface)');
colorbar;
clim([-40, 0]);      % Anything under 5dB is considered cutoff.
view(-90,90);
xline(40, 'g--', 'LineWidth', 2);
yline(10, 'r--', 'LineWidth', 2);
% view(-70, 40);

% --- Bottom plot for Phase (using unwrapped data) ---
subplot(2,1,2);
colormap turbo
% Swap axes and transpose data
surf(thrustArray, f, phaseC_unwrapped_deg', 'EdgeColor', 'none'); 
shading interp;
xlabel('Thrust [%]');           % <-- Swapped
ylabel('Frequency [Hz]');      % <-- Swapped
zlabel('Phase [degrees]');
title('ANF Phase (3D Surface - Unwrapped)');
colorbar;
clim([-180 0]);
view(-90,90);
% view(-70, 40);
xline(40, 'g--', 'LineWidth', 2);
yline(10, 'r--', 'LineWidth', 2);

% Link the camera angles so they rotate together
linkaxes([subplot(2,1,1), subplot(2,1,2)], 'xy');
sgtitle('Adaptive Notch Filter Sequence [ANF-S] Design for ASTRAv2');