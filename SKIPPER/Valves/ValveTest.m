%% superTADPOLE Sine Wave Test Script (Mode Comparison)
clear; clc; close all;

% Simulation Setup
tspan = [0 4];          % Simulate for 4 seconds 
x0 = [0; 0; 0; 0; 0];      % Initial states: [MotorPos, MotorVel, ValvePos, ValveVel]

% Define the Target Trajectory (0.5Hz, 20 deg amplitude)
Freq = 2;                % Hz
Amplitude = 10 * (pi / 180); % Convert degrees to radians
AngleTargetFCN = @(t) Amplitude * sin(2 * pi * Freq * t);

% ODE Solver Options
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5); 

% Run INTERNAL Mode
disp('Simulating Internal Mode...');
ode_wrapper_int = @(t, x) ValveDynamics(x, AngleTargetFCN(t), "Internal");
[t_int, x_int] = ode15s(ode_wrapper_int, tspan, x0, options);

% Run CLOSED Mode
disp('Simulating Closed Cascaded Mode...');
ode_wrapper_closed = @(t, x) ValveDynamics(x, AngleTargetFCN(t), "Closed");
[t_closed, x_closed] = ode15s(ode_wrapper_closed, tspan, x0, options);

disp('Simulations Complete.');

% Plot the Results
Target_deg = AngleTargetFCN(t_int) * (180 / pi);
ValvePos_Int_deg = x_int(:, 3) * (180 / pi);
ValvePos_Closed_deg = x_closed(:, 3) * (180 / pi);

figure('WindowStyle', 'docked', 'Name', 'Valve Tracking');
plot(t_int, Target_deg, '--k', 'LineWidth', 1.5); hold on;
plot(t_int, ValvePos_Int_deg, 'r', 'LineWidth', 1.5);
plot(t_closed, ValvePos_Closed_deg, 'b', 'LineWidth', 1.5);
ylim([-30 30]);

title('Valve Tracking: Internal vs. Closed Dual-Loop Controller');
xlabel('Time (s)');
ylabel('Angle (Degrees)');
legend('Target Trajectory', 'Internal Mode (Lags)', 'Closed Mode (Actively Tracks)', 'Location', 'best');
grid on;

MotorPos_Int_deg = x_int(:, 1) * (180 / pi) / 10;
MotorPos_Closed_deg = x_closed(:, 1) * (180 / pi) / 10;

figure('WindowStyle', 'docked', 'Name', 'Motor Tracking');
plot(t_int, Target_deg, '--k', 'LineWidth', 1.5); hold on;
plot(t_int, MotorPos_Int_deg, 'r', 'LineWidth', 1.5);
plot(t_closed, MotorPos_Closed_deg, 'b', 'LineWidth', 1.5);
ylim([-30 30]);

title('Motor Tracking (Geared Down): Internal vs. Closed Mode');
xlabel('Time (s)');
ylabel('Geared Motor Angle (Degrees)');
legend('Target Trajectory', 'Internal Mode Motor (Smooth)', 'Closed Mode Motor (Active Jumps)', 'Location', 'best');
grid on;