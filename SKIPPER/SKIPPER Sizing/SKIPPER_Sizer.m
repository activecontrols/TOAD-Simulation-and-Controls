%% Script to size thrust requierments for the upcoming PSP-AC Traditionally manufactured TCA
% Assumes Total OF of 1.
clear;
clc;
close all;

% Initialize arrays & parameters
MaxIter = 100;
GrossMass = zeros(MaxIter, 1);
PropMass = zeros(MaxIter, 1);
ThrustReq = zeros(MaxIter, 1);
TankHeights = zeros(MaxIter, 2);
TankMasses = zeros(MaxIter, 2);
GrossMass(1) = 205; % Wet Mass Guess [kg]
PropMass(1) = 41; % Total Propellant Mass Guess [kg]
ThrustReq(1) = 2446; % Engine Thrust Guess [N]

% Structures Mass (incl. all Avionics, Structures mass, COPVs, Fittings,
% and all things constant with thrust) Replace with better value later [kg]
StructMass = 157.45 - 5; % EXTREMELY SENSITIVE AND SUS
Rho_FU = 786; % IPA [kg/m^3]
Rho_OX = 1141.2; % Ox [kg/m^3]
chillin = 2; % Mass required for TCA chill-in [kg]
OF = 1; % Total OF (including 20% film)
TWR_Target = 1.35; 
Alpha = 0.3;
g0 = 9.8066;
FoS = 1.5;

%% Main Loop & Initial guessed from current TOAD Setup
for i = 1:MaxIter
    % Iteration step
    [PropMass(i), FlightTimes(i)] = SKIPPER_3DoF(ThrustReq(i), GrossMass(i));
    if PropMass(i) <= 0
        error(['Simulation failed at Iteration %d! \n' ...
               'KOMODO_3DoF could not find a valid landing trajectory for:\n' ...
               'Thrust = %.2f N (%.2f lbf) | Gross Mass = %.2f kg.\n' ...
               'Check your LQR weights or terminal conditions.'], ...
               i, ThrustReq(i), ThrustReq(i)/4.448, GrossMass(i));
    end
    PropMass(i) = PropMass(i) * FoS;
    [TankMasses(i, :), TankHeights(i, :)] = TankSizer(PropMass(i), Rho_FU, Rho_OX, OF, chillin);
    
    % New Gross Mass
    GrossMass(i + 1) = PropMass(i) + sum(TankMasses(i, :), 2) + StructMass;

    % Calculate thrust residual
    Residual = TWR_Target * GrossMass(i + 1) * g0 - ThrustReq(i);
    ThrustReq(i + 1) = ThrustReq(i) + Alpha * Residual;
    fprintf(['Iteration #%i Complete!\nResidual: %.2f lbf   ||   Thrust Requirement: ' ...
        '%.2f lbf \n'], i, Residual / 4.448, ThrustReq(i + 1) / 4.448);

    if abs(Residual) < 1
        LastIter = i;
        break;
    end
end

% Trim Arrays
% Trim Arrays to the last iteration
GrossMass = GrossMass(1:LastIter + 1);
PropMass = PropMass(1:LastIter);
ThrustReq = ThrustReq(1:LastIter + 1);
TankHeights = TankHeights(1:LastIter, :);
TankMasses = TankMasses(1:LastIter, :);

%% Plotting Section
% Create a docked figure
figure('Name', 'Sizing Loop Convergence', 'WindowStyle', 'docked');

% 1. Thrust Convergence (converted to lbf for readability)
subplot(3, 2, 1);
plot(1:LastIter+1, ThrustReq / 4.448, '-o', 'LineWidth', 1.5, 'Color', '#0072BD');
grid on;
xlabel('Iteration');
ylabel('Thrust (lbf)');
title('Thrust Requirement Convergence');

% 2. Mass Breakdown Convergence
subplot(3, 2, 2);
plot(1:LastIter+1, GrossMass * 2.205, '-k', 'LineWidth', 1.5, 'DisplayName', 'Gross Mass');
hold on;
plot(1:LastIter, PropMass * 2.205, '-o', 'LineWidth', 1.5, 'DisplayName', 'Propellant Mass');
plot(1:LastIter, sum(TankMasses, 2) * 2.205, '-s', 'LineWidth', 1.5, 'DisplayName', 'Total Tank Mass');
yline(StructMass * 2.205, '--', 'DisplayName', 'Fixed Struct Mass', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
grid on;
xlabel('Iteration');
ylabel('Mass (lbm)');
title('Vehicle Mass Sizing');
legend('Location', 'best');

% 3. Tank Geometry
subplot(3, 2, 3);
plot(1:LastIter, TankHeights(:, 1) * 12 * 3.281, '-o', 'LineWidth', 1.5, 'DisplayName', 'OX Cyl Height');
hold on;
plot(1:LastIter, TankHeights(:, 2) * 12 * 3.281, '-s', 'LineWidth', 1.5, 'DisplayName', 'FU Cyl Height');
grid on;
xlabel('Iteration');
ylabel('Cylindrical Height (in)');
title('Tank Dimension Sizing');
legend('Location', 'best');

% 4. Convergence Behavior (Step Size)
subplot(3, 2, 4);
% Calculate the absolute change in thrust between iterations to visualize decay
DeltaThrust = abs(diff(ThrustReq)) / 4.448; 
semilogy(1:LastIter, DeltaThrust, '-x', 'LineWidth', 1.5, 'Color', '#D95319');
grid on;
xlabel('Iteration');
ylabel('Absolute Step Size (lbf)');
title('Solver Convergence (Log Scale)');

% 5. Minimum Throttle Requirement
subplot(3, 2, [5, 6]); % Spans the bottom two columns
DryMass = sum(TankMasses, 2) + StructMass; 
MinThrottle = (0.9 .* DryMass .* g0) ./ ThrustReq(1:LastIter) .* 100; 

plot(1:LastIter, MinThrottle, '-d', 'LineWidth', 1.5, 'Color', '#7E2F8E', 'MarkerFaceColor', '#7E2F8E');
grid on;
xlabel('Iteration');
ylabel('Min Throttle (%)');
title('Minimum Throttle Requirement (Hover at 90% Dry Mass)');
% Add a text box to display the final converged value
finalThrottleStr = sprintf('Converged Min Throttle: %.1f%%', MinThrottle(end));
text(LastIter, MinThrottle(end), finalThrottleStr, 'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'FontSize', 10);

%% Tank Helper function
function [TankMass, TankHeight] = TankSizer(PropMass, Rho_FU, Rho_OX, OF, chillin)
    % Tank Radius & Wall Thickness
    Radius = 0.1524; % [m]
    Thickness = 0.00635; % [m]
    Ullage = 0.10; % Percent

    % Propellant Masses
    Mass_FU = PropMass / (1 + OF);
    Mass_OX = PropMass * OF / (1 + OF) + chillin; % NOT SURE WHY HAS 2 (CHILLIN?)

    % Propellant Volume Requirements
    Volume_FU = (1 + Ullage) * Mass_FU / Rho_FU;
    Volume_OX = (1 + Ullage) * Mass_OX / Rho_OX;

    % Predetermined Tank Masses
    Bulkhead_mass = 11.0509 / 2.205; % Mass of upper and lower bulkheads [kg]
    Vortex_mass = 0.0488 / 2.205; % Mass of tank vortex baffles
    Slosh_spacing = 2.8 / 12 / 3.281; % Spacing between each slosh baffle [m]
    slosh_mass = 1.478357 / 3 / 2.205; % Mass of individual slosh baffle [kg]

    % Tank Dimension (just cyllindrical for the TankHeight variables, we
    % then add endcap heights)
    Height_Endcap = Radius / sqrt(2);
    Volume_Endcaps = 4 / 3 * pi * Radius^2 * Height_Endcap;
    TankHeight_FU = (Volume_FU - Volume_Endcaps) / (pi * Radius^2);
    TankHeight_OX = (Volume_OX - Volume_Endcaps) / (pi * Radius^2);
    TankHeight(1) = TankHeight_OX;
    TankHeight(2) = TankHeight_FU;  

    % Tank Mass
    Al_Rho = 2700; % [kg/m^3]
    Ox_slosh_mass = slosh_mass * floor(TankHeight_OX / Slosh_spacing); % Mass for ox tank slosh baffles [kg]
    Fu_slosh_mass = slosh_mass * floor(TankHeight_FU / Slosh_spacing); % Mass for fuel tank slosh baffles [kg]
    Ox_cyl_mass = Al_Rho * TankHeight_OX * (pi * (Radius + Thickness)^2 - pi * Radius^2); % Mass of ox tank cylinder [kg]
    Fu_cyl_mass = Al_Rho * TankHeight_FU * (pi * (Radius + Thickness)^2 - pi * Radius^2); % Mass of ox tank cylinder [kg]

    % Masses
    TankMass(1) = Bulkhead_mass + Vortex_mass + Ox_slosh_mass + Ox_cyl_mass; % Oxygen Tank Mass [kg]
    TankMass(2) = Bulkhead_mass + Vortex_mass + Fu_slosh_mass + Fu_cyl_mass; % IPA Tank Mass [kg]
    % disp(floor(TankHeight_OX / Slosh_spacing))
    % disp(floor(TankHeight_FU / Slosh_spacing))
end
