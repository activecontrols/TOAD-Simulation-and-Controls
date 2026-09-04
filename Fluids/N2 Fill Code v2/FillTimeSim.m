%% COPV Fill Time Simulation for TOAD Vehicle
function Time2Settle = FillTimeSim(WindVel)
% clear; clc;
psi2Pa = 6894.75729;
Pa2psi = 1 / psi2Pa;
T_amb = 297.15;
% WindVel = 1.5; %m/s.

% Initialize Supply Tank
Tank.V      = 6 * 0.05;                    % m3 
Tank.P(1)   = 6000 * psi2Pa;               % Pa
Tank.T(1)   = 297.15;                      % K 
Tank.rho(1) = double(py.CoolProp.CoolProp.PropsSI('D', 'P', Tank.P(1), 'T', Tank.T(1), 'Nitrogen')); 
Tank.m(1)   = Tank.rho(1) * Tank.V;        % kg
Tank.h(1)   = double(py.CoolProp.CoolProp.PropsSI('Hmass', 'P', Tank.P(1), 'T', Tank.T(1), 'Nitrogen')); 
Tank.u(1)   = double(py.CoolProp.CoolProp.PropsSI('Umass', 'P', Tank.P(1), 'T', Tank.T(1), 'Nitrogen'));
Tank.U(1)   = Tank.u(1) * Tank.m(1);       % J

% Initialize Receiver COPV
COPV.V      = 0.036;                       % m3
COPV.P(1)   = 14.7 * psi2Pa;               % Pa 
COPV.T(1)   = 297.15;                      % K 
COPV.rho(1) = double(py.CoolProp.CoolProp.PropsSI('D', 'P', COPV.P(1), 'T', COPV.T(1), 'Nitrogen')); 
COPV.m(1)   = COPV.rho(1) * COPV.V;        % kg 
COPV.h(1)   = double(py.CoolProp.CoolProp.PropsSI('Hmass', 'P', COPV.P(1), 'T', COPV.T(1), 'Nitrogen')); 
COPV.u(1)   = double(py.CoolProp.CoolProp.PropsSI('Umass', 'P', COPV.P(1), 'T', COPV.T(1), 'Nitrogen'));
COPV.U(1)   = COPV.u(1) * COPV.m(1);       % J

COPV.TMax   = 320.00;                      % K
COPV_Temps  = [297.15, 297.15, 297.15];    % Aluminum, Carbon Fiber, and Fiber Glass temps

% Reg Properties
Reg.SetPress = 4500 * psi2Pa;              % Pa
Reg.SPE_Coef = 0.07;                       % psi drop per psi inlet drop 
Reg.CdA_MAX  = 2.4828e-6;                  % m2 

%% Simulation Setup
% Adaptive timestep controls
dt_min        = 0.01;     
dt_max_fill   = 5.0;      
dt_max_cool   = 5.0;      
dt_init       = 0.05;     
target_dmFrac = 0.10;     % max COPV mass change per ste
target_bandFrac = 0.10;   % fraction of regulator throttle band traversed per step

dT = dt_init;
isFilled = false;
i = 1;
Time(i) = 0;
dt_hist(1) = dT;
mdot_hist(1) = 0;
COPV_T_hist(1) = COPV_Temps(1);

while ~isFilled
    % Pause fill if temps exceed the level
    if COPV_Temps(1) < COPV.TMax
        % Calculate regulator massflow
        fillMode = 1; 
        gamma = double(py.CoolProp.CoolProp.PropsSI('isentropic_expansion_coefficient', ...
            'P', Tank.P(i), 'T', Tank.T(i), 'Nitrogen'));
        Tank_Rho = Tank.m(i) / Tank.V;
        mdot = RegulatorMdot(Reg, Tank.P(i), Tank_Rho, Tank.P(1), COPV.P(i), gamma);
        dt_limit = dt_max_fill;
    else
        fillMode = 0; 
        mdot = 0; 
        dt_limit = dt_max_cool;
    end

    %% Adaptive timestep
    if mdot > 0
        dt_mass_copv = target_dmFrac * COPV.m(i) / mdot;
        dt_mass_tank = target_dmFrac * Tank.m(i) / mdot;

        % Approximate dP/dt using P ~ const * m as a local scaling.
        % This is only a timestep limiter; the actual state is still obtained
        % from CoolProp after each accepted step.
        RegPress = Reg.SetPress + Reg.SPE_Coef * (Tank.P(1) - Tank.P(i));
        P_band = 0.03 * RegPress;
        dPdt_est = max(COPV.P(i) * mdot / COPV.m(i), eps);
        dt_reg_band = target_bandFrac * P_band / dPdt_est;

        dT = min([dt_limit, dt_mass_copv, dt_mass_tank, dt_reg_band]);
    else
        dT = dt_limit;
    end

    % Clamp to hard bounds.
    dT = max(dt_min, min(dT, dt_limit));
    
    % Heat Transfer out of the COPV, using a thermal network. Check all of
    % this as I am just copying from the old code and kinda flying blind
    [Q_gas, COPV_Temps] = COPVHeatTransfer(fillMode, dT, COPV.P(i), COPV.T(i), COPV_Temps, WindVel); 
    
    % State Integrations
    Tank.m(i+1) = Tank.m(i) - mdot * dT; 
    COPV.m(i+1) = COPV.m(i) + mdot * dT; 
    
    Tank.U(i+1) = Tank.U(i) - (mdot * Tank.h(i)) * dT; 
    COPV.U(i+1) = COPV.U(i) + (mdot * Tank.h(i)) * dT - Q_gas; 
    
    % Calculate intensive states so we can extract properties from
    % CoolProp
    Tank_rho = Tank.m(i+1) / Tank.V; 
    Tank_u   = Tank.U(i+1) / Tank.m(i+1); 
    
    COPV_rho = COPV.m(i+1) / COPV.V; 
    COPV_u   = COPV.U(i+1) / COPV.m(i+1); 
    
    % CoolProp Thermodynamic Lookups
    Tank.P(i+1) = double(py.CoolProp.CoolProp.PropsSI('P', 'Umass', Tank_u, 'Dmass', Tank_rho, 'Nitrogen')); 
    Tank.T(i+1) = double(py.CoolProp.CoolProp.PropsSI('T', 'Umass', Tank_u, 'Dmass', Tank_rho, 'Nitrogen')); 
    Tank.h(i+1) = double(py.CoolProp.CoolProp.PropsSI('Hmass', 'Umass', Tank_u, 'Dmass', Tank_rho, 'Nitrogen')); 
    
    COPV.P(i+1) = double(py.CoolProp.CoolProp.PropsSI('P', 'Umass', COPV_u, 'Dmass', COPV_rho, 'Nitrogen')); 
    COPV.T(i+1) = double(py.CoolProp.CoolProp.PropsSI('T', 'Umass', COPV_u, 'Dmass', COPV_rho, 'Nitrogen')); 
    COPV.h(i+1) = double(py.CoolProp.CoolProp.PropsSI('Hmass', 'Umass', COPV_u, 'Dmass', COPV_rho, 'Nitrogen')); 
    
    % Stop fill if the COPV press is at 98% of the set press and the temp
    % is at ambient
    isPressSettled = COPV.P(i+1) >= (Reg.SetPress * 0.995);
    isTempSettled  = COPV.T(i+1) <= (T_amb + 5);
    isFlowStopped  = mdot < 1e-2;
    
    if isPressSettled && isTempSettled && isFlowStopped
        isFilled = true; 
    end
    i = i + 1; 
    Time(i) = Time(i - 1) + dT;
    dt_hist(i) = dT;
    mdot_hist(i) = mdot;
    COPV_T_hist(i) = COPV_Temps(1);

    if mod(i, 100) == 0
        % fprintf('Time: %.1f min:    Fill Pct: %.2f%%    dt: %.3f s\n', ...
        %    Time(i) / 60, COPV.P(i) / Reg.SetPress * 100, dT);
    end
end

fprintf('Fill finished at t = %.1f s. Final COPV Pressure: %.1f psi\n', Time(end), COPV.P(end)*Pa2psi);

%% Supporting Regulator Function
function mdot = RegulatorMdot(Reg, P_tank, Tank_Rho, P_tank0, P_copv, gamma)
    % Press rise from inlet drop.
    RegPress = Reg.SetPress + Reg.SPE_Coef * (P_tank0 - P_tank);

    % Check this from actual flow curves from the reg
    P_band   = 0.03 * RegPress;
    
    % Throttling logic
    if P_copv >= RegPress
        Reg_CdA = 0;
    elseif P_copv <= (RegPress - P_band)
        Reg_CdA = Reg.CdA_MAX;
    else
        Reg_CdA = Reg.CdA_MAX * ((RegPress - P_copv) / P_band);
    end
    
    if Reg_CdA <= 0 || P_tank <= P_copv 
        mdot = 0; 
        return;
    end
    
    % Check for choked flow and use the equation that is needed
    r_crit = (2 / (gamma + 1))^(gamma / (gamma - 1)); 
    pr = P_copv / P_tank; 
    
    if pr <= r_crit 
        mdot = Reg_CdA * sqrt(gamma * P_tank * Tank_Rho) * ...
               (2 / (gamma + 1))^((gamma + 1) / (2 * (gamma - 1))); 
    else            
        mdot = Reg_CdA * sqrt(P_tank * Tank_Rho) * ...
               sqrt((2 * gamma / (gamma - 1)) * (pr^(2 / gamma) - pr^((gamma + 1) / gamma)));
    end
end
    
Time2Settle = Time(end);
end

% %% Plots & Post-Processing
% t_min = Time / 60;
% 
% % mdot is associated with each completed integration interval.
% t_flow = Time(1:end-1) / 60;
% mdot_plot = mdot_hist(2:end);
% 
% figure('Name', 'COPV Fill & Thermal Transient Analysis', 'Color', 'w');
% 
% % Pressure Transients
% subplot(2, 2, 1);
% plot(t_min, Tank.P * Pa2psi, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.1]); hold on;
% plot(t_min, COPV.P * Pa2psi, 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
% yline(Reg.SetPress * Pa2psi, '--r', 'Regulator Setpoint', 'LineWidth', 1.2);
% xlabel('Time (min)');
% ylabel('Pressure (psi)');
% title('Pressure History');
% grid on; grid minor;
% legend('Supply Tank', 'COPV Gas', 'Location', 'best');
% 
% % Temperature Transients
% subplot(2, 2, 2);
% plot(t_min, COPV_T_hist, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.1]); hold on;
% plot(t_min, Tank.T, 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
% yline(333.15, '--r', 'Absolute COPV Limit (60°C)', 'LineWidth', 1.2)
% yline(COPV.TMax, '--g', 'Target Max Temp', 'LineWidth', 1.2);
% yline(COPV.T(1), ':k', 'Ambient (24°C)', 'LineWidth', 1.2);
% xlabel('Time (min)');
% ylabel('Temperature (K)');
% title('Fluid Temperature Transients');
% grid on; grid minor;
% legend('COPV Gas', 'Supply Tank Gas', 'Location', 'best');
% 
% % Mass Flow Rate
% subplot(2, 2, 3);
% plot(t_flow, mdot_plot * 1000, 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
% xlabel('Time (min)');
% ylabel('Mass Flow Rate (g/s)');
% title('Regulator Delivery Rate (\dot{m})');
% grid on; grid minor;
% 
% % System Mass Distribution
% subplot(2, 2, 4);
% plot(t_min, COPV.m, 'LineWidth', 1.5, 'Color', [0 0.45 0.74]); hold on;
% plot(t_min, Tank.m, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.1]);
% xlabel('Time (min)');
% ylabel('Fluid Mass (kg)');
% title('Mass Tracking');
% grid on; grid minor;
% legend('COPV Gas Mass', 'Supply Tank Mass', 'Location', 'best');