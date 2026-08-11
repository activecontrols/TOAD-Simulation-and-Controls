%% This trajectory manager script shall take in the current time, and the
% name of a trajectory file and the current MET. It shall then find the
% corresponding gain matrix file and interpolate the gain matrices and
% state targets correspondingly. If MET > MaxT, hold final state, matrix
% and feedforward input. Outputs interpolated target state, feedforward
% input, and feedback gain matrix. 

function [X, U, P, LA, LT] = TrajectoryManager(t, constantsTOAD)
    %% Read the passed trajectory
    Time = constantsTOAD.Traj.Time;
    States = constantsTOAD.Traj.States;
    Inputs = constantsTOAD.Traj.Inputs;
    Feedback = constantsTOAD.Traj.FBCost;
    LA_Gain = constantsTOAD.Traj.LAGain;
    LT_Gain = constantsTOAD.Traj.LTGain;

    % Pre-allocate outputs immediately to lock fixed sizes for Simulink
    X = zeros(15, 1);
    U = zeros(4, 1);
    P = zeros(12, 12);
    LA = zeros(6, 3);
    LT = zeros(3, 1);
    
    %% Handle Boundary Conditions
    if t <= Time(1) || t == 0
        X(:) = States(1,:);
        U(:) = Inputs(1,:);
        P_temp = zeros(12, 12);
        P_temp(:) = Feedback(1, :, :); 
        P(:) = P_temp';

        LA_temp = zeros(6, 3);
        LA_temp(:) = LA_Gain(1, :, :);
        LA(:) = LA_temp;

        LT_temp = zeros(2, 1);
        LT_temp(:) = LT_Gain(1, :, :);
        LT(:) = LT_temp;

        return;
    elseif t >= Time(end)
        X(:) = States(end,:);
        m_lox_end = States(end, 14);
        m_ipa_end = States(end, 15);
        total_mass = constantsTOAD.m_dry + m_lox_end + m_ipa_end;
        
        U(:) = [0; 0; total_mass * constantsTOAD.g; 0];
        
        P_temp = zeros(12, 12);
        P_temp(:) = Feedback(end, :, :);
        P(:) = P_temp';

        LA_temp = zeros(6, 3);
        LA_temp(:) = LA_Gain(end, :, :);
        LA(:) = LA_temp;

        LT_temp = zeros(2, 1);
        LT_temp(:) = LT_Gain(end, :, :);
        LT(:) = LT_temp;

        return;
    end

    % Indexes (avoid find bc Simulink sucks)
    n_low = 1;
    for i = 1:(length(Time) - 1)
        if Time(i) <= t && Time(i+1) > t
            n_low = i;
            break;
        end
    end    
    n_high = n_low + 1;
    T_low = Time(n_low);     T_high = Time(n_high);
    dt = T_high - T_low;
    
    %% Interpolation
    % LESO
        LA_low_temp = zeros(6, 3);
        LA_low_temp(:) = LA_Gain(n_low, :, :);

        LA_up_temp = zeros(6, 3);
        LA_up_temp(:) = LA_Gain(n_high, :, :);

        LA(:) = LA_low_temp + (LA_up_temp - LA_low_temp) .* ((t - T_low) / dt);

        LT_low_temp = zeros(2, 1);
        LT_low_temp(:) = LT_Gain(n_low, :, :);

        LT_up_temp = zeros(2, 1);
        LT_up_temp(:) = LT_Gain(n_high, :, :);

        LT(:) = LT_low_temp + (LT_up_temp - LT_low_temp) .* ((t - T_low) / dt);
        
    % Gains
        K_low_temp = zeros(12, 12);
        K_low_temp(:) = Feedback(n_low, :, :);
        K_low = K_low_temp';

        K_up_temp = zeros(12, 12);
        K_up_temp(:) = Feedback(n_high, :, :);
        K_up = K_up_temp';
        
        P(:) = K_low + (K_up - K_low) .* ((t - T_low) / dt);
    
    % States
        X_low = States(n_low,:);
        X_up  = States(n_high,:);
        X(:) = X_low + (X_up - X_low) .* ((t - T_low) / dt);
    
    % Inputs
        U_low = Inputs(n_low,:);
        U_up  = Inputs(n_high,:);
        U(:) = U_low + (U_up - U_low) .* ((t - T_low) / dt);
end


