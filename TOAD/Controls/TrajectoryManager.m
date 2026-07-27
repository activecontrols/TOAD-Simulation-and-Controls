%% This trajectory manager script shall take in the current time, and the
% name of a trajectory file and the current MET. It shall then find the
% corresponding gain matrix file and interpolate the gain matrices and
% state targets correspondingly. If MET > MaxT, hold final state, matrix
% and feedforward input. Outputs interpolated target state, feedforward
% input, and feedback gain matrix. 

function [X, U, K] = TrajectoryManager(t, constantsTOAD)
    %% Read the passed trajectory
    Time = constantsTOAD.Traj.Time;
    States = constantsTOAD.Traj.States;
    Inputs = constantsTOAD.Traj.Inputs;
    Feedback = constantsTOAD.Traj.Feedback;

    % Pre-allocate outputs immediately to lock fixed sizes for Simulink
    X = zeros(15, 1);
    U = zeros(4, 1);
    K = zeros(4, 15);
    
    %% Handle Boundary Conditions
    if t <= Time(1) || t == 0
        X(:) = States(1,:);
        U(:) = Inputs(1,:);
        K_temp = zeros(15, 4);
        K_temp(:) = Feedback(1, :, :); 
        K(:) = K_temp';
        return;
    elseif t >= Time(end)
        X(:) = States(end,:);
        U(:) = Inputs(end,:);
        K_temp = zeros(15, 4);
        K_temp(:) = Feedback(end, :, :);
        K(:) = K_temp';
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
    % Gains
        K_low_temp = zeros(15, 4);
        K_low_temp(:) = Feedback(n_low, :, :);
        K_low = K_low_temp';

        K_up_temp = zeros(15, 4);
        K_up_temp(:) = Feedback(n_high, :, :);
        K_up = K_up_temp';
        
        K(:) = K_low + (K_up - K_low) .* ((t - T_low) / dt);
    
    % States
        X_low = States(n_low,:);
        X_up  = States(n_high,:);
        X(:) = X_low + (X_up - X_low) .* ((t - T_low) / dt);
    
    % Inputs
        U_low = Inputs(n_low,:);
        U_up  = Inputs(n_high,:);
        U(:) = U_low + (U_up - U_low) .* ((t - T_low) / dt);
end


