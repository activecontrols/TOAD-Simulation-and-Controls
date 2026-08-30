function [X, U, K, LA, LT] = TrajectoryManager(t, X_est, constantsTOAD)
    Time = constantsTOAD.Traj.Time;
    States = constantsTOAD.Traj.States;
    Inputs = constantsTOAD.Traj.Inputs;
    
    % FBGain maps to K_trans, FBCost maps to K_rot
    K_trans_Gain = constantsTOAD.Traj.KTGain;
    K_rot_Gain = constantsTOAD.Traj.KRGain; 
    LA_Gain = constantsTOAD.Traj.LAGain;
    LT_Gain = constantsTOAD.Traj.LTGain;

    X = zeros(15, 1);
    U = zeros(4, 1);
    K = zeros(6, 6);
    LA = zeros(6, 3);
    LT = zeros(9, 6);
    
    %% Handle Boundary Conditions
    if t <= Time(1) || t == 0
        
        X(:) = States(1,:);
        U(:) = Inputs(1,:);
        
        % Merge matrices
        K(1:3, 1:6) = K_trans_Gain(1, :, :);
        K(4:6, 1:6) = K_rot_Gain(1, :, :);
        
        LA(:) = LA_Gain(1, :, :);
        LT(:) = LT_Gain(1, :, :);
        return;
    elseif t >= Time(end)
        
        X(:) = States(end,:);
        X(1:4) = [1; 0; 0; 0];
        X(8:end) = 0;
        total_mass = constantsTOAD.m_dry + X_est(14,1) + X_est(15,1);
        U(:) = [0; 0; total_mass * constantsTOAD.g; 0];
        
        K(1:3, 1:6) = K_trans_Gain(end, :, :);
        K(4:6, 1:6) = K_rot_Gain(end, :, :);
        
        LA(:) = LA_Gain(end, :, :);
        LT(:) = LT_Gain(end, :, :);
        return;
    end
    
    n_low = 1;
    for i = 1:(length(Time) - 1)
        if Time(i) <= t && Time(i+1) > t
            n_low = i;
            break;
        end
    end    
    
    n_high = n_low + 1;
    T_low = Time(n_low);     
    T_high = Time(n_high);
    dt = T_high - T_low;
    ratio = (t - T_low) / dt;
    
    %% Interpolation
    LA_low = zeros(6,3); LA_low(:) = LA_Gain(n_low,:,:);
    LA_up  = zeros(6,3); LA_up(:)  = LA_Gain(n_high,:,:);
    LA(:) = LA_low + (LA_up - LA_low) .* ratio;

    LT_low = zeros(9,6); LT_low(:) = LT_Gain(n_low,:,:);
    LT_up  = zeros(9,6); LT_up(:)  = LT_Gain(n_high,:,:);
    LT(:) = LT_low + (LT_up - LT_low) .* ratio;
        
    K_trans_low = zeros(3,6); K_trans_low(:) = K_trans_Gain(n_low,:,:);
    K_trans_up  = zeros(3,6); K_trans_up(:)  = K_trans_Gain(n_high,:,:);
    K(1:3, 1:6) = K_trans_low + (K_trans_up - K_trans_low) .* ratio;

    K_rot_low = zeros(3,6); K_rot_low(:) = K_rot_Gain(n_low,:,:);
    K_rot_up  = zeros(3,6); K_rot_up(:)  = K_rot_Gain(n_high,:,:);
    K(4:6, 1:6) = K_rot_low + (K_rot_up - K_rot_low) .* ratio;
    
    X_low = States(n_low,:);
    X_up  = States(n_high,:);
    X(:) = X_low + (X_up - X_low) .* ratio;
    
    U_low = Inputs(n_low,:);
    U_up  = Inputs(n_high,:);
    U(:) = U_low + (U_up - U_low) .* ratio;
end