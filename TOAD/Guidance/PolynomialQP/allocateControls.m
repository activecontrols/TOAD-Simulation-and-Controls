function U = allocateControls(Wrench, Limits)
%% Algebraic Control Allocation Matching 6-DoF Plant Dynamics
    F_body   = Wrench.F_body;
    Tau_body = Wrench.Tau_body;
    CGz_vec  = Wrench.CGz;
    nPts     = size(F_body, 2);
    U        = zeros(4, nPts);

    for k = 1:nPts
        Fk   = F_body(:, k);
        tauk = Tau_body(:, k);
        l_k  = max(1e-3, CGz_vec(k));

        % Axial force from translational dynamics
        Fz = max(1e-3, Fk(3));

        if Limits.Vehicle == 0
            % Vehicle 0 (ASTRA): M_B = [0; 0; -l] x T_B + roll * t_dir
            % Coupling: tau_x = l*Fy + roll*(Fx/T), tau_y = -l*Fx + roll*(Fy/T), tau_z = roll*(Fz/T)
            % Initialize from pure moment arm, then 2-step exact fixed-point refinement for roll * t_dir
            Fx = -tauk(2) / l_k;
            Fy =  tauk(1) / l_k;
            T_mag = sqrt(Fx^2 + Fy^2 + Fz^2);
            roll_cmd = tauk(3) * (T_mag / Fz);
            for it = 1:2
                Fx = (-tauk(2) + roll_cmd * (Fy / T_mag)) / l_k;
                Fy = ( tauk(1) - roll_cmd * (Fx / T_mag)) / l_k;
                T_mag = sqrt(Fx^2 + Fy^2 + Fz^2);
                roll_cmd = tauk(3) * (T_mag / Fz);
            end
        else
            % Vehicle 1 (TOAD): M_B = [0; 0; -l] x T_B + [0; 0; roll]
            Fx = -tauk(2) / l_k;
            Fy =  tauk(1) / l_k;
            T_mag = sqrt(Fx^2 + Fy^2 + Fz^2);
            roll_cmd = tauk(3);
        end

        % Exact trigonometric gimbal inversion for T_B = T * [cos(th)*sin(ph); -sin(th); cos(th)*cos(ph)]
        sin_th = max(-1, min(1, -Fy / T_mag));
        theta  = asin(sin_th);
        phi    = atan2(Fx, Fz);

        U(:, k) = [theta; phi; T_mag; roll_cmd];
    end
end
