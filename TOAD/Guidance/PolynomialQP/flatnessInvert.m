function [X, Wrench, AccelCg] = flatnessInvert(QPTraj, t_vec, HDGRef, Constants, Limits)
%% Differential Flatness Inversion at Centre of Percussion
    nPts = numel(t_vec);
    P_pos   = PolyEval(QPTraj, t_vec, 0);
    P_vel   = PolyEval(QPTraj, t_vec, 1);
    P_acc   = PolyEval(QPTraj, t_vec, 2);
    P_jerk  = PolyEval(QPTraj, t_vec, 3);
    P_snap  = PolyEval(QPTraj, t_vec, 4);

    g = Limits.g;
    G_vec = [0; 0; g];
    rho = QPTraj.FlatOffset;

    % Heading reference vector in horizontal plane
    x_C = [HDGRef(1); HDGRef(2); 0];
    if norm(x_C) < 1e-9
        x_C = [1; 0; 0];
    else
        x_C = x_C / norm(x_C);
    end

    % First pass: Attitude, angular velocity, angular acceleration, and specific force
    q_s     = zeros(4, nPts);
    b1_s    = zeros(3, nPts);
    b2_s    = zeros(3, nPts);
    b3_s    = zeros(3, nPts);
    b3d_s   = zeros(3, nPts);
    b3dd_s  = zeros(3, nPts);
    w_s     = zeros(3, nPts);
    wd_s    = zeros(3, nPts);
    n_s     = zeros(1, nPts);

    for k = 1:nPts
        e  = P_acc(:, k) + G_vec;
        ed = P_jerk(:, k);
        edd = P_snap(:, k);

        n = norm(e);
        if n < 1e-6
            e = [0; 0; 1e-6];
            n = 1e-6;
        end
        u = e / n;
        nd = dot(u, ed);
        ud = (ed - u * nd) / n;
        ndd = dot(ud, ed) + dot(u, edd);
        udd = (edd - 2 * ud * nd - u * ndd) / n;

        % Body triad [b1, b2, b3] with b3 = u and b2 = norm(b3 x x_C)
        c  = cross(u, x_C);
        nc = norm(c);
        if nc < 1e-6
            x_C_alt = [0; 1; 0];
            c  = cross(u, x_C_alt);
            nc = norm(c);
        end
        b2 = c / nc;
        b1 = cross(b2, u);

        % First and second time derivatives of b2 and b1
        cd   = cross(ud, x_C);
        ncd  = dot(b2, cd);
        b2d  = (cd - b2 * ncd) / nc;
        cdd  = cross(udd, x_C);
        ncdd = dot(b2d, cd) + dot(b2, cdd);
        b2dd = (cdd - 2 * b2d * ncd - b2 * ncdd) / nc;

        b1d  = cross(b2d, u) + cross(b2, ud);
        b1dd = cross(b2dd, u) + 2 * cross(b2d, ud) + cross(b2, udd);

        % Body angular velocity w_B and angular acceleration wd_B
        wx = -dot(ud, b2);
        wy =  dot(ud, b1);
        wz =  dot(b1d, b2);

        wxd = -dot(udd, b2) - dot(ud, b2d);
        wyd =  dot(udd, b1) + dot(ud, b1d);
        wzd =  dot(b1dd, b2) + dot(b1d, b2d);

        % Quaternion from DCM C_IB = [b1, b2, b3] (DCM_Quat_Conversion expects C_IB)
        C_IB = [b1, b2, u];
        qk = DCM_Quat_Conversion(C_IB);
        if k > 1 && dot(q_s(:, k - 1), qk) < 0
            qk = -qk;
        end

        q_s(:, k)    = qk;
        b1_s(:, k)   = b1;
        b2_s(:, k)   = b2;
        b3_s(:, k)   = u;
        b3d_s(:, k)  = ud;
        b3dd_s(:, k) = udd;
        w_s(:, k)    = [wx; wy; wz];
        wd_s(:, k)   = [wxd; wyd; wzd];
        n_s(k)       = n;
    end

    % CG kinematics from P-space: r_cg = P - rho * b3
    r_s     = P_pos - rho * b3_s;
    v_s     = P_vel - rho * b3d_s;
    AccelCg = P_acc - rho * b3dd_s;

    % Second pass: Mass depletion and wrench synthesis
    ml_s = zeros(1, nPts);
    mi_s = zeros(1, nPts);
    if Limits.Vehicle == 1
        ml_s(1) = Limits.OxMass;
        mi_s(1) = Limits.FuMass;
    end

    F_body   = zeros(3, nPts);
    Tau_body = zeros(3, nPts);
    CGz_vec  = zeros(1, nPts);

    for k = 1:nPts
        m_k = Limits.m_dry + ml_s(k) + mi_s(k);
        if Limits.HasVehicle
            [J_tot, CGz] = ComputeJtot(ml_s(k), mi_s(k), Constants);
        else
            J_tot = Limits.J_init;
            CGz   = Limits.CGz_init;
        end
        CGz_vec(k) = CGz;

        % Rotational Euler dynamics: M_B = J_tot * wd + w x (J_tot * w)
        wk   = w_s(:, k);
        wdk  = wd_s(:, k);
        tauk = J_tot * wdk + cross(wk, J_tot * wk);
        Tau_body(:, k) = tauk;

        % Exact body force required to track CG acceleration: F_B = m * C_BI * (a_cg + G)
        C_BI = [b1_s(:, k), b2_s(:, k), b3_s(:, k)].';
        Fk   = m_k * (C_BI * (AccelCg(:, k) + G_vec));
        F_body(:, k) = Fk;

        % Forward integrate propellant depletion for Vehicle 1
        if k < nPts && Limits.Vehicle == 1
            dt_k = t_vec(k + 1) - t_vec(k);
            T_mag = norm(Fk);
            OF = Limits.OF;
            mdot_tot = (T_mag / Limits.MaxThrust) * Limits.MaxMdot;
            ml_s(k + 1) = max(0, ml_s(k) - dt_k * (OF / (1 + OF)) * mdot_tot);
            mi_s(k + 1) = max(0, mi_s(k) - dt_k * (1 / (1 + OF)) * mdot_tot);
        end
    end

    X = [q_s; r_s; v_s; w_s; ml_s; mi_s];
    Wrench = struct( ...
        'F_body',   F_body, ...
        'Tau_body', Tau_body, ...
        'CGz',      CGz_vec ...
    );
end
