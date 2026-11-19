function [mdot, isChoked, flowDir, Qdot] = CalculateLinkFlow(Link, P_up, P_down, UpProps, DownProps, dt)
% CALCULATELINKFLOW Computes mass flow rate and heat transfer for fluid & thermal links.
% Unified vectorized formulation:
% - Exact compressible gas dynamics (subsonic & choked) for gas-dominated flows
% - Incompressible hydraulic formulation with pipe friction for liquid-dominated flows
% - Dynamic dome regulator Supply Pressure Effect (SPE) and droop throttling
% - Smooth tanh check valve opening mask
% - Unified 'Thermal' link handling returning Qdot = G * (T_up - T_down)
%
% Calling Signatures:
%   Scalar:     [mdot, isChoked, flowDir, Qdot] = CalculateLinkFlow(Link, P_up, P_down, UpProps, DownProps, dt)
%   Vectorized: [mdot_vec, isChoked_vec, flowDir_vec, Qdot_vec] = CalculateLinkFlow(Links, P_up_vec, P_down_vec, UpProps, DownProps, dt)

    if nargin < 6 || isempty(dt), dt = 0.005; end

    N = numel(Link);
    P_up = reshape(P_up, [N, 1]);
    P_down = reshape(P_down, [N, 1]);

    % Fast-path: Native C++ MEX accelerator if available
    persistent hasMex;
    if isempty(hasMex)
        hasMex = (exist('CalculateLinkFlow_mex', 'file') == 3);
    end
    if hasMex
        [mdot, isChoked, flowDir, Qdot] = CalculateLinkFlow_mex(Link, P_up, P_down, UpProps, DownProps, dt);
        return;
    end

    % Pre-allocate outputs
    mdot     = zeros(N, 1);
    isChoked = false(N, 1);
    flowDir  = ones(N, 1);
    Qdot     = zeros(N, 1);

    % Extract link classification flags
    linkTypes = lower({Link.Type});
    isThermal   = strcmp(linkTypes, 'thermal')';
    isRegulator = strcmp(linkTypes, 'regulator')';
    isCheck     = strcmp(linkTypes, 'check')';
    isPipe      = strcmp(linkTypes, 'pipe')';

    %% 1. Handle Thermal Links (Qdot = G * (T_up - T_down), mdot = 0)
    if any(isThermal)
        thIdx = find(isThermal);
        T_up_th = [UpProps(thIdx).T]';
        T_down_th = [DownProps(thIdx).T]';
        
        G_vec = zeros(numel(thIdx), 1);
        for k = 1:numel(thIdx)
            idx = thIdx(k);
            if isfield(Link(idx), 'G') && ~isempty(Link(idx).G)
                G_vec(k) = Link(idx).G;
            elseif isfield(Link(idx), 'UA') && ~isempty(Link(idx).UA)
                G_vec(k) = Link(idx).UA;
            elseif isfield(Link(idx), 'Conductance') && ~isempty(Link(idx).Conductance)
                G_vec(k) = Link(idx).Conductance;
            end
        end
        Qdot(thIdx) = G_vec .* (T_up_th - T_down_th);
    end

    fluidIdx = find(~isThermal);
    if isempty(fluidIdx)
        return;
    end

    %% 2. Vectorized Fluid Link Hydraulics
    N_fl = numel(fluidIdx);
    P_u = P_up(fluidIdx);
    P_d = P_down(fluidIdx);
    
    DeltaP_Vec = P_u - P_d;
    flowDir(fluidIdx) = sign(DeltaP_Vec + 1e-12);

    % Determine flow direction and upwind properties
    isForward = (DeltaP_Vec >= 0);
    P_in  = P_u .* isForward + P_d .* (~isForward);
    P_out = P_d .* isForward + P_u .* (~isForward);

    T_u_all = [UpProps(fluidIdx).T]';
    T_d_all = [DownProps(fluidIdx).T]';
    T_in = T_u_all .* isForward + T_d_all .* (~isForward);
    T_in = max(T_in, 50.0);

    rho_u_all = [UpProps(fluidIdx).rho]';
    rho_d_all = [DownProps(fluidIdx).rho]';
    Rho_in = max(rho_u_all .* isForward + rho_d_all .* (~isForward), 0.01);

    gamma_u_all = [UpProps(fluidIdx).gamma]';
    gamma_d_all = [DownProps(fluidIdx).gamma]';
    Gamma_in = max(gamma_u_all .* isForward + gamma_d_all .* (~isForward), 1.05);

    % Determine active Cv, Area, and Type properties
    Cv_Vec = zeros(N_fl, 1);
    CdA_Vec = zeros(N_fl, 1);
    isReg_sub   = isRegulator(fluidIdx);
    isCheck_sub = isCheck(fluidIdx);
    isPipe_sub  = isPipe(fluidIdx);

    for k = 1:N_fl
        idx = fluidIdx(k);
        L = Link(idx);
        
        % Read Cv
        if isfield(L, 'Cv') && ~isempty(L.Cv)
            Cv_Vec(k) = L.Cv;
        elseif isfield(L, 'MaxCv') && ~isempty(L.MaxCv)
            st = 1.0;
            if isfield(L, 'State'), st = L.State; end
            Cv_Vec(k) = L.MaxCv * st;
        end

        % Read Area
        if isfield(L, 'A') && ~isempty(L.A) && L.A > 0
            cd_val = 0.70;
            if isfield(L, 'Cd') && ~isempty(L.Cd) && L.Cd > 0, cd_val = L.Cd; end
            CdA_Vec(k) = L.A * cd_val;
        end
    end

    % 2a. Dynamic Regulator Stroke & SPE Calculation
    if any(isReg_sub)
        regIdx = find(isReg_sub);
        for r = 1:numel(regIdx)
            k = regIdx(r);
            idx = fluidIdx(k);
            L = Link(idx);
            
            P_Set_val = 550 * 6894.757;
            if isfield(L, 'P_set'), P_Set_val = L.P_set; end
            
            SPE_val = 0.003;
            if isfield(L, 'SPE'), SPE_val = L.SPE; end
            
            Droop_val = 30 * 6894.757;
            if isfield(L, 'Droop'), Droop_val = L.Droop; end

            % Target pressure adjusted by Supply Pressure Effect (Unbalanced Poppet)
            P_Target = P_Set_val - SPE_val * P_u(k);
            
            % Pressure Error (How far below target are we?)
            P_Error = P_Target - P_d(k);
            
            % Calculate Stroke: 0 = closed (P_down >= P_Target), 1 = open (P_Error >= Droop)
            Norm_Err = P_Error / (Droop_val + 1e-6);
            Open_Fraction = max(0.0, min(1.0, Norm_Err));
            
            % Dynamically modulate regulator Cv
            max_cv = Cv_Vec(k);
            if max_cv <= 0 && isfield(L, 'MaxCv'), max_cv = L.MaxCv; end
            Cv_Vec(k) = max_cv * Open_Fraction;
        end
    end

    % Convert Cv to equivalent CdA for gas flow: CdA = Cv * 2.402e-5 * sqrt(0.5)
    for k = 1:N_fl
        if Cv_Vec(k) > 0
            CdA_Vec(k) = Cv_Vec(k) * 2.402e-5 * sqrt(0.5);
        end
    end

    % 2b. Flow evaluation: Gas (< 600 kg/m^3) vs Liquid (>= 600 kg/m^3)
    isGas = (Rho_in < 600.0);
    raw_mdot = zeros(N_fl, 1);
    choked_flags = false(N_fl, 1);

    if any(isGas)
        gIdx = find(isGas);
        P_i = P_in(gIdx);
        P_o = P_out(gIdx);
        T_i = T_in(gIdx);
        rho_i = Rho_in(gIdx);
        gam = Gamma_in(gIdx);
        cda = CdA_Vec(gIdx);

        R_spec = P_i ./ (rho_i .* T_i);
        PR = min(1.0, max(0.0, P_o ./ max(P_i, 1.0)));
        PR_crit = (2.0 ./ (gam + 1.0)) .^ (gam ./ (gam - 1.0));

        is_chk = (PR <= PR_crit);
        choked_flags(gIdx) = is_chk;

        % Choked flow rate
        choke_term = sqrt(gam ./ (R_spec .* T_i)) .* (2.0 ./ (gam + 1.0)) .^ ((gam + 1.0) ./ (2.0 * (gam - 1.0)));
        mdot_chk = cda .* P_i .* choke_term;

        % Subsonic compressible flow rate
        pr_pow1 = PR .^ (2.0 ./ gam);
        pr_pow2 = PR .^ ((gam + 1.0) ./ gam);
        sub_term = sqrt(max(0.0, 2.0 .* P_i .* rho_i .* (gam ./ (gam - 1.0)) .* (pr_pow1 - pr_pow2)));
        mdot_sub = cda .* sub_term;

        raw_mdot(gIdx) = is_chk .* mdot_chk + (~is_chk) .* mdot_sub;
    end

    if any(~isGas)
        lIdx = find(~isGas);
        dp_raw = abs(DeltaP_Vec(lIdx));
        rho_l  = Rho_in(lIdx);
        
        for p = 1:numel(lIdx)
            k = lIdx(p);
            idx = fluidIdx(k);
            L = Link(idx);
            if isPipe_sub(k)
                zeta = 20.0;
                if isfield(L, 'Zeta') && L.Zeta > 0, zeta = L.Zeta; end
                A_p = 5e-4;
                if isfield(L, 'A') && L.A > 0, A_p = L.A; end
                raw_mdot(k) = A_p * sqrt((2.0 * rho_l(p) * dp_raw(p)) / zeta);
            elseif Cv_Vec(k) > 0
                raw_mdot(k) = Cv_Vec(k) * 2.402e-5 * sqrt(rho_l(p) * dp_raw(p));
            else
                raw_mdot(k) = CdA_Vec(k) * sqrt(2.0 * rho_l(p) * dp_raw(p));
            end
        end
    end

    % 2c. Direction & Check Valve Masking
    sign_flow = 2.0 * isForward - 1.0;
    isOneWay = isCheck_sub | isReg_sub;
    for k = 1:N_fl
        idx = fluidIdx(k);
        if isfield(Link(idx), 'IsOneWay') && Link(idx).IsOneWay
            isOneWay(k) = true;
        end
    end

    Check_Open_Factor = 0.5 + 0.5 * tanh(50.0 * DeltaP_Vec);
    CheckMask = (~isOneWay) + (isOneWay .* Check_Open_Factor);
    FinalMassflow = sign_flow .* raw_mdot .* CheckMask;

    mdot(fluidIdx)     = FinalMassflow;
    isChoked(fluidIdx) = choked_flags;
end
