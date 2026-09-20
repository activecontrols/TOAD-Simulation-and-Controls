function [mdot, isChoked, flowDir] = CalculateLinkFlow(Link, P_up, P_down, UpProps, DownProps, dt)
% CALCULATELINKFLOW Computes quasi-steady mass flow rate for a single link
% (Valve, Orifice, Pipe, Check, or Regulator).
%
% Inputs:
%   Link      - Struct defining link properties (Type, Cv, A, Cd, Zeta, P_set, SPE, Droop, etc.)
%   P_up      - Upstream node pressure [Pa]
%   P_down    - Downstream node pressure [Pa]
%   UpProps   - Thermodynamic properties struct of upstream node
%   DownProps - Thermodynamic properties struct of downstream node
%   dt        - (Optional) Current timestep [s] for Courant limiter (default: 0.005)
%
% Outputs:
%   mdot      - Net mass flow rate [kg/s] (positive = Up -> Down)
%   isChoked  - Boolean indicating choked compressible flow
%   flowDir   - +1 for forward (Up -> Down), -1 for reverse (Down -> Up)

    if nargin < 6 || isempty(dt), dt = 0.005; end

    DeltaP = P_up - P_down;
    if abs(DeltaP) < 1.0 % Less than 1 Pa difference
        mdot = 0.0;
        isChoked = false;
        flowDir = 1;
        return;
    end

    if DeltaP >= 0
        flowDir = 1;
        P_in = P_up;
        P_out = P_down;
        Props_in = UpProps;
    else
        flowDir = -1;
        P_in = P_down;
        P_out = P_up;
        Props_in = DownProps;
    end

    rho_in = max(Props_in.rho, 0.01);
    gamma  = max(Props_in.gamma, 1.1);
    T_in   = max(Props_in.T, 50);

    % Determine if flow is gas-dominated or liquid-dominated
    isGasFlow = (rho_in < 200); % Gas or vapor density (< 200 kg/m^3)

    isChoked = false;
    linkType = lower(Link.Type);

    switch linkType
        case 'regulator'
            % Forward-flow only
            if flowDir == -1
                mdot = 0.0;
                return;
            end

            % Supply pressure effect (SPE) and droop throttling
            P_set_Pa = Link.P_set;
            SPE = 0.0;
            if isfield(Link, 'SPE'), SPE = Link.SPE; end

            % Supply pressure effect
            if isfield(Link, 'P_tank0')
                P_target = P_set_Pa + SPE * (Link.P_tank0 - P_up);
            elseif isfield(Link, 'P_up0')
                P_target = P_set_Pa + SPE * (Link.P_up0 - P_up);
            else
                P_target = P_set_Pa - SPE * (P_up - P_set_Pa);
            end

            if isfield(Link, 'P_band') && Link.P_band > 0
                Droop_Pa = Link.P_band;
            elseif isfield(Link, 'Droop') && Link.Droop > 0
                Droop_Pa = Link.Droop;
            else
                Droop_Pa = 0.03 * P_target;
            end
            Droop_Pa = max(Droop_Pa, 100);

            P_error  = P_target - P_down;
            openFrac = max(0.0, min(1.0, P_error / Droop_Pa));

            if openFrac <= 0.0 || P_up <= P_down
                mdot = 0.0;
                return;
            end

            % Base CdA from CdA_MAX or Cv
            if isfield(Link, 'CdA_MAX') && Link.CdA_MAX > 0
                CdA_base = Link.CdA_MAX;
            elseif isfield(Link, 'CdA') && Link.CdA > 0
                CdA_base = Link.CdA;
            else
                CdA_base = Link.Cv * 2.402e-5 * sqrt(1.0 / 2.0);
            end

            % Modulate by valve opening (Cv if used as a valve fraction, default 1.0)
            actuatorFrac = 1.0;
            if isfield(Link, 'Cv') && isfield(Link, 'CdA_MAX')
                actuatorFrac = max(0.0, min(1.0, Link.Cv));
            end

            CdA = CdA_base * openFrac * actuatorFrac;

            [mdot, isChoked] = calcGasOrificeFlow(CdA, P_up, P_down, rho_in, T_in, gamma);

        case 'check'
            % Forward flow only with cracking pressure
            P_crack = 0.0;
            if isfield(Link, 'P_crack'), P_crack = Link.P_crack; end

            if flowDir == -1 || DeltaP < P_crack
                mdot = 0.0;
                return;
            end

            effectiveDP = DeltaP - P_crack;
            if isfield(Link, 'Cv') && ~isempty(Link.Cv) && Link.Cv > 0
                mdot = Link.Cv * 2.402e-5 * sqrt(rho_in * effectiveDP);
            else
                CdA = Link.A * 0.6;
                mdot = CdA * sqrt(2 * rho_in * effectiveDP);
            end

        case 'orifice'
            A = Link.A;
            Cd = 0.7;
            if isfield(Link, 'Cd') && ~isempty(Link.Cd), Cd = Link.Cd; end
            CdA = Cd * A;

            % Injectors and nozzles are directional forward flow only
            isOneWay = contains(lower(Link.Name), 'inj') || contains(lower(Link.Name), 'nozzle');
            if isOneWay && flowDir == -1
                mdot = 0.0;
                isChoked = false;
                return;
            end

            if isGasFlow
                [mdot_mag, isChoked] = calcGasOrificeFlow(CdA, P_in, P_out, rho_in, T_in, gamma);
                mdot = flowDir * mdot_mag;
            else
                mdot = flowDir * CdA * sqrt(2 * rho_in * abs(DeltaP));
            end

        case {'solenoid', 'throttle', 'valve'}
            Cv = Link.Cv;
            if Cv <= 1e-6
                mdot = 0.0;
                return;
            end

            if isGasFlow
                % Convert Cv to equivalent CdA for compressible gas flow
                CdA = Cv * 2.402e-5 * sqrt(1.0 / 2.0);
                [mdot_mag, isChoked] = calcGasOrificeFlow(CdA, P_in, P_out, rho_in, T_in, gamma);
                mdot = flowDir * mdot_mag;
            else
                mdot = flowDir * Cv * 2.402e-5 * sqrt(rho_in * abs(DeltaP));
            end

        case 'pipe'
            % Frictional line pressure drop: DeltaP = zeta * (mdot^2) / (2 * rho * A^2)
            A = Link.A;
            zeta = max(Link.Zeta, 0.1);
            mdot = flowDir * A * sqrt((2 * rho_in * abs(DeltaP)) / zeta);

        otherwise
            % Default orifice
            CdA = 1e-6;
            if isfield(Link, 'A'), CdA = Link.A; end
            mdot = flowDir * CdA * sqrt(2 * rho_in * abs(DeltaP));
    end

    % Physical Courant mass-flux limiter: cannot drain > 50% of source mass per step
    if isfield(Props_in, 'm') && Props_in.m > 0
        maxFlow = 0.5 * Props_in.m / max(dt, 1e-4);
        mdot = sign(mdot) * min(abs(mdot), maxFlow);
    end
end

%% --- Helper: Compressible Gas Orifice Flow with Choking ---
function [mdot, isChoked] = calcGasOrificeFlow(CdA, P_in, P_out, rho_in, T_in, gamma)
    if P_out >= P_in
        mdot = 0.0;
        isChoked = false;
        return;
    end

    % Specific gas constant
    R_gas = P_in / (rho_in * T_in);
    PR = P_out / P_in;
    PR_crit = (2 / (gamma + 1))^(gamma / (gamma - 1));

    if PR <= PR_crit
        % Choked Flow
        isChoked = true;
        mdot = CdA * P_in * sqrt(gamma / (R_gas * T_in)) * ...
               (2 / (gamma + 1))^((gamma + 1) / (2 * (gamma - 1)));
    else
        % Subsonic Compressible Flow
        isChoked = false;
        mdot = CdA * sqrt(2 * P_in * rho_in * (gamma / (gamma - 1)) * ...
               (PR^(2 / gamma) - PR^((gamma + 1) / gamma)));
    end
end
