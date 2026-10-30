function out = FluidProperties(fluid, mode, in1, in2, in3, in4)
% FLUIDPROPERTIES Unified thermodynamic and transport property evaluation
% Supports Nitrogen, Oxygen, and Isopropanol (IPA).
%
% Calling Signatures:
%   props = FluidProperties(fluid, 'From_u_rho', u, rho)
%   props = FluidProperties(fluid, 'From_P_T', P, T)
%   val   = FluidProperties(fluid, prop_out, in1_name, in1_val, in2_name, in2_val)

    % Standardize fluid name
    fluidLower = lower(fluid);
    if strcmp(fluidLower, 'n2') || strcmp(fluidLower, 'nitrogen')
        fluidName = 'Nitrogen';
        isCoolProp = true;
    elseif strcmp(fluidLower, 'o2') || strcmp(fluidLower, 'lox') || strcmp(fluidLower, 'oxygen')
        fluidName = 'Oxygen';
        isCoolProp = true;
    elseif strcmp(fluidLower, 'ipa') || strcmp(fluidLower, 'fu') || strcmp(fluidLower, 'isopropanol')
        fluidName = 'IPA';
        isCoolProp = false;
    else
        error('FluidProperties: Unsupported fluid "%s"', fluid);
    end

    switch lower(mode)
        case 'from_u_rho'
            u   = in1;
            rho = in2;
            out = resolveFromURho(fluidName, isCoolProp, u, rho);

        case 'from_p_t'
            P = in1;
            T = in2;
            out = resolveFromPT(fluidName, isCoolProp, P, T);

        otherwise
            % Generic property query: FluidProperties(fluid, prop_out, in1_name, in1_val, in2_name, in2_val)
            prop_out = mode;
            in1_name = in1;
            in1_val  = in2;
            in2_name = in3;
            in2_val  = in4;

            if (strcmpi(in1_name, 'P') && strcmpi(in2_name, 'T'))
                P_eval = in1_val; T_eval = in2_val;
            elseif (strcmpi(in1_name, 'T') && strcmpi(in2_name, 'P'))
                T_eval = in1_val; P_eval = in2_val;
            else
                P_eval = 101325; T_eval = 293.15;
            end

            props_eval = resolveFromPT(fluidName, isCoolProp, P_eval, T_eval);
            switch upper(prop_out)
                case {'D', 'DMASS', 'RHO', 'DENSITY'}
                    out = props_eval.rho;
                case {'H', 'HMASS', 'ENTHALPY'}
                    out = props_eval.h;
                case {'U', 'UMASS', 'INTERNAL_ENERGY'}
                    out = props_eval.u;
                case {'CP', 'CPMASS'}
                    out = props_eval.cp;
                case {'V', 'VISCOSITY', 'MU', 'DYNAMIC_VISCOSITY'}
                    out = props_eval.mu;
                case {'L', 'CONDUCTIVITY', 'K', 'THERMAL_CONDUCTIVITY'}
                    out = props_eval.k;
                case {'P', 'PRESSURE'}
                    out = props_eval.P;
                case {'T', 'TEMPERATURE'}
                    out = props_eval.T;
                case {'GAMMA'}
                    out = props_eval.gamma;
                otherwise
                    out = 0.0;
            end
    end
end

%% --- Helper: Resolve All Properties from (u, rho) ---
function props = resolveFromURho(fluidName, isCoolProp, u, rho)
    persistent GI_P_n2 GI_T_n2 GI_gamma_n2
    persistent GI_P_ox GI_T_ox GI_gamma_ox
    persistent TablesLoaded

    % Initialize Interpolants from precomputed tables if not already loaded
    if isempty(TablesLoaded)
        TablesLoaded = false;
        try
            % Check relative paths for pre-generated tables
            matPaths = {
                fullfile(pwd, 'sandbox', 'experiments', 'Fluid Properties', 'FluidPropertyTables.mat'), ...
                fullfile(pwd, 'Fluid Properties', 'FluidPropertyTables.mat'), ...
                fullfile(fileparts(mfilename('fullpath')), '..', 'Fluid Properties', 'FluidPropertyTables.mat')
            };
            matFile = '';
            for mp = 1:length(matPaths)
                if exist(matPaths{mp}, 'file')
                    matFile = matPaths{mp};
                    break;
                end
            end

            if ~isempty(matFile)
                Tbl = load(matFile);
                GI_P_n2     = griddedInterpolant({Tbl.rho_vec_n2, Tbl.u_vec_n2}, Tbl.P_grid_n2, 'linear', 'nearest');
                GI_T_n2     = griddedInterpolant({Tbl.rho_vec_n2, Tbl.u_vec_n2}, Tbl.T_grid_n2, 'linear', 'nearest');
                GI_gamma_n2 = griddedInterpolant({Tbl.rho_vec_n2, Tbl.u_vec_n2}, Tbl.GAMMA_grid_n2, 'linear', 'nearest');

                GI_P_ox     = griddedInterpolant({Tbl.rho_vec_ox, Tbl.u_vec_ox}, Tbl.P_grid_ox, 'linear', 'nearest');
                GI_T_ox     = griddedInterpolant({Tbl.rho_vec_ox, Tbl.u_vec_ox}, Tbl.T_grid_ox, 'linear', 'nearest');
                GI_gamma_ox = griddedInterpolant({Tbl.rho_vec_ox, Tbl.u_vec_ox}, Tbl.GAMMA_grid_ox, 'linear', 'nearest');

                TablesLoaded = true;
            end
        catch
            TablesLoaded = false;
        end
    end

    % Safeguards
    rho = max(rho, 1e-4);

    if isCoolProp
        if TablesLoaded
            if strcmp(fluidName, 'Nitrogen')
                P     = max(GI_P_n2(rho, u), 1000);
                T     = max(GI_T_n2(rho, u), 60);
                h     = u + P / rho;
                gamma = GI_gamma_n2(rho, u);
                mu    = 1.78e-5 * ((T / 300)^0.72);
                k     = 0.026 * ((T / 300)^0.85);
                cp    = gamma * 296.8 / max(gamma - 1, 0.05);
            else % Oxygen
                if u < -50000 % Subcooled Liquid Oxygen
                    T     = max(60.0, 90.0 + (u + 1.25e5) / 1700.0);
                    rho_calc = max(1141.0 - 4.5 * (T - 90.0), 900.0);
                    P     = 101325;
                    h     = u + P / rho_calc;
                    gamma = 1.41;
                    mu    = 1.8e-4;
                    k     = 0.15;
                    cp    = 1700.0;
                else
                    P     = max(GI_P_ox(rho, u), 1000);
                    T     = max(GI_T_ox(rho, u), 60);
                    h     = u + P / rho;
                    gamma = GI_gamma_ox(rho, u);
                    mu    = 2.0e-5 * ((T / 300)^0.75);
                    k     = 0.026 * ((T / 300)^0.85);
                    cp    = gamma * 259.8 / max(gamma - 1, 0.05);
                end
            end
        else
            % Fallback to direct CoolProp PropsSI calls if tables are not available
            try
                P = double(py.CoolProp.CoolProp.PropsSI('P', 'Umass', u, 'Dmass', rho, fluidName));
                T = double(py.CoolProp.CoolProp.PropsSI('T', 'Umass', u, 'Dmass', rho, fluidName));
                h = double(py.CoolProp.CoolProp.PropsSI('Hmass', 'Umass', u, 'Dmass', rho, fluidName));
            catch
                P = max(101325, rho * 296.8 * 293.15);
                T = max(70, u / 1000 + 293.15);
                h = u + P / rho;
            end

            try
                cp = double(py.CoolProp.CoolProp.PropsSI('Cpmass', 'P', P, 'T', T, fluidName));
                cv = double(py.CoolProp.CoolProp.PropsSI('Cvmass', 'P', P, 'T', T, fluidName));
                gamma = max(1.1, min(1.67, cp / max(cv, 1e-3)));
            catch
                cp = 1040;
                gamma = 1.4;
            end

            try
                mu = double(py.CoolProp.CoolProp.PropsSI('V', 'P', P, 'T', T, fluidName));
            catch
                mu = 1.8e-5;
            end

            try
                k = double(py.CoolProp.CoolProp.PropsSI('L', 'P', P, 'T', T, fluidName));
            catch
                k = 0.026;
            end
        end

    else
        % IPA (Incompressible Liquid Evaluation)
        % Reference state: 293.15 K, 101325 Pa -> u ~ 0, rho ~ 786
        cp_ipa = 2600; % J/(kg-K) approximate average
        T = 293.15 + u / cp_ipa;
        T = max(190, min(400, T)); % Bounds

        % Pure analytical calibrated IPA liquid properties (zero interp2 overhead)
        dT = T - 293.15;
        cp = 2570.0 + 5.2 * dT;
        mu = 2.4e-3 * exp(-0.025 * max(-50, min(100, dT)));
        k  = 0.14;

        P = 101325; % Bulk thermodynamic pressure; hydraulic pressure set by network
        h = u + P / max(rho, 10);
        gamma = 1.15;
    end

    props.P     = P;
    props.T     = T;
    props.h     = h;
    props.u     = u;
    props.rho   = rho;
    props.cp    = cp;
    props.gamma = gamma;
    props.mu    = mu;
    props.k     = k;
end

%% --- Helper: Resolve All Properties from (P, T) ---
function props = resolveFromPT(fluidName, isCoolProp, P, T)
    P = max(P, 1000); % Positive pressure
    T = max(T, 50);   % Positive temperature

    if isCoolProp
        if strcmp(fluidName, 'Oxygen')
            try
                rho = double(py.CoolProp.CoolProp.PropsSI('D', 'P', P, 'T', T, 'Oxygen'));
                u   = double(py.CoolProp.CoolProp.PropsSI('Umass', 'P', P, 'T', T, 'Oxygen'));
                h   = double(py.CoolProp.CoolProp.PropsSI('Hmass', 'P', P, 'T', T, 'Oxygen'));
                cp  = double(py.CoolProp.CoolProp.PropsSI('Cpmass', 'P', P, 'T', T, 'Oxygen'));
                gamma = double(py.CoolProp.CoolProp.PropsSI('isentropic_expansion_coefficient', 'P', P, 'T', T, 'Oxygen'));
                mu  = double(py.CoolProp.CoolProp.PropsSI('V', 'P', P, 'T', T, 'Oxygen'));
                k   = double(py.CoolProp.CoolProp.PropsSI('L', 'P', P, 'T', T, 'Oxygen'));
            catch
                if T < 130
                    % Subcooled Liquid Oxygen (LOX)
                    rho = max(1141.0 - 4.5 * (T - 90.0), 900.0);
                    cp  = 1700.0;
                    u   = cp * (T - 90.0) - 1.3377e5;
                    h   = u + P / rho;
                    gamma = 1.41;
                    mu  = 1.8e-4;
                    k   = 0.15;
                else
                    % Gaseous Oxygen
                    R_ox = 259.8;
                    rho = P / (R_ox * T);
                    cp  = 918.0;
                    cv  = 658.0;
                    u   = cv * (T - 293.15) + 2.0e5;
                    h   = u + P / rho;
                    gamma = 1.395;
                    mu  = 2.0e-5;
                    k   = 0.026;
                end
            end
        else
            % Nitrogen (Gas / Supercritical)
            try
                rho = double(py.CoolProp.CoolProp.PropsSI('D', 'P', P, 'T', T, 'Nitrogen'));
                u   = double(py.CoolProp.CoolProp.PropsSI('Umass', 'P', P, 'T', T, 'Nitrogen'));
                h   = double(py.CoolProp.CoolProp.PropsSI('Hmass', 'P', P, 'T', T, 'Nitrogen'));
                cp  = double(py.CoolProp.CoolProp.PropsSI('Cpmass', 'P', P, 'T', T, 'Nitrogen'));
                gamma = double(py.CoolProp.CoolProp.PropsSI('isentropic_expansion_coefficient', 'P', P, 'T', T, 'Nitrogen'));
                mu  = double(py.CoolProp.CoolProp.PropsSI('V', 'P', P, 'T', T, 'Nitrogen'));
                k   = double(py.CoolProp.CoolProp.PropsSI('L', 'P', P, 'T', T, 'Nitrogen'));
            catch
                R_n2 = 296.8;
                % Real gas compressibility correction for high-pressure COPV
                Z = 1.0 + 0.00045 * (P / 1e5) * (293.15 / T);
                rho = P / (Z * R_n2 * T);
                cp  = 1040.0;
                cv  = 743.0;
                u   = cv * (T - 293.15) + 1.559e5;
                h   = u + P / rho;
                gamma = 1.40;
                mu  = 1.78e-5;
                k   = 0.026;
            end
        end
    else
        % IPA (Isopropanol - subcooled liquid analytical model)
        dT  = T - 293.15;
        rho = 786.0 - 0.85 * dT;
        cp  = 2570.0 + 5.2 * dT;
        mu  = 2.4e-3 * exp(-0.025 * max(-50, min(100, dT)));
        k   = 0.14;
        gamma = 1.15;
        u   = cp * dT;
        h   = u + P / rho;
    end

    props.P     = P;
    props.T     = T;
    props.h     = h;
    props.u     = u;
    props.rho   = rho;
    props.cp    = cp;
    props.gamma = gamma;
    props.mu    = mu;
    props.k     = k;
end

%% --- Helper: Evaluate Single Property for IPA ---
function val = evalIPAProperty(prop_out, in1_name, in1_val, in2_name, in2_val) %#ok<DEFNU>
    % Map inputs to P and T
    if (strcmpi(in1_name, 'P') && strcmpi(in2_name, 'T'))
        P = in1_val; T = in2_val;
    elseif (strcmpi(in1_name, 'T') && strcmpi(in2_name, 'P'))
        T = in1_val; P = in2_val;
    else
        P = 101325; T = 293.15;
    end

    dT  = T - 293.15;
    rho = 786.0 - 0.85 * dT;
    cp  = 2570.0 + 5.2 * dT;
    mu  = 2.4e-3 * exp(-0.025 * max(-50, min(100, dT)));
    k   = 0.14;

    switch upper(prop_out)
        case {'D', 'DMASS', 'RHO'}
            val = rho;
        case {'CP', 'CPMASS'}
            val = cp;
        case {'V', 'VISCOSITY'}
            val = mu;
        case {'L', 'CONDUCTIVITY'}
            val = k;
        case {'HMASS', 'H'}
            val = cp * dT + P / rho;
        case {'UMASS', 'U'}
            val = cp * dT;
        otherwise
            val = 0;
    end
end
