function out = FastFluidProps(fluid, mode, in1, in2, in3, in4)
% FASTFLUIDPROPS High-speed drop-in replacement for FluidProperties.m
% Uses persistent griddedInterpolant lookup tables pre-generated from CoolProp
% via generateFluidPropertyTables.m. Automatically falls back to direct CoolProp
% evaluation if tables are not yet generated.
%
% Supports Nitrogen, Oxygen, and Isopropanol (IPA) with 100% identical I/O to FluidProperties.
%
% Calling Signatures:
%   props = FastFluidProps(fluid, 'From_u_rho', u, rho)
%   props = FastFluidProps(fluid, 'From_P_T', P, T)
%   val   = FastFluidProps(fluid, prop_out, in1_name, in1_val, in2_name, in2_val)
%   FastFluidProps('reset') % Optional utility to force-reload tables from disk

    if nargin == 1 && strcmpi(fluid, 'reset')
        clear resolveFromURho resolveFromPT loadFluidTables;
        out = true;
        return;
    end

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
        error('FastFluidProps: Unsupported fluid "%s"', fluid);
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
            % Generic property query: FastFluidProps(fluid, prop_out, in1_name, in1_val, in2_name, in2_val)
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

%% =========================================================================
%  Helper: Centralized Singleton Loader for Precomputed Interpolants
%% =========================================================================
function DB = loadFluidTables()
    persistent CachedDB
    if ~isempty(CachedDB)
        DB = CachedDB;
        return;
    end

    DB.URhoLoaded = false;
    DB.PTLoaded   = false;
    DB.FullURho   = false;

    try
        thisDir = fileparts(mfilename('fullpath'));
        matPaths = {
            fullfile(thisDir, 'FluidPropertyTables.mat'), ...
            fullfile(pwd, 'FluidPropertyTables.mat'), ...
            fullfile(pwd, 'Fluid Properties', 'FluidPropertyTables.mat'), ...
            fullfile(pwd, 'sandbox', 'experiments', 'Fluid Properties', 'FluidPropertyTables.mat'), ...
            fullfile(thisDir, '..', 'Fluid Properties', 'FluidPropertyTables.mat')
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

            % 1. Check for Full (u, rho) CoolProp Tables (from generateFluidPropertyTables.m)
            if isfield(Tbl, 'URho_n2') && isfield(Tbl, 'URho_ox')
                DB.URho_n2 = buildInterpolantSet({Tbl.URho_n2.rho_vec, Tbl.URho_n2.u_vec}, Tbl.URho_n2, {'P','T','h','cp','gamma','mu','k'});
                DB.URho_ox = buildInterpolantSet({Tbl.URho_ox.rho_vec, Tbl.URho_ox.u_vec}, Tbl.URho_ox, {'P','T','h','cp','gamma','mu','k'});
                DB.URhoLoaded = true;
                DB.FullURho   = true;
            elseif isfield(Tbl, 'rho_vec_n2')
                % Legacy partial table fallback
                DB.URho_n2.P     = griddedInterpolant({Tbl.rho_vec_n2, Tbl.u_vec_n2}, Tbl.P_grid_n2, 'linear', 'nearest');
                DB.URho_n2.T     = griddedInterpolant({Tbl.rho_vec_n2, Tbl.u_vec_n2}, Tbl.T_grid_n2, 'linear', 'nearest');
                DB.URho_n2.gamma = griddedInterpolant({Tbl.rho_vec_n2, Tbl.u_vec_n2}, Tbl.GAMMA_grid_n2, 'linear', 'nearest');
                
                DB.URho_ox.P     = griddedInterpolant({Tbl.rho_vec_ox, Tbl.u_vec_ox}, Tbl.P_grid_ox, 'linear', 'nearest');
                DB.URho_ox.T     = griddedInterpolant({Tbl.rho_vec_ox, Tbl.u_vec_ox}, Tbl.T_grid_ox, 'linear', 'nearest');
                DB.URho_ox.gamma = griddedInterpolant({Tbl.rho_vec_ox, Tbl.u_vec_ox}, Tbl.GAMMA_grid_ox, 'linear', 'nearest');
                DB.URhoLoaded = true;
                DB.FullURho   = false;
            end

            % 2. Check for Full (P, T) CoolProp Tables (interpolated on log10(P) and T)
            if isfield(Tbl, 'PT_n2') && isfield(Tbl, 'PT_ox')
                DB.PT_n2 = buildInterpolantSet({Tbl.PT_n2.logP_vec, Tbl.PT_n2.T_vec}, Tbl.PT_n2, {'rho','u','h','cp','gamma','mu','k'});
                DB.PT_ox = buildInterpolantSet({Tbl.PT_ox.logP_vec, Tbl.PT_ox.T_vec}, Tbl.PT_ox, {'rho','u','h','cp','gamma','mu','k'});
                DB.PTLoaded = true;
            end
        end
    catch
        DB.URhoLoaded = false;
        DB.PTLoaded   = false;
    end

    CachedDB = DB;
end

function S = buildInterpolantSet(gridVecs, srcStruct, fieldList)
    S = struct();
    for k = 1:length(fieldList)
        fn = fieldList{k};
        S.(fn) = griddedInterpolant(gridVecs, srcStruct.(fn), 'linear', 'nearest');
    end
end

%% =========================================================================
%  Helper: Resolve All Properties from (u, rho)
%% =========================================================================
function props = resolveFromURho(fluidName, isCoolProp, u, rho)
    rho = max(rho, 1e-4); % Safeguard

    % Array broadcasting if inputs are mismatched sizes
    if isscalar(u) && ~isscalar(rho), u = u .* ones(size(rho)); end
    if isscalar(rho) && ~isscalar(u), rho = rho .* ones(size(u)); end

    if isCoolProp
        DB = loadFluidTables();
        if DB.URhoLoaded
            if strcmp(fluidName, 'Nitrogen')
                G = DB.URho_n2;
                R_spec = 296.8;
                mu_exp = 0.72; mu_ref = 1.78e-5;
            else
                G = DB.URho_ox;
                R_spec = 259.8;
                mu_exp = 0.75; mu_ref = 2.0e-5;
            end

            P     = max(G.P(rho, u), 1000);
            T     = max(G.T(rho, u), 60);
            gamma = G.gamma(rho, u);

            if DB.FullURho
                h  = G.h(rho, u);
                cp = G.cp(rho, u);
                mu = G.mu(rho, u);
                k  = G.k(rho, u);
            else
                h  = u + P ./ rho;
                mu = mu_ref .* ((T ./ 300).^mu_exp);
                k  = 0.026 .* ((T ./ 300).^0.85);
                cp = gamma .* R_spec ./ max(gamma - 1, 0.05);
            end
        else
            % Direct CoolProp fallback if .mat file is missing
            out_sz = size(u);
            P = zeros(out_sz); T = zeros(out_sz); h = zeros(out_sz);
            cp = zeros(out_sz); cv = zeros(out_sz); gamma = zeros(out_sz);
            mu = zeros(out_sz); k = zeros(out_sz);

            for i = 1:numel(u)
                P(i) = double(py.CoolProp.CoolProp.PropsSI('P', 'Umass', u(i), 'Dmass', rho(i), fluidName));
                T(i) = double(py.CoolProp.CoolProp.PropsSI('T', 'Umass', u(i), 'Dmass', rho(i), fluidName));
                h(i) = double(py.CoolProp.CoolProp.PropsSI('Hmass', 'Umass', u(i), 'Dmass', rho(i), fluidName));
                
                cp(i) = double(py.CoolProp.CoolProp.PropsSI('Cpmass', 'P', P(i), 'T', T(i), fluidName));
                cv(i) = double(py.CoolProp.CoolProp.PropsSI('Cvmass', 'P', P(i), 'T', T(i), fluidName));
                gamma(i) = max(1.1, min(1.67, cp(i) ./ max(cv(i), 1e-3)));
                
                mu(i) = double(py.CoolProp.CoolProp.PropsSI('V', 'P', P(i), 'T', T(i), fluidName));
                k(i)  = double(py.CoolProp.CoolProp.PropsSI('L', 'P', P(i), 'T', T(i), fluidName));
            end
        end
    else
        % IPA (Incompressible Liquid Evaluation - Natively Vectorized)
        cp_ipa = 2600; % J/(kg-K) approximate average
        T = 293.15 + u ./ cp_ipa;
        T = max(190, min(400, T)); % Bounds
        dT = T - 293.15;
        
        cp = 2570.0 + 5.2 .* dT;
        mu = 2.4e-3 .* exp(-0.025 .* max(-50, min(100, dT)));
        k  = 0.14;
        P = 101325 .* ones(size(T)); % Bulk thermodynamic pressure
        h = u + P ./ max(rho, 10);
        gamma = 1.15 .* ones(size(T));
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

%% =========================================================================
%  Helper: Resolve All Properties from (P, T)
%% =========================================================================
function props = resolveFromPT(fluidName, isCoolProp, P, T)
    P = max(P, 1000); % Positive pressure
    T = max(T, 50);   % Positive temperature

    % Array broadcasting if inputs are mismatched sizes
    if isscalar(P) && ~isscalar(T), P = P .* ones(size(T)); end
    if isscalar(T) && ~isscalar(P), T = T .* ones(size(P)); end

    if isCoolProp
        DB = loadFluidTables();
        if DB.PTLoaded
            if strcmp(fluidName, 'Nitrogen')
                G = DB.PT_n2;
            else
                G = DB.PT_ox;
            end
            logP = log10(P);
            rho   = G.rho(logP, T);
            u     = G.u(logP, T);
            h     = G.h(logP, T);
            cp    = G.cp(logP, T);
            gamma = G.gamma(logP, T);
            mu    = G.mu(logP, T);
            k     = G.k(logP, T);
        else
            % Direct CoolProp evaluation loop fallback
            out_sz = size(P);
            rho = zeros(out_sz); u = zeros(out_sz); h = zeros(out_sz);
            cp = zeros(out_sz); gamma = zeros(out_sz); mu = zeros(out_sz); k = zeros(out_sz);
            
            for i = 1:numel(P)
                rho(i)   = double(py.CoolProp.CoolProp.PropsSI('D', 'P', P(i), 'T', T(i), fluidName));
                u(i)     = double(py.CoolProp.CoolProp.PropsSI('Umass', 'P', P(i), 'T', T(i), fluidName));
                h(i)     = double(py.CoolProp.CoolProp.PropsSI('Hmass', 'P', P(i), 'T', T(i), fluidName));
                cp(i)    = double(py.CoolProp.CoolProp.PropsSI('Cpmass', 'P', P(i), 'T', T(i), fluidName));
                gamma(i) = double(py.CoolProp.CoolProp.PropsSI('isentropic_expansion_coefficient', 'P', P(i), 'T', T(i), fluidName));
                mu(i)    = double(py.CoolProp.CoolProp.PropsSI('V', 'P', P(i), 'T', T(i), fluidName));
                k(i)     = double(py.CoolProp.CoolProp.PropsSI('L', 'P', P(i), 'T', T(i), fluidName));
            end
        end
    else
        % IPA (Isopropanol - subcooled liquid analytical model - Natively Vectorized)
        dT  = T - 293.15;
        rho = 786.0 - 0.85 .* dT;
        cp  = 2570.0 + 5.2 .* dT;
        mu  = 2.4e-3 .* exp(-0.025 .* max(-50, min(100, dT)));
        k   = 0.14 .* ones(size(T));
        gamma = 1.15 .* ones(size(T));
        u   = cp .* dT;
        h   = u + P ./ rho;
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
