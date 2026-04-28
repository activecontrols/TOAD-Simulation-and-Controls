
function generate_IPA_rgp(outputFile, varargin)
%GENERATE_IPA_RGP Generate a CFX-style RGP file for liquid isopropyl alcohol.
%
% This generator is built around the user's IPA_properties(T,P) function.
% It writes the header/data organization to mirror a working CFX RGP sample:
%   $$$$HEADER
%   $$$<component>
%   1
%   $$PARAM
%   25
%   ...
%   $$$$DATA
%   $$$<component>
%   1
%   $$PARAM
%   25
%   ...
%   $$SUPER_TABLE
%   9
%   $TABLE_1 ... $TABLE_9
%
% IMPORTANT:
% 1) This script is intended for a subcooled-liquid IPA manifold case.
% 2) The supplied IPA_properties.m notes that it is currently reliable in
%    the liquid state and not intended as a full vapor-dome EOS.
% 3) For that reason, this script writes the 9 superheat tables and also a
%    standalone $$SAT_TABLE data block using liquid-side saturation values.
% 4) For robustness, the default T range extends above the user's operating
%    range so the boiling curve over 250-500 psia falls inside the RGP's
%    temperature envelope.
%
% Usage examples:
%   generate_IPA_rgp()
%   generate_IPA_rgp('IPA_liquid.rgp')
%   generate_IPA_rgp('IPA_liquid.rgp', 'TminK', 250, 'TmaxK', 500, ...
%                    'PminPsia', 250, 'PmaxPsia', 500, ...
%                    'nT', 251, 'nP', 51, 'componentKey', 'IPAL')
%
% Inputs are optional name/value pairs:
%   'TminK'        default 250
%   'TmaxK'        default 500   (recommended; operating window may be narrower)
%   'PminPsia'     default 250
%   'PmaxPsia'     default 500
%   'nT'           default 251
%   'nP'           default 51
%   'componentKey' default 'IPAL'
%   'materialName' default 'IPA liquid'
%   'databaseName' default 'User MATLAB generator'
%   'description'  default 'Generated from IPA_properties.m'
%   'bulkModulusPa' default 1.07e9
%
% Table mapping used here:
%   TABLE_1 = enthalpy              [J/kg]
%   TABLE_2 = speed of sound        [m/s]
%   TABLE_3 = specific volume       [m^3/kg]
%   TABLE_4 = cv                    [J/(kg*K)]
%   TABLE_5 = cp                    [J/(kg*K)]
%   TABLE_6 = dP/dv |_T             [Pa / (m^3/kg)]
%   TABLE_7 = entropy               [J/(kg*K)]
%   TABLE_8 = dynamic viscosity     [Pa*s]
%   TABLE_9 = thermal conductivity  [W/(m*K)]
%
% Assumptions used to complete tables not returned directly by IPA_properties:
%   - h(T,P) from numerical integration of cp(T,P) at each pressure column
%   - s(T,P) from numerical integration of cp(T,P)/T with a small liquid
%     pressure correction using alpha_p/rho
%   - w(T,P) from w = sqrt(K/rho) with constant bulk modulus K
%   - cv(T,P) from cp - T*alpha_p^2*K/rho
%   - dP/dv|T from constant-K liquid approximation: dP/dv|T = -K*rho
%
% This function writes fixed-width scientific notation with 5 values per line.
%
% Author: OpenAI helper
% Date: 2026-04-15

    if nargin < 1 || isempty(outputFile)
        outputFile = 'IPA_liquid_250to500K_250to500psia.rgp';
    end

    opts = parseInputs(varargin{:});

    % ---- Unit conversions / constants ----
    psi_to_Pa = 6894.757293168361;
    Pmin = opts.PminPsia * psi_to_Pa;
    Pmax = opts.PmaxPsia * psi_to_Pa;

    Tvec = linspace(opts.TminK, opts.TmaxK, opts.nT);
    Pvec = linspace(Pmin, Pmax, opts.nP);

    % Use the supplied function's lower bound as the triple-point proxy.
    % This is not a separate NIST triple-point derivation.
    T_triple = 185.2;
    P_triple = 0.0;

    % Critical constants for IPA
    T_critical = 508.3;       % K
    P_critical = 47.6e5;      % Pa  (47.6 bar)

    molecularWeight = 60.0950e-3; % kg/mol
    R_universal = 8.314462618;    % J/(mol*K)
    R_specific = R_universal / molecularWeight; % J/(kg*K)

    % ---- Evaluate properties on the T-P grid ----
    nT = numel(Tvec);
    nP = numel(Pvec);

    rho = zeros(nT, nP);
    cp  = zeros(nT, nP);
    mu  = zeros(nT, nP);
    ktc = zeros(nT, nP);
    stateGrid = strings(nT, nP);

    % Boiling curve versus pressure
    T_sat = zeros(1, nP);

    for j = 1:nP
        % boiling_temp is independent of input T inside the user's function,
        % so use a benign liquid-state temperature.
        [~, ~, ~, ~, boilT, ~] = IPA_properties(min(max(300, opts.TminK), opts.TmaxK), Pvec(j));
        T_sat(j) = boilT;

        for i = 1:nT
            [mu(i,j), ktc(i,j), rho(i,j), cp_kJkgK, ~, stateGrid(i,j)] = IPA_properties(Tvec(i), Pvec(j));
            cp(i,j) = cp_kJkgK * 1000.0; % convert kJ/kg-K -> J/kg-K
        end
    end

    if any(stateGrid(:) == "boiling")
        warning(['Some T-P grid points were labeled "boiling" by IPA_properties. ', ...
                 'For a pure subcooled-liquid RGP, reduce TmaxK or increase pressure range.']);
    end

    % ---- Derived properties needed by the 9-table format ----
    % alpha_p = -(1/rho) * (drho/dT)|P
    alpha_p = zeros(size(rho));
    for j = 1:nP
        drho_dT = gradient(rho(:,j), Tvec);
        alpha_p(:,j) = -drho_dT ./ rho(:,j);
    end

    Kbulk = opts.bulkModulusPa;      % Pa
    w  = sqrt(Kbulk ./ rho);         % m/s
    v  = 1.0 ./ rho;                 % m^3/kg
    cv = cp - (Tvec(:) .* (alpha_p.^2) .* Kbulk) ./ rho; % J/kg-K
    cv = min(cv, cp);
    cv = max(cv, 1.0);               % keep positive and sane

    dPdv_T = -Kbulk .* rho;          % Pa / (m^3/kg)

    % Enthalpy and entropy references:
    % set h = 0 and s = 0 at Tmin for each pressure column.
    h = zeros(nT, nP);
    s = zeros(nT, nP);
    Pref = Pvec(1);

    for j = 1:nP
        h(:,j) = cumtrapz(Tvec, cp(:,j));

        s_thermal = cumtrapz(Tvec, cp(:,j) ./ Tvec(:));
        % small liquid pressure correction relative to first pressure column
        s(:,j) = s_thermal - ((alpha_p(:,j) ./ rho(:,j)) .* (Pvec(j) - Pref));
    end

    % ---- Saturation-curve side arrays required at the end of each table ----
    % Evaluate slightly on the liquid side of the saturation curve.
    nSat = nP;
    T_sat_eval = max(opts.TminK, min(T_sat - 0.5, opts.TmaxK));

    rho_sat = zeros(1, nSat);
    cp_sat  = zeros(1, nSat);
    mu_sat  = zeros(1, nSat);
    k_sat   = zeros(1, nSat);

    for j = 1:nSat
        [mu_sat(j), k_sat(j), rho_sat(j), cp_kJkgK, ~, satState] = IPA_properties(T_sat_eval(j), Pvec(j));
        cp_sat(j) = cp_kJkgK * 1000.0;

        if satState == "boiling"
            warning('Saturation-side evaluation at j=%d is at/above boiling according to IPA_properties.', j);
        end
    end

    % alpha_p along saturation-side evaluation
    drhoSat_dT = gradient(rho_sat, T_sat_eval);
    alpha_sat = -drhoSat_dT ./ rho_sat;

    w_sat  = sqrt(Kbulk ./ rho_sat);
    v_sat  = 1.0 ./ rho_sat;
    cv_sat = cp_sat - (T_sat_eval .* (alpha_sat.^2) .* Kbulk) ./ rho_sat;
    cv_sat = min(cv_sat, cp_sat);
    cv_sat = max(cv_sat, 1.0);
    dPdv_sat = -Kbulk .* rho_sat;

    % simple saturation-side h and s built by integrating along T_sat_eval
    h_sat = cumtrapz(T_sat_eval, cp_sat);
    s_sat = cumtrapz(T_sat_eval, cp_sat ./ T_sat_eval) - ((alpha_sat ./ rho_sat) .* (Pvec - Pref));

    % ---- Validation ----
    assert(all(isfinite(rho(:))) && all(rho(:) > 0), 'Non-finite or non-positive density detected.');
    assert(all(isfinite(cp(:)))  && all(cp(:)  > 0), 'Non-finite or non-positive cp detected.');
    assert(all(isfinite(mu(:)))  && all(mu(:)  > 0), 'Non-finite or non-positive viscosity detected.');
    assert(all(isfinite(ktc(:))) && all(ktc(:) > 0), 'Non-finite or non-positive conductivity detected.');
    assert(all(diff(Tvec) > 0), 'Temperature grid must be strictly increasing.');
    assert(all(diff(Pvec) > 0), 'Pressure grid must be strictly increasing.');

    % ---- Write the RGP file ----
    fid = fopen(outputFile, 'w');
    if fid < 0
        error('Could not open output file: %s', outputFile);
    end

    cleaner = onCleanup(@() fclose(fid));

    % Header block
    fprintf(fid, '$$$$HEADER\n');
    writeParamBlock(fid, opts, R_specific, Pvec, Tvec, T_sat, ...
                    T_critical, P_critical, T_triple, P_triple, nT, nP, nSat);

    % Data block
    fprintf(fid, '$$$$DATA\n');
    writeParamBlock(fid, opts, R_specific, Pvec, Tvec, T_sat, ...
                    T_critical, P_critical, T_triple, P_triple, nT, nP, nSat);

    fprintf(fid, '$$SUPER_TABLE\n');
    fprintf(fid, '9\n');

    writeOneTable(fid, 1, Tvec, Pvec, h,      T_sat, h_sat);
    writeOneTable(fid, 2, Tvec, Pvec, w,      T_sat, w_sat);
    writeOneTable(fid, 3, Tvec, Pvec, v,      T_sat, v_sat);
    writeOneTable(fid, 4, Tvec, Pvec, cv,     T_sat, cv_sat);
    writeOneTable(fid, 5, Tvec, Pvec, cp,     T_sat, cp_sat);
    writeOneTable(fid, 6, Tvec, Pvec, dPdv_T, T_sat, dPdv_sat);
    writeOneTable(fid, 7, Tvec, Pvec, s,      T_sat, s_sat);
    writeOneTable(fid, 8, Tvec, Pvec, mu,     T_sat, mu_sat);
    writeOneTable(fid, 9, Tvec, Pvec, ktc,    T_sat, k_sat);

    % Standalone saturation table block. CFX docs require valid saturation
    % curves for clipping, and the historical RGP generator also writes a
    % $$SAT_TABLE section after the super tables.
    writeSatTable(fid, Pvec, T_sat, h_sat, cp_sat, rho_sat, dPdv_sat, s_sat, cv_sat, w_sat, mu_sat, k_sat);

    fprintf('Wrote RGP file: %s\n', outputFile);
    fprintf('Temperature range written: %.3f K to %.3f K\n', Tvec(1), Tvec(end));
    fprintf('Pressure range written: %.3f psia to %.3f psia\n', opts.PminPsia, opts.PmaxPsia);
    fprintf('Boiling curve over pressure range: %.3f K to %.3f K\n', min(T_sat), max(T_sat));
end

function opts = parseInputs(varargin)
    p = inputParser;
    addParameter(p, 'TminK', 250, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'TmaxK', 500, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'PminPsia', 250, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'PmaxPsia', 500, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'nT', 251, @(x) isnumeric(x) && isscalar(x) && x >= 2);
    addParameter(p, 'nP', 51, @(x) isnumeric(x) && isscalar(x) && x >= 2);
    addParameter(p, 'componentKey', 'IPAL', @(x) ischar(x) || isstring(x));
    addParameter(p, 'materialName', 'IPA liquid', @(x) ischar(x) || isstring(x));
    addParameter(p, 'databaseName', 'User MATLAB generator', @(x) ischar(x) || isstring(x));
    addParameter(p, 'description', 'Generated from IPA_properties.m', @(x) ischar(x) || isstring(x));
    addParameter(p, 'bulkModulusPa', 1.07e9, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, varargin{:});
    opts = p.Results;

    if opts.TmaxK <= opts.TminK
        error('TmaxK must be greater than TminK.');
    end
    if opts.PmaxPsia <= opts.PminPsia
        error('PmaxPsia must be greater than PminPsia.');
    end

    opts.componentKey = char(opts.componentKey);
    opts.materialName = char(opts.materialName);
    opts.databaseName = char(opts.databaseName);
    opts.description  = char(opts.description);
end

function writeParamBlock(fid, opts, R_specific, Pvec, Tvec, T_sat, ...
                         T_critical, P_critical, T_triple, P_triple, nT, nP, nSat)
    fprintf(fid, '$$$%s\n', opts.componentKey);
    fprintf(fid, '1\n');               % ignored by CFX in known examples
    fprintf(fid, '$$PARAM\n');
    fprintf(fid, '25\n');

    fprintf(fid, 'DESCRIPTION\n');
    fprintf(fid, '%s\n', opts.description);

    fprintf(fid, 'NAME\n');
    fprintf(fid, '%s\n', opts.materialName);

    fprintf(fid, 'DATABASE\n');
    fprintf(fid, '%s\n', opts.databaseName);

    fprintf(fid, 'UNITS\n');
    fprintf(fid, '1\n');

    fprintf(fid, 'PMIN_SUPERHEAT\n');
    writeScalar(fid, Pvec(1));

    fprintf(fid, 'PMAX_SUPERHEAT\n');
    writeScalar(fid, Pvec(end));

    fprintf(fid, 'TMIN_SUPERHEAT\n');
    writeScalar(fid, Tvec(1));

    fprintf(fid, 'TMAX_SUPERHEAT\n');
    writeScalar(fid, Tvec(end));

    fprintf(fid, 'TMIN_SATURATION\n');
    writeScalar(fid, min(T_sat));

    fprintf(fid, 'TMAX_SATURATION\n');
    writeScalar(fid, max(T_sat));

    fprintf(fid, 'P_CRITICAL\n');
    writeScalar(fid, P_critical);

    fprintf(fid, 'P_TRIPLE\n');
    writeScalar(fid, P_triple);

    fprintf(fid, 'T_CRITICAL\n');
    writeScalar(fid, T_critical);

    fprintf(fid, 'T_TRIPLE\n');
    writeScalar(fid, T_triple);

    fprintf(fid, 'GAS_CONSTANT\n');
    writeScalar(fid, R_specific);

    for k = 1:9
        fprintf(fid, 'TABLE_%d\n', k);
        fprintf(fid, '%d %d\n', nT, nP);
    end

    fprintf(fid, 'SAT_TABLE\n');
    fprintf(fid, '%d %d %d\n', nSat, 4, 9);
end

function writeOneTable(fid, tableNo, Tvec, Pvec, A, T_sat, A_sat)
    fprintf(fid, '$TABLE_%d\n', tableNo);
    fprintf(fid, '%d %d\n', numel(Tvec), numel(Pvec));

    writeVector(fid, Tvec);
    writeVector(fid, Pvec);

    % Write in j-major ordering to mimic the known Fortran writer:
    % ((A(i,j), i=1,nt), j=1,np)
    writeMatrixJMajor(fid, A);

    writeVector(fid, T_sat);
    writeVector(fid, A_sat);
end

function writeSatTable(fid, Pvec, T_sat, h_sat, cp_sat, rho_sat, dPdv_sat, s_sat, cv_sat, w_sat, mu_sat, k_sat)
    % Historical TASCflow RGP saturation table format:
    %   psat, tsat, blank1, blank2,
    %   liquid group (H,Cp,Rho,dPdrho|T,S,Cv,w,eta,tcx),
    %   vapor  group (same 9 arrays).
    %
    % This IPA property source is liquid-focused, so the vapor block is
    % duplicated from the liquid-side saturation values as a parser-safe
    % fallback for single-phase liquid applications.

    fprintf(fid, '$$SAT_TABLE\n');
    fprintf(fid, '%d %d %d\n', numel(Pvec), 4, 9);

    blank = zeros(size(Pvec));
    dPdrho_sat = -dPdv_sat .* (rho_sat.^2);

    writeVector(fid, Pvec);
    writeVector(fid, T_sat);
    writeVector(fid, blank);
    writeVector(fid, blank);

    % Liquid group
    writeVector(fid, h_sat);
    writeVector(fid, cp_sat);
    writeVector(fid, rho_sat);
    writeVector(fid, dPdrho_sat);
    writeVector(fid, s_sat);
    writeVector(fid, cv_sat);
    writeVector(fid, w_sat);
    writeVector(fid, mu_sat);
    writeVector(fid, k_sat);

    % Vapor group fallback
    writeVector(fid, h_sat);
    writeVector(fid, cp_sat);
    writeVector(fid, rho_sat);
    writeVector(fid, dPdrho_sat);
    writeVector(fid, s_sat);
    writeVector(fid, cv_sat);
    writeVector(fid, w_sat);
    writeVector(fid, mu_sat);
    writeVector(fid, k_sat);
end

function writeScalar(fid, x)
    fprintf(fid, '%.10E\n', x);
end

function writeVector(fid, x)
    x = x(:).';
    n = numel(x);
    perLine = 5;
    idx = 1;
    while idx <= n
        j2 = min(idx + perLine - 1, n);
        for j = idx:j2
            fprintf(fid, '%17.7E', x(j));
        end
        fprintf(fid, '\n');
        idx = j2 + 1;
    end
end

function writeMatrixJMajor(fid, A)
    [nT, nP] = size(A);
    flat = zeros(1, nT*nP);
    k = 1;
    for j = 1:nP
        for i = 1:nT
            flat(k) = A(i,j);
            k = k + 1;
        end
    end
    writeVector(fid, flat);
end
