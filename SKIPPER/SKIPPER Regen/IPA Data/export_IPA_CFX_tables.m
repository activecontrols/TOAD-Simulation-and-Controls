function export_IPA_CFX_tables(varargin)
%EXPORT_IPA_CFX_TABLES Export IPA property tables from IPA_properties.m for CFX CEL use.
%
% This script samples the user's IPA_properties(T,P) model on a structured
% temperature-pressure grid and writes CSV files that can be used to create
% table/interpolation functions in ANSYS CFX.
%
% Default range is chosen to stay inside the strongest native density/cp
% table range of IPA_properties.m (270-320 K) and inside the requested
% pressure band (250-500 psia).
%
% Outputs written to outputDir:
%   IPA_rho_long.csv      - Long-form [T_K, P_Pa, rho_kg_m3]
%   IPA_mu_long.csv       - Long-form [T_K, P_Pa, mu_Pa_s]
%   IPA_k_long.csv        - Long-form [T_K, P_Pa, k_W_m_K]
%   IPA_cp_long.csv       - Long-form [T_K, P_Pa, cp_J_kg_K]
%   IPA_boiling_long.csv  - Long-form [T_K, P_Pa, boiling_T_K]
%   IPA_state_long.csv    - Long-form [T_K, P_Pa, state, is_boiling]
%
%   IPA_rho_table.csv     - Matrix style, first row = pressure [Pa], first
%                           column = temperature [K]
%   IPA_mu_table.csv
%   IPA_k_table.csv
%   IPA_cp_table.csv
%
%   CFX_CEL_expressions.txt - Suggested CEL expressions to paste into CFX
%   IPA_export_summary.txt  - Grid and validity summary
%
% Example:
%   export_IPA_CFX_tables();
%
%   export_IPA_CFX_tables('outputDir','C:/temp/ipa_cfx', ...
%       'TminK',270,'TmaxK',320,'TstepK',2, ...
%       'PminPsia',250,'PmaxPsia',500,'nP',26);
%
% Notes:
% - IPA_properties.m expects T [K] and P [Pa].
% - IPA_properties.m returns cp in kJ/kg-K, so this script converts cp to
%   J/kg-K before export because that is what CFX expects for SI units.
% - This exporter is intended for liquid-property use in CFX CEL; it is not
%   a replacement for a full real-gas/two-phase property model.

    opts = parseInputs(varargin{:});
    ensurePropertyFunctionAvailable();
    if ~exist(opts.outputDir, 'dir')
        mkdir(opts.outputDir);
    end

    PvecPa = buildPressureVector(opts);
    TvecK  = buildTemperatureVector(opts);

    nT = numel(TvecK);
    nP = numel(PvecPa);

    rho = nan(nT, nP);
    mu  = nan(nT, nP);
    k   = nan(nT, nP);
    cp  = nan(nT, nP);
    boilT = nan(nT, nP);
    isBoiling = false(nT, nP);
    stateText = strings(nT, nP);

    for i = 1:nT
        Ti = TvecK(i);
        for j = 1:nP
            Pj = PvecPa(j);
            [mu_ij, k_ij, rho_ij, cp_ij, boil_ij, state_ij] = IPA_properties(Ti, Pj);
            mu(i,j) = mu_ij;
            k(i,j) = k_ij;
            rho(i,j) = rho_ij;
            cp(i,j) = cp_ij * 1000.0; % convert kJ/kg-K -> J/kg-K
            boilT(i,j) = boil_ij;
            stateText(i,j) = string(state_ij);
            isBoiling(i,j) = strcmpi(strtrim(char(string(state_ij))), 'boiling');
        end
    end

    validateGrid(rho, mu, k, cp, boilT, TvecK, PvecPa, isBoiling, opts);

    writeLongCsv(fullfile(opts.outputDir, 'IPA_rho_long.csv'), ...
        {'T_K','P_Pa','rho_kg_m3'}, TvecK, PvecPa, rho);
    writeLongCsv(fullfile(opts.outputDir, 'IPA_mu_long.csv'), ...
        {'T_K','P_Pa','mu_Pa_s'}, TvecK, PvecPa, mu);
    writeLongCsv(fullfile(opts.outputDir, 'IPA_k_long.csv'), ...
        {'T_K','P_Pa','k_W_m_K'}, TvecK, PvecPa, k);
    writeLongCsv(fullfile(opts.outputDir, 'IPA_cp_long.csv'), ...
        {'T_K','P_Pa','cp_J_kg_K'}, TvecK, PvecPa, cp);
    writeLongCsv(fullfile(opts.outputDir, 'IPA_boiling_long.csv'), ...
        {'T_K','P_Pa','boiling_T_K'}, TvecK, PvecPa, boilT);
    writeStateCsv(fullfile(opts.outputDir, 'IPA_state_long.csv'), TvecK, PvecPa, stateText, isBoiling);

    writeMatrixCsv(fullfile(opts.outputDir, 'IPA_rho_table.csv'), 'rho_kg_m3', TvecK, PvecPa, rho);
    writeMatrixCsv(fullfile(opts.outputDir, 'IPA_mu_table.csv'), 'mu_Pa_s', TvecK, PvecPa, mu);
    writeMatrixCsv(fullfile(opts.outputDir, 'IPA_k_table.csv'), 'k_W_m_K', TvecK, PvecPa, k);
    writeMatrixCsv(fullfile(opts.outputDir, 'IPA_cp_table.csv'), 'cp_J_kg_K', TvecK, PvecPa, cp);

    writeCELText(fullfile(opts.outputDir, 'CFX_CEL_expressions.txt'), opts);
    writeSummary(fullfile(opts.outputDir, 'IPA_export_summary.txt'), opts, TvecK, PvecPa, ...
        rho, mu, k, cp, boilT, isBoiling);

    fprintf('\nExport complete. Files written to:\n  %s\n', opts.outputDir);
    fprintf('Suggested first use in CFX:\n');
    fprintf('  1) Create/import table or interpolation functions from the CSV files.\n');
    fprintf('  2) Use the CEL expressions listed in CFX_CEL_expressions.txt.\n');
    fprintf('  3) Start with density only, then add mu, k, and cp.\n\n');
end

function opts = parseInputs(varargin)
    p = inputParser;
    p.FunctionName = mfilename;

    addParameter(p, 'outputDir', fullfile(pwd, 'IPA_CFX_Tables'), @(x) ischar(x) || isstring(x));
    addParameter(p, 'TminK', 270, @(x) validateattributes(x, {'numeric'}, {'scalar','real'}));
    addParameter(p, 'TmaxK', 320, @(x) validateattributes(x, {'numeric'}, {'scalar','real'}));
    addParameter(p, 'TstepK', 2, @(x) validateattributes(x, {'numeric'}, {'scalar','real','positive'}));
    addParameter(p, 'PminPsia', 250, @(x) validateattributes(x, {'numeric'}, {'scalar','real','positive'}));
    addParameter(p, 'PmaxPsia', 500, @(x) validateattributes(x, {'numeric'}, {'scalar','real','positive'}));
    addParameter(p, 'nP', 26, @(x) validateattributes(x, {'numeric'}, {'scalar','integer','>=',2}));
    addParameter(p, 'allowBoilingPoints', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'allowNonfinite', false, @(x) islogical(x) || isnumeric(x));

    parse(p, varargin{:});
    opts = p.Results;
    opts.outputDir = char(opts.outputDir);

    if opts.TmaxK <= opts.TminK
        error('TmaxK must be greater than TminK.');
    end
    if opts.PmaxPsia <= opts.PminPsia
        error('PmaxPsia must be greater than PminPsia.');
    end
end

function ensurePropertyFunctionAvailable()
    if exist('IPA_properties', 'file') ~= 2
        error(['Could not find IPA_properties.m on the MATLAB path. ', ...
               'Put IPA_properties.m in the current folder or add its folder to the path.']);
    end
end

function PvecPa = buildPressureVector(opts)
    psia_to_pa = 6894.757293168;
    PvecPa = linspace(opts.PminPsia * psia_to_pa, opts.PmaxPsia * psia_to_pa, opts.nP);
end

function TvecK = buildTemperatureVector(opts)
    TvecK = opts.TminK : opts.TstepK : opts.TmaxK;
    if abs(TvecK(end) - opts.TmaxK) > 1e-12
        TvecK = [TvecK, opts.TmaxK]; %#ok<AGROW>
    end
end

function validateGrid(rho, mu, k, cp, boilT, TvecK, PvecPa, isBoiling, opts)
    allArrays = {rho, mu, k, cp, boilT};
    names = {'rho','mu','k','cp','boiling temperature'};

    for idx = 1:numel(allArrays)
        A = allArrays{idx};
        badFinite = ~isfinite(A);
        if any(badFinite(:)) && ~opts.allowNonfinite
            [r,c] = find(badFinite, 1, 'first');
            error(['Non-finite %s detected at T = %.6g K, P = %.6g Pa. ', ...
                   'Narrow the range or inspect IPA_properties.m.'], ...
                   names{idx}, TvecK(r), PvecPa(c));
        end
    end

    if any(rho(:) <= 0)
        [r,c] = find(rho <= 0, 1, 'first');
        error('Non-positive density detected at T = %.6g K, P = %.6g Pa.', TvecK(r), PvecPa(c));
    end
    if any(mu(:) <= 0)
        [r,c] = find(mu <= 0, 1, 'first');
        error('Non-positive viscosity detected at T = %.6g K, P = %.6g Pa.', TvecK(r), PvecPa(c));
    end
    if any(k(:) <= 0)
        [r,c] = find(k <= 0, 1, 'first');
        error('Non-positive thermal conductivity detected at T = %.6g K, P = %.6g Pa.', TvecK(r), PvecPa(c));
    end
    if any(cp(:) <= 0)
        [r,c] = find(cp <= 0, 1, 'first');
        error('Non-positive specific heat detected at T = %.6g K, P = %.6g Pa.', TvecK(r), PvecPa(c));
    end

    if any(isBoiling(:)) && ~opts.allowBoilingPoints
        [r,c] = find(isBoiling, 1, 'first');
        error(['A sampled point is at/above the IPA boiling line according to IPA_properties.m ', ...
               '(T = %.6g K, P = %.6g Pa, T_boil = %.6g K). ', ...
               'Reduce TmaxK, increase pressure range, or set allowBoilingPoints=true if intentional.'], ...
               TvecK(r), PvecPa(c), boilT(r,c));
    end
end

function writeLongCsv(filename, header, TvecK, PvecPa, values)
    fid = fopen(filename, 'w');
    assert(fid >= 0, 'Could not open %s for writing.', filename);
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, '%s,%s,%s\n', header{1}, header{2}, header{3});
    for i = 1:numel(TvecK)
        for j = 1:numel(PvecPa)
            fprintf(fid, '%.15g,%.15g,%.15g\n', TvecK(i), PvecPa(j), values(i,j));
        end
    end
end

function writeStateCsv(filename, TvecK, PvecPa, stateText, isBoiling)
    fid = fopen(filename, 'w');
    assert(fid >= 0, 'Could not open %s for writing.', filename);
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'T_K,P_Pa,state,is_boiling\n');
    for i = 1:numel(TvecK)
        for j = 1:numel(PvecPa)
            fprintf(fid, '%.15g,%.15g,%s,%d\n', TvecK(i), PvecPa(j), char(stateText(i,j)), isBoiling(i,j));
        end
    end
end

function writeMatrixCsv(filename, valueLabel, TvecK, PvecPa, values)
    fid = fopen(filename, 'w');
    assert(fid >= 0, 'Could not open %s for writing.', filename);
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'T_K_vs_P_Pa__%s', valueLabel);
    for j = 1:numel(PvecPa)
        fprintf(fid, ',%.15g', PvecPa(j));
    end
    fprintf(fid, '\n');

    for i = 1:numel(TvecK)
        fprintf(fid, '%.15g', TvecK(i));
        for j = 1:numel(PvecPa)
            fprintf(fid, ',%.15g', values(i,j));
        end
        fprintf(fid, '\n');
    end
end

function writeCELText(filename, opts)
    psia_to_pa = 6894.757293168;
    PminPa = opts.PminPsia * psia_to_pa;
    PmaxPa = opts.PmaxPsia * psia_to_pa;

    fid = fopen(filename, 'w');
    assert(fid >= 0, 'Could not open %s for writing.', filename);
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'Suggested CFX CEL expressions for the exported IPA tables\n');
    fprintf(fid, '========================================================\n\n');
    fprintf(fid, '1) Create/import user functions from the CSV files in CFX.\n');
    fprintf(fid, '   Use argument units [K], [Pa] and result units matching each file.\n');
    fprintf(fid, '   Example function names below assume you create these functions:\n');
    fprintf(fid, '     rho_tab(T,P), mu_tab(T,P), k_tab(T,P), cp_tab(T,P)\n\n');
    fprintf(fid, '2) Define clipping expressions first:\n\n');
    fprintf(fid, 'Tclip = max(%.15g [K], min(%.15g [K], T))\n', opts.TminK, opts.TmaxK);
    fprintf(fid, 'Pclip = max(%.15g [Pa], min(%.15g [Pa], p))\n\n', PminPa, PmaxPa);
    fprintf(fid, '3) Then define material properties:\n\n');
    fprintf(fid, 'rho_IPA = rho_tab(Tclip, Pclip)\n');
    fprintf(fid, 'mu_IPA  = mu_tab(Tclip, Pclip)\n');
    fprintf(fid, 'k_IPA   = k_tab(Tclip, Pclip)\n');
    fprintf(fid, 'cp_IPA  = cp_tab(Tclip, Pclip)\n\n');
    fprintf(fid, '4) Assign in the material definition:\n');
    fprintf(fid, '   Density               = rho_IPA\n');
    fprintf(fid, '   Dynamic Viscosity     = mu_IPA\n');
    fprintf(fid, '   Thermal Conductivity  = k_IPA\n');
    fprintf(fid, '   Specific Heat Capacity= cp_IPA\n\n');
    fprintf(fid, 'Start by assigning rho only, solve, then add mu, then k, then cp.\n');
end

function writeSummary(filename, opts, TvecK, PvecPa, rho, mu, k, cp, boilT, isBoiling)
    psia_to_pa = 6894.757293168;
    fid = fopen(filename, 'w');
    assert(fid >= 0, 'Could not open %s for writing.', filename);
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'IPA CFX table export summary\n');
    fprintf(fid, '============================\n\n');
    fprintf(fid, 'Source property function: IPA_properties.m\n');
    fprintf(fid, 'Temperature range: %.6g to %.6g K\n', TvecK(1), TvecK(end));
    fprintf(fid, 'Pressure range   : %.6g to %.6g Pa\n', PvecPa(1), PvecPa(end));
    fprintf(fid, 'Pressure range   : %.6g to %.6g psia\n', PvecPa(1)/psia_to_pa, PvecPa(end)/psia_to_pa);
    fprintf(fid, 'Temperature step : %.6g K\n', opts.TstepK);
    fprintf(fid, 'Pressure points  : %d\n\n', numel(PvecPa));

    fprintf(fid, 'rho range [kg/m^3]      : %.6g to %.6g\n', min(rho(:)), max(rho(:)));
    fprintf(fid, 'mu range [Pa*s]         : %.6g to %.6g\n', min(mu(:)), max(mu(:)));
    fprintf(fid, 'k range [W/m-K]         : %.6g to %.6g\n', min(k(:)), max(k(:)));
    fprintf(fid, 'cp range [J/kg-K]       : %.6g to %.6g\n', min(cp(:)), max(cp(:)));
    fprintf(fid, 'boiling T range [K]     : %.6g to %.6g\n', min(boilT(:)), max(boilT(:)));
    fprintf(fid, 'boiling points sampled? : %s\n', string(any(isBoiling(:))));
end
