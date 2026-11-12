function [tbl, csvPath, trajStem] = ExportTrajectoryCSV(InitialGuess, constants6DoF, varargin)
%% Export 20-Column Simulation CSV and Generate TV-LQI Tracking Gains
    p = inputParser;
    addParameter(p, 'MissionName', "Waypoint");
    addParameter(p, 'Version', 1);
    addParameter(p, 'SaveDir', "");
    addParameter(p, 'Filename', "");
    addParameter(p, 'RunSaveGains', true);
    parse(p, varargin{:});
    opts = p.Results;

    % 1. Resolve vehicle name and output file stem
    isTOAD = isfield(constants6DoF, 'Vehicle') && constants6DoF.Vehicle == 1;
    if isTOAD
        vehName = "TOAD";
    else
        vehName = "ASTRA";
    end

    if strlength(string(opts.Filename)) > 0
        trajStem = string(opts.Filename);
    else
        trajStem = string(sprintf('%s_%s_v%03d', vehName, string(opts.MissionName), round(opts.Version)));
    end

    if strlength(string(opts.SaveDir)) > 0
        trajDir = char(opts.SaveDir);
    else
        trajDir = fullfile(pwd, 'Guidance', 'Trajectories');
    end
    if ~exist(trajDir, 'dir')
        mkdir(trajDir);
    end
    csvPath = fullfile(trajDir, char(trajStem + ".csv"));

    % 2. Build 20-column table matching LoadTOADSim.m and SaveGains.m
    ts = InitialGuess.Time(:);
    Xr = InitialGuess.X;
    Ur = InitialGuess.U;
    Nn = numel(ts);
    Nc = size(Ur, 2);

    if Nc == Nn - 1
        U_pad = [Ur, Ur(:, end)];
    elseif Nc == Nn
        U_pad = Ur;
    else
        error('ExportTrajectoryCSV:DimMismatch', ...
            'Control matrix columns (%d) must be Nn-1 (%d) or Nn (%d).', Nc, Nn - 1, Nn);
    end

    tbl = table( ...
        ts, ...
        Xr(1, :)', Xr(2, :)', Xr(3, :)', Xr(4, :)', ...
        Xr(5, :)', Xr(6, :)', Xr(7, :)', ...
        Xr(8, :)', Xr(9, :)', Xr(10, :)', ...
        Xr(11, :)', Xr(12, :)', Xr(13, :)', ...
        Xr(14, :)', Xr(15, :)', ...
        U_pad(1, :)', U_pad(2, :)', U_pad(3, :)', U_pad(4, :)', ...
        'VariableNames', { ...
            'Time', ...
            'QuatW', 'QuatX', 'QuatY', 'QuatZ', ...
            'PosX', 'PosY', 'PosZ', ...
            'VelX', 'VelY', 'VelZ', ...
            'AngRateX', 'AngRateY', 'AngRateZ', ...
            'MassLox', 'MassFuel', ...
            'GimbalTheta', 'GimbalPhi', 'ThrustMag', 'RollCmd' ...
        } ...
    );

    writetable(tbl, csvPath);
    fprintf('  [ExportTrajectoryCSV] Wrote %d samples to: %s\n', Nn, csvPath);

    % 3. Populate default TV-LQI / LESO tuning if not yet present in constants6DoF
    if opts.RunSaveGains
        constants6DoF = ensureControllerTuning(constants6DoF, isTOAD);
        gainsDir = fullfile(trajDir, 'Gains');
        if ~exist(gainsDir, 'dir')
            mkdir(gainsDir);
        end

        % SaveGains reads "Guidance\Trajectories\<trajStem>.csv" relative to pwd
        defaultTrajFile = fullfile(pwd, 'Guidance', 'Trajectories', char(trajStem + ".csv"));
        if strcmpi(csvPath, defaultTrajFile)
            fprintf('  [ExportTrajectoryCSV] Running SaveGains for %s...\n', trajStem);
            SaveGains(trajStem, constants6DoF);
        else
            % Custom SaveDir (e.g., sandbox testing): run GainGenerator directly and write to gainsDir
            fprintf('  [ExportTrajectoryCSV] Running GainGenerator for %s in %s...\n', trajStem, gainsDir);
            [K_trans_List, K_rot_List, LA, LT] = GainGenerator(Xr, Ur, ts', constants6DoF);
            writematrix(K_trans_List, fullfile(gainsDir, char("K_trans_" + trajStem + ".csv")));
            writematrix(K_rot_List,   fullfile(gainsDir, char("K_rot_"   + trajStem + ".csv")));
            writematrix(LA,           fullfile(gainsDir, char("LA_"      + trajStem + ".csv")));
            writematrix(LT,           fullfile(gainsDir, char("LT_"      + trajStem + ".csv")));
        end
        fprintf('  [ExportTrajectoryCSV] TV-LQI & LESO gains saved to: %s\n', gainsDir);
    end
end

function c = ensureControllerTuning(c, isTOAD)
%% Populate Default LQR & LESO Tuning Matrices Matching LoadTOADSim.m
    if ~isfield(c, 'Q_trans')
        if isTOAD
            max_x_trans = 1.3 * ones(1, 6);
        else
            max_x_trans = 1.25 * ones(1, 6);
        end
        c.Q_trans = diag(1 ./ max_x_trans.^2);
    end
    if ~isfield(c, 'R_trans')
        max_a_trans = 1.2;
        c.R_trans = eye(3) .* (1 / max_a_trans^2);
    end
    if ~isfield(c, 'OmegaThr')
        c.OmegaThr = 2.2;
    end
    if ~isfield(c, 'Q_rot') || ~isfield(c, 'R_rot')
        if isTOAD
            max_x_rot = [0.12, 0.12, 0.12, 0.22, 0.22, 0.21];
            c.R_rot = diag([35, 35, 1 / 4^2]);
        else
            max_x_rot = [0.15, 0.15, 0.12, 0.75, 0.75, 0.7];
            c.R_rot = diag([70, 70, 1 / 0.10^2]);
        end
        c.Q_rot = diag(1 ./ max_x_rot.^2);
    end
    if ~isfield(c, 'OmegaAtt')
        c.OmegaAtt = 2.2;
    end
end
