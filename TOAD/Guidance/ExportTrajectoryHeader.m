function headerPath = ExportTrajectoryHeader(trajNames, varargin)
%% Export MATLAB/CSV trajectories to C++ header
%
% Takes one or more trajectory CSV filenames, reads  20-column state, input, gains
%  generates C++ header file formatted with static const arrays and Eigen accessors\
%
% Optional Parameters:
%   'OutputFile'    - Output header filename or path (default: /cpp/ASTRA_Trajectory_Data.h)
%   'TrajDir'       - Trajectory CSV directory (default: <this_dir>/Trajectories)
%   'GainsDir'      - Gains CSV directory (default: <TrajDir>/Gains)
%   'DryMass'       - Vehicle dry mass in kg (default: 1.169 for ASTRA, 141.52 for TOAD)
%   'Gravity'       - Gravity
%   'Namespace'     - C++ namespace (default: 'AstraGuidance')

    p = inputParser;
    addRequired(p, 'trajNames');
    addParameter(p, 'OutputFile', '');
    addParameter(p, 'TrajDir', '');
    addParameter(p, 'GainsDir', '');
    addParameter(p, 'DryMass', 1.1690);
    addParameter(p, 'Gravity', 9.8015);
    addParameter(p, 'Namespace', 'AstraGuidance');
    addParameter(p, 'UseProgmem', false);
    parse(p, trajNames, varargin{:});
    opts = p.Results;

    % Normalize trajNames into a cell array of strings
    if ischar(trajNames) || isstring(trajNames)
        trajList = {char(trajNames)};
    elseif iscell(trajNames)
        trajList = cellfun(@char, trajNames, 'UniformOutput', false);
    else
        error('ExportTrajectoryHeader:InvalidInput', 'trajNames must be string, char, or cell array of strings.');
    end

    thisDir = fileparts(mfilename('fullpath'));
    if isempty(thisDir)
        thisDir = pwd;
    end

    % Resolve paths
    if isempty(opts.TrajDir)
        trajDir = fullfile(thisDir, 'Trajectories');
    else
        trajDir = opts.TrajDir;
    end

    if isempty(opts.GainsDir)
        gainsDir = fullfile(trajDir, 'Gains');
    else
        gainsDir = opts.GainsDir;
    end

    if isempty(opts.OutputFile)
        headerPath = fullfile(thisDir, '..', 'cpp', 'ASTRA_Trajectory_Data.h');
    else
        headerPath = opts.OutputFile;
        if ~isstring(headerPath) && ~ischar(headerPath)
            error('ExportTrajectoryHeader:InvalidOutputFile', 'OutputFile must be a string or path.');
        end
        if ~fileparts(headerPath)
            headerPath = fullfile(thisDir, '..', 'cpp', headerPath);
        end
    end

    headerDir = fileparts(headerPath);
    if ~exist(headerDir, 'dir') && ~isempty(headerDir)
        mkdir(headerDir);
    end

    numTrajs = numel(trajList);
    trajectories = cell(numTrajs, 1);

    
    fprintf('Target Header: %s\n', headerPath);
    fprintf('Trajectories to process: %d\n', numTrajs);


    for tIdx = 1:numTrajs
        rawName = trajList{tIdx};
        [~, stem, ~] = fileparts(rawName);

        % Resolve CSV file
        csvFile = fullfile(trajDir, [stem, '.csv']);
        if ~exist(csvFile, 'file')
            csvFile = fullfile(trajDir, rawName);
            if ~exist(csvFile, 'file')
                error('ExportTrajectoryHeader:FileNotFound', 'Trajectory CSV not found: %s', csvFile);
            end
        end

        fprintf('[%d/%d] Loading trajectory: %s\n', tIdx, numTrajs, stem);
        data = readmatrix(csvFile);

        % Validate columns: 1 time 15 state 4 input
        [N, numCols] = size(data);
        if numCols < 20
            error('ExportTrajectoryHeader:BadDimensions', 'Expected at least 20 columns in %s (got %d).', csvFile, numCols);
        end

        traj.stem = stem;
        traj.ident = regexprep(stem, '[^a-zA-Z0-9_]', '_');
        traj.numPoints = N;
        traj.time = data(:, 1);
        traj.states = data(:, 2:16);
        traj.inputs = data(:, 17:20);
        traj.dtNom = mean(diff(traj.time));

        % Load Gains
        traj.kTransFile = fullfile(gainsDir, ['K_trans_', stem, '.csv']);
        traj.kRotFile   = fullfile(gainsDir, ['K_rot_', stem, '.csv']);
        traj.laFile     = fullfile(gainsDir, ['LA_', stem, '.csv']);
        traj.ltFile     = fullfile(gainsDir, ['LT_', stem, '.csv']);

        traj.hasGains = exist(traj.kTransFile, 'file') && exist(traj.kRotFile, 'file') && exist(traj.laFile, 'file')     && exist(traj.ltFile, 'file');

        if traj.hasGains
            kTransMat = readmatrix(traj.kTransFile);
            kRotMat   = readmatrix(traj.kRotFile);
            laMat     = readmatrix(traj.laFile);
            ltMat     = readmatrix(traj.ltFile);

            % Verify rows match N
            if size(kTransMat, 1) ~= N || size(kRotMat, 1) ~= N || ...
               size(laMat, 1) ~= N     || size(ltMat, 1) ~= N
                warning('ExportTrajectoryHeader:GainRowMismatch', 'Gain CSV row count does not match trajectory row count (%d vs %d). Clamping/padding.', size(kTransMat, 1), N);
                kTransMat = padOrClampRows(kTransMat, N, 18);
                kRotMat   = padOrClampRows(kRotMat, N, 18);
                laMat     = padOrClampRows(laMat, N, 18);
                ltMat     = padOrClampRows(ltMat, N, 54);
            end

            traj.kTrans = kTransMat;
            traj.kRot   = kRotMat;
            traj.la     = laMat;
            traj.lt     = ltMat;
            fprintf('Found gain matrices: K_trans (3x6), K_rot (3x6), LA (6x3), LT (9x6)\n');
        else
            warning('ExportTrajectoryHeader:MissingGains', 'Gains not found for %s in %s. Emitting zeros for gain tables.', stem, gainsDir);
            traj.kTrans = zeros(N, 18);
            traj.kRot   = zeros(N, 18);
            traj.la     = zeros(N, 18);
            traj.lt     = zeros(N, 54);
        end

        trajectories{tIdx} = traj;
    end

    %% Write C++ Header File
    fid = fopen(headerPath, 'w');
    if fid == -1
        error('ExportTrajectoryHeader:FileWriteError', 'Cannot open %s for writing.', headerPath);
    end

    c = onCleanup(@() fclose(fid));

    progmemAttr = '';
    if isfield(opts, 'UseProgmem') && opts.UseProgmem
        progmemAttr = ' PROGMEM';
    end

    fprintf(fid, '// AUTOMATICALLY GENERATED FILE DO NOT EDIT MANUALLY\n');
    fprintf(fid, '// Generated by ExportTrajectoryHeader.m on %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));


    fprintf(fid, '#pragma once\n\n');
    fprintf(fid, '#include "matlab_funcs.h"\n');
    fprintf(fid, '#include <cmath>\n\n');

    fprintf(fid, 'namespace %s {\n\n', opts.Namespace);

    % Emit Trajectory IDs enum
    fprintf(fid, 'enum TrajectoryID {\n');
    for i = 1:numTrajs
        if i == 1
            fprintf(fid, '    TRAJ_%s = 0,\n', upper(trajectories{i}.ident));
        else
            fprintf(fid, '    TRAJ_%s,\n', upper(trajectories{i}.ident));
        end
    end
    fprintf(fid, '    TRAJ_COUNT = %d\n', numTrajs);
    fprintf(fid, '};\n\n');

    % Metadata struct declaration
    fprintf(fid, 'struct TrajectoryMetadata {\n');
    fprintf(fid, '    const char* name;\n');
    fprintf(fid, '    int num_points;\n');
    fprintf(fid, '    float t_start;\n');
    fprintf(fid, '    float t_end;\n');
    fprintf(fid, '    float dt_nominal;\n');
    fprintf(fid, '    float m_dry;\n');
    fprintf(fid, '    float g;\n');
    fprintf(fid, '};\n\n');

    % Trajectory Dataset struct
    fprintf(fid, 'struct TrajectoryDataset {\n');
    fprintf(fid, '    TrajectoryMetadata meta;\n');
    fprintf(fid, '    const float* time;\n');
    fprintf(fid, '    const float* states; // N x 15 (q(4), pos(3), vel(3), omega(3), mass(2))\n');
    fprintf(fid, '    const float* inputs; // N x 4  (pitch, yaw, thrust, roll)\n');
    fprintf(fid, '    const float* k_trans; // N x 18\n');
    fprintf(fid, '    const float* k_rot; // N x 18\n');
    fprintf(fid, '    const float* la; // N x 18\n');
    fprintf(fid, '    const float* lt; // N x 54\n');
    fprintf(fid, '};\n\n');

    % Emit data tables for each trajectory
    for i = 1:numTrajs
        t = trajectories{i};
        N = t.numPoints;
        prefix = ['TRAJ_', upper(t.ident)];

        fprintf(fid, '// Trajectory: %s (%d points, duration: %.3f s)\n', t.stem, N, t.time(end) - t.time(1));
        fprintf(fid, 'constexpr int %s_NUM_POINTS = %d;\n\n', prefix, N);

        % Time array
        fprintf(fid, 'static const float %s_TIME[%d]%s = {\n', prefix, N, progmemAttr);
        writeFloatArray(fid, t.time, 8, '    ');
        fprintf(fid, '};\n\n');

        % States array
        fprintf(fid, 'static const float %s_STATES[%d][15]%s = {\n', prefix, N, progmemAttr);
        for row = 1:N
            fprintf(fid, '    {');
            for col = 1:15
                if col < 15
                    fprintf(fid, '%s, ', formatFloat(t.states(row, col)));
                else
                    fprintf(fid, '%s', formatFloat(t.states(row, col)));
                end
            end
            if row < N
                fprintf(fid, '},\n');
            else
                fprintf(fid, '}\n');
            end
        end
        fprintf(fid, '};\n\n');

        % Inputs array
        fprintf(fid, 'static const float %s_INPUTS[%d][4]%s = {\n', prefix, N, progmemAttr);
        for row = 1:N
            fprintf(fid, '    {');
            for col = 1:4
                if col < 4
                    fprintf(fid, '%s, ', formatFloat(t.inputs(row, col)));
                else
                    fprintf(fid, '%s', formatFloat(t.inputs(row, col)));
                end
            end
            if row < N
                fprintf(fid, '},\n');
            else
                fprintf(fid, '}\n');
            end
        end
        fprintf(fid, '};\n\n');

        % K_trans
        fprintf(fid, 'static const float %s_K_TRANS[%d][18]%s = {\n', prefix, N, progmemAttr);
        writeMatrixArray(fid, t.kTrans, N, 18);
        fprintf(fid, '};\n\n');

        % K_rot
        fprintf(fid, 'static const float %s_K_ROT[%d][18]%s = {\n', prefix, N, progmemAttr);
        writeMatrixArray(fid, t.kRot, N, 18);
        fprintf(fid, '};\n\n');

        % LA
        fprintf(fid, 'static const float %s_LA[%d][18]%s = {\n', prefix, N, progmemAttr);
        writeMatrixArray(fid, t.la, N, 18);
        fprintf(fid, '};\n\n');

        % LT
        fprintf(fid, 'static const float %s_LT[%d][54]%s = {\n', prefix, N, progmemAttr);
        writeMatrixArray(fid, t.lt, N, 54);
        fprintf(fid, '};\n\n');
    end

    % Emit Trajectory Table Definitions
    fprintf(fid, '// Master Trajectory Datasets Table\n');
    fprintf(fid, 'static const TrajectoryDataset TRAJECTORY_REGISTRY[TRAJ_COUNT] = {\n');
    for i = 1:numTrajs
        t = trajectories{i};
        prefix = ['TRAJ_', upper(t.ident)];
        fprintf(fid, '    {\n');
        fprintf(fid, '        {"%s", %d, %.6gf, %.6gf, %.6gf, %.6gf, %.6gf},\n', ...
            t.stem, t.numPoints, single(t.time(1)), single(t.time(end)), single(t.dtNom), ...
            single(opts.DryMass), single(opts.Gravity));
        fprintf(fid, '        %s_TIME,\n', prefix);
        fprintf(fid, '        (const float*)%s_STATES,\n', prefix);
        fprintf(fid, '        (const float*)%s_INPUTS,\n', prefix);
        fprintf(fid, '        (const float*)%s_K_TRANS,\n', prefix);
        fprintf(fid, '        (const float*)%s_K_ROT,\n', prefix);
        fprintf(fid, '        (const float*)%s_LA,\n', prefix);
        fprintf(fid, '        (const float*)%s_LT\n', prefix);
        if i < numTrajs
            fprintf(fid, '    },\n');
        else
            fprintf(fid, '    }\n');
        end
    end
    fprintf(fid, '};\n\n');

    % Default active trajectory reference
    fprintf(fid, '// Default Active Trajectory (index 0: %s)\n', trajectories{1}.stem);
    fprintf(fid, 'inline const TrajectoryDataset& getDefaultTrajectory() {\n');
    fprintf(fid, '    return TRAJECTORY_REGISTRY[0];\n');
    fprintf(fid, '}\n\n');

    fprintf(fid, 'inline const TrajectoryDataset& getTrajectory(TrajectoryID id) {\n');
    fprintf(fid, '    if (id < 0 || id >= TRAJ_COUNT) {\n');
    fprintf(fid, '        return TRAJECTORY_REGISTRY[0];\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '    return TRAJECTORY_REGISTRY[id];\n');
    fprintf(fid, '}\n\n');

    fprintf(fid, '} // namespace %s\n', opts.Namespace);

    fprintf('\nSuccessfully generated header file: %s\n', headerPath);
end

%% Write 1D float array with line breaks
function writeFloatArray(fid, arr, perLine, indent)
    n = numel(arr);
    fprintf(fid, '%s', indent);
    for k = 1:n
        fprintf(fid, '%s', formatFloat(arr(k)));
        if k < n
            fprintf(fid, ', ');
        end
        if mod(k, perLine) == 0 && k < n
            fprintf(fid, '\n%s', indent);
        end
    end
end

%% Write 2D float array (N x M)
function writeMatrixArray(fid, mat, rows, cols)
    for r = 1:rows
        fprintf(fid, '    {');
        for c = 1:cols
            if c < cols
                fprintf(fid, '%s, ', formatFloat(mat(r, c)));
            else
                fprintf(fid, '%s', formatFloat(mat(r, c)));
            end
        end
        if r < rows
            fprintf(fid, '},\n');
        else
            fprintf(fid, '}\n');
        end
    end
end

%% Format single float for C++ literal (avoids invalid '0f', '1f')
function s = formatFloat(v)
    v = single(v);
    if v == 0
        s = '0.0f';
    else
        s = sprintf('%.8gf', v);
        if ~contains(s, '.') && ~contains(s, 'e') && ~contains(s, 'E')
            s = sprintf('%.1ff', v);
        end
    end
end

%% Pad or clamp rows if gain CSV length diverges
function out = padOrClampRows(mat, targetRows, targetCols)
    currRows = size(mat, 1);
    currCols = size(mat, 2);
    out = zeros(targetRows, targetCols);
    minR = min(currRows, targetRows);
    minC = min(currCols, targetCols);
    out(1:minR, 1:minC) = mat(1:minR, 1:minC);
    if currRows < targetRows
        % Repeat last row
        for r = (currRows + 1):targetRows
            out(r, :) = out(currRows, :);
        end
    end
end
