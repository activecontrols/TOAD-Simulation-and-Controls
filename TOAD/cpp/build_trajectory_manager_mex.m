% build_trajectory_manager_mex.m
% Compiles the C++ TrajectoryManager into a MEX function for verification in MATLAB

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)
    thisDir = pwd;
end

eigenDir = fullfile(thisDir, 'Extras', 'eigen-3.4.0', 'eigen-3.4.0');
if ~isfile(fullfile(eigenDir, 'Eigen', 'Dense'))
    eigenDir = fullfile(thisDir, 'Extras', 'eigen-3.4.0');
end

% Ensure ASTRA_Trajectory_Data.h exists
headerFile = fullfile(thisDir, 'ASTRA_Trajectory_Data.h');
if ~isfile(headerFile)
    fprintf('ASTRA_Trajectory_Data.h not found. Generating via ExportTrajectoryHeader...\n');
    guidanceDir = fullfile(thisDir, '..', 'Guidance');
    addpath(guidanceDir);
    ExportTrajectoryHeader('ASTRA_Circle_v003', 'OutputFile', headerFile);
end

mex('-D_USE_MATH_DEFINES', ...
    ['-I' fullfile(thisDir, 'Extras')], ...
    ['-I' eigenDir], ...
    ['-I' thisDir], ...
    fullfile(thisDir, 'trajectory_manager_mex.cpp'), ...
    fullfile(thisDir, 'TrajectoryManager.cpp'));

disp('TrajectoryManager MEX build finished successfully.');
