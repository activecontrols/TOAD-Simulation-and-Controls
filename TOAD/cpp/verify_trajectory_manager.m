%% verify_trajectory_manager.m
% Verification harness asserting 1:1 input/output numerical parity between:
%   - MATLAB implementation:  TOAD/Guidance/TrajectoryManager.m
%   - C++ implementation:     TOAD/cpp/TrajectoryManager.cpp (via trajectory_manager_mex)
%
% Tests boundary conditions (t < 0, t = 0, t > t_end), interior knot points,
% sub-interval points, and quaternion normalization.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)
    thisDir = pwd;
end
addpath(fullfile(thisDir, '..', 'Guidance'));
addpath(thisDir);

%% Build MEX if not present or out of date
mexName = ['trajectory_manager_mex.', mexext];
if ~exist(fullfile(thisDir, mexName), 'file')
    fprintf('[Step 1] Compiling trajectory_manager_mex...\n');
    run(fullfile(thisDir, 'build_trajectory_manager_mex.m'));
else
    fprintf('[Step 1] trajectory_manager_mex binary found.\n');
end

%% Initialize both implementations
fprintf('[Step 2] Initializing MATLAB & C++ Trajectory Managers...\n');
trajectory_manager_mex('init');
trajectory_manager_mex('reset');

% Load trajectory metadata for test sweep
csvPath = fullfile(thisDir, '..', 'Guidance', 'Trajectories', 'ASTRA_Circle_v003.csv');
if ~exist(csvPath, 'file')
    error('verify_trajectory_manager:MissingCSV', 'Trajectory CSV not found: %s', csvPath);
end
trajTable = readmatrix(csvPath);
t_knots = trajTable(:, 1);
t_start = t_knots(1);
t_end   = t_knots(end);

m_dry = 1.1690;
g_acc = 9.8015;

%% Construct Comprehensive Evaluation Time Sweep
t_fine = linspace(t_start, t_end, 2000)';
t_mid  = 0.5 * (t_knots(1:end-1) + t_knots(2:end));
t_test = sort(unique([
    -2.0; -0.001; 0.0;
    t_knots;
    t_mid;
    t_fine;
    t_end; t_end + 0.1; t_end + 5.0; t_end + 20.0
    ]));

numTests = numel(t_test);
fprintf('[Step 3] Running evaluation across %d test time points...\n', numTests);

% Mock estimated state
X_est_mock = zeros(15, 1);
X_est_mock(1)  = 0.999;
X_est_mock(2)  = 0.01;
X_est_mock(3)  = -0.02;
X_est_mock(4)  = 0.005;
X_est_mock(5:7) = [0.1; -0.2; 1.5];
X_est_mock(8:10) = [0.05; -0.02; 0.8];
X_est_mock(11:13) = [0.01; 0.02; -0.01];
X_est_mock(14) = 0.05;
X_est_mock(15) = 0.05;

%% Execute Sweep & Accumulate Discrepancies
err_X  = zeros(numTests, 1);
err_U  = zeros(numTests, 1);
err_KT = zeros(numTests, 1);
err_KR = zeros(numTests, 1);
err_LA = zeros(numTests, 1);
err_LT = zeros(numTests, 1);

tic;
for k = 1:numTests
    t_curr = t_test(k);

    % MATLAB Reference Evaluation
    [X_m, U_m, KT_m, KR_m, LA_m, LT_m] = TrajectoryManager(t_curr, X_est_mock);

    % C++ Flight Evaluation
    [X_c, U_c, KT_c, KR_c, LA_c, LT_c] = trajectory_manager_mex('step', t_curr, X_est_mock, m_dry, g_acc);

    % Compute Maximum Absolute Errors
    err_X(k)  = max(abs(X_m(:)  - X_c(:)));
    err_U(k)  = max(abs(U_m(:)  - U_c(:)));
    err_KT(k) = max(abs(KT_m(:) - KT_c(:)));
    err_KR(k) = max(abs(KR_m(:) - KR_c(:)));
    err_LA(k) = max(abs(LA_m(:) - LA_c(:)));
    err_LT(k) = max(abs(LT_m(:) - LT_c(:)));
end
tElapsed = toc;

max_err_X  = max(err_X);
max_err_U  = max(err_U);
max_err_KT = max(err_KT);
max_err_KR = max(err_KR);
max_err_LA = max(err_LA);
max_err_LT = max(err_LT);

%% Print Results & Assert Strict Tolerance
fprintf('\n--- PARITY VERIFICATION SUMMARY ---\n');
fprintf('Evaluated %d time steps in %.3f s (%.1f microseconds/eval)\n', ...
    numTests, tElapsed, (tElapsed / numTests) * 1e6);
fprintf('Max |X_matlab  - X_cpp|:       %.3e\n', max_err_X);
fprintf('Max |U_matlab  - U_cpp|:       %.3e\n', max_err_U);
fprintf('Max |KT_matlab - KT_cpp|:      %.3e\n', max_err_KT);
fprintf('Max |KR_matlab - KR_cpp|:      %.3e\n', max_err_KR);
fprintf('Max |LA_matlab - LA_cpp|:      %.3e\n', max_err_LA);
fprintf('Max |LT_matlab - LT_cpp|:      %.3e\n', max_err_LT);


TOL = 1e-5;
assert(max_err_X  < TOL, sprintf('VERIFICATION FAILED: State vector error (%.3e) exceeds tolerance (%.1e).', max_err_X, TOL));
assert(max_err_U  < TOL, sprintf('VERIFICATION FAILED: Control input error (%.3e) exceeds tolerance (%.1e).', max_err_U, TOL));
assert(max_err_KT < TOL, sprintf('VERIFICATION FAILED: K_trans gain error (%.3e) exceeds tolerance (%.1e).', max_err_KT, TOL));
assert(max_err_KR < TOL, sprintf('VERIFICATION FAILED: K_rot gain error (%.3e) exceeds tolerance (%.1e).', max_err_KR, TOL));
assert(max_err_LA < TOL, sprintf('VERIFICATION FAILED: LA gain error (%.3e) exceeds tolerance (%.1e).', max_err_LA, TOL));
assert(max_err_LT < TOL, sprintf('VERIFICATION FAILED: LT gain error (%.3e) exceeds tolerance (%.1e).', max_err_LT, TOL));

fprintf('\n>>> SUCCESS: MATLAB and C++ TrajectoryManager outputs match 1:1 within single-precision tolerance (%.1e)! <<<\n\n', TOL);
