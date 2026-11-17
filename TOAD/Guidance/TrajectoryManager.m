function [X, U, varargout] = TrajectoryManager(t, X_est, constantsTOAD)
%
% evb clap for claude documentation
%
% Interpolates nominal state targets, feedforward actuator commands, and time-varying
% feedback (TV-LQI) and disturbance observer (LESO) gains across flight time.
%
% Inputs:
%   t              - Current flight elapsed time (MET in seconds)
%   X_est          - Current estimated state (15x1 or 19x1 float/double vector)
%   constantsTOAD  - Vehicle constants struct containing constantsTOAD.Traj
%                    (If omitted or empty, loads default ASTRA_Circle_v003)
%
% Outputs (Supports both 5-output legacy Simulink API and 6-output flight API):
%   [X, U, K, LA, LT]                 - Legacy Simulink API (K is 6x6: K_trans; K_rot)
%   [X, U, K_trans, K_rot, LA, LT]   - Flight Software 1:1 API with C++
%
% Mathematical Flow (Identical to C++ TrajectoryManager.cpp):
%   1. Boundary checks:
%      - If t <= Time(1): Return initial point (t=0) values.
%      - If t >= Time(end): Hold terminal position, upright attitude [1,0,0,0],
%        zero velocity/rates, hover thrust (m_total * g), and terminal gains.
%   2. Interval localization:
%      - Find n_low such that Time(n_low) <= t < Time(n_high).
%   3. Linear interpolation parameter:
%      - ratio = (t - Time(n_low)) / (Time(n_high) - Time(n_low))
%   4. Linear state & control interpolation:
%      - X = X_low + ratio * (X_high - X_low)
%      - U = U_low + ratio * (U_high - U_low)
%   5. Quaternion normalization (Mandatory):
%      - q = X(1:4) / norm(X(1:4))
%   6. Gain interpolation:
%      - K_trans = K_trans_low + ratio * (K_trans_high - K_trans_low)
%      - K_rot   = K_rot_low   + ratio * (K_rot_high   - K_rot_low)
%      - LA      = LA_low      + ratio * (LA_high      - LA_low)
%      - LT      = LT_low      + ratio * (LT_high      - LT_low)

persistent cachedTraj

%% Resolve Trajectory Data
if nargin < 3 || isempty(constantsTOAD) || ~isfield(constantsTOAD, 'Traj')
    if isempty(cachedTraj)
        defaultCsv = fullfile(fileparts(mfilename('fullpath')), 'Trajectories', 'ASTRA_Circle_v003.csv');
        if exist(defaultCsv, 'file')
            rawData = readmatrix(defaultCsv);
            cachedTraj.Time   = rawData(:, 1);
            cachedTraj.States = rawData(:, 2:16);
            cachedTraj.Inputs = rawData(:, 17:20);

            % Load Gains
            gainsDir = fullfile(fileparts(mfilename('fullpath')), 'Trajectories', 'Gains');
            ktFile = fullfile(gainsDir, 'K_trans_ASTRA_Circle_v003.csv');
            krFile = fullfile(gainsDir, 'K_rot_ASTRA_Circle_v003.csv');
            laFile = fullfile(gainsDir, 'LA_ASTRA_Circle_v003.csv');
            ltFile = fullfile(gainsDir, 'LT_ASTRA_Circle_v003.csv');

            if exist(ktFile, 'file')
                cachedTraj.KTGain = reshape(readmatrix(ktFile), [], 3, 6);
                cachedTraj.KRGain = reshape(readmatrix(krFile), [], 3, 6);
                cachedTraj.LAGain = reshape(readmatrix(laFile), [], 6, 3);
                cachedTraj.LTGain = reshape(readmatrix(ltFile), [], 9, 6);
            else
                N = size(rawData, 1);
                cachedTraj.KTGain = zeros(N, 3, 6);
                cachedTraj.KRGain = zeros(N, 3, 6);
                cachedTraj.LAGain = zeros(N, 6, 3);
                cachedTraj.LTGain = zeros(N, 9, 6);
            end
            cachedTraj.m_dry = 1.1690;
            cachedTraj.g = 9.8015;
        else
            error('TrajectoryManager:NoData', 'No trajectory data supplied and default file not found.');
        end
    end
    Time   = cachedTraj.Time;
    States = cachedTraj.States;
    Inputs = cachedTraj.Inputs;
    KTGain = cachedTraj.KTGain;
    KRGain = cachedTraj.KRGain;
    LAGain = cachedTraj.LAGain;
    LTGain = cachedTraj.LTGain;
    m_dry  = cachedTraj.m_dry;
    g_acc  = cachedTraj.g;
else
    Time   = constantsTOAD.Traj.Time;
    States = constantsTOAD.Traj.States;
    Inputs = constantsTOAD.Traj.Inputs;
    KTGain = constantsTOAD.Traj.KTGain;
    KRGain = constantsTOAD.Traj.KRGain;
    LAGain = constantsTOAD.Traj.LAGain;
    LTGain = constantsTOAD.Traj.LTGain;
    if isfield(constantsTOAD, 'm_dry')
        m_dry = constantsTOAD.m_dry;
    else
        m_dry = 1.1690;
    end
    if isfield(constantsTOAD, 'g')
        g_acc = constantsTOAD.g;
    else
        g_acc = 9.8015;
    end
end

if nargin < 2 || isempty(X_est)
    X_est = zeros(15, 1);
    X_est(1) = 1.0;
end

% Preallocate outputs
X = zeros(15, 1);
U = zeros(4, 1);
K_trans = zeros(3, 6);
K_rot   = zeros(3, 6);
LA      = zeros(6, 3);
LT      = zeros(9, 6);

N = length(Time);

%% Boundary Condition Handling
if t <= Time(1)
    X(:) = States(1, :);
    U(:) = Inputs(1, :);
    K_trans(:, :) = squeeze(KTGain(1, :, :));
    K_rot(:, :)   = squeeze(KRGain(1, :, :));
    LA(:, :)      = squeeze(LAGain(1, :, :));
    LT(:, :)      = squeeze(LTGain(1, :, :));
    varargout = packOutputs(K_trans, K_rot, LA, LT, nargout);
    return;
elseif t >= Time(end)
    X(:) = States(end, :);
    % Command upright attitude at terminal hold
    X(1:4) = [1.0; 0.0; 0.0; 0.0];
    % Zero out velocities and body rates
    X(8:13) = 0.0;

    % Calculate total mass (dry mass + remaining fuel & oxidizer)
    if numel(X_est) >= 15
        m_prop = X_est(14) + X_est(15);
    elseif numel(X_est) >= 2
        m_prop = X(14) + X(15);
    else
        m_prop = 0.0;
    end
    total_mass = m_dry + m_prop;

    % Hover feedforward command
    U(:) = [0.0; 0.0; total_mass * g_acc; 0.0];

    K_trans(:, :) = squeeze(KTGain(end, :, :));
    K_rot(:, :)   = squeeze(KRGain(end, :, :));
    LA(:, :)      = squeeze(LAGain(end, :, :));
    LT(:, :)      = squeeze(LTGain(end, :, :));
    varargout = packOutputs(K_trans, K_rot, LA, LT, nargout);
    return;
end

%% Interval Search
n_low = 1;
for i = 1:(N - 1)
    if Time(i) <= t && Time(i + 1) > t
        n_low = i;
        break;
    end
end
n_high = n_low + 1;

%% Interpolation Parameter
T_low  = Time(n_low);
T_high = Time(n_high);
dt     = T_high - T_low;
ratio  = (t - T_low) / dt;

%% State and Control Interpolation
X_low = States(n_low, :)';
X_up  = States(n_high, :)';
X = X_low + (X_up - X_low) .* ratio;

% Quaternion normalization
q_norm = norm(X(1:4));
if q_norm > 1e-6
    X(1:4) = X(1:4) ./ q_norm;
end

U_low = Inputs(n_low, :)';
U_up  = Inputs(n_high, :)';
U = U_low + (U_up - U_low) .* ratio;

%% TV-LQI and LESO Gain Interpolation
KT_low = squeeze(KTGain(n_low, :, :));
KT_up  = squeeze(KTGain(n_high, :, :));
K_trans = KT_low + (KT_up - KT_low) .* ratio;

KR_low = squeeze(KRGain(n_low, :, :));
KR_up  = squeeze(KRGain(n_high, :, :));
K_rot = KR_low + (KR_up - KR_low) .* ratio;

LA_low = squeeze(LAGain(n_low, :, :));
LA_up  = squeeze(LAGain(n_high, :, :));
LA = LA_low + (LA_up - LA_low) .* ratio;

LT_low = squeeze(LTGain(n_low, :, :));
LT_up  = squeeze(LTGain(n_high, :, :));
LT = LT_low + (LT_up - LT_low) .* ratio;

varargout = packOutputs(K_trans, K_rot, LA, LT, nargout);
end

%% Helper to pack outputs according to nargout
function vout = packOutputs(K_trans, K_rot, LA, LT, nOut)
if nOut <= 2
    vout = {};
elseif nOut <= 5
    % simulink format
    K = zeros(6, 6);
    K(1:3, 1:6) = K_trans;
    K(4:6, 1:6) = K_rot;
    raw = {K, LA, LT};
    vout = raw(1:min(numel(raw), nOut - 2));
else
    % cpp format
    raw = {K_trans, K_rot, LA, LT};
    vout = raw(1:min(numel(raw), nOut - 2));
end
end
