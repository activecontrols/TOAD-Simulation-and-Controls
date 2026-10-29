function [X, U, K, LA, LT] = TrajectoryManager(t, X_est, constantsTOAD)
%% TRAJECTORYMANAGER  Time-indexed reference + gain interpolation, with in-flight abort.
%   Nominal: interpolates constantsTOAD.Traj (Time, States, Inputs, gains).
%   Abort:   at t >= ABORT the abort trajectory is solved once (PolyQP, Mode="Abort"),
%            stored in persistents, and the SAME interpolation below runs on it.
%
%   X_est ordering (15x1): [q(4); r(3); v(3); w(3); ml; mi]

NN    = 61;   % abort trajectory nodes = Opts.N + 1 (fixed so persistents have fixed size)
ABORT = 0;    % [s] sim time at which abort fires; 0 = disabled

persistent abortActive abortTried
persistent AbTime AbStates AbInputs AbKT AbKR AbLA AbLT
if isempty(abortActive)
    abortActive = false;
    abortTried  = false;
    AbTime   = zeros(NN, 1);
    AbStates = zeros(NN, 15);
    AbInputs = zeros(NN, 4);
    AbKT = zeros(NN, 3, 6);
    AbKR = zeros(NN, 3, 6);
    AbLA = zeros(NN, 6, 3);
    AbLT = zeros(NN, 9, 6);
    ends

%% Abort: one-shot solve, then switch trajectory source
calcABORT = (ABORT > 0) && (t >= ABORT) && ~abortTried;
if calcABORT
    abortTried = true;                      % never re-solve, even if infeasible

    mL = X_est(14, 1); mI = X_est(15, 1);
    Ab = struct('r', X_est(5:7, 1), 'v', X_est(8:10, 1), 'a', [0; 0; 0], ...
                'q', X_est(1:4, 1), 'w', X_est(11:13, 1), ...
                'm', constantsTOAD.m_dry + mL + mI, 'ml', mL, 'mi', mI, ...
                'r_lz', constantsTOAD.Traj.States(end, 5:7).');   % ASSUMPTION: last nominal node = LZ
    Opts = struct('Mode', "Abort", 'Vehicle', 1, 'N', NN - 1);

    [IG, ~] = PolyQP(Ab, [], constantsTOAD, Opts);

    if ~IG.Abort.Feasible
        warning('TrajectoryManager:AbortInfeasible', ...
            'Abort trajectory infeasible (thrust/att metric %.3g, corridor %.3g).', ...
            IG.Abort.ThrustAttMetric, IG.Abort.CorridorViol);
    end

    % Trajectory -> same layout as constantsTOAD.Traj (rows = time samples)
    AbTime   = t + IG.Time(:);                    % absolute time
    AbStates = IG.X.';                            % NN x 15
    AbInputs = [IG.U, IG.U(:, end)].';            % U has N cols; pad to NN rows

    % Abort has no scheduled gains: freeze the nominal gains at the abort instant
    n0 = find(constantsTOAD.Traj.Time <= t, 1, 'last');
    if isempty(n0), n0 = 1; end
    KT0 = reshape(constantsTOAD.Traj.KTGain(n0, :, :), 1, 3, 6);
    KR0 = reshape(constantsTOAD.Traj.KRGain(n0, :, :), 1, 3, 6);
    LA0 = reshape(constantsTOAD.Traj.LAGain(n0, :, :), 1, 6, 3);
    LT0 = reshape(constantsTOAD.Traj.LTGain(n0, :, :), 1, 9, 6);
    for k = 1:NN
        AbKT(k, :, :) = KT0;
        AbKR(k, :, :) = KR0;
        AbLA(k, :, :) = LA0;
        AbLT(k, :, :) = LT0;
    end
    abortActive = true;
end

%% Select trajectory source
if abortActive
    Time = AbTime;  States = AbStates;  Inputs = AbInputs;
    K_trans_Gain = AbKT;  K_rot_Gain = AbKR;  LA_Gain = AbLA;  LT_Gain = AbLT;
else
    Time = constantsTOAD.Traj.Time;
    States = constantsTOAD.Traj.States;
    Inputs = constantsTOAD.Traj.Inputs;
    % FBGain maps to K_trans, FBCost maps to K_rot
    K_trans_Gain = constantsTOAD.Traj.KTGain;
    K_rot_Gain = constantsTOAD.Traj.KRGain;
    LA_Gain = constantsTOAD.Traj.LAGain;
    LT_Gain = constantsTOAD.Traj.LTGain;
end

X = zeros(15, 1);
U = zeros(4, 1);
K = zeros(6, 6);
LA = zeros(6, 3);
LT = zeros(9, 6);

%% Handle Boundary Conditions
if t <= Time(1) || t == 0
    X(:) = States(1, :);
    U(:) = Inputs(1, :);
    K(1:3, 1:6) = K_trans_Gain(1, :, :);
    K(4:6, 1:6) = K_rot_Gain(1, :, :);
    LA(:) = LA_Gain(1, :, :);
    LT(:) = LT_Gain(1, :, :);
    return;
elseif t >= Time(end)
    X(:) = States(end, :);
    X(1:4) = [1; 0; 0; 0];
    X(8:end) = 0;
    total_mass = constantsTOAD.m_dry + X_est(14, 1) + X_est(15, 1);
    U(:) = [0; 0; total_mass * constantsTOAD.g; 0];
    K(1:3, 1:6) = K_trans_Gain(end, :, :);
    K(4:6, 1:6) = K_rot_Gain(end, :, :);
    LA(:) = LA_Gain(end, :, :);
    LT(:) = LT_Gain(end, :, :);
    return;
end

n_low = 1;
for i = 1:(length(Time) - 1)
    if Time(i) <= t && Time(i+1) > t
        n_low = i;
        break;
    end
end
n_high = n_low + 1;
T_low = Time(n_low);
T_high = Time(n_high);
dt = T_high - T_low;
ratio = (t - T_low) / dt;

%% Interpolation
LA_low = zeros(6,3); LA_low(:) = LA_Gain(n_low,:,:);
LA_up  = zeros(6,3); LA_up(:)  = LA_Gain(n_high,:,:);
LA(:) = LA_low + (LA_up - LA_low) .* ratio;

LT_low = zeros(9,6); LT_low(:) = LT_Gain(n_low,:,:);
LT_up  = zeros(9,6); LT_up(:)  = LT_Gain(n_high,:,:);
LT(:) = LT_low + (LT_up - LT_low) .* ratio;

K_trans_low = zeros(3,6); K_trans_low(:) = K_trans_Gain(n_low,:,:);
K_trans_up  = zeros(3,6); K_trans_up(:)  = K_trans_Gain(n_high,:,:);
K(1:3, 1:6) = K_trans_low + (K_trans_up - K_trans_low) .* ratio;

K_rot_low = zeros(3,6); K_rot_low(:) = K_rot_Gain(n_low,:,:);
K_rot_up  = zeros(3,6); K_rot_up(:)  = K_rot_Gain(n_high,:,:);
K(4:6, 1:6) = K_rot_low + (K_rot_up - K_rot_low) .* ratio;

X_low = States(n_low, :);
X_up  = States(n_high, :);
X(:) = X_low + (X_up - X_low) .* ratio;

% Quaternion normalization to eliminate norm shrinkage and DCM defect
q_norm = norm(X(1:4));
if q_norm > 1e-6
    X(1:4) = X(1:4) / q_norm;
end

U_low = Inputs(n_low, :);
U_up  = Inputs(n_high, :);
U(:) = U_low + (U_up - U_low) .* ratio;

end