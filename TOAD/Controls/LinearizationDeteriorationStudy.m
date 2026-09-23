%% LinearizationDeteriorationStudy.m
% Quantifies how the upright linearization (and the LQR gain designed from
% it) deteriorates as attitude departs from upright.
%
% Reference : JacobianX / JacobianU at upright quaternion, U = [0;0;m*g;0]
% Comparison: JacobianX / JacobianU at tilted attitudes (same X otherwise)
%
% Uses the same reduced-state mapping and rotational/translational
% partitioning as TOAD_TVLQI.m, and SolveLQR.m for gain design.
%
% Reduced state (12): [theta(3); pos(3); vel(3); omega(3)]
% Inputs (4)        : [gimbal1; gimbal2; thrust; roll torque]
% Rotational loop   : states [1:3, 10:12], inputs [1, 2, 4]
%


%% ------------------------- USER SETUP ---------------------------------
% constantsTOAD must exist in the workspace, or create it here.
% e.g. constantsTOAD = TOAD_Constants();
if ~exist('constants6DoF', 'var')
    error(['constants6DoF not found. Load/construct it before running ' ...
           '(needs m_dry, m_wet, g, MaxThrust).']);
end

X_fuel  = [0; 0];            % X(14:15) fuel/prop masses used for the study
dT      = 0.01;              % discretization step used for spectral radius

angles_deg = 0:5:180;        % tilt sweep
axes_list  = { [1;0;0], [0;1;0], [0;0;1], [1;1;0]/sqrt(2) };
axes_names = {'Roll (x)', 'Pitch (y)', 'Yaw (z)', 'Diagonal (x+y)'};

% Thrust at tilted points. false -> keep thrust = m*g (isolates the effect
% of attitude alone, which is what the upright reference uses).
% true -> m*g/cos(tilt), clamped, to mimic a hover-holding thrust.
compensateThrust = true;

%% ----------------------------------------------------------------------

m_tot = constants6DoF.m_dry + sum(X_fuel);
g     = constants6DoF.g;
rotIdx = [1:3, 10:12];
rotIn  = [1, 2, 4];

% Base state at upright
X0 = zeros(15,1);
X0(1:4) = [1; 0; 0; 0];
X0(14:15) = X_fuel;
U0 = [0; 0; m_tot * g; 0];

% Translational model: fixed analytic double integrator (matches RicattiRecursion).
A_trans_c = [zeros(3,3), eye(3,3); zeros(3,3), zeros(3,3)];
B_trans_c = [zeros(3,3); eye(3,3)];
M_c_trans     = [A_trans_c, B_trans_c; zeros(3,6), zeros(3,3)];
M_d_trans     = expm(M_c_trans * dT);
A_trans_d_ref = M_d_trans(1:6, 1:6);
B_trans_d_ref = M_d_trans(1:6, 7:9);


% Rotational model: attitude+rate rows/cols; B uses theta/phi/roll cols only (thrust excluded).
%% Upright reference linearization
[Ar0, Br0] = reduceLin(X0, U0);
A_rot_c = Ar0(rotIdx, rotIdx);
B_rot_c = Br0(rotIdx, rotIn);
M_c_rot     = [A_rot_c, B_rot_c; zeros(3,6), zeros(3,3)];
M_d_rot     = expm(M_c_rot * dT);
A_rot_d_ref = M_d_rot(1:6, 1:6);
B_rot_d_ref = M_d_rot(1:6, 7:9);

Q_trans = constants6DoF.Q_trans; R_trans = constants6DoF.R_trans;
Q_rot   = constants6DoF.Q_rot;   R_rot   = constants6DoF.R_rot;

P_trans_ref = idare(A_trans_d_ref, B_trans_d_ref, Q_trans, R_trans);
P_rot_ref   = idare(A_rot_d_ref,   B_rot_d_ref,   Q_rot,   R_rot);
K_trans_ref = (B_trans_d_ref'*P_trans_ref*B_trans_d_ref + R_trans) \ (B_trans_d_ref'*P_trans_ref*A_trans_d_ref);
%K_rot_ref   = (B_rot_d_ref'*P_rot_ref*B_rot_d_ref     + R_rot)   \ (B_rot_d_ref'*P_rot_ref*A_rot_d_ref);
K_rot_ref = SolveLQR(A_rot_d_ref, B_rot_d_ref, Q_rot, R_rot);

SR0     = max(abs(eig(A_rot_d_ref - B_rot_d_ref * K_rot_ref)));
SigMax0 = max(real(log(eig(A_rot_d_ref - B_rot_d_ref*K_rot_ref)) / dT));
fprintf('Upright reference: rot spec. radius = %.5f, max Re(eig) = %.4f\n\n', ...
        SR0, SigMax0);

%% Sweep
nA = numel(angles_deg);
nX = numel(axes_list);

dA_rot   = nan(nA, nX);   % relative change in rot-subsystem A
dB_rot   = nan(nA, nX);   % relative change in rot-subsystem B
coupling = nan(nA, nX);   % ||dvel/dtheta|| block (thrust-direction sensitivity)
SR_rot   = nan(nA, nX);   % closed-loop spectral radius, fixed upright K
Sig_rot  = nan(nA, nX);   % closed-loop max Re(eig), fixed upright K
dK_rot   = nan(nA, nX);   % ||K(re-solved at tilt) - K0|| / ||K0||
dK_trans = nan(nA, nX);   % ||K_trans - K0|| / ||K0||

for j = 1:nX
    ax = axes_list{j};
    for i = 1:nA
        a = deg2rad(angles_deg(i));
        X = X0;
        X(1:4) = [cos(a/2); sin(a/2) * ax];

        if i==floor(nA/2)
            disp(X)
        end

        if compensateThrust
            Tn = min(m_tot * g / max(cos(a), 0.2), constants6DoF.MaxThrust);
        else
            Tn = m_tot * g;
        end
        U = [0; 0; Tn; 0];

        [Ar, Br] = reduceLin(X, U);
        Arot = Ar(rotIdx, rotIdx);
        Brot = Br(rotIdx, rotIn);

        dA_rot(i,j)  = norm(Arot - A_rot_c, 'fro') / norm(A_rot_c, 'fro');
        dB_rot(i,j)  = norm(Brot - B_rot_c, 'fro') / norm(B_rot_c, 'fro');
        
        coupling(i,j) = norm(Ar(7:9, 1:3), 'fro');

        % Fixed upright gain on the tilted linearization
        [Ad, Bd] = c2dExpm(Arot, Brot, dT);
        SR_rot(i,j)  = max(abs(eig(Ad - Bd * K_rot_ref)));
        Sig_rot(i,j) = max(real(eig(Arot - Brot * K_rot_ref)));

        % Gain that LQR would pick if re-linearized at this attitude
        try
            
 
        P_rot   = idare(Ad,   Bd,   Q_rot,   R_rot);
        K_rot   = (Bd'*P_rot*Bd    + R_rot)   \ (Bd'*P_rot*Ad);

            dK_rot(i,j) = norm(K_rot - K_rot_ref, 'fro') / norm(K_rot_ref, 'fro');
        catch
            dK_rot(i,j) = NaN;
        end
    end
end

%% Summary: first tilt at which the upright design fails
fprintf('%-16s | first tilt (deg) with SpecRad >= 1 | with Re(eig) >= 0\n', 'Axis');
fprintf('%s\n', repmat('-', 1, 70));
for j = 1:nX
    iSR  = find(SR_rot(:,j)  >= 1, 1, 'first');
    iSig = find(Sig_rot(:,j) >= 0, 1, 'first');
    fprintf('%-16s | %-29s | %s\n', axes_names{j}, strOrNone(angles_deg, iSR), ...
            strOrNone(angles_deg, iSig));
end

for j = 1:nX
    T = table(angles_deg(:), dA_rot(:,j), dB_rot(:,j), SR_rot(:,j), dK_rot(:,j), ...
        'VariableNames', {'Tilt_deg','dA_rot','dB_rot','SpecRad','dK'});
    fprintf('\nDetail for %s axis:\n', axes_names{j});
    disp(T);

end


%% Plots
cols = lines(nX);
figure('Name', 'Linearization deterioration', 'Position', [100 100 1200 800]);

subplot(2,3,1); hold on; grid on;
for j = 1:nX, plot(angles_deg, dA_rot(:,j), 'Color', cols(j,:), 'LineWidth', 1.5); end
title('Rot A: ||A(\theta)-A_0|| / ||A_0||'); xlabel('Tilt [deg]'); ylabel('Relative change');

subplot(2,3,2); hold on; grid on;
for j = 1:nX, plot(angles_deg, dB_rot(:,j), 'Color', cols(j,:), 'LineWidth', 1.5); end
title('Rot B: ||B(\theta)-B_0|| / ||B_0||'); xlabel('Tilt [deg]');

subplot(2,3,3); hold on; grid on;
for j = 1:nX, plot(angles_deg, coupling(:,j), 'Color', cols(j,:), 'LineWidth', 1.5); end
title('||\partial v̇ / \partial\theta|| (thrust coupling)'); xlabel('Tilt [deg]');

subplot(2,3,4); hold on; grid on;
for j = 1:nX, plot(angles_deg, SR_rot(:,j), 'Color', cols(j,:), 'LineWidth', 1.5); end
yline(1, 'r--', 'Unstable'); yline(SR0, 'k:', 'Upright');
title('Closed-loop spectral radius (fixed upright K)');
xlabel('Tilt [deg]'); ylabel('max |\lambda_d|');

subplot(2,3,5); hold on; grid on;
for j = 1:nX, plot(angles_deg, Sig_rot(:,j), 'Color', cols(j,:), 'LineWidth', 1.5); end
yline(0, 'r--', 'Unstable');
title('Closed-loop max Re(\lambda) (fixed upright K)');
xlabel('Tilt [deg]'); ylabel('[1/s]');

subplot(2,3,6); hold on; grid on;
for j = 1:nX, plot(angles_deg, dK_rot(:,j), 'Color', cols(j,:), 'LineWidth', 1.5); end
title('Gain drift: ||K_{LQR}(\theta)-K_0|| / ||K_0||'); xlabel('Tilt [deg]');

legend(axes_names, 'Location', 'best');

%% ------------------------- Local functions ----------------------------
function [Ar, Br] = reduceLin(X, U)
    % Same reduction as TOAD_TVLQI: 15-state -> 12-state attitude-error form
    A = JacobianX(X, U);
    B = JacobianU(X, U);
    Tm = zeros(15, 12);
    Tm(1:4, 1:3)   = 0.5 * XiMat(X(1:4));
    Tm(5:13, 4:12) = eye(9);
    Ar = pinv(Tm) * A * Tm;
    Br = pinv(Tm) * B;
end

function [Ad, Bd] = c2dExpm(A, B, dT)
    n = size(A, 1);
    m = size(B, 2);
    M = expm([A, B; zeros(m, n + m)] * dT);
    Ad = M(1:n, 1:n);
    Bd = M(1:n, n+1:n+m);
end

function s = strOrNone(angles, idx)
    if isempty(idx), s = 'never (within sweep)';
    else,            s = sprintf('%d', angles(idx));
    end
end
