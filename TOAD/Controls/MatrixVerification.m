%% Test Script for TVLQR Gain Generation 
% Call and load in parameters
if ~exist("constantsTOAD")
    LoadTOADSim;
end

% Hand tuning for Q for now
a_weights = ones(12,1);
a_weights = a_weights / norm(a_weights);
max_x = [0.23, 0.23, 0.04, 4, 4, 4, 5, 5, 5, 0.72, 0.72, 1];
Q = eye(12) .* a_weights ./ max_x.^2;
R = diag([3.5, 3.5, 1/2400^2, 5]);

K_List=ReadGains("GainMatrixTrajectory1.csv");
m = readmatrix(".\Guidance\Trajectories\Trajectory1.csv");
dT = m(2,1)-m(1,1);
x = m(:, 2:16)';
u = m(:, 17:20)';

[A_lin, B_lin] = generateAB(x,u, constantsTOAD);
SanityCheck = lqrd(A_lin, B_lin, Q, R, dT);

%% Plots
N = size(K_List, 1);
error_norm = zeros(N, 1);
K_converged = zeros(size(SanityCheck, 1));

% Calculate the error between the recursion and DARE at every time step
for n = 1:N
    K_t = squeeze(K_List(n, :, :))'; 
    
    % Frobenius norm of the difference matrix
    error_norm(n) = norm(K_t - SanityCheck, 'fro'); 
end

K_converged = squeeze(K_List(1, :, :))';

% Command Window Summary
disp('--- LQR Recursion Diagnostics ---');
disp(['Horizon Length (N): ', num2str(N), ' steps']);
disp(['Initial Error (at terminal step N): ', num2str(error_norm(N))]);
disp(['Final Error (at step 1, fully converged): ', num2str(error_norm(1))]);
disp('---------------------------------');

% Visualization Plots
figure('Name', 'Riccati Recursion Convergence', 'WindowStyle', 'docked');

% Convergence Error over Time
subplot(1, 2, 1);
semilogy(N:-1:1, flipud(error_norm), 'LineWidth', 2, 'Color', 'b');
set(gca, 'XDir', 'reverse');
grid on;
title('Convergence to DARE (Frobenius Norm)');
xlabel('Backward Recursion Step (N \rightarrow 1)');
ylabel('|| K_t - K_{dare} ||_F');

% Trajectories of the Largest Gain Components
subplot(1, 2, 2);
hold on; grid on;

% Find the indices of the 5 largest elements in the SanityCheck matrix 
[~, sorted_idx] = sort(abs(SanityCheck(:)), 'descend');
top_5_idx = sorted_idx(1:min(5, length(sorted_idx)));

% Extract and plot these specific components over time
for i = 1:length(top_5_idx)
    [row, col] = ind2sub(size(SanityCheck), top_5_idx(i));
    component_trajectory = squeeze(K_List(:, col, row)); 
    
    plot(N:-1:1, flipud(component_trajectory), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('K_{%d,%d}', row, col));
end

% Plot the DARE targets as dashed lines
for i = 1:length(top_5_idx)
    [row, col] = ind2sub(size(SanityCheck), top_5_idx(i));
    dare_val = SanityCheck(row, col);
    yline(dare_val, '--k', 'HandleVisibility', 'off'); % Target line
end

set(gca, 'XDir', 'reverse');
title('Gain Component Trajectories');
xlabel('Backward Recursion Step (N \rightarrow 1)');
ylabel('Gain Value');
legend('Location', 'best');
hold off;

function [A_lin,B_lin] = generateAB(x, u, constantsTOAD)

% state variables, x
syms q0 q1 q2 q3
syms r1 r2 r3
syms v1 v2 v3
syms omega1 omega2 omega3
syms m_lox m_ipa
% control for u
syms theta phi thrust roll

% state x_n-1
syms q01 q11 q21 q31
syms r11 r21 r31
syms v11 v21 v31
syms omega11 omega21 omega31
syms m_lox1 m_ipa1

q = [q0;q1;q2;q3];
r = [r1;r2;r3];
v = [v1;v2;v3];
m = [m_lox; m_ipa];
omegaB = [omega1; omega2; omega3];

xn = [q;r;v; omegaB; m]; % x_n

% U
un = [theta;phi;thrust;roll];

% X_n-1
qn1 = [q01;q11;q21;q31];
rn1 = [r11;r21;r31];
vn1 = [v11;v21;v31];
mn1 = [m_lox1; m_ipa1];
omegaBn1 = [omega11;omega21;omega31];
xn1 = [qn1;rn1;vn1;omegaBn1;mn1]; % x_n-1

%% Dynamics
TB = thrust * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
m_dry = constantsTOAD.m_dry;
m = m_dry + m_lox + m_ipa;

% Propellant Fill height
OxFluidHeight = (m_lox / constantsTOAD.OxMass) * constantsTOAD.OxHeight * 0.9;
FuFluidHeight = (m_ipa / constantsTOAD.FuMass) * constantsTOAD.FuHeight * 0.9;

% Propellant inertias
J_xx = 1/12 * m_lox * (3 * constantsTOAD.OxRadius^2 + OxFluidHeight^2);
J_zz = 1/2 * m_lox * constantsTOAD.OxRadius^2;
J_lox = [J_xx, 0, 0;
         0, J_xx, 0;
         0,    0, J_zz];

J_xx = 1/12 * m_ipa * (3 * constantsTOAD.FuRadius^2 + FuFluidHeight^2);
J_zz = 1/2 * m_ipa * constantsTOAD.FuRadius^2;
J_ipa = [J_xx, 0, 0;
         0, J_xx, 0;
         0,    0, J_zz];

% Fluid fills & location of CoM
OxFluidLocation = constantsTOAD.Ox_Z + OxFluidHeight / 2;
FuFluidLocation = constantsTOAD.Fu_Z + FuFluidHeight / 2;
CGz = (constantsTOAD.m_dry * constantsTOAD.rTB + m_lox * OxFluidLocation + m_ipa * FuFluidLocation) / m;

% Distances to CG
d_dry = constantsTOAD.rTB - CGz;
d_lox = OxFluidLocation - CGz;
d_ipa = FuFluidLocation - CGz;

% Shifted Inertias
J_dry = constantsTOAD.J + m_dry * diag([d_dry^2, d_dry^2, 0]);
J_lox = J_lox + m_lox * diag([d_lox^2, d_lox^2, 0]);
J_ipa = J_ipa + m_ipa * diag([d_ipa^2, d_ipa^2, 0]);
% Bring everything to instantaneous CG (w.r.t Engine attachment frame)
J_tot = J_dry + J_lox + J_ipa;

% Angular dynamics
MB = zetaCross([0; 0; -CGz])*TB + (TB * roll)/thrust;
M = [q(1) -q(2) -q(3) -q(4);
     q(2)  q(1) -q(4)  q(3);
     q(3)  q(4)  q(1) -q(2);
     q(4) -q(3)  q(2)  q(1)];
qdot = 0.5 * M * [0; omegaB];
omegaBdot = J_tot \ (MB - zetaCross(omegaB) * J_tot * omegaB);

% Linear Dynamics
C_BI = quatRot(q);
C_IB = C_BI.';
FI = C_IB * TB + [0; 0;-m*constantsTOAD.g];

rdot = [v1;v2;v3];
vdot = FI/m;

% Mass Dynamics
mdot_lox = -thrust / constantsTOAD.MaxThrust * constantsTOAD.OF / (1 + constantsTOAD.OF) * (constantsTOAD.MaxMdot);
mdot_ipa = -thrust / constantsTOAD.MaxThrust * 1 / (1 + constantsTOAD.OF) * (constantsTOAD.MaxMdot);

xdot = [qdot;rdot;vdot; omegaBdot;mdot_lox;mdot_ipa];

% Ignore the constrained state
T = [zeros(1,12); eye(12);zeros(2,12)];
T(1:4,1:3) = [zeros(1,3); eye(3)];

A = jacobian(xdot, xn);
B = jacobian(xdot, un);

A = pinv(T) * A * T;
B = pinv(T) * B;

% Turn the jacobians into matlab functions for evaluation speed
A_fcn = matlabFunction(A, 'Vars', {xn, xn1, un});
B_fcn = matlabFunction(B, 'Vars', {xn, xn1, un});

% Extract states
x_n   = x(:, end);
x_n1  = x(:, end-1);
u_n1  = u(:, end-1);

% Jacobian Eval
A_lin = A_fcn(x_n, x_n1, u_n1);
B_lin = B_fcn(x_n, x_n1, u_n1);

end
