%% Test Script for TVLQR Gain Generation 
% Call and load in parameters
if ~exist("constantsTOAD")
    LoadTOADSim;
end
N = 100;
dT = 0.5;

% Base states
x_0 = [1; zeros(12,1); constantsTOAD.OxMass; constantsTOAD.FuMass];
u_0 = [0, 0, constantsTOAD.m_wet * constantsTOAD.g, 0]';

% Repeat 
List.States = repmat(x_0, [1, N]);
List.Inputs = repmat(u_0, [1, N]);

% Test Call
[K_List, SanityCheck] = GainGenerator(List.States, List.Inputs, dT, constantsTOAD);
K_List(1, :, :) = K_List(2, :, :);

%% Plots
N = size(K_List, 1);
error_norm = zeros(N, 1);
K_converged = zeros(size(SanityCheck));

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
