%% 1. Setup and Preallocation
num_sims = length(out);
total_error_RMS = zeros(num_sims, 1);
max_error       = zeros(num_sims, 1);
is_stable       = false(num_sims, 1);

% Extract scalar magnitudes from your input arrays for plotting
J_mag = zeros(num_sims, 1);
Lever_mag = zeros(num_sims, 1);

% Set threshold in DEGREES
stability_threshold = 200; 

%% 2. Loop Through All 100 Runs
for i = 1:num_sims
    % Calculate input magnitudes for the current run
    J_mag(i) = norm(J_d_vals{i}); 
    Lever_mag(i) = norm(TB_d_vals{i});
    
    % Extract the quaternion error data
    error_data = out(i).('AttErrorControl').Data; 
    
    % --- QUATERNION TO ANGLE CONVERSION ---
    % Assuming your Simulink block outputs Scalar-First quaternions [qw, qx, qy, qz]
    % (If your block uses Scalar-Last [qx, qy, qz, qw], change the 1 to a 4 below)
    qw = error_data(1, :, :); 
    
    % Calculate the rotational error magnitude in radians, then convert to degrees
    % Ensure qw stays within [-1, 1] to prevent complex numbers from floating point errors
    qw = max(min(qw, 1), -1); 
    error_angle_rad = 2 * acos(abs(qw));
    
    % This is your new error magnitude array in degrees
    error_mag = rad2deg(error_angle_rad); 
    
    % Calculate Total RMS Error and Max Error based on the rotation angle
    total_error_RMS(i) = sqrt(mean(error_mag.^2, 'omitnan'));
    max_error(i) = max(error_mag, [], 'omitnan');
    
    % Check Stability
    if max_error(i) < stability_threshold && ~isnan(max_error(i))
        is_stable(i) = true;
    else
        is_stable(i) = false;
    end
end

%% 3. Print Summary
disp('-----------------------------------');
disp(['Analysis Complete: ', num2str(sum(is_stable)), ' / ', num2str(num_sims), ' runs stable.']);
disp('-----------------------------------');

%% 4. Plot the 3D Results
figure('Name', 'Controller Monte Carlo Analysis', 'Color', 'w');
hold on; grid on;

% Plot stable runs (Blue Dots)
scatter3(J_mag(is_stable), Lever_mag(is_stable), total_error_RMS(is_stable), ...
    75, 'b', 'filled', 'MarkerEdgeColor', 'k');

% Plot unstable runs (Red X's)
scatter3(J_mag(~is_stable), Lever_mag(~is_stable), total_error_RMS(~is_stable), ...
    100, 'r', 'x', 'LineWidth', 2);

% Formatting
xlabel('Moment of Inertia (Norm)');
ylabel('Lever Arm (Norm)');
zlabel('Total RMS Error (Degrees)');
title('Controller Robustness: MoI vs Lever Arm');
legend('Stable', 'Unstable', 'Location', 'best');
view(45, 30);