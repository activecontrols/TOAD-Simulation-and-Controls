%% Helper Function to Print Formatted Eigen Matrices
function print_eigen_matrix(var_name, M)
[rows, cols] = size(M);
fprintf("%s << \n", var_name);
for r = 1:rows
    fprintf("    ");
    for c = 1:cols
        val = M(r, c);
        % Clean up tiny floating-point numerical noise (< 1e-12) to exact 0
        if abs(val) < 1e-10
            val_str = "0";
        else
            val_str = sprintf("%.6e", val);
        end

        % Check if last element of the entire matrix
        if r == rows && c == cols
            fprintf("%s; \n", val_str, r);
        else
            fprintf("%s, ", val_str);
        end
    end
    if r < rows
        fprintf("\n");
    end
end
fprintf("\n");
end

% Average Gains from Trajectory Series
K_trans_nominal = squeeze(mean(constants6DoF.Traj.KTGain, 1));
K_rot_nominal   = squeeze(mean(constants6DoF.Traj.KRGain, 1));
L_Thrust_nominal = squeeze(mean(constants6DoF.Traj.LTGain, 1));
L_Att_nominal    = squeeze(mean(constants6DoF.Traj.LAGain, 1));

% Print C++ Formatted Matrices for controller_and_estimator.cpp
print_eigen_matrix("constantsASTRA.K_trans", K_trans_nominal);
print_eigen_matrix("constantsASTRA.K_rot",   K_rot_nominal);
print_eigen_matrix("constantsASTRA.L_Thrust", L_Thrust_nominal);
print_eigen_matrix("constantsASTRA.L_Att",    L_Att_nominal);