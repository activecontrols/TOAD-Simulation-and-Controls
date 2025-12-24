function J = NumericalJacobian(Func, X)
    % Calculates Jacobian J via finite differences
    % J(:,i) = (F(X + dx) - F(X)) / dx
    
    epsilon = 1e-8; % Perturbation size
    N = length(X);
    F0 = Func(X);
    J = zeros(N, N);
    
    for i = 1:N
        X_perturb = X;
        % Handle magnitude differences (perturb large pressures and small massflows appropriately)
        perturbation = epsilon * max(abs(X(i)), 1); 
        X_perturb(i) = X(i) + perturbation;
        
        F_new = Func(X_perturb);
        J(:, i) = (F_new - F0) / perturbation;
    end
end