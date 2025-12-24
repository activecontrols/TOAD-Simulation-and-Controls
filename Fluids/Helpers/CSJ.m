function J = CSJ(Func, X0)
    % Computes Jacobian using Complex Step Differentiation
    % Accuracy is approx 1e-16 (Machine Precision)
    % Requires Func to be written using only complex-analytic operations
    % (e.g. no abs(), max(), min(), <, >) inside the smooth path.
    
    h = 1e-100; % Imaginary step size
    N = length(X0);
    J = zeros(N, N);
    
    % We do not need F0 for the derivative, only the perturbed steps
    for i = 1:N
        X_perturb = X0;
        X_perturb(i) = X0(i) + 1i * h;
        
        F_new = Func(X_perturb);
        
        % The derivative is simply the imaginary part divided by h
        J(:, i) = imag(F_new) / h;
    end
end