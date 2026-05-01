function [mu, sigma] = PredictGP(xstar, X_train, L, alpha, Theta)

    % Extract params
    lengthScales = Theta(1:end-2)';
    signalVar    = Theta(end-1);
 
    % K_star - N_train x N_test
    K_star = MaternKernel(X_train, xstar, lengthScales, signalVar);
 
    % Posterior mean
    mu = K_star' * alpha;
 
    % Posterior stdev
    v     = L \ K_star;
    var_v = signalVar^2 - sum(v .^ 2, 1)';
    sigma = sqrt(max(var_v, 0));
end