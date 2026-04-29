function [mu, sigma] = PredictGP(xstar, X_train, L, alpha, Theta)
    % Unpack parameters
    lengthScales = Theta(1:end-2)';
    signalVar = Theta(end-1);

    % Distance between training data and new point
    K_star = MaternKernel(X_train, xstar, lengthScales, signalVar);

    % Predicted mean
    mu = K_star' * alpha;

    % Predicted variance    
    v = L \ K_star;
    sigma = signalVar^2 - v' * v;
    sigma = sqrt(max(sigma, 0));
end