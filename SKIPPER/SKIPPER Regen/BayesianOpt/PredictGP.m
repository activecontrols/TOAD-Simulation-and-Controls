function [mu, sigma] = PredictGP(xstar, X_train, y_train, Theta)
    % Predicts a normal distribution given the training data. Used to
    % select next evaluation spot
    % Unpack parameters
    lengthScales = Theta(1:8)';
    signalVar = Theta(9);
    noiseVar = Theta(10);
    NumSamples = size(X_train, 1);

    K_train = MaternKernel(X_train, X_train, lengthScales, signalVar);
    K_train = K_train + noiseVar * eye(NumSamples);
    L = chol(K_train, 'lower');

    % Distance between training data and new point
    K_star = MaternKernel(X_train, xstar, lengthScales, signalVar);
    alpha = L' \ (L \ y_train);

    % Predicted mean
    mu = K_star' * alpha;

    % Predicted variance    
    v = L \ K_star;
    sigma = signalVar^2 - v' * v;
    sigma = sqrt(max(sigma, 0));
end