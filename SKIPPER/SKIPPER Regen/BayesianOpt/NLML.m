function NegLog = NLML(logTheta, Geometries, Objective)

    % Unpack theta and exponentiate
    Theta = exp(logTheta);
    lengthScales = Theta(1:end-2)';
    signalVar = Theta(end-1);
    noiseVar = Theta(end);
    NumSamples = size(Geometries, 1);

    % Kernel
    K_mat = MaternKernel(Geometries, Geometries, lengthScales, signalVar);
    K = K_mat + noiseVar * eye(NumSamples);

    % Decomposition
    try
        L = chol(K, 'lower');
    catch
        % fminunc rescue: Add a tiny jitter to the diagonal
        jitter = 1e-6 * max(diag(K));
        try
            L = chol(K + jitter * eye(NumSamples), 'lower');
        catch
            NegLog = 1e6; 
            return;
        end
    end
    
    % Compute the Negative Log Marginal Likelihood
    alpha = L' \ (L \ Objective);
    NegLog = 0.5 * Objective' * alpha + sum(log(diag(L))) + 0.5 * NumSamples * log(2*pi);

end