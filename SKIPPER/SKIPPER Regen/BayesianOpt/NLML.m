function NegLog = NLML(logTheta, Geometries, Objective)

    % Unpack theta and exponentiate
    Theta = exp(logTheta);
    lengthScales = Theta(1:end-2)';
    signalVar = Theta(end-1);
    noiseVar = Theta(end);
    NumSamples = size(Geometries, 1);

    % Kernel
    K_mat = MaternKernel(Geometries, Geometries, lengthScales, signalVar);
    
    % Force jitter immediately to prevent RCOND during the '\' division
    jitter = 1e-5 * max(diag(K_mat));
    K = K_mat + (noiseVar + jitter) * eye(NumSamples);
    
    % Decomposition without the try/catch rescue
    try
        L = chol(K, 'lower');
    catch
        NegLog = 1e6;
        return;
    end
    
    % Compute the Negative Log Marginal Likelihood
    alpha = L' \ (L \ Objective);
    NegLog = 0.5 * Objective' * alpha + sum(log(diag(L))) + 0.5 * NumSamples * log(2*pi);

end