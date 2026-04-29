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
        NegLog = realmax; % Massive penalty if factorization fails
        return;
    end
    % Compute the Negative Log Marginal Likelihood
    alpha = L' \ (L \ Objective);
    NegLog = 0.5 * Objective' * alpha + sum(log(diag(L))) + 0.5 * NumSamples * log(2*pi);


end