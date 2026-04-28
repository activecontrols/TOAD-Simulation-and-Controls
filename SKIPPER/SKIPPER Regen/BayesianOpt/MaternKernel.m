function K = MaternKernel(X1, X2, lengthScales, signalVar)
    % Scale the inputs by the length scales
    G1 = X1 ./ lengthScales;
    G2 = X2 ./ lengthScales;
    
    % Implicit expansion for pairwise distance
    % G1 becomes (N1 x 1 x D), G2 becomes (1 x N2 x D)
    G1_res = reshape(G1, size(G1, 1), 1, size(G1, 2));
    G2_res = reshape(G2, 1, size(G2, 1), size(G2, 2));
    
    % Broadcasting creates the (N1 x N2 x D) difference matrix
    DiffMatrix = G1_res - G2_res;
    r = sqrt(sum(DiffMatrix.^2, 3));
    
    % Matern 5/2 equation
    K = signalVar^2 * (1 + sqrt(5) * r + (5/3) * r.^2) .* exp(-sqrt(5) * r);
end