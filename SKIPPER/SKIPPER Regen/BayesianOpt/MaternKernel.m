function K = MaternKernel(X1, X2, lengthScales, signalVar)
     
    G1 = X1 ./ lengthScales;   % N1 x D
 
    if isequal(X1, X2)
        % Symmetric Case
        % One matrix multiply, reuse diagonal as squared norms
        G1G1 = G1 * G1';                       % N1 x N1
        sq   = diag(G1G1);                     % N1 x 1
        D2   = sq + sq' - 2 * G1G1;            % N1 x N1
    else
        % General case
        G2  = X2 ./ lengthScales;              % N2 x D
        sq1 = sum(G1 .^ 2, 2);                 % N1 x 1
        sq2 = sum(G2 .^ 2, 2);                 % N2 x 1
        D2  = sq1 + sq2' - 2 * (G1 * G2');     % N1 x N2
    end
 
    % Clamp tiny negatives that arise from floating-point cancellation
    r   = sqrt(max(D2, 0));
 
    % Matern Kernel
    sr5 = sqrt(5) * r;
    K   = signalVar^2 * (1 + sr5 + (5/3) * r .^ 2) .* exp(-sr5);
end