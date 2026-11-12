function [A0, Q0, H0] = assembleQP(r)
%% Assemble Normalised-Time QP Matrices (Paper Eq. ??)
    nc = 2 * r;
    A0 = zeros(nc, nc);
    Q0 = zeros(nc, nc);

    % Endpoint derivative map A0 at tau = 0 and tau = 1
    for k = 0:(r - 1)
        A0(k + 1, k + 1) = factorial(k);
        for i = k:(nc - 1)
            A0(r + k + 1, i + 1) = factorial(i) / factorial(i - k);
        end
    end

    % Cost Hessian Q0 for integral of squared r-th derivative over [0, 1]
    for i = r:(nc - 1)
        fi = factorial(i) / factorial(i - r);
        for j = r:(nc - 1)
            fj = factorial(j) / factorial(j - r);
            Q0(i + 1, j + 1) = (fi * fj) / (i + j - 2 * r + 1);
        end
    end

    % Exact block inversion of A0 = [D0, 0; A10, A11] to preserve PSD Hessian H0
    % Only the bottom r rows of A0^-1 (corresponding to c_r .. c_{2r-1}) enter Q0
    A10 = A0((r + 1):nc, 1:r);
    A11 = A0((r + 1):nc, (r + 1):nc);
    D0_inv = diag(1 ./ diag(A0(1:r, 1:r)));

    B1 = A11 \ eye(r);
    B0 = -B1 * (A10 * D0_inv);
    B_bot = [B0, B1];

    Q_rr = Q0((r + 1):nc, (r + 1):nc);
    H0 = B_bot.' * Q_rr * B_bot;
    H0 = 0.5 * (H0 + H0.');
end
