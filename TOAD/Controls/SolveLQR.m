%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The following function solves an infinite time horizon LQR problem by
% computing the unique solution P to the Continious-Time Algebraic Ricatti
% Equation (CARE) using the Hamiltonian matrix. This solution is then used
% to compute the optimal gain matrix for controls.
%
% By: Pablo Plata   -   10/22/25
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K = SolveLQR(A, B, Q, R)

% Define the Hamiltonian matrix H (size 2n x 2n)
R_inv = inv(R);
H = [A  -B * R_inv * B';
    -Q  -A'];

% Find eigenvalues of the Hamiltonian
%   (eig(H) function returns a matrix of
%   eigenvectors of H and a diagonal matrix of eigenvalues of H.)
[EVec_H, EVal_H] = eig(H);
EVal_H = diag(EVal_H);
NegEigen = real(EVal_H) < 0;

% Denote the V matrix, containing the eigenvectors corresponding to stable
% eigenvalues. Then decompose into V1 (first n rows) and V2 (last n rows)
V = EVec_H(:, NegEigen);
V1 = V(1:size(A,1), :);
V2 = V(size(A,1)+1:end,:);

% Compute the solution P to the Algebraic Riccati Equation, and ensure real
% coefficients and symmetry
P = V2 / V1;
P = real(P);
P = (P + P') / 2;

% Compute the gain matrix for the stabilizing solution
K = R_inv * B' * P;



