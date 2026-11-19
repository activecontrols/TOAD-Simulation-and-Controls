% Kinematic mappings, replaces the T approximation for small angles.
% Derivation for the Xi matrix comes from q = q_nom * delq
function Xi = XiMat(q)
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    Xi = [-q1, -q2, -q3;
           q0, -q3,  q2;
           q3,  q0, -q1;
          -q2,  q1,  q0];
end