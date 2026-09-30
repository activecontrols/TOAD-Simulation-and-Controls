% Optimal Gain Matrix
function K = OptGain(A, B, R, P_t)
    K = (B.'*P_t*B+R) \ B.'*P_t*A;
end