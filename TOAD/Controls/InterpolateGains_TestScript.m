function [U, State_ERR] = TOAD_Control_FCN_VLQR(PosTarget, X, constantsTOAD, t, MaxVel, VelFF, HDGRef, GND, K_List)

% Time Counter
persistent lastT VelErrorI AttErrorI lastAttError lastAccelZ
if isempty(lastT)
    lastT = 0;
    VelErrorI = zeros(3,1);
    AttErrorI = zeros(3,1);
    lastAttError = zeros(3,1);
end
dT = t - lastT;
lastT = t;

% not done, do with this what you will pablo