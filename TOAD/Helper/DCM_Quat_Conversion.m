% Convert a DCM to a quaternion robustly
function q = DCM_Quat_Conversion(R)
% Extract trace
tr = R(1,1) + R(2,2) + R(3,3);

% Pre-allocate output
q = zeros(4,1);

% CASE 1: Standard Case (Trace > 0)
% This works for most orientations (small angles, < 90 deg tilt)
if tr > 0
    S = sqrt(tr + 1.0) * 2; % S = 4 * qw
    q(1) = 0.25 * S;        
    q(2) = (R(3,2) - R(2,3)) / S; 
    q(3) = (R(1,3) - R(3,1)) / S; 
    q(4) = (R(2,1) - R(1,2)) / S; 
    
% CASE 2: Singularity Avoidance (Trace <= 0)
% We must find the largest diagonal element to avoid dividing by zero.
else
    if (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
        % Column 1 (X) is dominant
        S = sqrt(1.0 + R(1,1) - R(2,2) - R(3,3)) * 2; 
        q(1) = (R(3,2) - R(2,3)) / S; 
        q(2) = 0.25 * S;        
        q(3) = (R(1,2) + R(2,1)) / S; 
        q(4) = (R(1,3) + R(3,1)) / S; 
        
    elseif (R(2,2) > R(3,3))
        % Column 2 (Y) is dominant
        S = sqrt(1.0 + R(2,2) - R(1,1) - R(3,3)) * 2; 
        q(1) = (R(1,3) - R(3,1)) / S; 
        q(2) = (R(1,2) + R(2,1)) / S; 
        q(3) = 0.25 * S;        
        q(4) = (R(2,3) + R(3,2)) / S; 
        
    else
        % Column 3 (Z) is dominant
        S = sqrt(1.0 + R(3,3) - R(1,1) - R(2,2)) * 2;
        q(1) = (R(2,1) - R(1,2)) / S; 
        q(2) = (R(1,3) + R(3,1)) / S; 
        q(3) = (R(2,3) + R(3,2)) / S; 
        q(4) = 0.25 * S;                 
    end
end

% Enforce normalization to correct numerical drift
q = q / norm(q);