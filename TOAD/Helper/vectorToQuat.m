function q = vectorToQuat(v)
%% VECTORTOQUAT  Unit quaternion rotating body Z-axis [0; 0; 1] to target vector v.
%
%   q = vectorToQuat(v) computes the unit quaternion [w; x; y; z] aligning
%   the body z-axis e3 = [0; 0; 1] with the 3x1 target vector v.
%
%   Inputs:
%       v - 3x1 or 1x3 target vector
%
%   Outputs:
%       q - 4x1 unit quaternion [w; x; y; z]

    v = v(:);
    nv = norm(v);
    if nv < 1e-9
        q = [1; 0; 0; 0];
        return;
    end
    b3 = v / nv;
    e3 = [0; 0; 1];
    d = dot(e3, b3);
    if d > 0.999999
        q = [1; 0; 0; 0];
    elseif d < -0.999999
        q = [0; 1; 0; 0];
    else
        ax = cross(e3, b3);
        q = [1 + d; ax(1); ax(2); ax(3)];
        q = q / norm(q);
    end
end
