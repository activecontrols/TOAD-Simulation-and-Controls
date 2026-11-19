function u = inputfcn(K, x, t, thrust_max, AscentVEL)
    
    maxU = [thrust_max;0.122];
    minU = [thrust_max * 0.5;-0.122];
    MaxDeltaThrottle = 0.6 / 1.2;  %throttle change per sec

    % Input Saturation
    if nargin() <= 4
        u = -K*ref_generator(x, t);
    else
         u = -K*ref_generator(x, t, AscentVEL);
    end
    u = max(minU, u);
    u = min(maxU, u);

    persistent prevInput 
    persistent prevTime
    if isempty(prevInput)
        prevInput = [u(1); u(2)];
        prevTime = 0;
    end
    dt = t - prevTime;

    % Thrust Rate Limiter
    MaxAllowThrottle = prevInput(1) + MaxDeltaThrottle * dt * thrust_max;
    MinAllowThrottle = prevInput(1) - MaxDeltaThrottle * dt * thrust_max;

    u(1) = max(u(1), MinAllowThrottle);
    u(1) = min(u(1), MaxAllowThrottle);

end