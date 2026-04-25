function u = inputfcn(K, x, t)
    
    thrust_max = 2446.52;   %N
    maxU = [thrust_max;0.122];
    minU = [thrust_max * 0.4;-0.122];
    MaxDeltaThrottle = 0.6 / 1.2;  %throttle change per sec

    % Input Saturation
    u = -K*ref_generator(x, t);
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