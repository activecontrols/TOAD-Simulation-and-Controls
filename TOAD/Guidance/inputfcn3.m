function u = inputfcn3(u, t) 
    u = u';
    thrust_max = 1.5 * 9.8;   %N
    maxU = [pi/24, pi/24, thrust_max, thrust_max * 10];
    minU = [-pi/24, -pi/24, thrust_max * 0.4, -thrust_max * 10];
    MaxDeltaThrottle = 2;  %throttle change per sec
    MaxDeltaAngle = pi/2;       %rad/s

    % Input Saturation
    u = max(minU, u);
    u = min(maxU, u);


    % Avoid rate limiter due to testing Actuator Models

    % persistent prevInput 
    % persistent prevTime
    % if isempty(prevInput)
    %     prevInput = u;
    %     prevTime = 0;
    % end
    % dt = t - prevTime;
    % prevTime = t;
    % 
    % % Input Rate Limiter
    % MaxAllowThrottle = prevInput(3) + MaxDeltaThrottle * dt * thrust_max;
    % MinAllowThrottle = prevInput(3) - MaxDeltaThrottle * dt * thrust_max;
    % 
    % MaxAllowGimbal = prevInput(1:2) + MaxDeltaAngle * dt;
    % MinAllowGimbal = prevInput(1:2) - MaxDeltaAngle * dt;
    % 
    % % Output
    % u(3) = max(u(3), MinAllowThrottle);
    % u(3) = min(u(3), MaxAllowThrottle);
    % 
    % u(1:2) = max(MinAllowGimbal, u(1:2));
    % u(1:2) = min(MaxAllowGimbal, u(1:2));
    % prevInput = u;
end