function wp = TOADWaypoint(pos, options)
%% TOADWaypoint Builder
    arguments
        pos                    (3,1) double
        options.VelBC          (3,1) double = [nan; nan; nan]
        options.AccelBC        (3,1) double = [nan; nan; nan]
        options.JerkBC         (3,1) double = [nan; nan; nan]
        options.HoldTime       (1,1) double = nan
        options.MaxVel         (3,1) double = [4; 4; 5]
        options.MaxAccel       (3,1) double = [inf; inf; inf]
        options.StraightLine   (1,1) logical = false
        options.HDGRef         (3,1) double = [1; 0; 0]
        options.Label          (1,:) char   = ''
        options.B3BC           (3,1) double = [nan; nan; nan] % FIX-F03: thrust axis b3 at this knot for the P-space BC (NaN: b3 = e3)
        options.B3DotBC        (3,1) double = [nan; nan; nan] % FIX-F03: db3/dt at this knot for the P-space velocity BC (NaN: not applied)
    end

    % Fixed-length 9-character label for codegen
    FIXED_LEN = 9;
    currentLen = length(options.Label);
    if currentLen < FIXED_LEN
        finalLabel = [options.Label, repmat(' ', 1, FIXED_LEN - currentLen)];
    elseif currentLen > FIXED_LEN
        finalLabel = options.Label(1:FIXED_LEN);
    else
        finalLabel = options.Label;
    end

    % Rest knot check
    isPassAndGo = isnan(options.HoldTime);

    % Resolve rest knot derivative defaults
    velBC   = options.VelBC;
    accelBC = options.AccelBC;
    jerkBC  = options.JerkBC;
    if ~isPassAndGo
        for i = 1:3
            if isnan(velBC(i)),   velBC(i)   = 0; end
            if isnan(accelBC(i)), accelBC(i) = 0; end
            if isnan(jerkBC(i)),  jerkBC(i)  = 0; end
        end
    end

    % Construct fixed-layout waypoint struct
    wp = struct(...
        'PosBC',        pos, ...
        'VelBC',        velBC, ...
        'AccelBC',      accelBC, ...
        'JerkBC',       jerkBC, ...
        'HoldTime',     options.HoldTime, ...
        'MaxVel',       options.MaxVel, ...
        'MaxAccel',     options.MaxAccel, ...
        'StraightLine', options.StraightLine, ...
        'HDGRef',       options.HDGRef, ...
        'Label',        finalLabel, ...
        'IsPassAndGo',  isPassAndGo, ...
        'B3BC',         options.B3BC, ... % FIX-F03: optional knot attitude for the P-space BC
        'B3DotBC',      options.B3DotBC ... % FIX-F03: optional knot attitude rate for the P-space BC
    );
end
