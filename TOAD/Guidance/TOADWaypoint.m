function wp = TOADWaypoint(pos, options)
    arguments
        pos                 (3,1) double
        options.MaxVel      (3,1) double = [3; 3; 6]
        options.VelFF       (3,1) double = [0; 0; 0]
        options.HoldTime    (1,1) double = nan
        options.MaxTime     (1,1) double = inf
        options.PosTol      (1,1) double = 1.0
        options.Label       (1,:) char   = ''
    end
    
    % Pre-calculate the fixed-length label into a local variable
    FIXED_LEN = 9;
    currentLen = length(options.Label);
    
    if currentLen < FIXED_LEN
        % Pad with spaces
        finalLabel = [options.Label, repmat(' ', 1, FIXED_LEN - currentLen)];
    elseif currentLen > FIXED_LEN
        % Truncate
        finalLabel = options.Label(1:FIXED_LEN);
    else
        finalLabel = options.Label;
    end
    
    % Pre-calculate the boolean flag from the inputs
    calcIsPassAndGo  = isnan(options.HoldTime);
    
    % Build the entire struct in one single assignment.
    % This guarantees the C-compiler sees the exact memory layout instantly.
    wp = struct(...
        'Position',     pos, ...
        'MaxVel',       options.MaxVel, ...
        'VelFF',        options.VelFF, ...
        'HoldTime',     options.HoldTime, ...
        'MaxTime',      options.MaxTime, ...
        'PosTol',       options.PosTol, ...
        'Label',        finalLabel, ...
        'IsPassAndGo',  calcIsPassAndGo ...
    );
end