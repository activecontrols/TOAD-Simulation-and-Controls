function ref = ref_generator(x, t)
    
    MaxAscentSpeed = 4.57;      %m/s
    MaxDescentSpeed = -4.57;    %m/s
    MaxLatSpeed = 2;            %m/s
    MaxAngularVel = 5 * pi/180; %rad/s
    HoldTimeReqs = [5, 0.2, 0.2];    % Time needed to hold at each checkpoint

    persistent timeFlag
    persistent i
    persistent timeCounter
    persistent prevTime

    if isempty(timeFlag)
        timeFlag = 999;
        i = 1;
        timeCounter = 0;
        prevTime = 0;
    end
    dt = t - prevTime;
    prevTime = t;

    TargetX1 = [5, 5, 5];
    TargetX2 = [50, 0, 0];
    TargetX5 = [0, 0, 0];

    X1Gain = 0.23;
    X1Error = TargetX1(i) - x(1);
    TargetX3 = X1Gain * X1Error;
    TargetX3 = max(TargetX3, -MaxLatSpeed);
    TargetX3 = min(TargetX3, MaxLatSpeed);

    X2Gain = 0.7;
    X2Error = TargetX2(i) - x(2);
    TargetX4 = X2Gain * X2Error;
    TargetX4 = max(TargetX4, MaxDescentSpeed);
    TargetX4 = min(TargetX4, MaxAscentSpeed);

    X5Gain = 0;
    X5Error = TargetX5(i) - x(5);
    TargetX6 = X5Gain * X5Error;
    TargetX6 = max(TargetX6, -MaxAngularVel);
    TargetX6 = min(TargetX6, MaxAngularVel);
    TargetVec = [TargetX1(i); TargetX2(i); TargetX3; TargetX4; TargetX5(i); TargetX6; 0];

    ref = x - TargetVec;

    if abs(ref(2,1)) < 2
        timeCounter = timeCounter + dt;
        % fprintf("Time flag #%d: %.2f\n", i, t);
    % elseif abs(ref(2,1)) >= 2
    %     timeFlag = 999;
    end

    if timeCounter > HoldTimeReqs(i) && i < size(HoldTimeReqs,2)
        i = i + 1;
        timeCounter = 0;
        % fprintf("Hold time complete, chaing PoF.\n");
    end