function Output = ExpMovingAvg(Input, Last, Alpha)
    % Basic Exponential Moving Average Implementation to pre-process
    % measurements.

    % Filter difference equation
    Output = Alpha * Input + (1 - Alpha) * Last;
end
