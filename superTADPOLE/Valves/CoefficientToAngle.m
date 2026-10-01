function angle = CoefficientToAngle(valve_coefficient)
    
cv_list = [0, 0.01, 0.06, 0.21, 0.46, 0.83, 1.35, 2.04, 2.93]/(1.316e6);
angle_list = [90*0.2, 90*0.3, 90*0.4, 90*0.5, 90*0.6, 90*0.7, 90*0.8, 90*0.9, 90*1];
if valve_coefficient >= 2.93/(1.316e6)
    angle = 90;
elseif valve_coefficient < 0
    angle = 0;
else
    angle = interp1(cv_list, angle_list, valve_coefficient*1.165);
end
end