function valve_coefficient = AngleToCoefficient(angle)
    
cv_list = [0, 0.01, 0.06, 0.21, 0.46, 0.83, 1.35, 2.04, 2.93]/(1.316e6);
angle_list = [90*0.2, 90*0.3, 90*0.4, 90*0.5, 90*0.6, 90*0.7, 90*0.8, 90*0.9, 90*1];
if angle >= 90
    valve_coefficient = 2.93;
elseif angle <= 6 || angle == 1/0
    valve_coefficient = 0;
else
    valve_coefficient = interp1(angle_list, cv_list, angle);
end
end