%% WindSensitivity.m

% Define parameters for wind sensitivity analysis
windSpeeds = (0:1:14) / 2.23694;
fillTime = zeros(size(windSpeeds));

for i = 1:(length(windSpeeds))
    fillTime(i) = FillTimeSim(windSpeeds(i)) / 60;
    fprintf('%.2f%% Done!\n', i/length(windSpeeds) * 100);
end

%% Plot
figure;
plot(windSpeeds * 2.23694, fillTime, 'LineWidth', 2); grid on;
xline(7.75, 'g--', 'Average Windspeed');
xline(2.91, 'r--', '95th percentile Windspeed');
xlabel('Wind Speed (mph)');
ylabel('Fill Time (minutes)');
title('Wind Sensitivity Analysis');