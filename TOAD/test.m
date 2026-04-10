% close all
% 
% x = lognrnd(log(10^-6), .1, 1000000, 1);
% for i = 1:length(x)
%     x(i) = LogNormal(10^-6, 1.2);
% end
% mean(x)
% mode(x)
% 
% histogram(x);
% xline(mean(x), color='r', LineWidth=2)
% xlim([0 10^-5])
% 
% function samples = LogNormal(target_mode, sigma)
%     % Calculate the underlying normal mean (mu) based on the target mode
%     mu_normal = log(target_mode) + (sigma^2);
% 
%     % Generate using standard normal (randn), scale, shift, and exponentiate
%     samples = exp(mu_normal + sigma * randn());
% end

x = rand(10, 1);
y = rand(10, 1);
z = rand(10, 1);

figure
scatter(x, y, 'r'); hold on
scatter(x, z, 'b');