function K_inp = InterpolateGains(t, K_List, t_list)
%GAINS Based on checkpoints reached, find gain matrix for a given point
%   Linear interpolation between gain matrices

n_lower = find(t_list >= t, 1);
n_upper = find(t_list <= t, 1, "last");

K_low = K_List(n_lower,:);
K_up = K_List(n_upper,:);
K_inp = ((K_up - K_low) ./ (n_upper - n_lower) .* (t-t_list(n_lower))) + K_low;

end