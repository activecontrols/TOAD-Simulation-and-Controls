function K_inp = InterpolateGains(t, K_List, t_max)
%GAINS Based on checkpoints reached, find gain matrix for a given point
%   Linear interpolation between gain matrices

n_lower = floor(t/t_max*size(K_List,1));
n_upper = ceil(t/t_max*size(K_List,1));

K_low = K_List(n_lower,:);
K_up = K_List(n_upper,:);
K_inp = ((K_up - K_low) ./ (n_upper - n_lower) .* (t/t_max*size(K_List,1)-n_lower)) + K_low;

end