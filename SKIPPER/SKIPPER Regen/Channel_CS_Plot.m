function Channel_CS_Plot(N,t_w,h_c,w_c,R,figureNum,place) % plots the channel cross section for visualization

if figureNum < 1 % just creates a new figure 
    figure('Name', "Cross Section at" + place)
else             % rewrites over a specific figure for continuous updating
    figure(figureNum)
    clf(figureNum)
end

% --- Parameters - left over from before it was a function 
for j = 1:3
w_fin = 2*(R(j)+t_w(j)).*pi./N - w_c(j); %fin width

% --- Derived geometry ---
r_outer = R(j) + t_w(j) + h_c(j);   % outer wall radius
r_inner = R(j) + t_w(j);  % inner channel wall radius

% Convert arc length to angular width
theta_c = w_c(j) ./ r_outer;       % angular width of a channel
theta_step = 2.*pi ./ N;         % total angle per channel+fin
theta_f = theta_step - theta_c;% angular width of one fin

% --- Plot setup ---
subplot(1,length(t_w),j)
hold on; axis equal; box on;
xlabel('x [m]'); ylabel('y [m]');
title("Cross Section at" + place(j));


% --- Draw channels and fins ---
for i = 1:N
    % --- Channel (curved wedge) ---
    theta_start = (i-1) * theta_step;          % start of this channel
    theta_end   = theta_start + theta_c;       % end of channel

    % define the curved region of the channel
    theta_channel = linspace(theta_start, theta_end, 50);
    x_outer = r_outer * cos(theta_channel);
    y_outer = r_outer * sin(theta_channel);
    x_inner = r_inner * cos(fliplr(theta_channel));
    y_inner = r_inner * sin(fliplr(theta_channel));
    
    % fill curved annular patch (blue)
    fill([x_outer x_inner], [y_outer y_inner], 'b', 'EdgeColor', 'none');

    % --- Fin (curved wedge) ---
    theta_start_f = theta_end;
    theta_end_f   = theta_start + theta_step;

    theta_fin = linspace(theta_start_f, theta_end_f, 50);
    x_outer_f = r_outer * cos(theta_fin);
    y_outer_f = r_outer * sin(theta_fin);
    x_inner_f = r_inner * cos(fliplr(theta_fin));
    y_inner_f = r_inner * sin(fliplr(theta_fin));

    % fill curved annular patch (gray)
    fill([x_outer_f x_inner_f], [y_outer_f y_inner_f], [0.8 0.8 0.8], 'EdgeColor', 'none');
end

% --- Draw inner and outer walls ---
theta = linspace(0, 2*pi, 600);
fill([R(j)*cos(theta) r_inner*cos(fliplr(theta))], [R(j)*sin(theta) r_inner*sin(fliplr(theta))], [.8 .8 .8], 'EdgeColor', 'none'); %inner wall
fill([r_outer*cos(theta) (r_outer+.001)*cos(fliplr(theta))], [r_outer*sin(theta) (r_outer+.001)*sin(fliplr(theta))], [.8 .8 .8], 'EdgeColor', 'none'); %outer wall
plot(R(j)*cos(theta), R(j)*sin(theta), 'k', 'LineWidth', 2);  % chamber wall
plot(r_inner*cos(theta), r_inner*sin(theta), 'k', 'LineWidth', 2);  % inner channel wall
plot(r_outer*cos(theta), r_outer*sin(theta), 'k', 'LineWidth', 2);  % outer channel wall
plot((r_outer+.001)*cos(theta), (r_outer+0.001)*sin(theta), 'k', 'LineWidth', 2);  % outer channel wall
ylim([-r_outer-.001, r_outer+.001])
xlim([-r_outer-.001, r_outer+.001])

% Add text box with dimensions
text(-0.01,0,sprintf("Number of Channels: %d \n Channel Width: %.1f mm \nChannel Height: %.1f mm \n          Fin Width: %.1f mm \n Wall Thickness: %.1f mm",N,w_c(j)*1000,h_c(j)*1000,w_fin*1000,t_w(j)*1000));
end


