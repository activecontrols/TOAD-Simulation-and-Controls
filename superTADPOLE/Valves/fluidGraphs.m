u_o = 1.231e-06;
u_f = u_o;
u = [u_o; u_f];
init_cond = [0.4; 0.4; 1e+05; 3.447e+6; 3.447e+6];
[t, y] = ode45(@(t, y) fluidMechanics(t, y, u), [0, 5], init_cond);


% Graphs

% Font size
set(0, 'DefaultAxesFontSize', 18); 
set(0, 'DefaultTextFontSize', 18);

% figure (1)
% plot(t, y, 'LineWidth', 1.5); % Plots all 5 columns automatically
% title('Fluid Mechanics System State Variables');
% xlabel('Time [s]');
% ylabel('Values');
% legend('Oxygen Flow', 'Fuel Flow', 'Chamber Pressure', 'O2 Pressure', 'Fuel Pressure');
% grid on;
% 
% figure (2)
% plot(t, y(:,1), 'b', 'LineWidth', 1.5)
% xlabel('Time [s]');
% ylabel('Mass flow [kg/s]');
% title('Oxygen Mass Flow Rate');
% grid on;
% 
% figure (3)
% plot(t, y(:,2), 'r', 'LineWidth', 1.5)
% xlabel('Time [s]');
% ylabel('Mass flow [kg/s]');
% title('Fuel Mass Flow Rate');
% grid on;
% 
% figure (4)
% plot(t, y(:,3) / 1e6, 'k', 'LineWidth', 1.5)
% xlabel('Time [s]');
% ylabel('Chamber Pressure [MPa]');
% title('Chamber Pressure');
% grid on;
% 
% figure (5)
% plot(t, y(:,4) / 1e6, 'b--', 'LineWidth', 1.5)
% xlabel('Time [s]');
% ylabel('Oxygen Line Pressure [MPa]');
% title('Oxygen Line Pressure');
% grid on;
% 
% figure (6)
% plot(t, y(:,5) / 1e6, 'r--', 'LineWidth', 1.5)
% xlabel('Time [s]');
% ylabel('Fuel Line Pressure [MPa]');
% title('Fuel Line Pressure');
% grid on;
% 
