clc; 
clear all; 
close all;
u = convertUnits;
[num, text, raw] = xlsread('ipadata.xlsx'); 
size = size(text); 
Datamatrix = zeros(size(1),14);
for i = 1:size(1)
    string = text(i,1); 
    splitString = split(string);
    numbersarray = transpose(str2double(splitString));
    Datamatrix(i,1:length(numbersarray)) = numbersarray;
end
%Cells = num2cell(num)
%writecell(Cells,'ipadata.xlsx','Sheet',3,'Range','A2:G38')
a = "done"
%coolant = 'methanol'; %coolant definition
alpha = 5.4574*10^(-9);
coolant = 'ethanol';
T_prop = linspace(273.15, 500, 100);
P_l1 = 14.7 * u.PSI2PA;
P_l2 = 100 * u.PSI2PA;
P_l3 = 200 * u.PSI2PA;
P_l4 = 300 * u.PSI2PA;
P_l5 = 400 * u.PSI2PA;
P_l6 = 1000 * u.PSI2PA;
for i = 1:length(T_prop) 
    mu_1(i) = py.CoolProp.CoolProp.PropsSI('V','T', T_prop(i), 'P', P_l1, coolant); % viscosity of bulk coolant [Pa-s]
    mu_2(i) = py.CoolProp.CoolProp.PropsSI('V','T', T_prop(i), 'P', P_l2, coolant); % viscosity of bulk coolant [Pa-s]
    mu_3(i) = py.CoolProp.CoolProp.PropsSI('V','T', T_prop(i), 'P', P_l3, coolant); % viscosity of bulk coolant [Pa-s]
    mu_4(i) = py.CoolProp.CoolProp.PropsSI('V','T', T_prop(i), 'P', P_l4, coolant); % viscosity of bulk coolant [Pa-s]
    mu_5(i) = py.CoolProp.CoolProp.PropsSI('V','T', T_prop(i), 'P', P_l5, coolant); % viscosity of bulk coolant [Pa-s]
    mu_6(i) = py.CoolProp.CoolProp.PropsSI('V','T', T_prop(i), 'P', P_l6, coolant); % viscosity of bulk coolant [Pa-s]
    k_1(i) = py.CoolProp.CoolProp.PropsSI('L', 'T', T_prop(i), 'P', P_l1, coolant);
    k_2(i) = py.CoolProp.CoolProp.PropsSI('L', 'T', T_prop(i), 'P', P_l2, coolant);
    k_3(i) = py.CoolProp.CoolProp.PropsSI('L', 'T', T_prop(i), 'P', P_l3, coolant);
    k_4(i) = py.CoolProp.CoolProp.PropsSI('L', 'T', T_prop(i), 'P', P_l4, coolant);
    k_5(i) = py.CoolProp.CoolProp.PropsSI('L', 'T', T_prop(i), 'P', P_l5, coolant);
    k_6(i) = py.CoolProp.CoolProp.PropsSI('L', 'T', T_prop(i), 'P', P_l6, coolant);
    IPA_fit(i) = 4.1936*exp(-.026*(T_prop(i)-273.15));
    mu_fit(i) = 1.6666*exp(-.018*(T_prop(i)-273.15));
    mu_fit2(i) = mu_fit(i)*exp(alpha*(P_l2-14.7));
    IPA_fit2(i) = 4.5054*exp(-.031*(T_prop(i)-273.15));
    if (T_prop(i) - 273.15) < 36.5 
        IPA_fit3(i) = 4.5054*exp(-.031*(T_prop(i)-273.15));
        
    else 
        IPA_fit3(i) = 3.3724*exp(-.023*(T_prop(i)-273.15));
    end
    IPA_fit4(i) = IPA_fit3(i)*exp(alpha*(P_l2-14.7));
    IPA_fit5(i) = IPA_fit3(i)*exp(alpha*(P_l6-14.7));
     
    [u1(i), k1(i), state1(i), bp1(i)] = IPA_transportprop(T_prop(i),P_l6)
    [u2(i),k2(i), state2(i),bp2(i)] = IPA_transportprop(T_prop(i),P_l1)
    
end
figure(1)
plot(T_prop, mu_1,T_prop, mu_2,T_prop, mu_3,T_prop, mu_4,T_prop, mu_5,T_prop, mu_6,T_prop, mu_fit./1000,T_prop, mu_fit2./1000,T_prop,IPA_fit4./1000, T_prop, IPA_fit5./1000, T_prop, IPA_fit3./1000)
legend("14.7","100","200","300","400","500","fit14.7","fit100","IPA100","IPA600", "IPA14.7")
x = [P_l1,P_l2,P_l3,P_l4,P_l5,P_l6];
y = [log(mu_1(1)),log(mu_2(1)),log(mu_3(1)),log(mu_4(1)),log(mu_5(1)),log(mu_6(1))];

figure(2)
plot(x,y)
a = (y(6) - y(1))/(x(6)-x(1)) 

figure(3)
plot(T_prop,u1)

figure(4)
plot(T_prop,k_1,T_prop,k_2,T_prop,k_3,T_prop,k_4,T_prop,k_5,T_prop,k_6,T_prop, k1,T_prop,k2)
legend("14.7","100","200","300","400","1000")

figure(5)
x2 = [P_l1,P_l2,P_l3,P_l4,P_l5,P_l6];
y2 = [k_1(1),k_2(1),k_3(1),k_4(1),k_5(1),k_6(1)];
X = x2 - 101325;
Y = (y2./.1689  ) -1; 


plot(X,Y)

% 
% if (T_search - 273.15) < 36.5
%     visc_temp = 4.5054*exp(-.031*(T_search-273.15));
% else
%     visc_temp = 3.3724*exp(-.023*(T_search-273.15));
% end
% 
% visc_tp =  visc_temp*exp(alpha*(P_search-14.7));
% visc_final = visc_tp /1000