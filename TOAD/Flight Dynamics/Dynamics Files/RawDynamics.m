function [xdot_kin,xdot_mass,J_tot,netTau] = RawDynamics(in1,in2,in3,in4,MaxMdot_d)
%RawDynamics
%    [XDOT_KIN,XDOT_MASS,J_tot,netTau] = RawDynamics(IN1,IN2,IN3,IN4,MaxMdot_d)


J_d1_1 = in3(1);
J_d1_2 = in3(4);
J_d1_3 = in3(7);
J_d2_1 = in3(2);
J_d2_2 = in3(5);
J_d2_3 = in3(8);
J_d3_1 = in3(3);
J_d3_2 = in3(6);
J_d3_3 = in3(9);
TB_d3 = in4(3,:);
m_ipa = in1(15,:);
m_lox = in1(14,:);
omega1 = in1(11,:);
omega2 = in1(12,:);
omega3 = in1(13,:);
phi = in2(2,:);
q0 = in1(1,:);
q1 = in1(2,:);
q2 = in1(3,:);
q3 = in1(4,:);
roll = in2(4,:);
theta = in2(1,:);
thrust = in2(3,:);
v1 = in1(8,:);
v2 = in1(9,:);
v3 = in1(10,:);
t2 = cos(phi);
t3 = cos(theta);
t4 = sin(phi);
t5 = sin(theta);
t6 = q1.^2;
t7 = q2.^2;
t8 = q3.^2;
t9 = q0.*q1.*2.0;
t10 = q0.*q2.*2.0;
t11 = q0.*q3.*2.0;
t12 = q1.*q2.*2.0;
t13 = q1.*q3.*2.0;
t14 = q2.*q3.*2.0;
t18 = m_ipa./2.0;
t19 = m_ipa./4.0;
t20 = m_lox./2.0;
t21 = m_lox./4.0;
t22 = m_ipa+m_lox+5.1e+1./4.0e+1;
t24 = m_ipa+m_lox+3.315e-1;
t15 = t6.*2.0;
t16 = t7.*2.0;
t17 = t8.*2.0;
t23 = 1.0./t22;
mt1 = [omega1.*q1.*(-1.0./2.0)-(omega2.*q2)./2.0-(omega3.*q3)./2.0;(omega1.*q0)./2.0-(omega2.*q3)./2.0+(omega3.*q2)./2.0;(omega2.*q0)./2.0+(omega1.*q3)./2.0-(omega3.*q1)./2.0;omega1.*q2.*(-1.0./2.0)+(omega2.*q1)./2.0+(omega3.*q0)./2.0;v1;v2;v3;t23.*(t5.*thrust.*(t11-t12)+t2.*t3.*thrust.*(t10+t13)-t3.*t4.*thrust.*(t16+t17-1.0));t23.*(t5.*thrust.*(t15+t17-1.0)+t3.*t4.*thrust.*(t11+t12)-t2.*t3.*thrust.*(t9-t14))];
mt2 = [-t23.*(m_ipa.*9.801450000000001+m_lox.*9.801450000000001+t5.*thrust.*(t9+t14)+t2.*t3.*thrust.*(t15+t16-1.0)+t3.*t4.*thrust.*(t10-t13)+1.249684875e+1)];
xdot_kin = [mt1;mt2];
if nargout > 1
    xdot_mass = [0.0;0.0];
end
if nargout > 2
    t25 = J_d3_3+t18+t20+1.0./5.0e+1;
    t26 = t23.*t24;
    t27 = -t26;
    t28 = t26-1.0;
    t29 = t28.^2;
    t32 = TB_d3+t27+1.3e+1./5.0e+1;
    t30 = m_ipa.*t29;
    t31 = m_lox.*t29;
    t33 = t32.^2;
    t34 = t33.*(5.1e+1./4.0e+1);
    t35 = J_d1_1+t19+t21+t30+t31+t34+6.7e+1./1.0e+3;
    t36 = J_d2_2+t19+t21+t30+t31+t34+6.7e+1./1.0e+3;
    J_tot = reshape([t35,J_d2_1,J_d3_1,J_d1_2,t36,J_d3_2,J_d1_3,J_d2_3,t25],[3,3]);
end
if nargout > 3
    netTau = [omega1.*(J_d2_1.*omega3-J_d3_1.*omega2)+omega3.*(J_d2_3.*omega3-omega2.*t25)-omega2.*(J_d3_2.*omega2-omega3.*t36)+roll.*t3.*t4+t5.*t27.*thrust;-roll.*t5-omega2.*(J_d1_2.*omega3-J_d3_2.*omega1)-omega3.*(J_d1_3.*omega3-omega1.*t25)+omega1.*(J_d3_1.*omega1-omega3.*t35)+t3.*t4.*t27.*thrust;omega3.*(J_d1_3.*omega2-J_d2_3.*omega1)+omega2.*(J_d1_2.*omega2-omega1.*t36)-omega1.*(J_d2_1.*omega1-omega2.*t35)+roll.*t2.*t3];
end
end
