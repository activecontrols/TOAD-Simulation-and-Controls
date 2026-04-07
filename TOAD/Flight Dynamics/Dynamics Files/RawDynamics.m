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
TB_d1 = in4(1,:);
TB_d2 = in4(2,:);
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
t6 = m_ipa.^2;
t7 = m_lox.^2;
t8 = q1.^2;
t9 = q2.^2;
t10 = q3.^2;
t11 = q0.*q1.*2.0;
t12 = q0.*q2.*2.0;
t13 = q0.*q3.*2.0;
t14 = q1.*q2.*2.0;
t15 = q1.*q3.*2.0;
t16 = q2.*q3.*2.0;
t20 = m_ipa+m_lox+2.65e+2./2.0;
t22 = MaxMdot_d+1.3204;
t23 = m_ipa.*1.125568468923699e-2;
t24 = m_lox.*8.164100096246391e-3;
t25 = m_ipa.*1.0658e-2;
t26 = m_lox.*1.0658e-2;
t17 = t8.*2.0;
t18 = t9.*2.0;
t19 = t10.*2.0;
t21 = 1.0./t20;
mt1 = [omega1.*q1.*(-1.0./2.0)-(omega2.*q2)./2.0-(omega3.*q3)./2.0;(omega1.*q0)./2.0-(omega2.*q3)./2.0+(omega3.*q2)./2.0;(omega2.*q0)./2.0+(omega1.*q3)./2.0-(omega3.*q1)./2.0;omega1.*q2.*(-1.0./2.0)+(omega2.*q1)./2.0+(omega3.*q0)./2.0;v1;v2;v3;t21.*(t5.*thrust.*(t13-t14)+t2.*t3.*thrust.*(t12+t15)-t3.*t4.*thrust.*(t18+t19-1.0));t21.*(t5.*thrust.*(t17+t19-1.0)+t3.*t4.*thrust.*(t13+t14)-t2.*t3.*thrust.*(t11-t16))];
mt2 = [-t21.*(m_ipa.*9.801450000000001+m_lox.*9.801450000000001+t5.*thrust.*(t11+t16)+t2.*t3.*thrust.*(t17+t18-1.0)+t3.*t4.*thrust.*(t12-t15)+1.298692125e+3)];
xdot_kin = [mt1;mt2];
if nargout > 1
    t27 = t23+2.7e+1./2.0e+1;
    t28 = t24+1.7e+1./2.0e+1;
    t31 = t22.*thrust.*2.043719242025408e-4;
    t32 = t6.*5.067617512940958e-4;
    t34 = t7.*2.666101215261213e-4;
    t35 = J_d3_3+t25+t26+1.5e+1;
    t29 = m_ipa.*t27;
    t30 = m_lox.*t28;
    t33 = -t31;
    xdot_mass = [t33;t33];
end
if nargout > 2
    t36 = t32+6.3948e-2;
    t37 = t34+6.3948e-2;
    t38 = (m_ipa.*t36)./1.2e+1;
    t39 = (m_lox.*t37)./1.2e+1;
    t40 = t29+t30+5.83e+2./4.0;
    t41 = t21.*t40;
    t42 = -t41;
    t43 = TB_d3+t42;
    t47 = t27+t42;
    t48 = t28+t42;
    t44 = t43+1.1e+1./1.0e+1;
    t49 = t47.^2;
    t50 = t48.^2;
    t45 = t44.^2;
    t51 = m_ipa.*t49;
    t52 = m_lox.*t50;
    t46 = t45.*(2.65e+2./2.0);
    t53 = J_d1_1+t38+t39+t46+t51+t52+5.0e+1;
    t54 = J_d2_2+t38+t39+t46+t51+t52+5.0e+1;
    J_tot = reshape([t53,J_d2_1,J_d3_1,J_d1_2,t54,J_d3_2,J_d1_3,J_d2_3,t35],[3,3]);
end
if nargout > 3
    netTau = [omega1.*(J_d2_1.*omega3-J_d3_1.*omega2)+omega3.*(J_d2_3.*omega3-omega2.*t35)-omega2.*(J_d3_2.*omega2-omega3.*t54)+roll.*t3.*t4+t5.*t43.*thrust+TB_d2.*t2.*t3.*thrust;-roll.*t5-omega2.*(J_d1_2.*omega3-J_d3_2.*omega1)-omega3.*(J_d1_3.*omega3-omega1.*t35)+omega1.*(J_d3_1.*omega1-omega3.*t53)-TB_d1.*t2.*t3.*thrust+t3.*t4.*t43.*thrust;omega3.*(J_d1_3.*omega2-J_d2_3.*omega1)+omega2.*(J_d1_2.*omega2-omega1.*t54)-omega1.*(J_d2_1.*omega1-omega2.*t53)-TB_d1.*t5.*thrust+roll.*t2.*t3-TB_d2.*t3.*t4.*thrust];
end
end
