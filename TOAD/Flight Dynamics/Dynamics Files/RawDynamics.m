function [xdot_kin,xdot_mass,J_tot,netTau] = RawDynamics(in1,in2)
%RawDynamics
%    [XDOT_KIN,XDOT_MASS,J_tot,netTau] = RawDynamics(IN1,IN2)


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
t22 = m_ipa.*1.125568468923699e-2;
t23 = m_lox.*8.164100096246391e-3;
t24 = m_ipa.*1.0658e-2;
t25 = m_lox.*1.0658e-2;
t26 = thrust.*2.698526887170348e-4;
t17 = t8.*2.0;
t18 = t9.*2.0;
t19 = t10.*2.0;
t21 = 1.0./t20;
mt1 = [omega1.*q1.*(-1.0./2.0)-(omega2.*q2)./2.0-(omega3.*q3)./2.0;(omega1.*q0)./2.0-(omega2.*q3)./2.0+(omega3.*q2)./2.0;(omega2.*q0)./2.0+(omega1.*q3)./2.0-(omega3.*q1)./2.0;omega1.*q2.*(-1.0./2.0)+(omega2.*q1)./2.0+(omega3.*q0)./2.0;v1;v2;v3;t21.*(t5.*thrust.*(t13-t14)+t2.*t3.*thrust.*(t12+t15)-t3.*t4.*thrust.*(t18+t19-1.0));t21.*(t5.*thrust.*(t17+t19-1.0)+t3.*t4.*thrust.*(t13+t14)-t2.*t3.*thrust.*(t11-t16))];
mt2 = [-t21.*(m_ipa.*9.801450000000001+m_lox.*9.801450000000001+t5.*thrust.*(t11+t16)+t2.*t3.*thrust.*(t17+t18-1.0)+t3.*t4.*thrust.*(t12-t15)+1.298692125e+3)];
xdot_kin = [mt1;mt2];
if nargout > 1
    t27 = -t26;
    xdot_mass = [t27;t27];
end
if nargout > 2
    t28 = t22+2.7e+1./2.0e+1;
    t29 = t23+1.7e+1./2.0e+1;
    t32 = t6.*5.067617512940958e-4;
    t33 = t7.*2.666101215261213e-4;
    t34 = t24+t25+1.5e+1;
    t30 = m_ipa.*t28;
    t31 = m_lox.*t29;
    t35 = t32+6.3948e-2;
    t36 = t33+6.3948e-2;
    t37 = (m_ipa.*t35)./1.2e+1;
    t38 = (m_lox.*t36)./1.2e+1;
    t39 = t30+t31+5.83e+2./4.0;
    t40 = t21.*t39;
    t41 = -t40;
    t42 = t40-1.1e+1./1.0e+1;
    t43 = t42.^2;
    t45 = t28+t41;
    t46 = t29+t41;
    t44 = t43.*(2.65e+2./2.0);
    t47 = t45.^2;
    t48 = t46.^2;
    t49 = m_ipa.*t47;
    t50 = m_lox.*t48;
    t51 = t37+t38+t44+t49+t50+5.0e+1;
    J_tot = reshape([t51,0.0,0.0,0.0,t51,0.0,0.0,0.0,t34],[3,3]);
end
if nargout > 3
    netTau = [-omega2.*omega3.*t34+omega2.*omega3.*t51+roll.*t3.*t4+t5.*t41.*thrust;-roll.*t5+omega1.*omega3.*t34-omega1.*omega3.*t51+t3.*t4.*t41.*thrust;roll.*t2.*t3];
end
end
