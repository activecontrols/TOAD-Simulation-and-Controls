function [J_tot, CGz] = ComputeJtot(m_lox, m_ipa, constantsTOAD)

    m_dry = constantsTOAD.m_dry;
    m = m_dry + m_lox + m_ipa;

    % Propellant fill height
    OxFluidHeight = (m_lox / constantsTOAD.OxMass) * constantsTOAD.OxHeight * 0.9;
    FuFluidHeight = (m_ipa / constantsTOAD.FuMass) * constantsTOAD.FuHeight * 0.9;

    % Propellant inertias (about their own fluid-column centroid)
    J_xx = 1/12 * m_lox * (3 * constantsTOAD.OxRadius^2 + OxFluidHeight^2);
    J_zz = 1/2  * m_lox * constantsTOAD.OxRadius^2;
    J_lox = diag([J_xx, J_xx, J_zz]);

    J_xx = 1/12 * m_ipa * (3 * constantsTOAD.FuRadius^2 + FuFluidHeight^2);
    J_zz = 1/2  * m_ipa * constantsTOAD.FuRadius^2;
    J_ipa = diag([J_xx, J_xx, J_zz]);

    % Fluid fill location & instantaneous CG (engine-attachment frame)
    OxFluidLocation = constantsTOAD.Ox_Z + OxFluidHeight / 2;
    FuFluidLocation = constantsTOAD.Fu_Z + FuFluidHeight / 2;
    CGz = (m_dry * constantsTOAD.rTB + m_lox * OxFluidLocation + m_ipa * FuFluidLocation) / m;

    % Parallel-axis shift of each component to the instantaneous CG
    d_dry = constantsTOAD.rTB - CGz;
    d_lox = OxFluidLocation - CGz;
    d_ipa = FuFluidLocation - CGz;

    J_dry = constantsTOAD.J + m_dry * diag([d_dry^2, d_dry^2, 0]);
    J_lox = J_lox + m_lox * diag([d_lox^2, d_lox^2, 0]);
    J_ipa = J_ipa + m_ipa * diag([d_ipa^2, d_ipa^2, 0]);

    % Bring everything to the instantaneous CG
    J_tot = J_dry + J_lox + J_ipa;
end