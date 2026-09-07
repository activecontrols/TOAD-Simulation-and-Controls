function [J_tot, CGz] = ComputeJtot(m_lox, m_ipa, constants6DoF)

    m_dry = constants6DoF.m_dry;
    m = m_dry + m_lox + m_ipa;

    % Propellant fill height
    if constants6DoF.Vehicle == "ASTRAv2"
        OxFluidHeight = 0;
        FuFluidHeight = 0;
    else
        OxFluidHeight = (m_lox / constants6DoF.OxMass) * constants6DoF.OxHeight * 0.9;
        FuFluidHeight = (m_ipa / constants6DoF.FuMass) * constants6DoF.FuHeight * 0.9;
    end

    % Propellant inertias (about their own fluid-column centroid)
    J_xx = 1/12 * m_lox * (3 * constants6DoF.OxRadius^2 + OxFluidHeight^2);
    J_zz = 1/2  * m_lox * constants6DoF.OxRadius^2;
    J_lox = diag([J_xx, J_xx, J_zz]);

    J_xx = 1/12 * m_ipa * (3 * constants6DoF.FuRadius^2 + FuFluidHeight^2);
    J_zz = 1/2  * m_ipa * constants6DoF.FuRadius^2;
    J_ipa = diag([J_xx, J_xx, J_zz]);

    % Fluid fill location & instantaneous CG (engine-attachment frame)
    OxFluidLocation = constants6DoF.Ox_Z + OxFluidHeight / 2;
    FuFluidLocation = constants6DoF.Fu_Z + FuFluidHeight / 2;
    CGz = (m_dry * constants6DoF.rTB + m_lox * OxFluidLocation + m_ipa * FuFluidLocation) / m;

    % Parallel-axis shift of each component to the instantaneous CG
    d_dry = constants6DoF.rTB - CGz;
    d_lox = OxFluidLocation - CGz;
    d_ipa = FuFluidLocation - CGz;

    J_dry = constants6DoF.J + m_dry * diag([d_dry^2, d_dry^2, 0]);
    J_lox = J_lox + m_lox * diag([d_lox^2, d_lox^2, 0]);
    J_ipa = J_ipa + m_ipa * diag([d_ipa^2, d_ipa^2, 0]);

    % Bring everything to the instantaneous CG
    J_tot = J_dry + J_lox + J_ipa;
end