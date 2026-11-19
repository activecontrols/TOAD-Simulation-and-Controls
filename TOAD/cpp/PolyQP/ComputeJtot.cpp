// Transpiled from TOAD/Helper/ComputeJtot.m @ 5ac6bd23
#include "PolyQP.h"
#include "matlab_compat.h"

namespace PolyQP {

PolyQPStatus ComputeJtot(Real m_lox, Real m_ipa, const PolyQPConstants &constants6DoF, Mat3 &J_tot, Real &CGz) {

    const Real m_dry = constants6DoF.m_dry;
    const Real m = m_dry + m_lox + m_ipa;

    // Propellant fill height
    Real OxFluidHeight;
    Real FuFluidHeight;
    if (constants6DoF.Vehicle == 0) {
        OxFluidHeight = 0;
        FuFluidHeight = 0;
    } else {
        OxFluidHeight = (m_lox / constants6DoF.OxMass) * constants6DoF.OxHeight * Real(0.9);
        FuFluidHeight = (m_ipa / constants6DoF.FuMass) * constants6DoF.FuHeight * Real(0.9);
    }

    // Propellant inertias (about their own fluid-column centroid)
    Real J_xx = Real(1) / Real(12) * m_lox * (3 * std::pow(constants6DoF.OxRadius, Real(2)) + std::pow(OxFluidHeight, Real(2)));
    Real J_zz = Real(1) / Real(2) * m_lox * std::pow(constants6DoF.OxRadius, Real(2));
    Mat3 J_lox = Vec3(J_xx, J_xx, J_zz).asDiagonal();

    J_xx = Real(1) / Real(12) * m_ipa * (3 * std::pow(constants6DoF.FuRadius, Real(2)) + std::pow(FuFluidHeight, Real(2)));
    J_zz = Real(1) / Real(2) * m_ipa * std::pow(constants6DoF.FuRadius, Real(2));
    Mat3 J_ipa = Vec3(J_xx, J_xx, J_zz).asDiagonal();

    // Fluid fill location & instantaneous CG (engine-attachment frame)
    const Real OxFluidLocation = constants6DoF.Ox_Z + OxFluidHeight / 2;
    const Real FuFluidLocation = constants6DoF.Fu_Z + FuFluidHeight / 2;
    CGz = (m_dry * constants6DoF.rTB + m_lox * OxFluidLocation + m_ipa * FuFluidLocation) / m;

    // Parallel-axis shift of each component to the instantaneous CG
    const Real d_dry = constants6DoF.rTB - CGz;
    const Real d_lox = OxFluidLocation - CGz;
    const Real d_ipa = FuFluidLocation - CGz;

    const Mat3 diag_dry = Vec3(std::pow(d_dry, Real(2)), std::pow(d_dry, Real(2)), Real(0)).asDiagonal();
    const Mat3 J_dry = constants6DoF.J + m_dry * diag_dry;
    const Mat3 diag_lox = Vec3(std::pow(d_lox, Real(2)), std::pow(d_lox, Real(2)), Real(0)).asDiagonal();
    J_lox = J_lox + m_lox * diag_lox;
    const Mat3 diag_ipa = Vec3(std::pow(d_ipa, Real(2)), std::pow(d_ipa, Real(2)), Real(0)).asDiagonal();
    J_ipa = J_ipa + m_ipa * diag_ipa;

    // Bring everything to the instantaneous CG
    J_tot = J_dry + J_lox + J_ipa;
    return PolyQPStatus::Ok;
}

}
