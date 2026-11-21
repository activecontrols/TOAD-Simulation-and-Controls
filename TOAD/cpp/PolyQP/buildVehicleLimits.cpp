// Transpiled from TOAD/Guidance/PolynomialQP/buildVehicleLimits.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

PolyQPStatus buildVehicleLimits(LimitsT &Limits) {
    // Build Vehicle Limits and Flatness Parameters
    // TRANSPILE-FLAG F38 (superseded: ASTRA carries no propellant, FIX-F29): the MATLAB Opts argument (InitOxMass/InitFuMass) and Constants argument are not carried; the constants are ASTRA_PolyQP_Constants.h.
    // Fixed planner constants
    const Real ThrustMargin = Real(0.10);
    const Real GimbalMargin = Real(0.10);
    const Real MaxGimbalAngle = deg2rad(10);
    const Real MaxGimbalRate = deg2rad(30);
    const Real MaxRollRate = 5;
    const Real MaxTilt = deg2rad(30);
    const Real MinSegTime = Real(0.2);
    const Real DtSample = Real(0.1);

    // Vehicle selection
    // ASTRA-ONLY: the vehicle 0 (ASTRA) values; the vehicle 1 (TOAD) branch was removed.
    const Real MT = ASTRA::MaxThrust;
    const Real MaxThrustRate = 30;
    const Real minThrottle = Real(0.25);
    const Real m_lox0 = 0;
    const Real m_ipa0 = 0;
    const Vec3 defAccel = Vec3(Real(2.5), Real(2.5), Real(3.5));

    // Initial mass, inertia, and center of percussion offset
    const Real m_init = ASTRA::m_dry + m_lox0 + m_ipa0;
    Mat3 J_init;
    Real CGz_init;
    ComputeJtot(m_lox0, m_ipa0, J_init, CGz_init);
    const Real rho = M1(J_init, 1, 1) / (m_init * CGz_init);

    // Effective actuator bands with controller margins
    const Real T_lo = (minThrottle + ThrustMargin) * MT;
    const Real T_hi = (Real(1.00) - ThrustMargin) * MT;
    const Real gimbalLim = (Real(1.00) - GimbalMargin) * MaxGimbalAngle;
    const Real rollLim = (Real(1.00) - ThrustMargin) * MaxRollRate;

    Limits.HasVehicle = true;
    Limits.g = ASTRA::g;
    Limits.m_init = m_init;
    Limits.m_dry = ASTRA::m_dry;
    Limits.MaxThrust = MT;
    Limits.rTB = ASTRA::rTB;
    Limits.CGz_init = CGz_init;
    Limits.J_init = J_init;
    Limits.Rho = rho;
    Limits.T_lo = T_lo;
    Limits.T_hi = T_hi;
    Limits.GimbalLimit = gimbalLim;
    Limits.GimbalRateLim = MaxGimbalRate;
    Limits.ThrustRateLim = MaxThrustRate;
    Limits.RollLimit = rollLim;
    Limits.RollRateLim = MaxRollRate;
    Limits.MaxTilt = MaxTilt;
    Limits.DefaultMaxAccel = defAccel;
    Limits.MinSegTime = MinSegTime;
    Limits.DtSample = DtSample;
    return PolyQPStatus::Ok;
}

}
