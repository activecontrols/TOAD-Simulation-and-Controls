// Transpiled from TOAD/Guidance/PolynomialQP/buildVehicleLimits.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

static_assert(PolyQPCapacity::m == ASTRA::m_dry && PolyQPCapacity::g == ASTRA::g && PolyQPCapacity::MaxThrust == ASTRA::MaxThrust,
              "STATIC: the grid capacity derivation in PolyQP_config.h must use the ASTRA constants");

PolyQPStatus buildVehicleLimits(LimitsT &Limits) {
    // Build Vehicle Limits and Flatness Parameters
    // TRANSPILE-FLAG F38 (superseded: ASTRA carries no propellant, FIX-F29): the MATLAB Opts argument (InitOxMass/InitFuMass) and Constants argument are not carried; the constants are ASTRA_PolyQP_Constants.h.
    // Fixed planner constants
    const double ThrustMargin = double(0.10);
    const double GimbalMargin = double(0.10);
    const double MaxGimbalAngle = deg2rad(10);
    const double MaxGimbalRate = deg2rad(30);
    const double MaxRollRate = 5;
    const double MaxTilt = deg2rad(30);
    const double MinSegTime = double(0.2);
    const double DtSample = double(0.1);

    // Vehicle selection
    // ASTRA-ONLY: the vehicle 0 (ASTRA) values; the vehicle 1 (TOAD) branch was removed.
    const double MT = ASTRA::MaxThrust;
    const double MaxThrustRate = 30;
    const double minThrottle = double(0.25);
    const double m_lox0 = 0;
    const double m_ipa0 = 0;
    const Vector3d defAccel = Vector3d(double(2.5), double(2.5), double(3.5));

    // Initial mass, inertia, and center of percussion offset
    const double m_init = ASTRA::m_dry + m_lox0 + m_ipa0;
    Matrix3_3d J_init;
    double CGz_init;
    ComputeJtot(m_lox0, m_ipa0, J_init, CGz_init);
    const double rho = M1(J_init, 1, 1) / (m_init * CGz_init);

    // Effective actuator bands with controller margins
    const double T_lo = (minThrottle + ThrustMargin) * MT;
    const double T_hi = (double(1.00) - ThrustMargin) * MT;
    const double gimbalLim = (double(1.00) - GimbalMargin) * MaxGimbalAngle;
    const double rollLim = (double(1.00) - ThrustMargin) * MaxRollRate;

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
