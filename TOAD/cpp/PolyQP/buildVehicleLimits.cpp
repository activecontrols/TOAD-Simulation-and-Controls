// Transpiled from TOAD/Guidance/PolynomialQP/buildVehicleLimits.m @ 5ac6bd23
#include "PolyQP.h"
#include "matlab_compat.h"

namespace PolyQP {

PolyQPStatus buildVehicleLimits(const PolyQPConstants &Constants, LimitsT &Limits, const PolyQPOpts &Opts) {
    // Build Vehicle Limits and Flatness Parameters
    // TRANSPILE-FLAG F38: the optional MATLAB Opts (nargin >= 2, isfield, isempty) is the default argument Opts = PolyQPOpts() placed last; a missing, empty or non-finite InitOxMass/InitFuMass is the NaN default (full tank).
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
    const Real veh = Constants.Vehicle;
    const Real MT = Constants.MaxThrust;
    Real MaxThrustRate;
    Real minThrottle;
    Real m_lox0;
    Real m_ipa0;
    Vec3 defAccel;
    if (veh == 1) {
        MaxThrustRate = 1500;
        minThrottle = Real(0.50);
        m_lox0 = Constants.OxMass;
        m_ipa0 = Constants.FuMass;
        // FIX-F29: start from the estimated propellant when given (Constants.OxMass/FuMass stay the full-tank reference for ComputeJtot)
        if (std::isfinite(Opts.InitOxMass)) {
            m_lox0 = Opts.InitOxMass;
        }
        if (std::isfinite(Opts.InitFuMass)) {
            m_ipa0 = Opts.InitFuMass;
        }
        defAccel = Vec3(Real(2.0), Real(2.0), Real(3.0));
    } else {
        MaxThrustRate = 30;
        minThrottle = Real(0.25);
        m_lox0 = 0;
        m_ipa0 = 0;
        defAccel = Vec3(Real(2.5), Real(2.5), Real(3.5));
    }

    // Initial mass, inertia, and center of percussion offset
    const Real m_init = Constants.m_dry + m_lox0 + m_ipa0;
    Mat3 J_init;
    Real CGz_init;
    POLYQP_TRY(ComputeJtot(m_lox0, m_ipa0, Constants, J_init, CGz_init));
    const Real rho = M1(J_init, 1, 1) / (m_init * CGz_init);

    // Effective actuator bands with controller margins
    const Real T_lo = (minThrottle + ThrustMargin) * MT;
    const Real T_hi = (Real(1.00) - ThrustMargin) * MT;
    const Real gimbalLim = (Real(1.00) - GimbalMargin) * MaxGimbalAngle;
    const Real rollLim = (Real(1.00) - ThrustMargin) * MaxRollRate;

    Limits.HasVehicle = true;
    Limits.Vehicle = veh;
    Limits.g = Constants.g;
    Limits.m_init = m_init;
    Limits.m_dry = Constants.m_dry;
    Limits.OxMass = m_lox0;
    Limits.FuMass = m_ipa0;
    Limits.MaxThrust = MT;
    Limits.MaxMdot = Constants.MaxMdot;
    Limits.OF = Constants.OF;
    Limits.rTB = Constants.rTB;
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
