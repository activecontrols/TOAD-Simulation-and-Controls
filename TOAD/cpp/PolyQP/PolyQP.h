// PolyQP transpile support header (public declarations) for TOAD/Guidance/PolynomialQP @ 5ac6bd23
// ASTRA only: the vehicle constants are the constexpr values of ASTRA_PolyQP_Constants.h (the MATLAB Constants /
// constants6DoF argument is not carried) and the polynomial order is the constexpr POLYQP_R = 5 (the MATLAB r / Opts.R
// argument is not carried).
#pragma once

#include "PolyQP_types.h"

namespace PolyQP {

PolyQPStatus PolyQP(const TOADWaypointArray &Waypoints, InitialGuessT &InitialGuess, QPTrajT &QPTraj, ReportT &Report);

PolyQPStatus resolveWaypoints(const TOADWaypointArray &Waypoints, const LimitsT &Limits, PlanT &Plan);

PolyQPStatus assembleQP(MatNC &A0, MatNC &Q0, MatNC &H0);

PolyQPStatus solveLeg(const Fix3 &FixLeg, const Bool3 &FixLegFixed, const RowSeg &SegTimes, const MatNC &A0, const MatNC &H0, Coef3 &CoefPhys, Fix3 &Deriv, Real &Cost);

PolyQPStatus allocateTimes(const Fix3 &FixLeg, const Bool3 &FixLegFixed, const RowSeg &SegTimesInit, const MatNC &A0, const MatNC &H0, Real MinSegTime, RowSeg &SegTimes, Real &Cost);

PolyQPStatus evaluateConstraints(const QPTrajT &QPTraj, const PlanT &Plan, const LimitsT &Limits, EvalOutT &EvalOut, bool &Feasible);

PolyQPStatus flatnessInvert(const QPTrajT &QPTraj, const RowGrid &t_vec, const Vec3 &HDGRef, const LimitsT &Limits, Mat15Grid &X, WrenchT &Wrench,
                            Mat3Grid &AccelCg);

PolyQPStatus allocateControls(const WrenchT &Wrench, const LimitsT &Limits, Mat4Grid &U);

PolyQPStatus sampleOutputs(const QPTrajT &QPTraj, const PlanT &Plan, const LimitsT &Limits, InitialGuessT &InitialGuess);

PolyQPStatus PolyEval(const QPTrajT &QPTraj, const RowGrid &t, Mat3Grid &Out, int k = 0);

PolyQPStatus buildVehicleLimits(LimitsT &Limits);

// NO-NAN: posFixed marks the fixed position components (MATLAB: a NaN position component is free).
TOADWaypointT TOADWaypoint(const Vec3 &pos, const TOADWaypointOptions &options = TOADWaypointOptions(), const Vector3b &posFixed = Vector3b::Constant(true));

}
