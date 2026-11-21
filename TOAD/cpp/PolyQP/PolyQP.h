// PolyQP transpile support header (public declarations) for TOAD/Guidance/PolynomialQP @ 5ac6bd23
#pragma once

#include "PolyQP_types.h"

namespace PolyQP {

PolyQPStatus PolyQP(const TOADWaypointArray &Waypoints, const PolyQPConstants &Constants, InitialGuessT &InitialGuess, QPTrajT &QPTraj, ReportT &Report,
                    const PolyQPOpts &Opts = PolyQPOpts());

PolyQPStatus resolveWaypoints(const TOADWaypointArray &Waypoints, int r, const LimitsT &Limits, PlanT &Plan);

PolyQPStatus assembleQP(int r, MatNC &A0, MatNC &Q0, MatNC &H0);

PolyQPStatus solveLeg(const Fix3 &FixLeg, const RowSeg &SegTimes, const MatNC &A0, const MatNC &H0, Coef3 &CoefPhys, Fix3 &Deriv, Real &Cost);

PolyQPStatus allocateTimes(const Fix3 &FixLeg, const RowSeg &SegTimesInit, const MatNC &A0, const MatNC &H0, Real MinSegTime, RowSeg &SegTimes, Real &Cost);

PolyQPStatus evaluateConstraints(const QPTrajT &QPTraj, const PlanT &Plan, const PolyQPConstants &Constants, const LimitsT &Limits, EvalOutT &EvalOut,
                                 bool &Feasible);

PolyQPStatus flatnessInvert(const QPTrajT &QPTraj, const RowGrid &t_vec, const Vec3 &HDGRef, const PolyQPConstants &Constants, const LimitsT &Limits,
                            Mat15Grid &X, WrenchT &Wrench, Mat3Grid &AccelCg);

PolyQPStatus allocateControls(const WrenchT &Wrench, const LimitsT &Limits, Mat4Grid &U);

PolyQPStatus sampleOutputs(const QPTrajT &QPTraj, const PlanT &Plan, const PolyQPConstants &Constants, const LimitsT &Limits, InitialGuessT &InitialGuess);

PolyQPStatus PolyEval(const QPTrajT &QPTraj, const RowGrid &t, Mat3Grid &Out, int k = 0);

PolyQPStatus buildVehicleLimits(const PolyQPConstants &Constants, LimitsT &Limits, const PolyQPOpts &Opts = PolyQPOpts());

TOADWaypointT TOADWaypoint(const Vec3 &pos, const TOADWaypointOptions &options = TOADWaypointOptions());




}
