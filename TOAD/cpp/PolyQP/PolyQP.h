// PolyQP transpile support header (public declarations) for TOAD/Guidance/PolynomialQP @ 5ac6bd23
// ASTRA only: the vehicle constants are the constexpr values of ASTRA_PolyQP_Constants.h (the MATLAB Constants /
// constants6DoF argument is not carried) and the polynomial order is the constexpr POLYQP_R = 5 (the MATLAB r / Opts.R
// argument is not carried). STATIC: every array is fixed-size; grid functions take the number of points explicitly.
// STREAM: the output nodes go to a NodeSink (PolyQP_types.h) instead of InitialGuess arrays.
#pragma once

#include "PolyQP_types.h"

namespace PolyQP {

PolyQPStatus PolyQP(const TOADWaypointArray &Waypoints, const NodeSink &Sink, InitialGuessT &InitialGuess, QPTrajT &QPTraj, ReportT &Report);

PolyQPStatus resolveWaypoints(const TOADWaypointArray &Waypoints, const LimitsT &Limits, PlanT &Plan);

PolyQPStatus assembleQP(Matrix10_10d &A0, Matrix10_10d &Q0, Matrix10_10d &H0);

PolyQPStatus solveLeg(const Array3_5_Kd &FixLeg, const Array3_5_Kb &FixLegFixed, const Matrix1_Sd &SegTimes, const Matrix10_10d &A0, const Matrix10_10d &H0,
                      Array10_S_3d &CoefPhys, Array3_5_Kd &Deriv, double &Cost);

PolyQPStatus allocateTimes(const Array3_5_Kd &FixLeg, const Array3_5_Kb &FixLegFixed, const Matrix1_Sd &SegTimesInit, const Matrix10_10d &A0, const Matrix10_10d &H0,
                           double MinSegTime, Matrix1_Sd &SegTimes, double &Cost);

PolyQPStatus evaluateConstraints(const QPTrajT &QPTraj, const PlanT &Plan, const LimitsT &Limits, EvalOutT &EvalOut, bool &Feasible);

// STREAM: one time sample of flatnessInvert. qPrev carries the previous sample's quaternion for the sign continuity
// (hasPrevQ = false for the first sample); it is updated to this sample's quaternion.
PolyQPStatus flatnessInvertPoint(const QPTrajT &QPTraj, double t, const Vector3d &HDGRef, const LimitsT &Limits, bool hasPrevQ, Vector4d &qPrev, FlatPointT &Out);

PolyQPStatus flatnessInvert(const QPTrajT &QPTraj, const Matrix1_Nd &t_vec, int nPts, const Vector3d &HDGRef, const LimitsT &Limits, Matrix15_Nd &X, WrenchT &Wrench,
                            Matrix3_Nd &AccelCg);

// STREAM: allocateControls for one sample.
Vector4d allocateControlsPoint(const Vector3d &Fk, const Vector3d &tauk, double CGz, const LimitsT &Limits);

PolyQPStatus allocateControls(const WrenchT &Wrench, int nPts, const LimitsT &Limits, Matrix4_Nd &U);

PolyQPStatus sampleOutputs(const QPTrajT &QPTraj, const PlanT &Plan, const LimitsT &Limits, const NodeSink &Sink, InitialGuessT &InitialGuess);

// STREAM: the k-th derivative of the piecewise polynomial at one time.
Vector3d PolyEvalPoint(const QPTrajT &QPTraj, double t, int k = 0);

PolyQPStatus PolyEval(const QPTrajT &QPTraj, const Matrix1_Nd &t, int nPts, Matrix3_Nd &Out, int k = 0);

PolyQPStatus buildVehicleLimits(LimitsT &Limits);

// NO-NAN: posFixed marks the fixed position components (MATLAB: a NaN position component is free).
TOADWaypointT TOADWaypoint(const Vector3d &pos, const TOADWaypointOptions &options = TOADWaypointOptions(), const Vector3b &posFixed = Vector3b::Constant(true));

}
