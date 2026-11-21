// PolyQP transpile support header (types and structs) for TOAD/Guidance/PolynomialQP @ 5ac6bd23
#pragma once

#include "ASTRA_PolyQP_Constants.h"
#include "PolyQP_config.h"

#include "matlab_funcs.h"

#include <cstdint>
#include <limits>

#if defined(POLYQP_ERROR_MESSAGES)
#include <cstdio>
#endif

namespace PolyQP {

// STATIC: no type aliases here; every array type is a fully fixed-size typedef in cpp/matlab_funcs.h.

enum class PolyQPStatus : int {
    Ok = 0,
    InvalidWaypointCount,
    InvalidOrder,
    InvalidHDGRef,
    InconsistentHDGRef,
    InvalidLimits,
    IllPosedQP,
    SingularRPP,
    CapacityExceeded
};

// NO-NAN: MATLAB marks a free boundary condition with NaN. Here every value-or-NaN field has a parallel "Fixed" mask
// (true = the MATLAB value is not NaN); a free entry holds 0 in the value field. A NaN HoldTime (pass-and-go) is
// IsPassAndGo = true with HoldTime = 0.
struct TOADWaypointT {
    Vector3d PosBC;
    Vector3b PosFixed;
    Vector3d VelBC;
    Vector3b VelFixed;
    Vector3d AccelBC;
    Vector3b AccelFixed;
    Vector3d JerkBC;
    Vector3b JerkFixed;
    double HoldTime;
    Vector3d MaxVel;
    Vector3d MaxAccel;
    bool StraightLine;
    Vector3d HDGRef;
    char Label[10];
    bool IsPassAndGo;
    Vector3d B3BC;    // FIX-F03
    Vector3b B3Fixed;
    Vector3d B3DotBC; // FIX-F03
    Vector3b B3DotFixed;
};

// The MATLAB TOADWaypoint name-value defaults; NO-NAN: the NaN defaults are "not fixed" (value 0).
struct TOADWaypointOptions {
    Vector3d VelBC = Vector3d::Zero();
    Vector3b VelFixed = Vector3b::Constant(false);
    Vector3d AccelBC = Vector3d::Zero();
    Vector3b AccelFixed = Vector3b::Constant(false);
    Vector3d JerkBC = Vector3d::Zero();
    Vector3b JerkFixed = Vector3b::Constant(false);
    double HoldTime = 0;
    bool HoldTimeGiven = false; // NO-NAN: MATLAB default HoldTime = NaN (pass-and-go)
    Vector3d MaxVel = Vector3d(double(4), double(4), double(5));
    Vector3d MaxAccel = Vector3d::Constant(std::numeric_limits<double>::infinity());
    bool StraightLine = false;
    Vector3d HDGRef = Vector3d(double(1), double(0), double(0));
    const char *Label = "";
    Vector3d B3BC = Vector3d::Zero();    // FIX-F03
    Vector3b B3Fixed = Vector3b::Constant(false);
    Vector3d B3DotBC = Vector3d::Zero(); // FIX-F03
    Vector3b B3DotFixed = Vector3b::Constant(false);
};

struct TOADWaypointArray {
    TOADWaypointT data[POLYQP_MAX_WP];
    int count = 0;

    int numel() const { return count; }
    TOADWaypointT &operator()(int i) {
        POLYQP_ASSERT(i >= 1 && i <= count);
        return data[i - 1];
    }
    const TOADWaypointT &operator()(int i) const {
        POLYQP_ASSERT(i >= 1 && i <= count);
        return data[i - 1];
    }
};

struct LimitsT {
    bool HasVehicle;
    double g;
    double m_init;
    double m_dry;
    double MaxThrust;
    double rTB;
    double CGz_init;
    Matrix3_3d J_init;
    double Rho;
    double T_lo;
    double T_hi;
    double GimbalLimit;
    double GimbalRateLim;
    double ThrustRateLim;
    double RollLimit;
    double RollRateLim;
    double MaxTilt;
    Vector3d DefaultMaxAccel;
    double MinSegTime;
    double DtSample;
};

struct PlanT {
    int R;
    int KnotCount;
    int SegCount;
    Array3_5_Kd FixCg;
    Array3_5_Kb FixCgFixed; // NO-NAN: ~isnan(FixCg)
    Array1_Kb KnotIsRest;
    Matrix1_Sd SegT;
    Array1_Sb IsHold;
    Array1_Sb IsCruise;
    Matrix3_Sd SegMaxVel;
    Matrix3_Sd SegMaxAccel;
    int LegCount;
    Array1_Si LegStartSeg;
    Array1_Si LegEndSeg;
    Vector3d HDGRef;
    Matrix3_Kd B3BC;    // FIX-F03
    Array3_Kb B3Fixed;    // NO-NAN: ~isnan(B3BC)
    Matrix3_Kd B3DotBC; // FIX-F03
    Array3_Kb B3DotFixed; // NO-NAN: ~isnan(B3DotBC)
};

struct QPTrajT {
    int R;
    int Count;
    Matrix1_Sd T;
    Matrix1_Kd KeyTimes;
    Array10_S_3d Coef;
    Array1_Sb IsHold;
    Array1_Sb IsCruise;
    Array3_6_Kd Fix;
    Array3_6_Kb FixFixed; // NO-NAN: ~isnan(Fix)
    double Cost;
    double T_total;
    double FlatOffset;
};

struct EvalOutT {
    Matrix9_Sd SegUtil;
    Vector9d MaxUtilPerType;
    double MaxOverallUtil;
    Matrix1_Sd SegDilate;
};

struct ReportT {
    bool Feasible;
    double MaxOverallUtil;
    double UtilVel;
    double UtilAccel;
    double UtilThrustHi;
    double UtilThrustLo;
    double UtilGimbal;
    double UtilGimbalRate;
    double UtilThrustRate;
    double UtilRoll;
    double UtilTilt;
    int OuterIterations;
    double Cost;
    Matrix1_Sd SegTimes;
    double T_total;
};

// STREAM: MATLAB InitialGuess.Time / X / U are not stored. sampleOutputs (and the abort fallback) hand each output node to a
// NodeSink as it is evaluated, so the samples go straight into the caller's arrays (the flight trajectory arrays, or a MATLAB
// array in the verification MEX). InitialGuessT keeps the scalars.
struct InitialGuessT {
    int Nn;         // number of output nodes, numel(InitialGuess.Time)
    double T_total; // InitialGuess.T_total
    double dt;      // every entry of InitialGuess.dt_row
};

// STREAM: receives output node k (1-based) of an abort trajectory: time t, the 15-state X, and the 4 inputs U. For k <= Nn - 1,
// U = InitialGuess.U(:, k); for the last node, U = InitialGuess.U(:, end), the duplicated last row of the abort block's
// Inputs = [U'; U(:, end)']. capacity is the number of nodes the sink can hold; a longer trajectory is CapacityExceeded.
struct NodeSink {
    void *ctx;
    int capacity;
    void (*put)(void *ctx, int k, double t, const Vector15d &X, const Vector4d &U);
};

// Grid outputs of flatnessInvert / allocateControls (kept for verification; the flight path streams per point).
struct WrenchT {
    Matrix3_Nd F_body;
    Matrix3_Nd Tau_body;
    Matrix1_Nd CGz;
};

// STREAM: one sample of flatnessInvert (the per-point function the grid version loops over).
struct FlatPointT {
    Vector15d X;      // X(:, k) = [q; r_cg; v_cg; w_B; m_lox; m_ipa]
    Vector3d F_body;  // Wrench.F_body(:, k)
    Vector3d Tau_body; // Wrench.Tau_body(:, k)
    double CGz;       // Wrench.CGz(k)
    Vector3d AccelCg; // AccelCg(:, k)
};

// STATIC: every solver buffer that is not a small local lives in this one workspace (one object, PolyQP.cpp), instead of
// function-local statics spread over the files. The stages run sequentially and never re-enter one another's members:
// PolyQP -> resolveWaypoints, solveAllLegs (-> allocateTimes -> solveLeg, solveLeg), evaluateConstraints, sampleOutputs.
// With the grids streamed (evaluateConstraints, sampleOutputs), no grid-sized buffer remains, so nothing large is left to
// overlap. The solver is therefore not reentrant: call it from one task only, never from an ISR (TRANSPILE_FLAGS.md).
struct PolyQPWorkspace {
    struct SolveLegT {
        MatrixNV_NVd Rmat;
        Array10_10_Sd Hseg;
        MatrixNV_NVd Rns;
        VectorNVd fx;
        VectorNVb fxFixed;
        VectorNVd dv;
        VectorNVd dF;
        MatrixF_Fd R_PP;
        MatrixF_NVd R_PF;
        MatrixF_Fd Tz;
        MatrixF_Fd Rns_PP;
        MatrixF_Fd TzT_Rns_PP;
        MatrixF_Fd R_PP_T;
        MatrixF_NVd Rns_PF;
        VectorFd d_diag;
        VectorFd s_jac;
        MatrixF_Fd s_jac_s_jacT;
        MatrixF_Fd R_PP_eq;
        MatrixF_Fd R_PP_eq_T;
        VectorFd R_PF_dF;
        VectorFd rhs_eq;
        Eigen::PartialPivLU<MatrixF_Fd> lu_rc;
        Eigen::LLT<MatrixF_Fd> llt;
        Eigen::PartialPivLU<MatrixF_Fd> lu;
        VectorFd R_PP_eq_rhs_eq;
        VectorFd z;
        VectorFd Tz_z;
    } solveLeg;
    struct SolveAllLegsT {
        Array3_5_Kd FixP;
        Array3_5_Kb FixPFixed;
        Array10_S_3d CoefAll;
        Array3_5_Kd DerivAll;
        Array3_5_Kb DerivAllFixed;
        Array3_5_Kd FixLegAT; // allocateTimes branch
        Array3_5_Kb FixLegFixedAT;
        Array3_5_Kd FixLeg;
        Array3_5_Kb FixLegFixed;
        Array10_S_3d CoefLeg;
        Array3_5_Kd DerivLeg;
        Array3_6_Kd FixOut;
        Array3_6_Kb FixOutFixed;
    } solveAllLegs;
    struct AllocateTimesT {
        Array10_S_3d CoefIgnored;
        Array3_5_Kd DerivIgnored;
        Array3_5_Kd Deriv;
    } allocateTimes;
    struct ResolveWaypointsT {
        Array3_5_Kd FixCg;
        Array3_5_Kb FixCgFixed;
    } resolveWaypoints;
    struct PolyQPT {
        PlanT Plan;
        EvalOutT EvalOut;
        Array3_5_Kd dLeg;
        Array3_5_Kb dLegFixed;
    } PolyQP;
    struct AbortGuidanceT {
        TOADWaypointArray Waypoints;
        QPTrajT QPTraj;
    } AbortGuidance;
};
PolyQPWorkspace &workspace();

#if defined(POLYQP_ERROR_MESSAGES)
// Verification only: the padded R_PP_eq of the last solveLeg call per axis, for the padded-vs-unpadded rcond check.
struct PolyQPRcondRecord {
    double rc[3];
    int nP[3];
    MatrixF_Fd R_PP_eq[3];
};
inline PolyQPRcondRecord &lastRcond() {
    static PolyQPRcondRecord rec{};
    return rec;
}
#endif

#if defined(POLYQP_ERROR_MESSAGES)
struct PolyQPErrorRecord {
    char identifier[48];
    char message[256];
};
inline PolyQPErrorRecord &lastError() {
    static PolyQPErrorRecord rec{};
    return rec;
}
#endif

inline bool repmat(const TOADWaypointT &A, int m, int n, TOADWaypointArray &out) {
    if (m * n > POLYQP_MAX_WP) return false;
    out.count = m * n;
    for (int q = 1; q <= m * n; ++q) out(q) = A;
    return true;
}

}

#if defined(POLYQP_ERROR_MESSAGES)
#define POLYQP_ERROR(ID, ...)                                                                                                                                  \
    do {                                                                                                                                                       \
        std::snprintf(::PolyQP::lastError().identifier, sizeof(::PolyQP::lastError().identifier), "%s", "PolyQP:" #ID);                                       \
        std::snprintf(::PolyQP::lastError().message, sizeof(::PolyQP::lastError().message), __VA_ARGS__);                                                      \
        return ::PolyQP::PolyQPStatus::ID;                                                                                                                     \
    } while (0)
#else
#define POLYQP_ERROR(ID, ...) return ::PolyQP::PolyQPStatus::ID
#endif

#define POLYQP_CAPACITY(cond)                                                                                                                                  \
    do {                                                                                                                                                       \
        if (!(cond)) {                                                                                                                                         \
            POLYQP_ERROR(CapacityExceeded, "PolyQP capacity exceeded: %s", #cond);                                                                             \
        }                                                                                                                                                      \
    } while (0)
