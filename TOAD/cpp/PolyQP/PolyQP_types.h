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

using Real = double;

template <int MR, int MC> using MatX = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, MR, MC>;
template <int MR> using VecX = Eigen::Matrix<Real, Eigen::Dynamic, 1, Eigen::ColMajor, MR, 1>;
template <int MC> using RowX = Eigen::Matrix<Real, 1, Eigen::Dynamic, Eigen::RowMajor, 1, MC>;
template <int MC> using BoolRowX = Eigen::Array<bool, 1, Eigen::Dynamic, Eigen::RowMajor, 1, MC>;
template <int MC> using IntRowX = Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor, 1, MC>;

using Vec3 = Vector3d;
using Vec4 = Vector4d;
using Mat3 = Matrix3_3d;

using MatR = MatX<POLYQP_MAX_R, POLYQP_MAX_R>;
using VecR = VecX<POLYQP_MAX_R>;
using MatNC = MatX<POLYQP_MAX_NC, POLYQP_MAX_NC>;
using VecNC = VecX<POLYQP_MAX_NC>;
using MatRNC = MatX<POLYQP_MAX_R, POLYQP_MAX_NC>;
using MatNV = MatX<POLYQP_MAX_NV, POLYQP_MAX_NV>;
using VecNV = VecX<POLYQP_MAX_NV>;

using RowSeg = RowX<POLYQP_MAX_SEG>;
using RowKnot = RowX<POLYQP_MAX_KNOTS>;
using BoolSeg = BoolRowX<POLYQP_MAX_SEG>;
using BoolKnot = BoolRowX<POLYQP_MAX_KNOTS>;
using IntSeg = IntRowX<POLYQP_MAX_SEG>;
using IntGrid = IntRowX<POLYQP_MAX_NODES>;
using Mat3Seg = MatX<3, POLYQP_MAX_SEG>;
using Mat3Knot = MatX<3, POLYQP_MAX_KNOTS>;
using Mat9Seg = MatX<9, POLYQP_MAX_SEG>;
using Vec9 = VecX<9>;
using VecSeg = VecX<POLYQP_MAX_SEG>;

using RowGrid = RowX<POLYQP_MAX_NODES>;
using Mat3Grid = MatX<3, POLYQP_MAX_NODES>;
using Mat4Grid = MatX<4, POLYQP_MAX_NODES>;
using Mat15Grid = MatX<15, POLYQP_MAX_NODES>;
using BoolGrid = BoolRowX<POLYQP_MAX_NODES>;

using Fix3 = Array3<Real, 3, POLYQP_MAX_R, POLYQP_MAX_KNOTS>;
using FixOut3 = Array3<Real, 3, 6, POLYQP_MAX_KNOTS>;
using Coef3 = Array3<Real, POLYQP_MAX_NC, POLYQP_MAX_SEG, 3>;
using Bool3 = Array3<bool, 3, POLYQP_MAX_R, POLYQP_MAX_KNOTS>;

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

struct TOADWaypointT {
    Vec3 PosBC;
    Vec3 VelBC;
    Vec3 AccelBC;
    Vec3 JerkBC;
    Real HoldTime;
    Vec3 MaxVel;
    Vec3 MaxAccel;
    bool StraightLine;
    Vec3 HDGRef;
    char Label[10];
    bool IsPassAndGo;
    Vec3 B3BC;    // FIX-F03
    Vec3 B3DotBC; // FIX-F03
};

struct TOADWaypointOptions {
    Vec3 VelBC = Vec3::Constant(std::numeric_limits<Real>::quiet_NaN());
    Vec3 AccelBC = Vec3::Constant(std::numeric_limits<Real>::quiet_NaN());
    Vec3 JerkBC = Vec3::Constant(std::numeric_limits<Real>::quiet_NaN());
    Real HoldTime = std::numeric_limits<Real>::quiet_NaN();
    Vec3 MaxVel = Vec3(Real(4), Real(4), Real(5));
    Vec3 MaxAccel = Vec3::Constant(std::numeric_limits<Real>::infinity());
    bool StraightLine = false;
    Vec3 HDGRef = Vec3(Real(1), Real(0), Real(0));
    const char *Label = "";
    Vec3 B3BC = Vec3::Constant(std::numeric_limits<Real>::quiet_NaN());    // FIX-F03
    Vec3 B3DotBC = Vec3::Constant(std::numeric_limits<Real>::quiet_NaN()); // FIX-F03
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
    Real g;
    Real m_init;
    Real m_dry;
    Real MaxThrust;
    Real rTB;
    Real CGz_init;
    Mat3 J_init;
    Real Rho;
    Real T_lo;
    Real T_hi;
    Real GimbalLimit;
    Real GimbalRateLim;
    Real ThrustRateLim;
    Real RollLimit;
    Real RollRateLim;
    Real MaxTilt;
    Vec3 DefaultMaxAccel;
    Real MinSegTime;
    Real DtSample;
};

struct PlanT {
    int R;
    int KnotCount;
    int SegCount;
    Fix3 FixCg;
    BoolKnot KnotIsRest;
    RowSeg SegT;
    BoolSeg IsHold;
    BoolSeg IsCruise;
    Mat3Seg SegMaxVel;
    Mat3Seg SegMaxAccel;
    int LegCount;
    IntSeg LegStartSeg;
    IntSeg LegEndSeg;
    Vec3 HDGRef;
    Mat3Knot B3BC;    // FIX-F03
    Mat3Knot B3DotBC; // FIX-F03
};

struct QPTrajT {
    int R;
    int Count;
    RowSeg T;
    RowKnot KeyTimes;
    Coef3 Coef;
    BoolSeg IsHold;
    BoolSeg IsCruise;
    FixOut3 Fix;
    Real Cost;
    Real T_total;
    Real FlatOffset;
};

struct EvalOutT {
    Mat9Seg SegUtil;
    Vec9 MaxUtilPerType;
    Real MaxOverallUtil;
    RowSeg SegDilate;
};

struct ReportT {
    bool Feasible;
    Real MaxOverallUtil;
    Real UtilVel;
    Real UtilAccel;
    Real UtilThrustHi;
    Real UtilThrustLo;
    Real UtilGimbal;
    Real UtilGimbalRate;
    Real UtilThrustRate;
    Real UtilRoll;
    Real UtilTilt;
    int OuterIterations;
    Real Cost;
    RowSeg SegTimes;
    Real T_total;
};

struct InitialGuessT {
    RowGrid Time;
    Mat15Grid X;
    Mat4Grid U;
    Real T_total;
    RowGrid dt_row;
};

struct WrenchT {
    Mat3Grid F_body;
    Mat3Grid Tau_body;
    RowGrid CGz;
};

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
