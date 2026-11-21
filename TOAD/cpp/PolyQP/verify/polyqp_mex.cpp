// Verification only (never part of the flight build): MATLAB MEX dispatcher exposing every public PolyQP function.
// Usage: [out..., status, msg] = polyqp_mex('<FunctionName>', args...)
// status is '' on success or the MATLAB error identifier (e.g. 'PolyQP:SingularRPP'); msg is the formatted message.
//
// The boundary converts between MATLAB and the ASTRA-only, NaN-free, fixed-size C++:
//  - NO-NAN: a MATLAB NaN boundary condition is (value 0, fixed = false) in C++, and back.
//  - STATIC: the C++ arrays have the abort capacities of PolyQP_config.h (polyqp_mex) or the larger verification capacities
//    (polyqp_mex_wide). MATLAB's own structs keep its capacities (MAX_SEG = 40, MAX_KNOTS = 41, MAX_R = 6), so inputs are read
//    up to the C++ capacity (entries beyond it must hold MATLAB's default, else a capacity error) and outputs are padded back.
//  - STREAM: InitialGuess and the abort trajectory arrive through a NodeSink and are assembled into MATLAB arrays here.
#include "mex.h"

#include "AbortGuidance.h"
#include "PolyQP.h"

#include <cmath>
#include <cstring>
#include <string>
#include <vector>

using namespace PolyQP;

namespace {

constexpr int MATLAB_MAX_SEG = 40;   // resolveWaypoints.m / PolyQP.m
constexpr int MATLAB_MAX_KNOTS = 41; // MAX_SEG + 1
constexpr int MATLAB_MAX_R = 6;      // resolveWaypoints.m FixCg = nan(3, MAX_R, MAX_KNOTS)
constexpr int MATLAB_MAX_NODES = 4000; // sampleOutputs.m: Nc <= 3999

void boundsHookImpl(const char *file, int line) { mexErrMsgIdAndTxt("polyqp_mex:bounds", "PolyQP bounds check failed at %s:%d", file, line); }

bool sameOrBothNaN(double a, double b) { return a == b || (std::isnan(a) && std::isnan(b)); }

// ------------------------------------------------------------------ readers
const mxArray *getField(const mxArray *s, const char *name, mwIndex i = 0) {
    if (!mxIsStruct(s)) mexErrMsgIdAndTxt("polyqp_mex:type", "Expected a struct while reading '%s'", name);
    const mxArray *f = mxGetField(s, i, name);
    if (f == nullptr) mexErrMsgIdAndTxt("polyqp_mex:field", "Missing struct field '%s'", name);
    return f;
}

double elem(const mxArray *a, size_t i) {
    if (mxIsLogical(a)) return mxGetLogicals(a)[i] ? 1.0 : 0.0;
    if (mxIsDouble(a) && !mxIsComplex(a)) return mxGetPr(a)[i];
    mexErrMsgIdAndTxt("polyqp_mex:type", "Expected a real double or logical array");
    return 0.0;
}

double num(const mxArray *a, const char *what) {
    if (a == nullptr || mxGetNumberOfElements(a) != 1) mexErrMsgIdAndTxt("polyqp_mex:type", "'%s' must be a scalar", what);
    return elem(a, 0);
}
double numField(const mxArray *s, const char *name, mwIndex i = 0) { return num(getField(s, name, i), name); }

// STATIC: a fixed-size matrix of R rows and capacity C columns from a MATLAB R x n array: columns beyond C must equal fill.
template <typename MatT> void readCols(const mxArray *a, MatT &M, double fill, const char *what) {
    const int rows = static_cast<int>(mxGetM(a));
    const int cols = static_cast<int>(mxGetN(a));
    if (rows != MatT::RowsAtCompileTime) mexErrMsgIdAndTxt("polyqp_mex:size", "'%s' has %d rows, expected %d", what, rows, static_cast<int>(MatT::RowsAtCompileTime));
    using S = typename MatT::Scalar;
    M.setConstant(static_cast<S>(fill));
    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i < rows; ++i) {
            const double v = elem(a, static_cast<size_t>(i + j * rows));
            if (j < MatT::ColsAtCompileTime) {
                M(i, j) = static_cast<S>(v);
            } else if (!sameOrBothNaN(v, fill)) {
                mexErrMsgIdAndTxt("polyqp_mex:capacity", "'%s' column %d is used, beyond the C++ capacity %d", what, j + 1, static_cast<int>(MatT::ColsAtCompileTime));
            }
        }
    }
}
template <typename MatT> void readColsField(const mxArray *s, const char *name, MatT &M, double fill) { readCols(getField(s, name), M, fill, name); }

// A fixed-size vector (3, 4, 15, ...) from a MATLAB vector of exactly that many elements
template <typename VecT> void readVecN(const mxArray *a, VecT &v, const char *what) {
    const int n = static_cast<int>(mxGetNumberOfElements(a));
    if (n != VecT::SizeAtCompileTime) mexErrMsgIdAndTxt("polyqp_mex:size", "'%s' must have %d elements", what, static_cast<int>(VecT::SizeAtCompileTime));
    for (int i = 0; i < n; ++i) v(i) = static_cast<typename VecT::Scalar>(elem(a, static_cast<size_t>(i)));
}
Vector3d readVec3(const mxArray *a, const char *what) {
    Vector3d v;
    readVecN(a, v, what);
    return v;
}

// A fixed-capacity row (1 x C) and its count from a MATLAB vector of n <= C elements
template <typename RowT> int readRowCount(const mxArray *a, RowT &v, const char *what) {
    const int n = static_cast<int>(mxGetNumberOfElements(a));
    if (n > RowT::SizeAtCompileTime) mexErrMsgIdAndTxt("polyqp_mex:capacity", "'%s' has %d elements, capacity %d", what, n, static_cast<int>(RowT::SizeAtCompileTime));
    v.setZero();
    for (int i = 0; i < n; ++i) v(i) = static_cast<typename RowT::Scalar>(elem(a, static_cast<size_t>(i)));
    return n;
}

// STATIC: an Array3 from a MATLAB array whose dimensions may exceed the C++ capacity; the excess must equal fill.
template <typename A3> void readArray3Cap(const mxArray *a, A3 &A, double fill, const char *what) {
    const mwSize nd = mxGetNumberOfDimensions(a);
    const mwSize *dims = mxGetDimensions(a);
    if (nd > 3) mexErrMsgIdAndTxt("polyqp_mex:size", "'%s' must have at most 3 dimensions", what);
    const int n1 = static_cast<int>(dims[0]);
    const int n2 = static_cast<int>(dims[1]);
    const int n3 = nd == 3 ? static_cast<int>(dims[2]) : 1;
    if (n1 > A3::Cap1) mexErrMsgIdAndTxt("polyqp_mex:capacity", "'%s' is %dx%dx%d, capacity %dx%dx%d", what, n1, n2, n3, A3::Cap1, A3::Cap2, A3::Cap3);
    using T = typename std::remove_reference<decltype(A.data[0])>::type;
    A.resize(n1, std::min(n2, int(A3::Cap2)), std::min(n3, int(A3::Cap3)));
    for (int k = 1; k <= n3; ++k) {
        for (int j = 1; j <= n2; ++j) {
            for (int i = 1; i <= n1; ++i) {
                const double v = elem(a, static_cast<size_t>((i - 1) + n1 * ((j - 1) + n2 * (k - 1))));
                if (j <= A3::Cap2 && k <= A3::Cap3) {
                    A(i, j, k) = static_cast<T>(v);
                } else if (!sameOrBothNaN(v, fill)) {
                    mexErrMsgIdAndTxt("polyqp_mex:capacity", "'%s'(%d, %d, %d) is used, beyond the C++ capacity %dx%dx%d", what, i, j, k, A3::Cap1, A3::Cap2, A3::Cap3);
                }
            }
        }
    }
}

// NO-NAN boundary: the C++ solver carries free boundary conditions as (value 0, fixed = false); MATLAB marks them NaN.
void readVec3Mask(const mxArray *a, const char *what, Vector3d &v, Vector3b &fixed) {
    v = readVec3(a, what);
    for (int i = 0; i < 3; ++i) {
        fixed(i) = !std::isnan(v(i));
        if (!fixed(i)) v(i) = 0;
    }
}

template <typename A3, typename B3> void readArray3Mask(const mxArray *a, A3 &A, B3 &F, const char *what) {
    readArray3Cap(a, A, mxGetNaN(), what);
    F.resize(A.n1, A.n2, A.n3);
    for (int q = 0; q < A.numel(); ++q) {
        F.data[q] = !std::isnan(static_cast<double>(A.data[q]));
        if (!F.data[q]) A.data[q] = 0;
    }
}

void readKnotMask(const mxArray *a, const char *what, Matrix3_Kd &M, Array3_Kb &F) {
    readCols(a, M, mxGetNaN(), what);
    for (int j = 0; j < M.cols(); ++j)
        for (int i = 0; i < 3; ++i) {
            F(i, j) = !std::isnan(M(i, j));
            if (!F(i, j)) M(i, j) = 0;
        }
}

// MATLAB resolveWaypoints.m sizes Plan.FixCg with its own capacity MAX_R = 6 (3 x 6 x MAX_KNOTS); the C++ stores the r = 5
// columns and its own knot capacity. Column 6 and the knots beyond the C++ capacity must be NaN (unused).
void readFixCg(const mxArray *a, Array3_5_Kd &F, Array3_5_Kb &Ffixed) {
    const mwSize nd = mxGetNumberOfDimensions(a);
    const mwSize *dims = mxGetDimensions(a);
    const int n1 = static_cast<int>(dims[0]);
    const int n2 = static_cast<int>(dims[1]);
    const int n3 = nd == 3 ? static_cast<int>(dims[2]) : 1;
    if (n1 != 3 || (n2 != POLYQP_R && n2 != MATLAB_MAX_R)) mexErrMsgIdAndTxt("polyqp_mex:capacity", "'FixCg' is %dx%dx%d", n1, n2, n3);
    const int K = std::min(n3, int(Array3_5_Kd::Cap3));
    F.resize(3, POLYQP_R, K);
    Ffixed.resize(3, POLYQP_R, K);
    for (int k = 1; k <= n3; ++k) {
        for (int j = 1; j <= n2; ++j) {
            for (int i = 1; i <= 3; ++i) {
                const double v = elem(a, static_cast<size_t>((i - 1) + 3 * ((j - 1) + n2 * (k - 1))));
                if (j <= POLYQP_R && k <= K) {
                    Ffixed(i, j, k) = !std::isnan(v);
                    F(i, j, k) = Ffixed(i, j, k) ? static_cast<double>(v) : double(0);
                } else if (!std::isnan(v)) {
                    mexErrMsgIdAndTxt("polyqp_mex:capacity", "'FixCg'(%d, %d, %d) is used, beyond the C++ storage 3x%dx%d", i, j, k, POLYQP_R, K);
                }
            }
        }
    }
}

// ASTRA only: the solver's constants are compiled in (ASTRA_PolyQP_Constants.h). The MATLAB constants struct is still
// passed by the verification scripts; it must be LoadTOADParams(0), so check every value the solver uses.
void checkConstants(const mxArray *s) {
    Matrix3_3d J;
    readCols(getField(s, "J"), J, 0.0, "J");
    const Matrix3_3d Jc = Vector3d(ASTRA::J_xx, ASTRA::J_yy, ASTRA::J_zz).asDiagonal();
    if (numField(s, "Vehicle") != 0.0 || numField(s, "m_dry") != ASTRA::m_dry || numField(s, "g") != ASTRA::g || numField(s, "rTB") != ASTRA::rTB ||
        numField(s, "MaxThrust") != ASTRA::MaxThrust || J != Jc)
        mexErrMsgIdAndTxt("polyqp_mex:constants", "The constants struct is not LoadTOADParams(0): the C++ solver is ASTRA-only with compiled-in constants");
}

void checkOrder(int r) {
    if (r != POLYQP_R) mexErrMsgIdAndTxt("polyqp_mex:order", "The C++ solver hardcodes r = %d (got r = %d)", POLYQP_R, r);
}

LimitsT readLimits(const mxArray *s) {
    LimitsT L;
    L.HasVehicle = numField(s, "HasVehicle") != 0.0;
    L.g = numField(s, "g");
    L.m_init = numField(s, "m_init");
    L.m_dry = numField(s, "m_dry");
    L.MaxThrust = numField(s, "MaxThrust");
    L.rTB = numField(s, "rTB");
    L.CGz_init = numField(s, "CGz_init");
    readCols(getField(s, "J_init"), L.J_init, 0.0, "J_init");
    L.Rho = numField(s, "Rho");
    L.T_lo = numField(s, "T_lo");
    L.T_hi = numField(s, "T_hi");
    L.GimbalLimit = numField(s, "GimbalLimit");
    L.GimbalRateLim = numField(s, "GimbalRateLim");
    L.ThrustRateLim = numField(s, "ThrustRateLim");
    L.RollLimit = numField(s, "RollLimit");
    L.RollRateLim = numField(s, "RollRateLim");
    L.MaxTilt = numField(s, "MaxTilt");
    L.DefaultMaxAccel = readVec3(getField(s, "DefaultMaxAccel"), "DefaultMaxAccel");
    L.MinSegTime = numField(s, "MinSegTime");
    L.DtSample = numField(s, "DtSample");
    return L;
}

void readPlan(const mxArray *s, PlanT &P) {
    P.R = static_cast<int>(numField(s, "R"));
    checkOrder(P.R);
    P.KnotCount = static_cast<int>(numField(s, "KnotCount"));
    P.SegCount = static_cast<int>(numField(s, "SegCount"));
    if (P.KnotCount > POLYQP_MAX_KNOTS || P.SegCount > POLYQP_MAX_SEG) mexErrMsgIdAndTxt("polyqp_mex:capacity", "Plan has %d knots / %d segments", P.KnotCount, P.SegCount);
    readFixCg(getField(s, "FixCg"), P.FixCg, P.FixCgFixed);
    readColsField(s, "KnotIsRest", P.KnotIsRest, 0.0);
    readColsField(s, "SegT", P.SegT, 0.0);
    readColsField(s, "IsHold", P.IsHold, 0.0);
    readColsField(s, "IsCruise", P.IsCruise, 0.0);
    readColsField(s, "SegMaxVel", P.SegMaxVel, mxGetInf());
    readColsField(s, "SegMaxAccel", P.SegMaxAccel, mxGetInf());
    P.LegCount = static_cast<int>(numField(s, "LegCount"));
    readColsField(s, "LegStartSeg", P.LegStartSeg, 0.0);
    readColsField(s, "LegEndSeg", P.LegEndSeg, 0.0);
    P.HDGRef = readVec3(getField(s, "HDGRef"), "HDGRef");
    readKnotMask(getField(s, "B3BC"), "B3BC", P.B3BC, P.B3Fixed);
    readKnotMask(getField(s, "B3DotBC"), "B3DotBC", P.B3DotBC, P.B3DotFixed);
}

void readQPTraj(const mxArray *s, QPTrajT &Q) {
    Q.R = static_cast<int>(numField(s, "R"));
    checkOrder(Q.R);
    Q.Count = static_cast<int>(numField(s, "Count"));
    if (Q.Count > POLYQP_MAX_SEG) mexErrMsgIdAndTxt("polyqp_mex:capacity", "QPTraj has %d segments", Q.Count);
    readColsField(s, "T", Q.T, 0.0);
    readColsField(s, "KeyTimes", Q.KeyTimes, 0.0);
    readArray3Cap(getField(s, "Coef"), Q.Coef, 0.0, "Coef");
    readColsField(s, "IsHold", Q.IsHold, 0.0);
    readColsField(s, "IsCruise", Q.IsCruise, 0.0);
    readArray3Mask(getField(s, "Fix"), Q.Fix, Q.FixFixed, "Fix");
    Q.Cost = numField(s, "Cost");
    Q.T_total = numField(s, "T_total");
    Q.FlatOffset = numField(s, "FlatOffset");
}

int readGridRow(const mxArray *a, Matrix1_Nd &t, const char *what) { return readRowCount(a, t, what); }

int readWrench(const mxArray *s, WrenchT &W) {
    const int n = static_cast<int>(mxGetN(getField(s, "F_body")));
    if (n > POLYQP_MAX_NODES) mexErrMsgIdAndTxt("polyqp_mex:capacity", "Wrench has %d samples, capacity %d", n, POLYQP_MAX_NODES);
    readColsField(s, "F_body", W.F_body, 0.0);
    readColsField(s, "Tau_body", W.Tau_body, 0.0);
    readColsField(s, "CGz", W.CGz, 0.0);
    return n;
}

void readWaypoints(const mxArray *s, TOADWaypointArray &W) {
    const size_t n = mxGetNumberOfElements(s);
    if (static_cast<int>(n) > POLYQP_MAX_WP) mexErrMsgIdAndTxt("polyqp_mex:capacity", "Waypoint count %d exceeds POLYQP_MAX_WP", static_cast<int>(n));
    W.count = static_cast<int>(n);
    for (size_t i = 0; i < n; ++i) {
        TOADWaypointT &wp = W.data[i];
        const mwIndex k = static_cast<mwIndex>(i);
        readVec3Mask(getField(s, "PosBC", k), "PosBC", wp.PosBC, wp.PosFixed);
        readVec3Mask(getField(s, "VelBC", k), "VelBC", wp.VelBC, wp.VelFixed);
        readVec3Mask(getField(s, "AccelBC", k), "AccelBC", wp.AccelBC, wp.AccelFixed);
        readVec3Mask(getField(s, "JerkBC", k), "JerkBC", wp.JerkBC, wp.JerkFixed);
        wp.HoldTime = numField(s, "HoldTime", k);
        if (std::isnan(wp.HoldTime)) wp.HoldTime = 0; // NO-NAN: pass-and-go (IsPassAndGo below)
        wp.MaxVel = readVec3(getField(s, "MaxVel", k), "MaxVel");
        wp.MaxAccel = readVec3(getField(s, "MaxAccel", k), "MaxAccel");
        wp.StraightLine = numField(s, "StraightLine", k) != 0.0;
        wp.HDGRef = readVec3(getField(s, "HDGRef", k), "HDGRef");
        std::memset(wp.Label, 0, sizeof(wp.Label));
        mxGetString(getField(s, "Label", k), wp.Label, sizeof(wp.Label));
        wp.IsPassAndGo = numField(s, "IsPassAndGo", k) != 0.0;
        readVec3Mask(getField(s, "B3BC", k), "B3BC", wp.B3BC, wp.B3Fixed);
        readVec3Mask(getField(s, "B3DotBC", k), "B3DotBC", wp.B3DotBC, wp.B3DotFixed);
    }
}

// ------------------------------------------------------------------ writers
mxArray *scal(double v) { return mxCreateDoubleScalar(v); }
mxArray *lscal(bool v) { return mxCreateLogicalScalar(v); }

// rows x cols block of M (cols <= its capacity), padded with fill up to padCols columns
template <typename Derived> mxArray *toMxPad(const Eigen::DenseBase<Derived> &M, int rows, int cols, int padCols, double fill) {
    mxArray *a = mxCreateDoubleMatrix(static_cast<mwSize>(rows), static_cast<mwSize>(padCols), mxREAL);
    double *p = mxGetPr(a);
    for (int j = 0; j < padCols; ++j)
        for (int i = 0; i < rows; ++i) p[i + j * rows] = j < cols ? static_cast<double>(M.derived().coeff(i, j)) : fill;
    return a;
}
template <typename Derived> mxArray *toMx(const Eigen::DenseBase<Derived> &M) { return toMxPad(M, static_cast<int>(M.rows()), static_cast<int>(M.cols()), static_cast<int>(M.cols()), 0.0); }
template <typename Derived> mxArray *toMxCols(const Eigen::DenseBase<Derived> &M, int cols) { return toMxPad(M, static_cast<int>(M.rows()), cols, cols, 0.0); }
template <typename Derived> mxArray *toMxLogicalPad(const Eigen::DenseBase<Derived> &M, int cols, int padCols) {
    mxArray *a = mxCreateLogicalMatrix(static_cast<mwSize>(M.rows()), static_cast<mwSize>(padCols));
    mxLogical *p = mxGetLogicals(a);
    for (int j = 0; j < padCols; ++j)
        for (Eigen::Index i = 0; i < M.rows(); ++i) p[i + j * M.rows()] = (j < cols && M.derived().coeff(i, j)) ? 1 : 0;
    return a;
}
template <typename A3> mxArray *toMx3(const A3 &A) {
    const mwSize dims[3] = {static_cast<mwSize>(A.n1), static_cast<mwSize>(A.n2), static_cast<mwSize>(A.n3)};
    mxArray *a = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double *p = mxGetPr(a);
    for (int q = 0; q < A.numel(); ++q) p[q] = static_cast<double>(A.data[q]);
    return a;
}
// Array3 (with an optional mask: masked-out entries are NaN), padded to p2 x p3 with fill
template <typename A3> mxArray *toMx3Pad(const A3 &A, const bool *mask, int p2, int p3, double fill) {
    const mwSize dims[3] = {static_cast<mwSize>(A.n1), static_cast<mwSize>(p2), static_cast<mwSize>(p3)};
    mxArray *a = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double *p = mxGetPr(a);
    for (int k = 1; k <= p3; ++k)
        for (int j = 1; j <= p2; ++j)
            for (int i = 1; i <= A.n1; ++i) {
                double v = fill;
                if (j <= A.n2 && k <= A.n3) {
                    const int q = (i - 1) + A.n1 * ((j - 1) + A.n2 * (k - 1));
                    v = (mask == nullptr || mask[q]) ? static_cast<double>(A.data[q]) : mxGetNaN();
                }
                p[(i - 1) + A.n1 * ((j - 1) + p2 * (k - 1))] = v;
            }
    return a;
}

mxArray *writeVec3Mask(const Vector3d &v, const Vector3b &fixed) {
    mxArray *a = mxCreateDoubleMatrix(3, 1, mxREAL);
    double *p = mxGetPr(a);
    for (int i = 0; i < 3; ++i) p[i] = fixed(i) ? static_cast<double>(v(i)) : mxGetNaN();
    return a;
}

mxArray *writeKnotMask(const Matrix3_Kd &M, const Array3_Kb &F) {
    mxArray *a = mxCreateDoubleMatrix(3, MATLAB_MAX_KNOTS, mxREAL);
    double *p = mxGetPr(a);
    for (int j = 0; j < MATLAB_MAX_KNOTS; ++j)
        for (int i = 0; i < 3; ++i) p[i + 3 * j] = (j < M.cols() && F(i, j)) ? static_cast<double>(M(i, j)) : mxGetNaN();
    return a;
}

mxArray *writeFixCg(const Array3_5_Kd &F, const Array3_5_Kb &Ffixed) {
    const mwSize dims[3] = {3, static_cast<mwSize>(MATLAB_MAX_R), static_cast<mwSize>(MATLAB_MAX_KNOTS)};
    mxArray *a = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double *p = mxGetPr(a);
    for (int k = 1; k <= MATLAB_MAX_KNOTS; ++k)
        for (int j = 1; j <= MATLAB_MAX_R; ++j)
            for (int i = 1; i <= 3; ++i)
                p[(i - 1) + 3 * ((j - 1) + MATLAB_MAX_R * (k - 1))] = (j <= F.n2 && k <= F.n3 && Ffixed(i, j, k)) ? static_cast<double>(F(i, j, k)) : mxGetNaN();
    return a;
}

struct StructBuilder {
    std::vector<const char *> names;
    std::vector<mxArray *> vals;
    StructBuilder &add(const char *n, mxArray *v) {
        names.push_back(n);
        vals.push_back(v);
        return *this;
    }
    mxArray *build() {
        mxArray *s = mxCreateStructMatrix(1, 1, static_cast<int>(names.size()), names.data());
        for (size_t i = 0; i < names.size(); ++i) mxSetFieldByNumber(s, 0, static_cast<int>(i), vals[i]);
        return s;
    }
};

mxArray *writeLimits(const LimitsT &L) {
    return StructBuilder()
        .add("HasVehicle", lscal(L.HasVehicle))
        .add("g", scal(L.g))
        .add("m_init", scal(L.m_init))
        .add("m_dry", scal(L.m_dry))
        .add("MaxThrust", scal(L.MaxThrust))
        .add("rTB", scal(L.rTB))
        .add("CGz_init", scal(L.CGz_init))
        .add("J_init", toMx(L.J_init))
        .add("Rho", scal(L.Rho))
        .add("T_lo", scal(L.T_lo))
        .add("T_hi", scal(L.T_hi))
        .add("GimbalLimit", scal(L.GimbalLimit))
        .add("GimbalRateLim", scal(L.GimbalRateLim))
        .add("ThrustRateLim", scal(L.ThrustRateLim))
        .add("RollLimit", scal(L.RollLimit))
        .add("RollRateLim", scal(L.RollRateLim))
        .add("MaxTilt", scal(L.MaxTilt))
        .add("DefaultMaxAccel", toMx(L.DefaultMaxAccel))
        .add("MinSegTime", scal(L.MinSegTime))
        .add("DtSample", scal(L.DtSample))
        .build();
}

mxArray *writePlan(const PlanT &P) {
    const int S = POLYQP_MAX_SEG;
    const int K = POLYQP_MAX_KNOTS;
    return StructBuilder()
        .add("R", scal(P.R))
        .add("KnotCount", scal(P.KnotCount))
        .add("SegCount", scal(P.SegCount))
        .add("FixCg", writeFixCg(P.FixCg, P.FixCgFixed))
        .add("KnotIsRest", toMxLogicalPad(P.KnotIsRest, K, MATLAB_MAX_KNOTS))
        .add("SegT", toMxPad(P.SegT, 1, S, MATLAB_MAX_SEG, 0.0))
        .add("IsHold", toMxLogicalPad(P.IsHold, S, MATLAB_MAX_SEG))
        .add("IsCruise", toMxLogicalPad(P.IsCruise, S, MATLAB_MAX_SEG))
        .add("SegMaxVel", toMxPad(P.SegMaxVel, 3, S, MATLAB_MAX_SEG, mxGetInf()))
        .add("SegMaxAccel", toMxPad(P.SegMaxAccel, 3, S, MATLAB_MAX_SEG, mxGetInf()))
        .add("LegCount", scal(P.LegCount))
        .add("LegStartSeg", toMxPad(P.LegStartSeg, 1, S, MATLAB_MAX_SEG, 0.0))
        .add("LegEndSeg", toMxPad(P.LegEndSeg, 1, S, MATLAB_MAX_SEG, 0.0))
        .add("HDGRef", toMx(P.HDGRef))
        .add("B3BC", writeKnotMask(P.B3BC, P.B3Fixed))
        .add("B3DotBC", writeKnotMask(P.B3DotBC, P.B3DotFixed))
        .build();
}

mxArray *writeQPTraj(const QPTrajT &Q) {
    const int S = POLYQP_MAX_SEG;
    const int K = POLYQP_MAX_KNOTS;
    return StructBuilder()
        .add("R", scal(Q.R))
        .add("Count", scal(Q.Count))
        .add("T", toMxPad(Q.T, 1, S, MATLAB_MAX_SEG, 0.0))
        .add("KeyTimes", toMxPad(Q.KeyTimes, 1, K, MATLAB_MAX_KNOTS, 0.0))
        .add("Coef", toMx3Pad(Q.Coef, nullptr, MATLAB_MAX_SEG, 3, 0.0))
        .add("IsHold", toMxLogicalPad(Q.IsHold, S, MATLAB_MAX_SEG))
        .add("IsCruise", toMxLogicalPad(Q.IsCruise, S, MATLAB_MAX_SEG))
        .add("Fix", toMx3Pad(Q.Fix, Q.FixFixed.data, 6, MATLAB_MAX_KNOTS, mxGetNaN()))
        .add("Cost", scal(Q.Cost))
        .add("T_total", scal(Q.T_total))
        .add("FlatOffset", scal(Q.FlatOffset))
        .build();
}

mxArray *writeReport(const ReportT &R, int M) {
    return StructBuilder()
        .add("Feasible", lscal(R.Feasible))
        .add("MaxOverallUtil", scal(R.MaxOverallUtil))
        .add("UtilVel", scal(R.UtilVel))
        .add("UtilAccel", scal(R.UtilAccel))
        .add("UtilThrustHi", scal(R.UtilThrustHi))
        .add("UtilThrustLo", scal(R.UtilThrustLo))
        .add("UtilGimbal", scal(R.UtilGimbal))
        .add("UtilGimbalRate", scal(R.UtilGimbalRate))
        .add("UtilThrustRate", scal(R.UtilThrustRate))
        .add("UtilRoll", scal(R.UtilRoll))
        .add("UtilTilt", scal(R.UtilTilt))
        .add("OuterIterations", scal(R.OuterIterations))
        .add("Cost", scal(R.Cost))
        .add("SegTimes", toMxCols(R.SegTimes, M))
        .add("T_total", scal(R.T_total))
        .build();
}

mxArray *writeEvalOut(const EvalOutT &E, int M) {
    return StructBuilder()
        .add("SegUtil", toMxCols(E.SegUtil, M))
        .add("MaxUtilPerType", toMx(E.MaxUtilPerType))
        .add("MaxOverallUtil", scal(E.MaxOverallUtil))
        .add("SegDilate", toMxCols(E.SegDilate, M))
        .build();
}

mxArray *writeWrench(const WrenchT &W, int n) {
    return StructBuilder().add("F_body", toMxCols(W.F_body, n)).add("Tau_body", toMxCols(W.Tau_body, n)).add("CGz", toMxCols(W.CGz, n)).build();
}

mxArray *writeWaypoint(const TOADWaypointT &wp) {
    return StructBuilder()
        .add("PosBC", writeVec3Mask(wp.PosBC, wp.PosFixed))
        .add("VelBC", writeVec3Mask(wp.VelBC, wp.VelFixed))
        .add("AccelBC", writeVec3Mask(wp.AccelBC, wp.AccelFixed))
        .add("JerkBC", writeVec3Mask(wp.JerkBC, wp.JerkFixed))
        .add("HoldTime", scal(wp.IsPassAndGo ? mxGetNaN() : static_cast<double>(wp.HoldTime)))
        .add("MaxVel", toMx(wp.MaxVel))
        .add("MaxAccel", toMx(wp.MaxAccel))
        .add("StraightLine", lscal(wp.StraightLine))
        .add("HDGRef", toMx(wp.HDGRef))
        .add("Label", mxCreateString(wp.Label))
        .add("IsPassAndGo", lscal(wp.IsPassAndGo))
        .add("B3BC", writeVec3Mask(wp.B3BC, wp.B3Fixed))
        .add("B3DotBC", writeVec3Mask(wp.B3DotBC, wp.B3DotFixed))
        .build();
}

mxArray *emptyStruct() { return mxCreateStructMatrix(0, 0, 0, nullptr); }
mxArray *empty() { return mxCreateDoubleMatrix(0, 0, mxREAL); }

// ------------------------------------------------------------------ STREAM: the MEX's NodeSink
struct MexSink {
    double T[MATLAB_MAX_NODES];
    double X[MATLAB_MAX_NODES][15];
    double U[MATLAB_MAX_NODES][4];
    int n = 0;
    static void put(void *ctx, int k, double t, const Vector15d &Xk, const Vector4d &Uk) {
        MexSink &s = *static_cast<MexSink *>(ctx);
        if (k < 1 || k > MATLAB_MAX_NODES) mexErrMsgIdAndTxt("polyqp_mex:sink", "sink row %d out of range", k);
        s.T[k - 1] = t;
        for (int i = 0; i < 15; ++i) s.X[k - 1][i] = Xk(i);
        for (int i = 0; i < 4; ++i) s.U[k - 1][i] = Uk(i);
        if (k > s.n) s.n = k;
    }
    NodeSink sink(int capacity) {
        n = 0;
        return NodeSink{this, capacity, &MexSink::put};
    }
    // MATLAB InitialGuess: Time 1 x Nn, X 15 x Nn, U 4 x Nc, dt_row 1 x Nc
    mxArray *initialGuess(const InitialGuessT &G) const {
        const int Nn = G.Nn;
        const int Nc = Nn - 1;
        mxArray *Time = mxCreateDoubleMatrix(1, static_cast<mwSize>(Nn), mxREAL);
        mxArray *Xm = mxCreateDoubleMatrix(15, static_cast<mwSize>(Nn), mxREAL);
        mxArray *Um = mxCreateDoubleMatrix(4, static_cast<mwSize>(Nc), mxREAL);
        mxArray *dt = mxCreateDoubleMatrix(1, static_cast<mwSize>(Nc), mxREAL);
        for (int k = 0; k < Nn; ++k) {
            mxGetPr(Time)[k] = T[k];
            for (int i = 0; i < 15; ++i) mxGetPr(Xm)[i + 15 * k] = X[k][i];
            if (k < Nc) {
                for (int i = 0; i < 4; ++i) mxGetPr(Um)[i + 4 * k] = U[k][i];
                mxGetPr(dt)[k] = G.dt;
            }
        }
        return StructBuilder().add("Time", Time).add("X", Xm).add("U", Um).add("T_total", scal(G.T_total)).add("dt_row", dt).build();
    }
    // The abort block's Time (Nn x 1), States (Nn x 15), Inputs (Nn x 4)
    void trajectory(mxArray *&Time, mxArray *&States, mxArray *&Inputs) const {
        Time = mxCreateDoubleMatrix(static_cast<mwSize>(n), 1, mxREAL);
        States = mxCreateDoubleMatrix(static_cast<mwSize>(n), 15, mxREAL);
        Inputs = mxCreateDoubleMatrix(static_cast<mwSize>(n), 4, mxREAL);
        for (int k = 0; k < n; ++k) {
            mxGetPr(Time)[k] = T[k];
            for (int i = 0; i < 15; ++i) mxGetPr(States)[k + n * i] = X[k][i];
            for (int i = 0; i < 4; ++i) mxGetPr(Inputs)[k + n * i] = U[k][i];
        }
    }
};
MexSink g_sink;

// Verification of the identity-padded solve: rcond of the leading N x N block alone (fixed-size, N = 1..5)
template <int N> double rcondTopLeft(const MatrixF_Fd &A) {
    const Eigen::Matrix<double, N, N> B = A.template topLeftCorner<N, N>();
    Eigen::PartialPivLU<Eigen::Matrix<double, N, N>> lu(B);
    return lu.rcond();
}
double rcondUnpadded(const MatrixF_Fd &A, int n) {
    switch (n) {
    case 1: return rcondTopLeft<1>(A);
    case 2: return rcondTopLeft<2>(A);
    case 3: return rcondTopLeft<3>(A);
    case 4: return rcondTopLeft<4>(A);
    case 5: return rcondTopLeft<5>(A);
    default: return mxGetNaN();
    }
}

// ------------------------------------------------------------------ output plumbing
struct Outputs {
    std::vector<mxArray *> v;
    void push(mxArray *a) { v.push_back(a); }
    void finish(PolyQPStatus st, int nlhs, mxArray *plhs[]) {
        if (st == PolyQPStatus::Ok) {
            v.push_back(mxCreateString(""));
            v.push_back(mxCreateString(""));
        } else {
            v.push_back(mxCreateString(lastError().identifier));
            v.push_back(mxCreateString(lastError().message));
        }
        const int nOut = nlhs < 1 ? 1 : nlhs;
        if (nOut > static_cast<int>(v.size())) mexErrMsgIdAndTxt("polyqp_mex:nargout", "Too many output arguments (max %d)", static_cast<int>(v.size()));
        for (int i = 0; i < static_cast<int>(v.size()); ++i) {
            if (i < nOut) {
                plhs[i] = v[static_cast<size_t>(i)];
            } else {
                mxDestroyArray(v[static_cast<size_t>(i)]);
            }
        }
    }
};

void needArgs(int nrhs, int n, const char *usage) {
    if (nrhs != n + 1) mexErrMsgIdAndTxt("polyqp_mex:usage", "Usage: %s", usage);
}

std::string getString(const mxArray *a) {
    char buf[256];
    if (!mxIsChar(a)) mexErrMsgIdAndTxt("polyqp_mex:type", "Expected a char array");
    mxGetString(a, buf, sizeof(buf));
    return std::string(buf);
}

} // namespace

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    boundsHook() = &boundsHookImpl;
    if (nrhs < 1 || !mxIsChar(prhs[0])) mexErrMsgIdAndTxt("polyqp_mex:usage", "First argument must be a function name");
    const std::string cmd = getString(prhs[0]);
    lastError() = PolyQPErrorRecord{};
    Outputs out;

    if (cmd == "PolyQP") {
        if (nrhs < 3 || nrhs > 4) mexErrMsgIdAndTxt("polyqp_mex:usage", "Usage: [InitialGuess, QPTraj, Report, status, msg] = polyqp_mex('PolyQP', Waypoints, Constants[, Opts])");
        static TOADWaypointArray W;
        readWaypoints(prhs[1], W);
        checkConstants(prhs[2]);
        if (nrhs == 4 && !mxIsEmpty(prhs[3]) && mxIsStruct(prhs[3]) && mxGetField(prhs[3], 0, "R") != nullptr && !mxIsEmpty(mxGetField(prhs[3], 0, "R")))
            checkOrder(static_cast<int>(num(mxGetField(prhs[3], 0, "R"), "Opts.R")));
        InitialGuessT G{};
        static QPTrajT Q;
        static ReportT R;
        // The PolyQP check compares with MATLAB, which samples up to 4000 nodes; the flight node capacity is checked by AbortGuidance.
        const PolyQPStatus st = PolyQP::PolyQP(W, g_sink.sink(MATLAB_MAX_NODES), G, Q, R);
        const bool ok = st == PolyQPStatus::Ok;
        out.push(ok ? g_sink.initialGuess(G) : emptyStruct());
        out.push(ok ? writeQPTraj(Q) : emptyStruct());
        out.push(ok ? writeReport(R, Q.Count) : emptyStruct());
        out.finish(st, nlhs, plhs);
    } else if (cmd == "resolveWaypoints") {
        needArgs(nrhs, 3, "[Plan, status, msg] = polyqp_mex('resolveWaypoints', Waypoints, r, Limits)");
        static TOADWaypointArray W;
        readWaypoints(prhs[1], W);
        checkOrder(static_cast<int>(num(prhs[2], "r")));
        const LimitsT L = readLimits(prhs[3]);
        static PlanT P;
        const PolyQPStatus st = resolveWaypoints(W, L, P);
        out.push(st == PolyQPStatus::Ok ? writePlan(P) : emptyStruct());
        out.finish(st, nlhs, plhs);
    } else if (cmd == "assembleQP") {
        needArgs(nrhs, 1, "[A0, Q0, H0, status, msg] = polyqp_mex('assembleQP', r)");
        Matrix10_10d A0, Q0, H0;
        checkOrder(static_cast<int>(num(prhs[1], "r")));
        const PolyQPStatus st = assembleQP(A0, Q0, H0);
        out.push(toMx(A0));
        out.push(toMx(Q0));
        out.push(toMx(H0));
        out.finish(st, nlhs, plhs);
    } else if (cmd == "solveLeg" || cmd == "solveLegRcond") {
        needArgs(nrhs, 4, "[CoefPhys, Deriv, Cost, status, msg] = polyqp_mex('solveLeg', FixLeg, SegTimes, A0, H0)");
        static Array3_5_Kd FixLeg;
        static Array3_5_Kb FixLegFixed;
        readArray3Mask(prhs[1], FixLeg, FixLegFixed, "FixLeg");
        checkOrder(FixLeg.n2);
        if (static_cast<int>(mxGetDimensions(prhs[1])[mxGetNumberOfDimensions(prhs[1]) == 3 ? 2 : 0]) > POLYQP_MAX_KNOTS && mxGetNumberOfDimensions(prhs[1]) == 3)
            mexErrMsgIdAndTxt("polyqp_mex:capacity", "FixLeg has more knots than POLYQP_MAX_KNOTS = %d", POLYQP_MAX_KNOTS);
        Matrix1_Sd SegTimes;
        readRowCount(prhs[2], SegTimes, "SegTimes");
        Matrix10_10d A0, H0;
        readCols(prhs[3], A0, 0.0, "A0");
        readCols(prhs[4], H0, 0.0, "H0");
        static Array10_S_3d CoefPhys;
        static Array3_5_Kd Deriv;
        double Cost = 0;
        CoefPhys.resize(0, 0, 0);
        Deriv.resize(0, 0, 0);
        lastRcond() = PolyQPRcondRecord{};
        const PolyQPStatus st = solveLeg(FixLeg, FixLegFixed, SegTimes, A0, H0, CoefPhys, Deriv, Cost);
        const bool ok = st == PolyQPStatus::Ok;
        if (cmd == "solveLeg") {
            out.push(ok ? toMx3(CoefPhys) : empty());
            out.push(ok ? toMx3(Deriv) : empty());
            out.push(ok ? scal(Cost) : empty());
        } else {
            // [rcGuard, nFree, rcUnpadded, rcPadded, status, msg]: per axis, the rcond solveLeg's guard used, the free count, the rcond
            // of the unpadded nFree x nFree block computed here independently (NaN above 5), and Eigen's estimate on the padded system
            mxArray *rg = mxCreateDoubleMatrix(1, 3, mxREAL);
            mxArray *nf = mxCreateDoubleMatrix(1, 3, mxREAL);
            mxArray *ru = mxCreateDoubleMatrix(1, 3, mxREAL);
            mxArray *rp = mxCreateDoubleMatrix(1, 3, mxREAL);
            for (int d = 0; d < 3; ++d) {
                const int nP = lastRcond().nP[d];
                mxGetPr(nf)[d] = nP;
                mxGetPr(rg)[d] = nP > 0 ? lastRcond().rc[d] : mxGetNaN();
                mxGetPr(ru)[d] = nP > 0 ? rcondUnpadded(lastRcond().R_PP_eq[d], nP) : mxGetNaN();
                mxGetPr(rp)[d] = nP > 0 ? lastRcond().rcPadded[d] : mxGetNaN();
            }
            out.push(rg);
            out.push(nf);
            out.push(ru);
            out.push(rp);
        }
        out.finish(st, nlhs, plhs);
    } else if (cmd == "allocateTimes") {
        needArgs(nrhs, 5, "[SegTimes, Cost, status, msg] = polyqp_mex('allocateTimes', FixLeg, SegTimesInit, A0, H0, MinSegTime)");
        static Array3_5_Kd FixLeg;
        static Array3_5_Kb FixLegFixed;
        readArray3Mask(prhs[1], FixLeg, FixLegFixed, "FixLeg");
        checkOrder(FixLeg.n2);
        Matrix1_Sd SegTimesInit;
        const int M = readRowCount(prhs[2], SegTimesInit, "SegTimesInit");
        Matrix10_10d A0, H0;
        readCols(prhs[3], A0, 0.0, "A0");
        readCols(prhs[4], H0, 0.0, "H0");
        Matrix1_Sd SegTimes = Matrix1_Sd::Zero();
        double Cost = 0;
        const PolyQPStatus st = allocateTimes(FixLeg, FixLegFixed, SegTimesInit, A0, H0, num(prhs[5], "MinSegTime"), SegTimes, Cost);
        const bool ok = st == PolyQPStatus::Ok;
        out.push(ok ? toMxCols(SegTimes, M) : empty());
        out.push(ok ? scal(Cost) : empty());
        out.finish(st, nlhs, plhs);
    } else if (cmd == "evaluateConstraints") {
        needArgs(nrhs, 4, "[EvalOut, Feasible, status, msg] = polyqp_mex('evaluateConstraints', QPTraj, Plan, Constants, Limits)");
        static QPTrajT Q;
        readQPTraj(prhs[1], Q);
        static PlanT P;
        readPlan(prhs[2], P);
        checkConstants(prhs[3]);
        const LimitsT L = readLimits(prhs[4]);
        static EvalOutT E;
        bool Feasible = false;
        const PolyQPStatus st = evaluateConstraints(Q, P, L, E, Feasible);
        const bool ok = st == PolyQPStatus::Ok;
        out.push(ok ? writeEvalOut(E, Q.Count) : emptyStruct());
        out.push(lscal(Feasible));
        out.finish(st, nlhs, plhs);
    } else if (cmd == "flatnessInvert") {
        needArgs(nrhs, 5, "[X, Wrench, AccelCg, status, msg] = polyqp_mex('flatnessInvert', QPTraj, t_vec, HDGRef, Constants, Limits)");
        static QPTrajT Q;
        readQPTraj(prhs[1], Q);
        static Matrix1_Nd t_vec;
        const int n = readGridRow(prhs[2], t_vec, "t_vec");
        const Vector3d HDGRef = readVec3(prhs[3], "HDGRef");
        checkConstants(prhs[4]);
        const LimitsT L = readLimits(prhs[5]);
        static Matrix15_Nd X;
        static WrenchT Wr;
        static Matrix3_Nd AccelCg;
        const PolyQPStatus st = flatnessInvert(Q, t_vec, n, HDGRef, L, X, Wr, AccelCg);
        const bool ok = st == PolyQPStatus::Ok;
        out.push(ok ? toMxCols(X, n) : empty());
        out.push(ok ? writeWrench(Wr, n) : emptyStruct());
        out.push(ok ? toMxCols(AccelCg, n) : empty());
        out.finish(st, nlhs, plhs);
    } else if (cmd == "allocateControls") {
        needArgs(nrhs, 2, "[U, status, msg] = polyqp_mex('allocateControls', Wrench, Limits)");
        static WrenchT Wr;
        const int n = readWrench(prhs[1], Wr);
        const LimitsT L = readLimits(prhs[2]);
        static Matrix4_Nd U;
        const PolyQPStatus st = allocateControls(Wr, n, L, U);
        out.push(toMxCols(U, n));
        out.finish(st, nlhs, plhs);
    } else if (cmd == "sampleOutputs") {
        needArgs(nrhs, 4, "[InitialGuess, status, msg] = polyqp_mex('sampleOutputs', QPTraj, Plan, Constants, Limits)");
        static QPTrajT Q;
        readQPTraj(prhs[1], Q);
        static PlanT P;
        readPlan(prhs[2], P);
        checkConstants(prhs[3]);
        const LimitsT L = readLimits(prhs[4]);
        InitialGuessT G{};
        const PolyQPStatus st = sampleOutputs(Q, P, L, g_sink.sink(MATLAB_MAX_NODES), G);
        out.push(st == PolyQPStatus::Ok ? g_sink.initialGuess(G) : emptyStruct());
        out.finish(st, nlhs, plhs);
    } else if (cmd == "PolyEval") {
        if (nrhs < 3 || nrhs > 4) mexErrMsgIdAndTxt("polyqp_mex:usage", "Usage: [Out, status, msg] = polyqp_mex('PolyEval', QPTraj, t[, k])");
        static QPTrajT Q;
        readQPTraj(prhs[1], Q);
        static Matrix1_Nd t;
        const int n = readGridRow(prhs[2], t, "t");
        static Matrix3_Nd O;
        PolyQPStatus st;
        if (nrhs == 4) {
            st = PolyEval(Q, t, n, O, static_cast<int>(num(prhs[3], "k")));
        } else {
            st = PolyEval(Q, t, n, O);
        }
        out.push(toMxCols(O, n));
        out.finish(st, nlhs, plhs);
    } else if (cmd == "buildVehicleLimits") {
        if (nrhs < 2 || nrhs > 3) mexErrMsgIdAndTxt("polyqp_mex:usage", "Usage: [Limits, status, msg] = polyqp_mex('buildVehicleLimits', Constants[, Opts])");
        LimitsT L;
        checkConstants(prhs[1]);
        const PolyQPStatus st = buildVehicleLimits(L);
        out.push(writeLimits(L));
        out.finish(st, nlhs, plhs);
    } else if (cmd == "TOADWaypoint") {
        if (nrhs < 2 || (nrhs % 2) != 0) mexErrMsgIdAndTxt("polyqp_mex:usage", "Usage: [wp, status, msg] = polyqp_mex('TOADWaypoint', pos, 'Name', value, ...)");
        Vector3d pos;
        Vector3b posFixed;
        readVec3Mask(prhs[1], "pos", pos, posFixed);
        TOADWaypointOptions options;
        static char labelBuf[256];
        labelBuf[0] = '\0';
        for (int i = 2; i < nrhs; i += 2) {
            const std::string name = getString(prhs[i]);
            const mxArray *v = prhs[i + 1];
            if (name == "VelBC") readVec3Mask(v, "VelBC", options.VelBC, options.VelFixed);
            else if (name == "AccelBC") readVec3Mask(v, "AccelBC", options.AccelBC, options.AccelFixed);
            else if (name == "JerkBC") readVec3Mask(v, "JerkBC", options.JerkBC, options.JerkFixed);
            else if (name == "HoldTime") {
                const double h = num(v, "HoldTime");
                options.HoldTimeGiven = !std::isnan(h);
                options.HoldTime = options.HoldTimeGiven ? h : 0.0;
            } else if (name == "MaxVel") options.MaxVel = readVec3(v, "MaxVel");
            else if (name == "MaxAccel") options.MaxAccel = readVec3(v, "MaxAccel");
            else if (name == "StraightLine") options.StraightLine = num(v, "StraightLine") != 0.0;
            else if (name == "HDGRef") options.HDGRef = readVec3(v, "HDGRef");
            else if (name == "B3BC") readVec3Mask(v, "B3BC", options.B3BC, options.B3Fixed);
            else if (name == "B3DotBC") readVec3Mask(v, "B3DotBC", options.B3DotBC, options.B3DotFixed);
            else if (name == "Label") {
                mxGetString(v, labelBuf, sizeof(labelBuf));
                options.Label = labelBuf;
            } else mexErrMsgIdAndTxt("polyqp_mex:usage", "Unknown TOADWaypoint option '%s'", name.c_str());
        }
        out.push(writeWaypoint(TOADWaypoint(pos, options, posFixed)));
        out.finish(PolyQPStatus::Ok, nlhs, plhs);
    } else if (cmd == "ComputeJtot") {
        needArgs(nrhs, 3, "[J_tot, CGz, status, msg] = polyqp_mex('ComputeJtot', m_lox, m_ipa, constants6DoF)");
        Matrix3_3d J_tot;
        double CGz = 0;
        checkConstants(prhs[3]);
        ::ComputeJtot(num(prhs[1], "m_lox"), num(prhs[2], "m_ipa"), J_tot, CGz);
        out.push(toMx(J_tot));
        out.push(scal(CGz));
        out.finish(PolyQPStatus::Ok, nlhs, plhs);
    } else if (cmd == "DCM_Quat_Conversion") {
        needArgs(nrhs, 1, "[q, status, msg] = polyqp_mex('DCM_Quat_Conversion', R)");
        Matrix3_3d R;
        readCols(prhs[1], R, 0.0, "R");
        out.push(toMx(::DCM_Quat_Conversion(R)));
        out.finish(PolyQPStatus::Ok, nlhs, plhs);
    } else if (cmd == "AbortGuidance") {
        if (nrhs < 4 || nrhs > 8)
            mexErrMsgIdAndTxt("polyqp_mex:usage", "Usage: [Time, States, Inputs, tAbort, abortCalculated, Report, abortStatus, planStatus, msg] = "
                                                  "polyqp_mex('AbortGuidance', t, X_est, constantsTOAD, state[, ABORT[, u_cmd[, a_dist]]])");
        const double t = num(prhs[1], "t");
        Vector15d X_est;
        readVecN(prhs[2], X_est, "X_est");
        checkConstants(prhs[3]);
        AbortGuidanceState state;
        if (nrhs >= 5 && mxIsStruct(prhs[4]) && !mxIsEmpty(prhs[4])) {
            state.abortCalculated = numField(prhs[4], "abortCalculated") != 0.0;
            state.tAbort = numField(prhs[4], "tAbort");
            if (mxGetField(prhs[4], 0, "abortStatus") != nullptr) state.abortStatus = static_cast<AbortStatus>(static_cast<int>(numField(prhs[4], "abortStatus")));
        }
        static ReportT Rep;
        const bool before = state.abortCalculated;
        // STATIC: the flight node capacity, as the harness gives it (Abort::MAX_POINTS = POLYQP_MAX_NODES)
        const NodeSink sk = g_sink.sink(POLYQP_MAX_NODES);
        // FIX-F02: optional u_cmd (4) and a_dist (3); MATLAB defaults: hover thrust with the gimbal centred, no disturbance
        const double ABORT = nrhs >= 6 ? num(prhs[5], "ABORT") : 0.0;
        Vector4d u_cmd(0, 0, ASTRA::m_dry * ASTRA::g, 0);
        Vector3d a_dist = Vector3d::Zero();
        if (nrhs >= 7 && !mxIsEmpty(prhs[6])) readVecN(prhs[6], u_cmd, "u_cmd");
        if (nrhs >= 8 && !mxIsEmpty(prhs[7])) readVecN(prhs[7], a_dist, "a_dist");
        const AbortStatus ast = AbortGuidance(t, X_est, state, sk, Rep, ABORT, u_cmd, a_dist);
        const bool fired = state.abortCalculated && !before;
        mxArray *Time = empty();
        mxArray *States = empty();
        mxArray *Inputs = empty();
        if (fired) {
            mxDestroyArray(Time);
            mxDestroyArray(States);
            mxDestroyArray(Inputs);
            g_sink.trajectory(Time, States, Inputs);
        }
        out.push(Time);
        out.push(States);
        out.push(Inputs);
        out.push(scal(state.tAbort));
        out.push(lscal(state.abortCalculated));
        out.push(fired && state.planStatus == PolyQPStatus::Ok ? writeReport(Rep, workspace().AbortGuidance.QPTraj.Count) : emptyStruct());
        out.push(scal(static_cast<double>(static_cast<int>(ast))));
        out.finish(state.planStatus, nlhs, plhs);
    } else if (cmd == "quatRot") {
        needArgs(nrhs, 1, "[CIB, status, msg] = polyqp_mex('quatRot', q)");
        Vector4d qv;
        readVecN(prhs[1], qv, "q");
        out.push(toMx(::quatRot(qv)));
        out.finish(PolyQPStatus::Ok, nlhs, plhs);
    } else if (cmd == "abortFallback") {
        if (nrhs != 4) mexErrMsgIdAndTxt("polyqp_mex:usage", "Usage: [Time, States, Inputs, status, msg] = polyqp_mex('abortFallback', X_est, Constants, vCapZ)");
        Vector15d X_est;
        readVecN(prhs[1], X_est, "X_est");
        checkConstants(prhs[2]);
        abortFallback(X_est, num(prhs[3], "vCapZ"), g_sink.sink(MATLAB_MAX_NODES));
        mxArray *Time;
        mxArray *States;
        mxArray *Inputs;
        g_sink.trajectory(Time, States, Inputs);
        out.push(Time);
        out.push(States);
        out.push(Inputs);
        out.finish(PolyQPStatus::Ok, nlhs, plhs);
    } else if (cmd == "capacities") {
        // [MAX_WP, MAX_SEG, MAX_KNOTS, MAX_FREE, MAX_NODES]
        mxArray *c = mxCreateDoubleMatrix(1, 5, mxREAL);
        const double v[5] = {POLYQP_MAX_WP, POLYQP_MAX_SEG, POLYQP_MAX_KNOTS, POLYQP_MAX_FREE, POLYQP_MAX_NODES};
        for (int i = 0; i < 5; ++i) mxGetPr(c)[i] = v[i];
        out.push(c);
        out.finish(PolyQPStatus::Ok, nlhs, plhs);
    } else {
        mexErrMsgIdAndTxt("polyqp_mex:usage", "Unknown function '%s'", cmd.c_str());
    }
}
