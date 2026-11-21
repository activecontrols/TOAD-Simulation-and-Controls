// Verification only (never part of the flight build): MATLAB MEX dispatcher exposing every public PolyQP function.
// Usage: [out..., status, msg] = polyqp_mex('<FunctionName>', args...)
// status is '' on success or the MATLAB error identifier (e.g. 'PolyQP:SingularRPP'); msg is the formatted message.
#include "mex.h"

#include "AbortGuidance.h"
#include "PolyQP.h"

#include <cstring>
#include <string>
#include <vector>

using namespace PolyQP;

namespace {

void boundsHookImpl(const char *file, int line) { mexErrMsgIdAndTxt("polyqp_mex:bounds", "PolyQP bounds check failed at %s:%d", file, line); }

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

template <typename MatT> void readMat(const mxArray *a, MatT &M, const char *what) {
    const size_t rows = mxGetM(a);
    const size_t cols = mxGetN(a);
    const int maxR = MatT::MaxRowsAtCompileTime;
    const int maxC = MatT::MaxColsAtCompileTime;
    if (static_cast<int>(rows) > maxR || static_cast<int>(cols) > maxC)
        mexErrMsgIdAndTxt("polyqp_mex:capacity", "'%s' is %dx%d, capacity %dx%d", what, static_cast<int>(rows), static_cast<int>(cols), maxR, maxC);
    if ((MatT::RowsAtCompileTime != Eigen::Dynamic && static_cast<int>(rows) != MatT::RowsAtCompileTime) ||
        (MatT::ColsAtCompileTime != Eigen::Dynamic && static_cast<int>(cols) != MatT::ColsAtCompileTime))
        mexErrMsgIdAndTxt("polyqp_mex:size", "'%s' has the wrong size (%dx%d)", what, static_cast<int>(rows), static_cast<int>(cols));
    M.resize(static_cast<Eigen::Index>(rows), static_cast<Eigen::Index>(cols));
    for (size_t j = 0; j < cols; ++j)
        for (size_t i = 0; i < rows; ++i) M(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = static_cast<typename MatT::Scalar>(elem(a, i + j * rows));
}
template <typename MatT> void readMatField(const mxArray *s, const char *name, MatT &M, mwIndex i = 0) { readMat(getField(s, name, i), M, name); }

template <typename VecT> void readVecAny(const mxArray *a, VecT &v, const char *what) {
    const size_t n = mxGetNumberOfElements(a);
    if (static_cast<int>(n) > VecT::SizeAtCompileTime && VecT::SizeAtCompileTime != Eigen::Dynamic)
        mexErrMsgIdAndTxt("polyqp_mex:size", "'%s' must have %d elements", what, static_cast<int>(VecT::SizeAtCompileTime));
    if (VecT::SizeAtCompileTime != Eigen::Dynamic && static_cast<int>(n) != VecT::SizeAtCompileTime)
        mexErrMsgIdAndTxt("polyqp_mex:size", "'%s' must have %d elements", what, static_cast<int>(VecT::SizeAtCompileTime));
    if (static_cast<int>(n) > VecT::MaxSizeAtCompileTime) mexErrMsgIdAndTxt("polyqp_mex:capacity", "'%s' exceeds capacity", what);
    v.resize(static_cast<Eigen::Index>(n));
    for (size_t i = 0; i < n; ++i) v(static_cast<Eigen::Index>(i)) = static_cast<typename VecT::Scalar>(elem(a, i));
}

Vec3 readVec3(const mxArray *a, const char *what) {
    Vec3 v;
    readVecAny(a, v, what);
    return v;
}

template <typename A3> void readArray3(const mxArray *a, A3 &A, const char *what) {
    const mwSize nd = mxGetNumberOfDimensions(a);
    const mwSize *dims = mxGetDimensions(a);
    if (nd > 3) mexErrMsgIdAndTxt("polyqp_mex:size", "'%s' must have at most 3 dimensions", what);
    const int n1 = static_cast<int>(dims[0]);
    const int n2 = static_cast<int>(dims[1]);
    const int n3 = nd == 3 ? static_cast<int>(dims[2]) : 1;
    if (n1 > A3::Cap1 || n2 > A3::Cap2 || n3 > A3::Cap3)
        mexErrMsgIdAndTxt("polyqp_mex:capacity", "'%s' is %dx%dx%d, capacity %dx%dx%d", what, n1, n2, n3, A3::Cap1, A3::Cap2, A3::Cap3);
    A.resize(n1, n2, n3);
    for (int q = 0; q < n1 * n2 * n3; ++q) A.data[q] = static_cast<decltype(A.data[0] + 0)>(elem(a, static_cast<size_t>(q)));
}

template <typename BoolArr> void readBoolRow(const mxArray *a, BoolArr &b, const char *what) {
    const size_t n = mxGetNumberOfElements(a);
    if (static_cast<int>(n) > BoolArr::MaxSizeAtCompileTime) mexErrMsgIdAndTxt("polyqp_mex:capacity", "'%s' exceeds capacity", what);
    b.resize(static_cast<Eigen::Index>(n));
    for (size_t i = 0; i < n; ++i) b(static_cast<Eigen::Index>(i)) = elem(a, i) != 0.0;
}
template <typename IntArr> void readIntRow(const mxArray *a, IntArr &b, const char *what) {
    const size_t n = mxGetNumberOfElements(a);
    if (static_cast<int>(n) > IntArr::MaxSizeAtCompileTime) mexErrMsgIdAndTxt("polyqp_mex:capacity", "'%s' exceeds capacity", what);
    b.resize(static_cast<Eigen::Index>(n));
    for (size_t i = 0; i < n; ++i) b(static_cast<Eigen::Index>(i)) = static_cast<int>(elem(a, i));
}

// ASTRA only: the solver's constants are compiled in (ASTRA_PolyQP_Constants.h). The MATLAB constants struct is still
// passed by the verification scripts; it must be LoadTOADParams(0), so check every value the solver uses.
void checkConstants(const mxArray *s) {
    Mat3 J;
    readMatField(s, "J", J);
    const Mat3 Jc = Vec3(ASTRA::J_xx, ASTRA::J_yy, ASTRA::J_zz).asDiagonal();
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
    L.g = static_cast<Real>(numField(s, "g"));
    L.m_init = static_cast<Real>(numField(s, "m_init"));
    L.m_dry = static_cast<Real>(numField(s, "m_dry"));
    L.MaxThrust = static_cast<Real>(numField(s, "MaxThrust"));
    L.rTB = static_cast<Real>(numField(s, "rTB"));
    L.CGz_init = static_cast<Real>(numField(s, "CGz_init"));
    readMatField(s, "J_init", L.J_init);
    L.Rho = static_cast<Real>(numField(s, "Rho"));
    L.T_lo = static_cast<Real>(numField(s, "T_lo"));
    L.T_hi = static_cast<Real>(numField(s, "T_hi"));
    L.GimbalLimit = static_cast<Real>(numField(s, "GimbalLimit"));
    L.GimbalRateLim = static_cast<Real>(numField(s, "GimbalRateLim"));
    L.ThrustRateLim = static_cast<Real>(numField(s, "ThrustRateLim"));
    L.RollLimit = static_cast<Real>(numField(s, "RollLimit"));
    L.RollRateLim = static_cast<Real>(numField(s, "RollRateLim"));
    L.MaxTilt = static_cast<Real>(numField(s, "MaxTilt"));
    L.DefaultMaxAccel = readVec3(getField(s, "DefaultMaxAccel"), "DefaultMaxAccel");
    L.MinSegTime = static_cast<Real>(numField(s, "MinSegTime"));
    L.DtSample = static_cast<Real>(numField(s, "DtSample"));
    return L;
}

void readPlan(const mxArray *s, PlanT &P) {
    P.R = static_cast<int>(numField(s, "R"));
    P.KnotCount = static_cast<int>(numField(s, "KnotCount"));
    P.SegCount = static_cast<int>(numField(s, "SegCount"));
    readArray3(getField(s, "FixCg"), P.FixCg, "FixCg");
    readBoolRow(getField(s, "KnotIsRest"), P.KnotIsRest, "KnotIsRest");
    readMatField(s, "SegT", P.SegT);
    readBoolRow(getField(s, "IsHold"), P.IsHold, "IsHold");
    readBoolRow(getField(s, "IsCruise"), P.IsCruise, "IsCruise");
    readMatField(s, "SegMaxVel", P.SegMaxVel);
    readMatField(s, "SegMaxAccel", P.SegMaxAccel);
    P.LegCount = static_cast<int>(numField(s, "LegCount"));
    readIntRow(getField(s, "LegStartSeg"), P.LegStartSeg, "LegStartSeg");
    readIntRow(getField(s, "LegEndSeg"), P.LegEndSeg, "LegEndSeg");
    P.HDGRef = readVec3(getField(s, "HDGRef"), "HDGRef");
    readMatField(s, "B3BC", P.B3BC);
    readMatField(s, "B3DotBC", P.B3DotBC);
}

void readQPTraj(const mxArray *s, QPTrajT &Q) {
    Q.R = static_cast<int>(numField(s, "R"));
    Q.Count = static_cast<int>(numField(s, "Count"));
    readMatField(s, "T", Q.T);
    readMatField(s, "KeyTimes", Q.KeyTimes);
    readArray3(getField(s, "Coef"), Q.Coef, "Coef");
    readBoolRow(getField(s, "IsHold"), Q.IsHold, "IsHold");
    readBoolRow(getField(s, "IsCruise"), Q.IsCruise, "IsCruise");
    readArray3(getField(s, "Fix"), Q.Fix, "Fix");
    Q.Cost = static_cast<Real>(numField(s, "Cost"));
    Q.T_total = static_cast<Real>(numField(s, "T_total"));
    Q.FlatOffset = static_cast<Real>(numField(s, "FlatOffset"));
}

void readWrench(const mxArray *s, WrenchT &W) {
    readMatField(s, "F_body", W.F_body);
    readMatField(s, "Tau_body", W.Tau_body);
    readMatField(s, "CGz", W.CGz);
}

void readWaypoints(const mxArray *s, TOADWaypointArray &W) {
    const size_t n = mxGetNumberOfElements(s);
    if (static_cast<int>(n) > POLYQP_MAX_WP) mexErrMsgIdAndTxt("polyqp_mex:capacity", "Waypoint count %d exceeds POLYQP_MAX_WP", static_cast<int>(n));
    W.count = static_cast<int>(n);
    for (size_t i = 0; i < n; ++i) {
        TOADWaypointT &wp = W.data[i];
        const mwIndex k = static_cast<mwIndex>(i);
        wp.PosBC = readVec3(getField(s, "PosBC", k), "PosBC");
        wp.VelBC = readVec3(getField(s, "VelBC", k), "VelBC");
        wp.AccelBC = readVec3(getField(s, "AccelBC", k), "AccelBC");
        wp.JerkBC = readVec3(getField(s, "JerkBC", k), "JerkBC");
        wp.HoldTime = static_cast<Real>(numField(s, "HoldTime", k));
        wp.MaxVel = readVec3(getField(s, "MaxVel", k), "MaxVel");
        wp.MaxAccel = readVec3(getField(s, "MaxAccel", k), "MaxAccel");
        wp.StraightLine = numField(s, "StraightLine", k) != 0.0;
        wp.HDGRef = readVec3(getField(s, "HDGRef", k), "HDGRef");
        std::memset(wp.Label, 0, sizeof(wp.Label));
        mxGetString(getField(s, "Label", k), wp.Label, sizeof(wp.Label));
        wp.IsPassAndGo = numField(s, "IsPassAndGo", k) != 0.0;
        wp.B3BC = readVec3(getField(s, "B3BC", k), "B3BC");
        wp.B3DotBC = readVec3(getField(s, "B3DotBC", k), "B3DotBC");
    }
}

// ------------------------------------------------------------------ writers
mxArray *scal(double v) { return mxCreateDoubleScalar(v); }
mxArray *lscal(bool v) { return mxCreateLogicalScalar(v); }

template <typename Derived> mxArray *toMx(const Eigen::DenseBase<Derived> &M) {
    mxArray *a = mxCreateDoubleMatrix(static_cast<mwSize>(M.rows()), static_cast<mwSize>(M.cols()), mxREAL);
    double *p = mxGetPr(a);
    for (Eigen::Index j = 0; j < M.cols(); ++j)
        for (Eigen::Index i = 0; i < M.rows(); ++i) p[i + j * M.rows()] = static_cast<double>(M.derived().coeff(i, j));
    return a;
}
template <typename Derived> mxArray *toMxLogical(const Eigen::DenseBase<Derived> &M) {
    mxArray *a = mxCreateLogicalMatrix(static_cast<mwSize>(M.rows()), static_cast<mwSize>(M.cols()));
    mxLogical *p = mxGetLogicals(a);
    for (Eigen::Index j = 0; j < M.cols(); ++j)
        for (Eigen::Index i = 0; i < M.rows(); ++i) p[i + j * M.rows()] = M.derived().coeff(i, j) ? 1 : 0;
    return a;
}
template <typename Derived> mxArray *toMxIntAsDouble(const Eigen::DenseBase<Derived> &M) {
    mxArray *a = mxCreateDoubleMatrix(static_cast<mwSize>(M.rows()), static_cast<mwSize>(M.cols()), mxREAL);
    double *p = mxGetPr(a);
    for (Eigen::Index j = 0; j < M.cols(); ++j)
        for (Eigen::Index i = 0; i < M.rows(); ++i) p[i + j * M.rows()] = static_cast<double>(M.derived().coeff(i, j));
    return a;
}
template <typename A3> mxArray *toMx3(const A3 &A) {
    const mwSize dims[3] = {static_cast<mwSize>(A.n1), static_cast<mwSize>(A.n2), static_cast<mwSize>(A.n3)};
    mxArray *a = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double *p = mxGetPr(a);
    for (int q = 0; q < A.numel(); ++q) p[q] = static_cast<double>(A.data[q]);
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
    return StructBuilder()
        .add("R", scal(P.R))
        .add("KnotCount", scal(P.KnotCount))
        .add("SegCount", scal(P.SegCount))
        .add("FixCg", toMx3(P.FixCg))
        .add("KnotIsRest", toMxLogical(P.KnotIsRest))
        .add("SegT", toMx(P.SegT))
        .add("IsHold", toMxLogical(P.IsHold))
        .add("IsCruise", toMxLogical(P.IsCruise))
        .add("SegMaxVel", toMx(P.SegMaxVel))
        .add("SegMaxAccel", toMx(P.SegMaxAccel))
        .add("LegCount", scal(P.LegCount))
        .add("LegStartSeg", toMxIntAsDouble(P.LegStartSeg))
        .add("LegEndSeg", toMxIntAsDouble(P.LegEndSeg))
        .add("HDGRef", toMx(P.HDGRef))
        .add("B3BC", toMx(P.B3BC))
        .add("B3DotBC", toMx(P.B3DotBC))
        .build();
}

mxArray *writeQPTraj(const QPTrajT &Q) {
    return StructBuilder()
        .add("R", scal(Q.R))
        .add("Count", scal(Q.Count))
        .add("T", toMx(Q.T))
        .add("KeyTimes", toMx(Q.KeyTimes))
        .add("Coef", toMx3(Q.Coef))
        .add("IsHold", toMxLogical(Q.IsHold))
        .add("IsCruise", toMxLogical(Q.IsCruise))
        .add("Fix", toMx3(Q.Fix))
        .add("Cost", scal(Q.Cost))
        .add("T_total", scal(Q.T_total))
        .add("FlatOffset", scal(Q.FlatOffset))
        .build();
}

mxArray *writeInitialGuess(const InitialGuessT &G) {
    return StructBuilder()
        .add("Time", toMx(G.Time))
        .add("X", toMx(G.X))
        .add("U", toMx(G.U))
        .add("T_total", scal(G.T_total))
        .add("dt_row", toMx(G.dt_row))
        .build();
}

mxArray *writeReport(const ReportT &R) {
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
        .add("SegTimes", toMx(R.SegTimes))
        .add("T_total", scal(R.T_total))
        .build();
}

mxArray *writeEvalOut(const EvalOutT &E) {
    return StructBuilder()
        .add("SegUtil", toMx(E.SegUtil))
        .add("MaxUtilPerType", toMx(E.MaxUtilPerType))
        .add("MaxOverallUtil", scal(E.MaxOverallUtil))
        .add("SegDilate", toMx(E.SegDilate))
        .build();
}

mxArray *writeWrench(const WrenchT &W) {
    return StructBuilder().add("F_body", toMx(W.F_body)).add("Tau_body", toMx(W.Tau_body)).add("CGz", toMx(W.CGz)).build();
}

mxArray *writeWaypoint(const TOADWaypointT &wp) {
    return StructBuilder()
        .add("PosBC", toMx(wp.PosBC))
        .add("VelBC", toMx(wp.VelBC))
        .add("AccelBC", toMx(wp.AccelBC))
        .add("JerkBC", toMx(wp.JerkBC))
        .add("HoldTime", scal(wp.HoldTime))
        .add("MaxVel", toMx(wp.MaxVel))
        .add("MaxAccel", toMx(wp.MaxAccel))
        .add("StraightLine", lscal(wp.StraightLine))
        .add("HDGRef", toMx(wp.HDGRef))
        .add("Label", mxCreateString(wp.Label))
        .add("IsPassAndGo", lscal(wp.IsPassAndGo))
        .add("B3BC", toMx(wp.B3BC))
        .add("B3DotBC", toMx(wp.B3DotBC))
        .build();
}

mxArray *emptyStruct() { return mxCreateStructMatrix(0, 0, 0, nullptr); }

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
        static InitialGuessT G;
        static QPTrajT Q;
        static ReportT R;
        const PolyQPStatus st = PolyQP::PolyQP(W, G, Q, R);
        const bool ok = st == PolyQPStatus::Ok;
        out.push(ok ? writeInitialGuess(G) : emptyStruct());
        out.push(ok ? writeQPTraj(Q) : emptyStruct());
        out.push(ok ? writeReport(R) : emptyStruct());
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
        MatNC A0, Q0, H0;
        checkOrder(static_cast<int>(num(prhs[1], "r")));
        const PolyQPStatus st = assembleQP(A0, Q0, H0);
        out.push(toMx(A0));
        out.push(toMx(Q0));
        out.push(toMx(H0));
        out.finish(st, nlhs, plhs);
    } else if (cmd == "solveLeg") {
        needArgs(nrhs, 4, "[CoefPhys, Deriv, Cost, status, msg] = polyqp_mex('solveLeg', FixLeg, SegTimes, A0, H0)");
        static Fix3 FixLeg;
        readArray3(prhs[1], FixLeg, "FixLeg");
        RowSeg SegTimes;
        readVecAny(prhs[2], SegTimes, "SegTimes");
        MatNC A0, H0;
        readMat(prhs[3], A0, "A0");
        readMat(prhs[4], H0, "H0");
        static Coef3 CoefPhys;
        static Fix3 Deriv;
        Real Cost = 0;
        CoefPhys.resize(0, 0, 0);
        Deriv.resize(0, 0, 0);
        const PolyQPStatus st = solveLeg(FixLeg, SegTimes, A0, H0, CoefPhys, Deriv, Cost);
        const bool ok = st == PolyQPStatus::Ok;
        out.push(ok ? toMx3(CoefPhys) : mxCreateDoubleMatrix(0, 0, mxREAL));
        out.push(ok ? toMx3(Deriv) : mxCreateDoubleMatrix(0, 0, mxREAL));
        out.push(ok ? scal(Cost) : mxCreateDoubleMatrix(0, 0, mxREAL));
        out.finish(st, nlhs, plhs);
    } else if (cmd == "allocateTimes") {
        needArgs(nrhs, 5, "[SegTimes, Cost, status, msg] = polyqp_mex('allocateTimes', FixLeg, SegTimesInit, A0, H0, MinSegTime)");
        static Fix3 FixLeg;
        readArray3(prhs[1], FixLeg, "FixLeg");
        RowSeg SegTimesInit;
        readVecAny(prhs[2], SegTimesInit, "SegTimesInit");
        MatNC A0, H0;
        readMat(prhs[3], A0, "A0");
        readMat(prhs[4], H0, "H0");
        RowSeg SegTimes;
        Real Cost = 0;
        const PolyQPStatus st = allocateTimes(FixLeg, SegTimesInit, A0, H0, static_cast<Real>(num(prhs[5], "MinSegTime")), SegTimes, Cost);
        const bool ok = st == PolyQPStatus::Ok;
        out.push(ok ? toMx(SegTimes) : mxCreateDoubleMatrix(0, 0, mxREAL));
        out.push(ok ? scal(Cost) : mxCreateDoubleMatrix(0, 0, mxREAL));
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
        out.push(ok ? writeEvalOut(E) : emptyStruct());
        out.push(lscal(Feasible));
        out.finish(st, nlhs, plhs);
    } else if (cmd == "flatnessInvert") {
        needArgs(nrhs, 5, "[X, Wrench, AccelCg, status, msg] = polyqp_mex('flatnessInvert', QPTraj, t_vec, HDGRef, Constants, Limits)");
        static QPTrajT Q;
        readQPTraj(prhs[1], Q);
        static RowGrid t_vec;
        readVecAny(prhs[2], t_vec, "t_vec");
        const Vec3 HDGRef = readVec3(prhs[3], "HDGRef");
        checkConstants(prhs[4]);
        const LimitsT L = readLimits(prhs[5]);
        static Mat15Grid X;
        static WrenchT Wr;
        static Mat3Grid AccelCg;
        const PolyQPStatus st = flatnessInvert(Q, t_vec, HDGRef, L, X, Wr, AccelCg);
        const bool ok = st == PolyQPStatus::Ok;
        out.push(ok ? toMx(X) : mxCreateDoubleMatrix(0, 0, mxREAL));
        out.push(ok ? writeWrench(Wr) : emptyStruct());
        out.push(ok ? toMx(AccelCg) : mxCreateDoubleMatrix(0, 0, mxREAL));
        out.finish(st, nlhs, plhs);
    } else if (cmd == "allocateControls") {
        needArgs(nrhs, 2, "[U, status, msg] = polyqp_mex('allocateControls', Wrench, Limits)");
        static WrenchT Wr;
        readWrench(prhs[1], Wr);
        const LimitsT L = readLimits(prhs[2]);
        static Mat4Grid U;
        const PolyQPStatus st = allocateControls(Wr, L, U);
        out.push(toMx(U));
        out.finish(st, nlhs, plhs);
    } else if (cmd == "sampleOutputs") {
        needArgs(nrhs, 4, "[InitialGuess, status, msg] = polyqp_mex('sampleOutputs', QPTraj, Plan, Constants, Limits)");
        static QPTrajT Q;
        readQPTraj(prhs[1], Q);
        static PlanT P;
        readPlan(prhs[2], P);
        checkConstants(prhs[3]);
        const LimitsT L = readLimits(prhs[4]);
        static InitialGuessT G;
        const PolyQPStatus st = sampleOutputs(Q, P, L, G);
        out.push(st == PolyQPStatus::Ok ? writeInitialGuess(G) : emptyStruct());
        out.finish(st, nlhs, plhs);
    } else if (cmd == "PolyEval") {
        if (nrhs < 3 || nrhs > 4) mexErrMsgIdAndTxt("polyqp_mex:usage", "Usage: [Out, status, msg] = polyqp_mex('PolyEval', QPTraj, t[, k])");
        static QPTrajT Q;
        readQPTraj(prhs[1], Q);
        static RowGrid t;
        readVecAny(prhs[2], t, "t");
        static Mat3Grid O;
        PolyQPStatus st;
        if (nrhs == 4) {
            st = PolyEval(Q, t, O, static_cast<int>(num(prhs[3], "k")));
        } else {
            st = PolyEval(Q, t, O);
        }
        out.push(toMx(O));
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
        const Vec3 pos = readVec3(prhs[1], "pos");
        TOADWaypointOptions options;
        static char labelBuf[256];
        labelBuf[0] = '\0';
        for (int i = 2; i < nrhs; i += 2) {
            const std::string name = getString(prhs[i]);
            const mxArray *v = prhs[i + 1];
            if (name == "VelBC") options.VelBC = readVec3(v, "VelBC");
            else if (name == "AccelBC") options.AccelBC = readVec3(v, "AccelBC");
            else if (name == "JerkBC") options.JerkBC = readVec3(v, "JerkBC");
            else if (name == "HoldTime") options.HoldTime = static_cast<Real>(num(v, "HoldTime"));
            else if (name == "MaxVel") options.MaxVel = readVec3(v, "MaxVel");
            else if (name == "MaxAccel") options.MaxAccel = readVec3(v, "MaxAccel");
            else if (name == "StraightLine") options.StraightLine = num(v, "StraightLine") != 0.0;
            else if (name == "HDGRef") options.HDGRef = readVec3(v, "HDGRef");
            else if (name == "B3BC") options.B3BC = readVec3(v, "B3BC");
            else if (name == "B3DotBC") options.B3DotBC = readVec3(v, "B3DotBC");
            else if (name == "Label") {
                mxGetString(v, labelBuf, sizeof(labelBuf));
                options.Label = labelBuf;
            } else mexErrMsgIdAndTxt("polyqp_mex:usage", "Unknown TOADWaypoint option '%s'", name.c_str());
        }
        out.push(writeWaypoint(TOADWaypoint(pos, options)));
        out.finish(PolyQPStatus::Ok, nlhs, plhs);
    } else if (cmd == "ComputeJtot") {
        needArgs(nrhs, 3, "[J_tot, CGz, status, msg] = polyqp_mex('ComputeJtot', m_lox, m_ipa, constants6DoF)");
        Mat3 J_tot;
        Real CGz = 0;
        checkConstants(prhs[3]);
        ::ComputeJtot(static_cast<Real>(num(prhs[1], "m_lox")), static_cast<Real>(num(prhs[2], "m_ipa")), J_tot, CGz);
        const PolyQPStatus st = PolyQPStatus::Ok;
        out.push(toMx(J_tot));
        out.push(scal(CGz));
        out.finish(st, nlhs, plhs);
    } else if (cmd == "DCM_Quat_Conversion") {
        needArgs(nrhs, 1, "[q, status, msg] = polyqp_mex('DCM_Quat_Conversion', R)");
        Mat3 R;
        readMat(prhs[1], R, "R");
        out.push(toMx(::DCM_Quat_Conversion(R)));
        out.finish(PolyQPStatus::Ok, nlhs, plhs);
    } else if (cmd == "AbortGuidance") {
        if (nrhs < 4 || nrhs > 7)
            mexErrMsgIdAndTxt("polyqp_mex:usage", "Usage: [Time, States, Inputs, tAbort, abortCalculated, Report, abortStatus, planStatus, msg] = "
                                                  "polyqp_mex('AbortGuidance', t, X_est, constantsTOAD, state[, ABORT[, ABORT_PIN_THRUST]])");
        const Real t = static_cast<Real>(num(prhs[1], "t"));
        AbortXEstT X_est;
        readVecAny(prhs[2], X_est, "X_est");
        if (X_est.size() != 15) mexErrMsgIdAndTxt("polyqp_mex:size", "X_est must have 15 elements");
        checkConstants(prhs[3]);
        AbortGuidanceState state;
        if (nrhs >= 5 && mxIsStruct(prhs[4]) && !mxIsEmpty(prhs[4])) {
            state.abortCalculated = numField(prhs[4], "abortCalculated") != 0.0;
            state.tAbort = static_cast<Real>(numField(prhs[4], "tAbort"));
            if (mxGetField(prhs[4], 0, "abortStatus") != nullptr) state.abortStatus = static_cast<AbortStatus>(static_cast<int>(numField(prhs[4], "abortStatus")));
        }
        static AbortTimeT Time;
        static AbortStatesT States;
        static AbortInputsT Inputs;
        Time.resize(0);
        States.resize(0, 15);
        Inputs.resize(0, 4);
        static ReportT Rep;
        const bool before = state.abortCalculated;
        AbortStatus ast;
        if (nrhs == 7) {
            ast = AbortGuidance(t, X_est, state, Time, States, Inputs, Rep, static_cast<Real>(num(prhs[5], "ABORT")), num(prhs[6], "ABORT_PIN_THRUST") != 0.0);
        } else if (nrhs == 6) {
            ast = AbortGuidance(t, X_est, state, Time, States, Inputs, Rep, static_cast<Real>(num(prhs[5], "ABORT")));
        } else {
            ast = AbortGuidance(t, X_est, state, Time, States, Inputs, Rep);
        }
        const bool fired = state.abortCalculated && !before;
        out.push(fired ? toMx(Time) : mxCreateDoubleMatrix(0, 0, mxREAL));
        out.push(fired ? toMx(States) : mxCreateDoubleMatrix(0, 0, mxREAL));
        out.push(fired ? toMx(Inputs) : mxCreateDoubleMatrix(0, 0, mxREAL));
        out.push(scal(state.tAbort));
        out.push(lscal(state.abortCalculated));
        out.push(fired && state.planStatus == PolyQPStatus::Ok ? writeReport(Rep) : emptyStruct());
        out.push(scal(static_cast<double>(static_cast<int>(ast))));
        out.finish(state.planStatus, nlhs, plhs);
    } else if (cmd == "quatRot") {
        needArgs(nrhs, 1, "[CIB, status, msg] = polyqp_mex('quatRot', q)");
        VecX<4> qv;
        readVecAny(prhs[1], qv, "q");
        if (qv.size() != 4) mexErrMsgIdAndTxt("polyqp_mex:size", "q must have 4 elements");
        out.push(toMx(::quatRot(Vec4(qv))));
        out.finish(PolyQPStatus::Ok, nlhs, plhs);
    } else if (cmd == "abortFallback") {
        if (nrhs != 4) mexErrMsgIdAndTxt("polyqp_mex:usage", "Usage: [Time, States, Inputs, status, msg] = polyqp_mex('abortFallback', X_est, Constants, vCapZ)");
        AbortXEstT X_est;
        readVecAny(prhs[1], X_est, "X_est");
        if (X_est.size() != 15) mexErrMsgIdAndTxt("polyqp_mex:size", "X_est must have 15 elements");
        checkConstants(prhs[2]);
        static AbortTimeT Time;
        static AbortStatesT States;
        static AbortInputsT Inputs;
        abortFallback(X_est, static_cast<Real>(num(prhs[3], "vCapZ")), Time, States, Inputs);
        out.push(toMx(Time));
        out.push(toMx(States));
        out.push(toMx(Inputs));
        out.finish(PolyQPStatus::Ok, nlhs, plhs);
    } else {
        mexErrMsgIdAndTxt("polyqp_mex:usage", "Unknown function '%s'", cmd.c_str());
    }
}
