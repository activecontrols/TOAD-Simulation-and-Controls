#pragma once
#ifndef MATLAB_FUNCS_H
#define MATLAB_FUNCS_H

#include "astra_structs.h"
#include <ArduinoEigenDense.h>

// Note - using matlab funcs based off of ASTRA simulation: 3a2af9e756f0981da411d93f59ee6b421d5fa234
// TODO - update matlab func commit

using Vector2 = Eigen::Matrix<float, 2, 1>;
using Vector3 = Eigen::Matrix<float, 3, 1>;
using Vector4 = Eigen::Matrix<float, 4, 1>;
using Vector6 = Eigen::Matrix<float, 6, 1>;
using Vector9 = Eigen::Matrix<float, 9, 1>;
using Vector12 = Eigen::Matrix<float, 12, 1>;
using Vector13 = Eigen::Matrix<float, 13, 1>;
using Vector15 = Eigen::Matrix<float, 15, 1>;
using Vector16 = Eigen::Matrix<float, 16, 1>;
using Vector18 = Eigen::Matrix<float, 18, 1>;
using Vector19 = Eigen::Matrix<float, 19, 1>;
using Matrix1_15 = Eigen::Matrix<float, 1, 15>;
using Matrix15_1 = Eigen::Matrix<float, 15, 1>;
using Matrix2_2 = Eigen::Matrix<float, 2, 2>;
using Matrix3_3 = Eigen::Matrix<float, 3, 3>;
using Matrix3_6 = Eigen::Matrix<float, 3, 6>;
using Matrix6_6 = Eigen::Matrix<float, 6, 6>;
using Matrix6_3 = Eigen::Matrix<float, 6, 3>;
using Matrix3_9 = Eigen::Matrix<float, 3, 9>;
using Matrix3_12 = Eigen::Matrix<float, 3, 12>;
using Matrix3_15 = Eigen::Matrix<float, 3, 15>;
using Matrix4_4 = Eigen::Matrix<float, 4, 4>;
using Matrix4_15 = Eigen::Matrix<float, 4, 15>;
using Matrix6_9 = Eigen::Matrix<float, 6, 9>;
using Matrix6_12 = Eigen::Matrix<float, 6, 12>;
using Matrix6_15 = Eigen::Matrix<float, 6, 15>;
using Matrix9_3 = Eigen::Matrix<float, 9, 3>;
using Matrix9_6 = Eigen::Matrix<float, 9, 6>;
using Matrix9_9 = Eigen::Matrix<float, 9, 9>;
using Matrix12_3 = Eigen::Matrix<float, 12, 3>;
using Matrix12_6 = Eigen::Matrix<float, 12, 6>;
using Matrix12_12 = Eigen::Matrix<float, 12, 12>;
using Matrix15_3 = Eigen::Matrix<float, 15, 3>;
using Matrix15_4 = Eigen::Matrix<float, 15, 4>;
using Matrix15_6 = Eigen::Matrix<float, 15, 6>;
using Matrix15_15 = Eigen::Matrix<float, 15, 15>;
using Matrix6_18 = Eigen::Matrix<float, 6, 18>;
using Matrix18_6 = Eigen::Matrix<float, 18, 6>;
using Matrix18_18 = Eigen::Matrix<float, 18, 18>;

typedef struct {
    // Kalman Filt
    Vector3 mag;
    Matrix18_18 Q;
    Matrix6_6 R;

    // Controls
    Matrix3_6 K_trans;
    Matrix3_6 K_rot;
    Matrix3_3 K_I;

    // LESOs
    Matrix9_6 L_Thrust;
    Matrix6_3 L_Att;

    // Vehicle Params
    float g;
    float m;
    float rTB;
    float MaxThrust;
    Vector3 J;

} constantsASTRA_t;

Matrix3_3 zetaCross(Vector3 zeta);
Matrix3_3 quatRot(Vector4 q);
Matrix4_4 HamiltonianProd(Vector4 q);
Vector19 GroundEstimator(Vector19 x_est, constantsASTRA_t constantsASTRA, Vector15 z, float dT, Matrix15_15 &P, bool new_mag_packet, bool new_gps_packet, Matrix3_3 gps_vel_covar,
                         Matrix3_3 gps_pos_covar);
Vector19 FlightEstimator(Vector19 x_est, constantsASTRA_t constantsASTRA, Vector15 z, float dT, Matrix9_9 &P, bool new_gps_packet, Matrix3_3 gps_vel_covar, Matrix3_3 gps_pos_covar);
Vector19 FlightEstimator2(Vector19 x_est, constantsASTRA_t constantsASTRA, Vector15 z, float dT, Matrix15_15 &P, bool new_mag_packet, bool new_gps_packet, Matrix3_3 gps_vel_covar,
                          Matrix3_3 gps_pos_covar);
Vector16 StateAUG(Vector13 XKF, Vector3 G);
Vector4 DCM_Quat_Conversion(Matrix3_3 R);
Vector9 DigitalNF(Vector9 IN, bool GND, float THRUST_PERC, float dT, Matrix9_4 &X, Matrix9_4 &Y);
Matrix18_18 Q_gen(float dt = 0.002f);
void ASTRAv2_Controller_reset();
Vector4 ASTRAv2_Controller(bool GND, const Vector15 &X_est, const Vector15 &X_trg, const Vector4 &U_ff, const constantsASTRA_t &constantsTOAD, const Vector6 &Dist_LESO, float dT,
                           Controller_Internals *cs);
Vector3 LESO_Position(bool GND, const Vector15 &X_est, const Vector4 &U_cmd, const Matrix9_6 &L_Thrust, const constantsASTRA_t &constants, float dT, Vector9 &xhat);
Vector3 LESO_Attitude(bool GND, const Vector15 &X_est, const Vector4 &U_cmd, const Matrix6_3 &L_Att, const constantsASTRA_t &constants, float dT, Vector6 &xhat);

// call like so: matrixExpPade6<Matrix9_9>
template <typename MatrixN_N> MatrixN_N matrixExpPade6(MatrixN_N A) {
    using Scalar = typename MatrixN_N::Scalar;
    const Scalar b[] = {1.0, 0.5, 0.12, 0.0183333333333, 0.00199275362319, 0.000160590438368, 0.00000939085239292};

    MatrixN_N A2 = A * A;
    MatrixN_N A4 = A2 * A2;
    MatrixN_N A6 = A4 * A2;

    MatrixN_N U = A * (b[1] * MatrixN_N::Identity(A.rows(), A.cols()) + b[3] * A2 + b[5] * A4);
    U += A6 * b[6];

    MatrixN_N V = b[0] * MatrixN_N::Identity(A.rows(), A.cols()) + b[2] * A2 + b[4] * A4 + b[6] * A6;

    MatrixN_N numer = V + U;
    MatrixN_N denom = V - U;

    return denom.inverse() * numer;
}

// ============================================================================================================
// Double-precision MATLAB-equivalent helpers used by the PolyQP abort solver (cpp/PolyQP/).
// Moved here from cpp/PolyQP/matlab_compat.h, DCM_Quat_Conversion.cpp and quatRot.cpp. Templates and inline
// functions live in this header; the non-template helpers are in matlab_helpers.cpp. Nothing above this line
// changed: the float helpers keep their signatures and behaviour.
// ============================================================================================================

#include <algorithm>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <type_traits>
#if defined(POLYQP_DEBUG_BOUNDS)
#include <cstdlib>
#endif

using Vector3d = Eigen::Matrix<double, 3, 1>;
using Vector4d = Eigen::Matrix<double, 4, 1>;
using Matrix3_3d = Eigen::Matrix<double, 3, 3>;

// Verification-only index checks (-DPOLYQP_DEBUG_BOUNDS); compiled out in flight builds.
#if defined(POLYQP_DEBUG_BOUNDS)
using BoundsHookT = void (*)(const char *file, int line);
inline BoundsHookT &boundsHook() {
    static BoundsHookT hook = nullptr;
    return hook;
}
[[noreturn]] inline void boundsFail(const char *file, int line) {
    if (boundsHook() != nullptr) boundsHook()(file, line);
    std::abort();
}
#define POLYQP_ASSERT(cond)                                                                                                                                    \
    do {                                                                                                                                                       \
        if (!(cond)) ::boundsFail(__FILE__, __LINE__);                                                                                                         \
    } while (0)
#else
#define POLYQP_ASSERT(cond)                                                                                                                                    \
    do {                                                                                                                                                       \
    } while (0)
#endif

// MATLAB error() -> status return: propagate any non-Ok status (an enum class with an Ok enumerator).
#define POLYQP_TRY(expr)                                                                                                                                       \
    do {                                                                                                                                                       \
        const auto polyqp_status_ = (expr);                                                                                                                    \
        if (polyqp_status_ != decltype(polyqp_status_)::Ok) return polyqp_status_;                                                                             \
    } while (0)

// Enabled only for double-precision Eigen arguments, so these overloads never take part in resolving a call with a
// float matrix or float expression: the float helpers above are selected exactly as before.
template <typename Derived> using EnableIfDouble = typename std::enable_if<std::is_same<typename Derived::Scalar, double>::value, int>::type;

inline double NaN() { return std::numeric_limits<double>::quiet_NaN(); }
inline double Inf() { return std::numeric_limits<double>::infinity(); }
inline double pi() { return 3.141592653589793238462643383279502884; }

// 3-D fixed-capacity array with MATLAB column-major layout and 1-based (i, j, k) access.
template <typename T, int D1, int D2, int D3> struct Array3 {
    T data[D1 * D2 * D3];
    int n1 = 0;
    int n2 = 0;
    int n3 = 0;

    static constexpr int Cap1 = D1;
    static constexpr int Cap2 = D2;
    static constexpr int Cap3 = D3;

    void resize(int a, int b, int c) {
        POLYQP_ASSERT(a >= 0 && a <= D1 && b >= 0 && b <= D2 && c >= 0 && c <= D3);
        n1 = a;
        n2 = b;
        n3 = c;
    }
    void setConstant(int a, int b, int c, T v) {
        resize(a, b, c);
        for (int q = 0; q < a * b * c; ++q) data[q] = v;
    }
    int numel() const { return n1 * n2 * n3; }
    int size(int dim) const { return dim == 1 ? n1 : (dim == 2 ? n2 : n3); }

    T &operator()(int i, int j, int k) {
        POLYQP_ASSERT(i >= 1 && i <= n1 && j >= 1 && j <= n2 && k >= 1 && k <= n3);
        return data[(i - 1) + n1 * ((j - 1) + n2 * (k - 1))];
    }
    const T &operator()(int i, int j, int k) const {
        POLYQP_ASSERT(i >= 1 && i <= n1 && j >= 1 && j <= n2 && k >= 1 && k <= n3);
        return data[(i - 1) + n1 * ((j - 1) + n2 * (k - 1))];
    }
};

// 1-based element access A(i, j) / v(i), MATLAB-style
template <typename Derived> inline typename Derived::Scalar &M1(Eigen::DenseBase<Derived> &A, int i, int j) {
    POLYQP_ASSERT(i >= 1 && i <= A.rows() && j >= 1 && j <= A.cols());
    return A.derived().coeffRef(i - 1, j - 1);
}
template <typename Derived> inline typename Derived::Scalar M1(const Eigen::DenseBase<Derived> &A, int i, int j) {
    POLYQP_ASSERT(i >= 1 && i <= A.rows() && j >= 1 && j <= A.cols());
    return A.derived().coeff(i - 1, j - 1);
}

template <typename Derived> inline typename Derived::Scalar &V1(Eigen::DenseBase<Derived> &v, int i) {
    POLYQP_ASSERT(i >= 1 && i <= v.size());
    return v.derived().coeffRef(i - 1);
}
template <typename Derived> inline typename Derived::Scalar V1(const Eigen::DenseBase<Derived> &v, int i) {
    POLYQP_ASSERT(i >= 1 && i <= v.size());
    return v.derived().coeff(i - 1);
}

// MATLAB a:b as an Eigen index sequence (empty if b < a)
inline auto R1(int a, int b) { return Eigen::seqN(Eigen::Index(a - 1), Eigen::Index(std::max(0, b - a + 1))); }

template <typename Derived> inline auto Row1(Eigen::DenseBase<Derived> &A, int i) {
    POLYQP_ASSERT(i >= 1 && i <= A.rows());
    return A.derived().row(i - 1);
}
template <typename Derived> inline auto Row1(const Eigen::DenseBase<Derived> &A, int i) {
    POLYQP_ASSERT(i >= 1 && i <= A.rows());
    return A.derived().row(i - 1);
}
template <typename Derived> inline auto Col1(Eigen::DenseBase<Derived> &A, int j) {
    POLYQP_ASSERT(j >= 1 && j <= A.cols());
    return A.derived().col(j - 1);
}
template <typename Derived> inline auto Col1(const Eigen::DenseBase<Derived> &A, int j) {
    POLYQP_ASSERT(j >= 1 && j <= A.cols());
    return A.derived().col(j - 1);
}

// Element-wise exp, log, sqrt, abs
template <typename Derived> inline auto mexp(const Eigen::MatrixBase<Derived> &a) {
    return a.unaryExpr([](double x) { return std::exp(x); });
}
template <typename Derived> inline auto mlog(const Eigen::MatrixBase<Derived> &a) {
    return a.unaryExpr([](double x) { return std::log(x); });
}
template <typename Derived> inline auto msqrt(const Eigen::MatrixBase<Derived> &a) {
    return a.unaryExpr([](double x) { return std::sqrt(x); });
}
template <typename Derived> inline auto mabs(const Eigen::MatrixBase<Derived> &a) {
    return a.unaryExpr([](double x) { return std::fabs(x); });
}

// MATLAB find(mask): 1-based linear indices of the true entries (column-major)
template <int N> struct FindResult {
    int idx[N];
    int n = 0;

    struct ZeroBased {
        const int *p;
        int n;
        Eigen::Index size() const { return n; }
        Eigen::Index operator[](Eigen::Index i) const { return p[i] - 1; }
    };

    int numel() const { return n; }
    bool isempty() const { return n == 0; }
    int operator()(int q) const {
        POLYQP_ASSERT(q >= 1 && q <= n);
        return idx[q - 1];
    }
    ZeroBased zb() const { return ZeroBased{idx, n}; }
};

template <int N, typename Derived> inline FindResult<N> find(const Eigen::DenseBase<Derived> &mask) {
    FindResult<N> out;
    int lin = 0;
    for (Eigen::Index j = 0; j < mask.cols(); ++j) {
        for (Eigen::Index i = 0; i < mask.rows(); ++i) {
            ++lin;
            if (mask.derived().coeff(i, j)) {
                POLYQP_ASSERT(out.n < N);
                out.idx[out.n] = lin;
                ++out.n;
            }
        }
    }
    return out;
}

// FIX-F07: find(idx == v), the 1-based position of the first element equal to v (0 if none)
template <int N> inline int findFirst(const FindResult<N> &idx, int v) {
    for (int q = 1; q <= idx.numel(); ++q) {
        if (idx(q) == v) return q;
    }
    return 0;
}

// MATLAB factorial(n) for 0 <= n <= 12 (exact in double)
inline double factorial(int n) {
    static constexpr double tbl[13] = {1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0, 39916800.0, 479001600.0};
    POLYQP_ASSERT(n >= 0 && n <= 12);
    return tbl[n];
}

inline double deg2rad(double angleInDegrees) { return (pi() / 180.0) * angleInDegrees; }

// MATLAB max/min: NaN-ignoring (fmax/fmin), element-wise for matrices
inline double mmax(double a, double b) { return std::fmax(a, b); }
inline double mmin(double a, double b) { return std::fmin(a, b); }

template <typename DA, typename DB> inline auto mmax(const Eigen::MatrixBase<DA> &a, const Eigen::MatrixBase<DB> &b) {
    return a.binaryExpr(b, [](double x, double y) { return std::fmax(x, y); });
}
template <typename DA, typename DB> inline auto mmin(const Eigen::MatrixBase<DA> &a, const Eigen::MatrixBase<DB> &b) {
    return a.binaryExpr(b, [](double x, double y) { return std::fmin(x, y); });
}
template <typename DA> inline auto mmax(const Eigen::MatrixBase<DA> &a, double b) {
    return a.unaryExpr([b](double x) { return std::fmax(x, b); });
}
template <typename DA> inline auto mmin(const Eigen::MatrixBase<DA> &a, double b) {
    return a.unaryExpr([b](double x) { return std::fmin(x, b); });
}

// MATLAB max(v) / min(v) of a vector (NaN-ignoring; NaN only if every element is NaN)
template <typename Derived> inline double maxv(const Eigen::DenseBase<Derived> &v) {
    POLYQP_ASSERT(v.size() >= 1 && (v.rows() == 1 || v.cols() == 1));
    double m = NaN();
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        const double x = v.derived().coeff(i);
        if (!std::isnan(x) && (std::isnan(m) || x > m)) m = x;
    }
    return m;
}
template <typename Derived> inline double minv(const Eigen::DenseBase<Derived> &v) {
    POLYQP_ASSERT(v.size() >= 1 && (v.rows() == 1 || v.cols() == 1));
    double m = NaN();
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        const double x = v.derived().coeff(i);
        if (!std::isnan(x) && (std::isnan(m) || x < m)) m = x;
    }
    return m;
}
inline double maxv(std::initializer_list<double> v) {
    POLYQP_ASSERT(v.size() >= 1);
    double m = NaN();
    for (const double x : v) {
        if (!std::isnan(x) && (std::isnan(m) || x > m)) m = x;
    }
    return m;
}
inline double minv(std::initializer_list<double> v) {
    POLYQP_ASSERT(v.size() >= 1);
    double m = NaN();
    for (const double x : v) {
        if (!std::isnan(x) && (std::isnan(m) || x < m)) m = x;
    }
    return m;
}

// MATLAB max(X, [], 1) / max(X, [], 2)
template <typename Derived, typename Out> inline void maxDim1(const Eigen::DenseBase<Derived> &X, Out &y) {
    y.resize(1, X.cols());
    for (Eigen::Index j = 0; j < X.cols(); ++j) {
        double m = NaN();
        for (Eigen::Index i = 0; i < X.rows(); ++i) {
            const double x = X.derived().coeff(i, j);
            if (!std::isnan(x) && (std::isnan(m) || x > m)) m = x;
        }
        y(0, j) = m;
    }
}
template <typename Derived, typename Out> inline void maxDim2(const Eigen::DenseBase<Derived> &X, Out &y) {
    y.resize(X.rows(), 1);
    for (Eigen::Index i = 0; i < X.rows(); ++i) {
        double m = NaN();
        for (Eigen::Index j = 0; j < X.cols(); ++j) {
            const double x = X.derived().coeff(i, j);
            if (!std::isnan(x) && (std::isnan(m) || x > m)) m = x;
        }
        y(i, 0) = m;
    }
}

// MATLAB sum(v), mean(v), cumsum(v) in left-to-right order
template <typename Derived> inline double sumv(const Eigen::DenseBase<Derived> &v) {
    POLYQP_ASSERT(v.rows() == 1 || v.cols() == 1);
    double s = 0.0;
    for (Eigen::Index i = 0; i < v.size(); ++i) s = s + v.derived().coeff(i);
    return s;
}
template <typename Derived> inline double mean(const Eigen::DenseBase<Derived> &v) { return sumv(v) / double(v.size()); }

template <typename Derived, typename Out> inline void cumsum(const Eigen::DenseBase<Derived> &v, Out &y) {
    POLYQP_ASSERT(v.rows() == 1 || v.cols() == 1);
    y.resize(v.rows(), v.cols());
    double s = 0.0;
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        s = s + v.derived().coeff(i);
        y(i) = s;
    }
}

// MATLAB diff(U, 1, 2)
template <typename Derived, typename Out> inline void diff_dim2(const Eigen::DenseBase<Derived> &U, Out &y) {
    const Eigen::Index nc = std::max<Eigen::Index>(0, U.cols() - 1);
    y.resize(U.rows(), nc);
    for (Eigen::Index j = 0; j < nc; ++j) {
        for (Eigen::Index i = 0; i < U.rows(); ++i) y(i, j) = U.derived().coeff(i, j + 1) - U.derived().coeff(i, j);
    }
}

// MATLAB dot, norm, cross (dot accumulates left to right)
template <typename DA, typename DB> inline double mdot(const Eigen::MatrixBase<DA> &a, const Eigen::MatrixBase<DB> &b) {
    POLYQP_ASSERT(a.size() == b.size());
    double s = 0.0;
    for (Eigen::Index i = 0; i < a.size(); ++i) s = s + a.derived().coeff(i) * b.derived().coeff(i);
    return s;
}
template <typename DA> inline double mnorm(const Eigen::MatrixBase<DA> &a) { return std::sqrt(mdot(a, a)); }
inline Vector3d mcross(const Vector3d &a, const Vector3d &b) {
    const double c1 = V1(a, 2) * V1(b, 3) - V1(a, 3) * V1(b, 2);
    const double c2 = V1(a, 3) * V1(b, 1) - V1(a, 1) * V1(b, 3);
    const double c3 = V1(a, 1) * V1(b, 2) - V1(a, 2) * V1(b, 1);
    return Vector3d(c1, c2, c3);
}

// MATLAB linspace(d1, d2, n), including MATLAB's symmetric and overflow branches
template <typename Out> inline void linspace(double d1, double d2, int n, Out &y) {
    const int n1 = n - 1;
    y.resize(1, std::max(0, n));
    if (d1 == -d2 && n > 2) {
        for (int i = 0; i < n; ++i) y(i) = double(-n1 + 2 * i) * (d2 / double(n1));
        y(0) = d1;
        y(n - 1) = d2;
        if (n1 % 2 == 0) {
            y(n1 / 2) = 0.0;
        }
    } else {
        const double c = (d2 - d1) * double(n1 - 1);
        if (std::isinf(c)) {
            if (std::isinf(d2 - d1)) {
                for (int i = 0; i < n; ++i) y(i) = d1 + (d2 / double(n1)) * double(i) - (d1 / double(n1)) * double(i);
            } else {
                for (int i = 0; i < n; ++i) y(i) = d1 + double(i) * ((d2 - d1) / double(n1));
            }
        } else {
            for (int i = 0; i < n; ++i) y(i) = d1 + double(i) * (d2 - d1) / double(n1);
        }
        if (n > 0) {
            if (n == 1) {
                y(0) = d2;
            } else if (d1 == d2) {
                for (int i = 0; i < n; ++i) y(i) = d1;
            } else {
                y(0) = d1;
                y(n - 1) = d2;
            }
        }
    }
}

// MATLAB A(:) of an Array3
template <typename A> inline auto colon(const A &a) {
    using T = typename std::remove_const<typename std::remove_reference<decltype(a.data[0])>::type>::type;
    return Eigen::Map<const Eigen::Array<T, Eigen::Dynamic, 1>>(a.data, a.numel());
}

// MATLAB dst = src(i1:i2, j1:j2, k1:k2)
template <typename Src, typename Dst> inline void sub3(const Src &src, int i1, int i2, int j1, int j2, int k1, int k2, Dst &dst) {
    const int a = std::max(0, i2 - i1 + 1);
    const int b = std::max(0, j2 - j1 + 1);
    const int c = std::max(0, k2 - k1 + 1);
    dst.resize(a, b, c);
    for (int k = 1; k <= c; ++k) {
        for (int j = 1; j <= b; ++j) {
            for (int i = 1; i <= a; ++i) dst(i, j, k) = src(i1 + i - 1, j1 + j - 1, k1 + k - 1);
        }
    }
}

// MATLAB dst(i1:i2, j1:j2, k1:k2) = src
template <typename Dst, typename Src> inline void setSub3(Dst &dst, int i1, int i2, int j1, int j2, int k1, int k2, const Src &src) {
    const int a = std::max(0, i2 - i1 + 1);
    const int b = std::max(0, j2 - j1 + 1);
    const int c = std::max(0, k2 - k1 + 1);
    POLYQP_ASSERT(src.n1 == a && src.n2 == b && src.n3 == c);
    for (int k = 1; k <= c; ++k) {
        for (int j = 1; j <= b; ++j) {
            for (int i = 1; i <= a; ++i) dst(i1 + i - 1, j1 + j - 1, k1 + k - 1) = src(i, j, k);
        }
    }
}

// MATLAB dst(i1:i2, j1:j2, k1:k2) = v
template <typename Dst, typename T> inline void fillSub3(Dst &dst, int i1, int i2, int j1, int j2, int k1, int k2, T v) {
    for (int k = k1; k <= k2; ++k) {
        for (int j = j1; j <= j2; ++j) {
            for (int i = i1; i <= i2; ++i) dst(i, j, k) = v;
        }
    }
}

// MATLAB dst(:, 1:r, k) = slab
template <typename Dst, typename Derived> inline void setSlab(Dst &dst, int k, int r, const Eigen::DenseBase<Derived> &slab) {
    POLYQP_ASSERT(slab.rows() == dst.n1 && slab.cols() == r);
    for (int j = 1; j <= r; ++j) {
        for (int i = 1; i <= dst.n1; ++i) dst(i, j, k) = slab.derived().coeff(i - 1, j - 1);
    }
}

// MATLAB slab = src(:, 1:r, k)
template <typename Src, typename Out> inline void getSlab(const Src &src, int k, int r, Out &slab) {
    slab.resize(src.n1, r);
    for (int j = 1; j <= r; ++j) {
        for (int i = 1; i <= src.n1; ++i) slab(i - 1, j - 1) = src(i, j, k);
    }
}

// MATLAB src(:, j, k) of a 3 x n2 x n3 array
template <typename Src> inline Vector3d getFiber3(const Src &src, int j, int k) {
    POLYQP_ASSERT(src.n1 == 3);
    return Vector3d(src(1, j, k), src(2, j, k), src(3, j, k));
}

// MATLAB v = reshape(src(d, :, k1:k2), [], 1)
template <typename Src, typename Out> inline void reshapeRowSlice(const Src &src, int d, int k1, int k2, Out &v) {
    const int r = src.n2;
    const int c = std::max(0, k2 - k1 + 1);
    v.resize(r * c, 1);
    for (int k = 1; k <= c; ++k) {
        for (int j = 1; j <= r; ++j) v((j - 1) + (k - 1) * r) = src(d, j, k1 + k - 1);
    }
}

// MATLAB dst(d, :, :) = reshape(v, 1, r, K)
template <typename Dst, typename Derived> inline void assignRowSlice(Dst &dst, int d, const Eigen::DenseBase<Derived> &v) {
    const int r = dst.n2;
    const int K = dst.n3;
    POLYQP_ASSERT(v.size() == r * K);
    for (int k = 1; k <= K; ++k) {
        for (int j = 1; j <= r; ++j) dst(d, j, k) = v.derived().coeff((j - 1) + (k - 1) * r);
    }
}

// MATLAB isnan(A), any(B, dim), all(B), and if-of-an-array for Array3
template <typename Src, typename BoolArr> inline void isnan3(const Src &src, BoolArr &out) {
    out.resize(src.n1, src.n2, src.n3);
    for (int q = 0; q < src.numel(); ++q) out.data[q] = std::isnan(src.data[q]);
}

template <typename BoolArr> inline void any3(const BoolArr &B, int dim, BoolArr &out) {
    const int a = (dim == 1) ? 1 : B.n1;
    const int b = (dim == 2) ? 1 : B.n2;
    const int c = (dim == 3) ? 1 : B.n3;
    out.setConstant(a, b, c, false);
    for (int k = 1; k <= B.n3; ++k) {
        for (int j = 1; j <= B.n2; ++j) {
            for (int i = 1; i <= B.n1; ++i) {
                if (B(i, j, k)) out(dim == 1 ? 1 : i, dim == 2 ? 1 : j, dim == 3 ? 1 : k) = true;
            }
        }
    }
}

template <typename BoolArr> inline void all3(const BoolArr &B, BoolArr &out) {
    int dim = 1;
    if (B.n1 != 1) {
        dim = 1;
    } else if (B.n2 != 1) {
        dim = 2;
    } else if (B.n3 != 1) {
        dim = 3;
    } else {
        out.setConstant(1, 1, 1, B(1, 1, 1));
        return;
    }
    const int a = (dim == 1) ? 1 : B.n1;
    const int b = (dim == 2) ? 1 : B.n2;
    const int c = (dim == 3) ? 1 : B.n3;
    out.setConstant(a, b, c, true);
    for (int k = 1; k <= B.n3; ++k) {
        for (int j = 1; j <= B.n2; ++j) {
            for (int i = 1; i <= B.n1; ++i) {
                if (!B(i, j, k)) out(dim == 1 ? 1 : i, dim == 2 ? 1 : j, dim == 3 ? 1 : k) = false;
            }
        }
    }
}

template <typename BoolArr> inline bool ifTrue(const BoolArr &B) {
    if (B.numel() == 0) return false;
    for (int q = 0; q < B.numel(); ++q) {
        if (!B.data[q]) return false;
    }
    return true;
}

// Transpiled from TOAD/Helper/quatRot.m @ 5ac6bd23 (double overload; the float quatRot above is unchanged)
template <typename Derived, EnableIfDouble<Derived> = 0> inline Matrix3_3d quatRot(const Eigen::MatrixBase<Derived> &q) {
    // TRANSPILE-FLAG F36: duplicates the float quatRot in cpp/matlab_helpers.cpp (double here); the two must stay in sync. Added for FIX-F03.
    const double q0 = V1(q, 1);
    const double q1 = V1(q, 2);
    const double q2 = V1(q, 3);
    const double q3 = V1(q, 4);
    Matrix3_3d CIB;
    CIB << 1 - 2 * (std::pow(q2, 2.0) + std::pow(q3, 2.0)), 2 * (q1 * q2 + q0 * q3), 2 * (q1 * q3 - q0 * q2),
        2 * (q1 * q2 - q0 * q3), 1 - 2 * (std::pow(q1, 2.0) + std::pow(q3, 2.0)), 2 * (q2 * q3 + q0 * q1),
        2 * (q1 * q3 + q0 * q2), 2 * (q2 * q3 - q0 * q1), 1 - 2 * (std::pow(q1, 2.0) + std::pow(q2, 2.0));
    return CIB;
}

// Transpiled from TOAD/Helper/DCM_Quat_Conversion.m @ 5ac6bd23 (double overload; the float DCM_Quat_Conversion above is unchanged)
// TRANSPILE-FLAG F15: duplicates the float DCM_Quat_Conversion in cpp/matlab_helpers.cpp; the two must stay in sync. Transpiled as-is per MATLAB.
// Convert a DCM to a quaternion robustly
template <typename Derived, EnableIfDouble<Derived> = 0> inline Vector4d DCM_Quat_Conversion(const Eigen::MatrixBase<Derived> &R) {
    // Extract trace
    const double tr = M1(R, 1, 1) + M1(R, 2, 2) + M1(R, 3, 3);

    // Pre-allocate output
    Vector4d q = Vector4d::Zero();

    // CASE 1: Standard Case (Trace > 0)
    // This works for most orientations (small angles, < 90 deg tilt)
    if (tr > 0) {
        const double S = std::sqrt(tr + 1.0) * 2; // S = 4 * qw
        V1(q, 1) = 0.25 * S;
        V1(q, 2) = (M1(R, 3, 2) - M1(R, 2, 3)) / S;
        V1(q, 3) = (M1(R, 1, 3) - M1(R, 3, 1)) / S;
        V1(q, 4) = (M1(R, 2, 1) - M1(R, 1, 2)) / S;

    }
    // CASE 2: Singularity Avoidance (Trace <= 0)
    // We must find the largest diagonal element to avoid dividing by zero.
    else {
        if ((M1(R, 1, 1) > M1(R, 2, 2)) && (M1(R, 1, 1) > M1(R, 3, 3))) {
            // Column 1 (X) is dominant
            const double S = std::sqrt(1.0 + M1(R, 1, 1) - M1(R, 2, 2) - M1(R, 3, 3)) * 2;
            V1(q, 1) = (M1(R, 3, 2) - M1(R, 2, 3)) / S;
            V1(q, 2) = 0.25 * S;
            V1(q, 3) = (M1(R, 1, 2) + M1(R, 2, 1)) / S;
            V1(q, 4) = (M1(R, 1, 3) + M1(R, 3, 1)) / S;

        } else if (M1(R, 2, 2) > M1(R, 3, 3)) {
            // Column 2 (Y) is dominant
            const double S = std::sqrt(1.0 + M1(R, 2, 2) - M1(R, 1, 1) - M1(R, 3, 3)) * 2;
            V1(q, 1) = (M1(R, 1, 3) - M1(R, 3, 1)) / S;
            V1(q, 2) = (M1(R, 1, 2) + M1(R, 2, 1)) / S;
            V1(q, 3) = 0.25 * S;
            V1(q, 4) = (M1(R, 2, 3) + M1(R, 3, 2)) / S;

        } else {
            // Column 3 (Z) is dominant
            const double S = std::sqrt(1.0 + M1(R, 3, 3) - M1(R, 1, 1) - M1(R, 2, 2)) * 2;
            V1(q, 1) = (M1(R, 2, 1) - M1(R, 1, 2)) / S;
            V1(q, 2) = (M1(R, 1, 3) + M1(R, 3, 1)) / S;
            V1(q, 3) = (M1(R, 2, 3) + M1(R, 3, 2)) / S;
            V1(q, 4) = 0.25 * S;
        }
    }

    // Enforce normalization to correct numerical drift
    q = q / mnorm(q);
    return q;
}

// Transpiled from TOAD/Helper/ComputeJtot.m @ 5ac6bd23 (defined in matlab_helpers.cpp)
namespace PolyQP {
struct PolyQPConstants;
}
void ComputeJtot(double m_lox, double m_ipa, const PolyQP::PolyQPConstants &constants6DoF, Matrix3_3d &J_tot, double &CGz);

#endif // MATLAB_FUNCS_H
