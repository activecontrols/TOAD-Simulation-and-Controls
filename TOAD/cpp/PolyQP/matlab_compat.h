// PolyQP transpile support header (MATLAB-semantics helpers) for TOAD/Guidance/PolynomialQP @ 5ac6bd23
#pragma once

// TRANSPILE-FLAG F26: __FINITE_MATH_ONLY__ is predefined as 0 by GCC/Clang in IEEE builds, so the guard tests its value, not only defined(). Transpiled as-is per MATLAB.
#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__)
#error "PolyQP requires IEEE NaN/Inf semantics: do not build with -ffast-math / -ffinite-math-only"
#endif

#include "PolyQP_types.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <type_traits>

#if defined(POLYQP_ERROR_MESSAGES)
#include <cstdio>
#endif

#if defined(POLYQP_ERROR_MESSAGES)
namespace PolyQP {
struct PolyQPErrorRecord {
    char identifier[48];
    char message[256];
};
inline PolyQPErrorRecord &lastError() {
    static PolyQPErrorRecord rec{};
    return rec;
}
}
#define POLYQP_ERROR(ID, ...)                                                                                                                                  \
    do {                                                                                                                                                       \
        std::snprintf(::PolyQP::lastError().identifier, sizeof(::PolyQP::lastError().identifier), "%s", "PolyQP:" #ID);                                       \
        std::snprintf(::PolyQP::lastError().message, sizeof(::PolyQP::lastError().message), __VA_ARGS__);                                                      \
        return ::PolyQP::PolyQPStatus::ID;                                                                                                                     \
    } while (0)
#else
#define POLYQP_ERROR(ID, ...) return ::PolyQP::PolyQPStatus::ID
#endif

#define POLYQP_TRY(expr)                                                                                                                                       \
    do {                                                                                                                                                       \
        const ::PolyQP::PolyQPStatus polyqp_status_ = (expr);                                                                                                  \
        if (polyqp_status_ != ::PolyQP::PolyQPStatus::Ok) return polyqp_status_;                                                                               \
    } while (0)

#define POLYQP_CAPACITY(cond)                                                                                                                                  \
    do {                                                                                                                                                       \
        if (!(cond)) {                                                                                                                                         \
            POLYQP_ERROR(CapacityExceeded, "PolyQP capacity exceeded: %s", #cond);                                                                             \
        }                                                                                                                                                      \
    } while (0)

namespace PolyQP {

inline Real NaN() { return std::numeric_limits<Real>::quiet_NaN(); }
inline Real Inf() { return std::numeric_limits<Real>::infinity(); }
inline Real pi() { return Real(3.141592653589793238462643383279502884); }

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

template <typename Derived> inline auto mexp(const Eigen::MatrixBase<Derived> &a) {
    return a.unaryExpr([](Real x) { return std::exp(x); });
}
template <typename Derived> inline auto mlog(const Eigen::MatrixBase<Derived> &a) {
    return a.unaryExpr([](Real x) { return std::log(x); });
}
template <typename Derived> inline auto msqrt(const Eigen::MatrixBase<Derived> &a) {
    return a.unaryExpr([](Real x) { return std::sqrt(x); });
}
template <typename Derived> inline auto mabs(const Eigen::MatrixBase<Derived> &a) {
    return a.unaryExpr([](Real x) { return std::fabs(x); });
}

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

inline Real factorial(int n) {
    static constexpr Real tbl[13] = {Real(1),      Real(1),       Real(2),        Real(6),         Real(24),         Real(120),       Real(720),
                                     Real(5040),   Real(40320),   Real(362880),   Real(3628800),   Real(39916800),   Real(479001600)};
    POLYQP_ASSERT(n >= 0 && n <= 12);
    return tbl[n];
}

inline Real deg2rad(Real angleInDegrees) { return (pi() / Real(180)) * angleInDegrees; }

inline Real mmax(Real a, Real b) { return std::fmax(a, b); }
inline Real mmin(Real a, Real b) { return std::fmin(a, b); }

template <typename DA, typename DB> inline auto mmax(const Eigen::MatrixBase<DA> &a, const Eigen::MatrixBase<DB> &b) {
    return a.binaryExpr(b, [](Real x, Real y) { return std::fmax(x, y); });
}
template <typename DA, typename DB> inline auto mmin(const Eigen::MatrixBase<DA> &a, const Eigen::MatrixBase<DB> &b) {
    return a.binaryExpr(b, [](Real x, Real y) { return std::fmin(x, y); });
}
template <typename DA> inline auto mmax(const Eigen::MatrixBase<DA> &a, Real b) {
    return a.unaryExpr([b](Real x) { return std::fmax(x, b); });
}
template <typename DA> inline auto mmin(const Eigen::MatrixBase<DA> &a, Real b) {
    return a.unaryExpr([b](Real x) { return std::fmin(x, b); });
}

template <typename Derived> inline Real maxv(const Eigen::DenseBase<Derived> &v) {
    POLYQP_ASSERT(v.size() >= 1 && (v.rows() == 1 || v.cols() == 1));
    Real m = NaN();
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        const Real x = v.derived().coeff(i);
        if (!std::isnan(x) && (std::isnan(m) || x > m)) m = x;
    }
    return m;
}
template <typename Derived> inline Real minv(const Eigen::DenseBase<Derived> &v) {
    POLYQP_ASSERT(v.size() >= 1 && (v.rows() == 1 || v.cols() == 1));
    Real m = NaN();
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        const Real x = v.derived().coeff(i);
        if (!std::isnan(x) && (std::isnan(m) || x < m)) m = x;
    }
    return m;
}
inline Real maxv(std::initializer_list<Real> v) {
    POLYQP_ASSERT(v.size() >= 1);
    Real m = NaN();
    for (const Real x : v) {
        if (!std::isnan(x) && (std::isnan(m) || x > m)) m = x;
    }
    return m;
}
inline Real minv(std::initializer_list<Real> v) {
    POLYQP_ASSERT(v.size() >= 1);
    Real m = NaN();
    for (const Real x : v) {
        if (!std::isnan(x) && (std::isnan(m) || x < m)) m = x;
    }
    return m;
}

template <typename Derived, typename Out> inline void maxDim1(const Eigen::DenseBase<Derived> &X, Out &y) {
    y.resize(1, X.cols());
    for (Eigen::Index j = 0; j < X.cols(); ++j) {
        Real m = NaN();
        for (Eigen::Index i = 0; i < X.rows(); ++i) {
            const Real x = X.derived().coeff(i, j);
            if (!std::isnan(x) && (std::isnan(m) || x > m)) m = x;
        }
        y(0, j) = m;
    }
}
template <typename Derived, typename Out> inline void maxDim2(const Eigen::DenseBase<Derived> &X, Out &y) {
    y.resize(X.rows(), 1);
    for (Eigen::Index i = 0; i < X.rows(); ++i) {
        Real m = NaN();
        for (Eigen::Index j = 0; j < X.cols(); ++j) {
            const Real x = X.derived().coeff(i, j);
            if (!std::isnan(x) && (std::isnan(m) || x > m)) m = x;
        }
        y(i, 0) = m;
    }
}

template <typename Derived> inline Real sumv(const Eigen::DenseBase<Derived> &v) {
    POLYQP_ASSERT(v.rows() == 1 || v.cols() == 1);
    Real s = Real(0);
    for (Eigen::Index i = 0; i < v.size(); ++i) s = s + v.derived().coeff(i);
    return s;
}
template <typename Derived> inline Real mean(const Eigen::DenseBase<Derived> &v) { return sumv(v) / Real(v.size()); }

template <typename Derived, typename Out> inline void cumsum(const Eigen::DenseBase<Derived> &v, Out &y) {
    POLYQP_ASSERT(v.rows() == 1 || v.cols() == 1);
    y.resize(v.rows(), v.cols());
    Real s = Real(0);
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        s = s + v.derived().coeff(i);
        y(i) = s;
    }
}

template <typename Derived, typename Out> inline void diff_dim2(const Eigen::DenseBase<Derived> &U, Out &y) {
    const Eigen::Index nc = std::max<Eigen::Index>(0, U.cols() - 1);
    y.resize(U.rows(), nc);
    for (Eigen::Index j = 0; j < nc; ++j) {
        for (Eigen::Index i = 0; i < U.rows(); ++i) y(i, j) = U.derived().coeff(i, j + 1) - U.derived().coeff(i, j);
    }
}

template <typename DA, typename DB> inline Real mdot(const Eigen::MatrixBase<DA> &a, const Eigen::MatrixBase<DB> &b) {
    POLYQP_ASSERT(a.size() == b.size());
    Real s = Real(0);
    for (Eigen::Index i = 0; i < a.size(); ++i) s = s + a.derived().coeff(i) * b.derived().coeff(i);
    return s;
}
template <typename DA> inline Real mnorm(const Eigen::MatrixBase<DA> &a) { return std::sqrt(mdot(a, a)); }
inline Vec3 mcross(const Vec3 &a, const Vec3 &b) {
    const Real c1 = V1(a, 2) * V1(b, 3) - V1(a, 3) * V1(b, 2);
    const Real c2 = V1(a, 3) * V1(b, 1) - V1(a, 1) * V1(b, 3);
    const Real c3 = V1(a, 1) * V1(b, 2) - V1(a, 2) * V1(b, 1);
    return Vec3(c1, c2, c3);
}

template <typename Out> inline void linspace(Real d1, Real d2, int n, Out &y) {
    const int n1 = n - 1;
    y.resize(1, std::max(0, n));
    if (d1 == -d2 && n > 2) {
        for (int i = 0; i < n; ++i) y(i) = Real(-n1 + 2 * i) * (d2 / Real(n1));
        y(0) = d1;
        y(n - 1) = d2;
        if (n1 % 2 == 0) {
            y(n1 / 2) = Real(0);
        }
    } else {
        const Real c = (d2 - d1) * Real(n1 - 1);
        if (std::isinf(c)) {
            if (std::isinf(d2 - d1)) {
                for (int i = 0; i < n; ++i) y(i) = d1 + (d2 / Real(n1)) * Real(i) - (d1 / Real(n1)) * Real(i);
            } else {
                for (int i = 0; i < n; ++i) y(i) = d1 + Real(i) * ((d2 - d1) / Real(n1));
            }
        } else {
            for (int i = 0; i < n; ++i) y(i) = d1 + Real(i) * (d2 - d1) / Real(n1);
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

template <typename A> inline auto colon(const A &a) {
    using T = typename std::remove_const<typename std::remove_reference<decltype(a.data[0])>::type>::type;
    return Eigen::Map<const Eigen::Array<T, Eigen::Dynamic, 1>>(a.data, a.numel());
}

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

template <typename Dst, typename T> inline void fillSub3(Dst &dst, int i1, int i2, int j1, int j2, int k1, int k2, T v) {
    for (int k = k1; k <= k2; ++k) {
        for (int j = j1; j <= j2; ++j) {
            for (int i = i1; i <= i2; ++i) dst(i, j, k) = v;
        }
    }
}

template <typename Dst, typename Derived> inline void setSlab(Dst &dst, int k, int r, const Eigen::DenseBase<Derived> &slab) {
    POLYQP_ASSERT(slab.rows() == dst.n1 && slab.cols() == r);
    for (int j = 1; j <= r; ++j) {
        for (int i = 1; i <= dst.n1; ++i) dst(i, j, k) = slab.derived().coeff(i - 1, j - 1);
    }
}

template <typename Src, typename Out> inline void getSlab(const Src &src, int k, int r, Out &slab) {
    slab.resize(src.n1, r);
    for (int j = 1; j <= r; ++j) {
        for (int i = 1; i <= src.n1; ++i) slab(i - 1, j - 1) = src(i, j, k);
    }
}

template <typename Src> inline Vec3 getFiber3(const Src &src, int j, int k) {
    POLYQP_ASSERT(src.n1 == 3);
    return Vec3(src(1, j, k), src(2, j, k), src(3, j, k));
}

template <typename Src, typename Out> inline void reshapeRowSlice(const Src &src, int d, int k1, int k2, Out &v) {
    const int r = src.n2;
    const int c = std::max(0, k2 - k1 + 1);
    v.resize(r * c, 1);
    for (int k = 1; k <= c; ++k) {
        for (int j = 1; j <= r; ++j) v((j - 1) + (k - 1) * r) = src(d, j, k1 + k - 1);
    }
}

template <typename Dst, typename Derived> inline void assignRowSlice(Dst &dst, int d, const Eigen::DenseBase<Derived> &v) {
    const int r = dst.n2;
    const int K = dst.n3;
    POLYQP_ASSERT(v.size() == r * K);
    for (int k = 1; k <= K; ++k) {
        for (int j = 1; j <= r; ++j) dst(d, j, k) = v.derived().coeff((j - 1) + (k - 1) * r);
    }
}

template <typename Src> inline void isnan3(const Src &src, Bool3 &out) {
    out.resize(src.n1, src.n2, src.n3);
    for (int q = 0; q < src.numel(); ++q) out.data[q] = std::isnan(src.data[q]);
}

inline void any3(const Bool3 &B, int dim, Bool3 &out) {
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

inline void all3(const Bool3 &B, Bool3 &out) {
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

inline bool ifTrue(const Bool3 &B) {
    if (B.numel() == 0) return false;
    for (int q = 0; q < B.numel(); ++q) {
        if (!B.data[q]) return false;
    }
    return true;
}

inline bool repmat(const TOADWaypointT &A, int m, int n, TOADWaypointArray &out) {
    if (m * n > POLYQP_MAX_WP) return false;
    out.count = m * n;
    for (int q = 1; q <= m * n; ++q) out(q) = A;
    return true;
}

}
