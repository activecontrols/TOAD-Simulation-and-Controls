// Transpiled from TOAD/Guidance/PolynomialQP/abortFallback.m @ abort-fixes
#include "AbortGuidance.h"

namespace PolyQP {

void abortFallback(const AbortXEstT &X_est, Real vCapZ, AbortTimeT &Time, AbortStatesT &States, AbortInputsT &Inputs) {
    // FIX-F08: Closed-Form Vertical Descent Fallback (cannot fail)
    // FIX-F08: hold the lateral position, descend rest-to-rest from z0 to 0 with the degree-9 (r = 5) profile, all derivatives zero at both ends
    const Real MinSegTime = Real(0.2);
    const Real DtSample = Real(0.1);
    const Real c_v = Real(315) / Real(128);

    const Real x0 = V1(X_est, 5);
    const Real y0 = V1(X_est, 6);
    const Real z0 = V1(X_est, 7);
    const Real m_lox = V1(X_est, 14);
    const Real m_ipa = V1(X_est, 15);
    // NO-NAN: the MATLAB replacement of non-finite estimates was removed (the estimate is always real)

    // FIX-F08: duration from the peak velocity c_v * |z0| / T of the rest-to-rest profile at the vertical cap
    const Real T = std::max(MinSegTime, c_v * std::abs(z0) / vCapZ);

    // FIX-F08: same output grid as sampleOutputs
    // TRANSPILE-FLAG F34: the 3999 cap is also clamped to POLYQP_MAX_GRID so the fallback cannot exceed a reduced flight capacity profile (identical to MATLAB at the default capacities).
    const int Nc = static_cast<int>(std::min(std::min(Real(3999), Real(POLYQP_MAX_GRID)), std::max(Real(1), std::ceil(T / DtSample))));
    const int Nn = Nc + 1;
    static RowGrid t_vec;
    linspace(0, T, Nn, t_vec);

    // FIX-F08: hover feedforward, upright attitude, zero rates, propellant held at the estimate
    const Real m0 = ASTRA::m_dry + m_lox + m_ipa;
    const Real T_hover = m0 * ASTRA::g;

    States.setZero(Nn, 15);
    for (int k = 1; k <= Nn; ++k) {
        const Real tau = V1(t_vec, k) / T;
        const Real s = Real(126) * std::pow(tau, Real(5)) - Real(420) * std::pow(tau, Real(6)) + Real(540) * std::pow(tau, Real(7)) -
                       Real(315) * std::pow(tau, Real(8)) + Real(70) * std::pow(tau, Real(9));
        const Real sd = Real(630) * std::pow(tau, Real(4)) - Real(2520) * std::pow(tau, Real(5)) + Real(3780) * std::pow(tau, Real(6)) -
                        Real(2520) * std::pow(tau, Real(7)) + Real(630) * std::pow(tau, Real(8));
        M1(States, k, 1) = 1;
        M1(States, k, 5) = x0;
        M1(States, k, 6) = y0;
        M1(States, k, 7) = z0 - z0 * s;
        M1(States, k, 10) = -z0 * sd / T;
        M1(States, k, 14) = m_lox;
        M1(States, k, 15) = m_ipa;
    }

    Time = t_vec(R1(1, Nn)).transpose();
    Inputs.setZero(Nn, 4);
    for (int k = 1; k <= Nn; ++k) {
        M1(Inputs, k, 3) = T_hover;
    }
}

}
