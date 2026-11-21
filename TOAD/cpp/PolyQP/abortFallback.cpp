// Transpiled from TOAD/Guidance/PolynomialQP/abortFallback.m @ abort-fixes
#include "AbortGuidance.h"

namespace PolyQP {

void abortFallback(const Vector15d &X_est, double vCapZ, const NodeSink &Sink) {
    // FIX-F08: Closed-Form Vertical Descent Fallback (cannot fail)
    // FIX-F08: hold the lateral position, descend rest-to-rest from z0 to 0 with the degree-9 (r = 5) profile, all derivatives zero at both ends
    const double MinSegTime = double(0.2);
    const double DtSample = double(0.1);
    const double c_v = double(315) / double(128);

    const double x0 = V1(X_est, 5);
    const double y0 = V1(X_est, 6);
    const double z0 = V1(X_est, 7);
    const double m_lox = V1(X_est, 14);
    const double m_ipa = V1(X_est, 15);
    // NO-NAN: the MATLAB replacement of non-finite estimates was removed (the estimate is always real)

    // FIX-F08: duration from the peak velocity c_v * |z0| / T of the rest-to-rest profile at the vertical cap
    const double T = std::max(MinSegTime, c_v * std::abs(z0) / vCapZ);

    // FIX-F08: same output grid as sampleOutputs
    // TRANSPILE-FLAG F34: the 3999 cap is also clamped to the sink capacity (at most POLYQP_MAX_NODES) so the fallback always fits the caller's arrays: a long descent is sampled more coarsely (identical to MATLAB when it fits).
    const int capNc = std::min(POLYQP_MAX_GRID, Sink.capacity - 1);
    POLYQP_ASSERT(capNc >= 1);
    const int Nc = static_cast<int>(std::min(std::min(double(3999), double(capNc)), std::max(double(1), std::ceil(T / DtSample))));
    const int Nn = Nc + 1;

    // FIX-F08: hover feedforward, upright attitude, zero rates, propellant held at the estimate
    const double m0 = ASTRA::m_dry + m_lox + m_ipa;
    const double T_hover = m0 * ASTRA::g;

    // STREAM: States(k, :), Time(k) and Inputs(k, :) are written row by row into Sink
    Vector15d States_k = Vector15d::Zero();
    Vector4d Inputs_k = Vector4d::Zero();
    for (int k = 1; k <= Nn; ++k) {
        const double t_k = linspaceAt(0, T, Nn, k - 1);
        const double tau = t_k / T;
        const double s = double(126) * std::pow(tau, double(5)) - double(420) * std::pow(tau, double(6)) + double(540) * std::pow(tau, double(7)) -
                         double(315) * std::pow(tau, double(8)) + double(70) * std::pow(tau, double(9));
        const double sd = double(630) * std::pow(tau, double(4)) - double(2520) * std::pow(tau, double(5)) + double(3780) * std::pow(tau, double(6)) -
                          double(2520) * std::pow(tau, double(7)) + double(630) * std::pow(tau, double(8));
        V1(States_k, 1) = 1;
        V1(States_k, 5) = x0;
        V1(States_k, 6) = y0;
        V1(States_k, 7) = z0 - z0 * s;
        V1(States_k, 10) = -z0 * sd / T;
        V1(States_k, 14) = m_lox;
        V1(States_k, 15) = m_ipa;
        V1(Inputs_k, 3) = T_hover;
        Sink.put(Sink.ctx, k, t_k, States_k, Inputs_k);
    }
}

}
