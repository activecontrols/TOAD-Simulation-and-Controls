// Transpiled from TOAD/Guidance/PolynomialQP/sampleOutputs.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

PolyQPStatus sampleOutputs(const QPTrajT &QPTraj, const PlanT &Plan, const LimitsT &Limits, const NodeSink &Sink, InitialGuessT &InitialGuess) {
    // Sample Trajectory on Uniform DtSample Grid
    // STREAM: each node is evaluated (flatnessInvertPoint, allocateControlsPoint) and handed to Sink as it is computed, so the
    // STREAM: samples go straight into the caller's arrays; nothing grid-sized is stored here.
    const double T_total = QPTraj.T_total;
    const double dt_target = Limits.DtSample;
    // TRANSPILE-FLAG F11: the 3999 cap silently truncates the output grid for long trajectories. Transpiled as-is per MATLAB.
    const int Nc = static_cast<int>(std::min(double(3999), std::max(double(1), std::ceil(T_total / dt_target))));
    const int Nn = Nc + 1;
    POLYQP_CAPACITY(Nn <= Sink.capacity); // STATIC: checked before anything is written, so a failed plan leaves the sink untouched
    const double dt = T_total / double(Nc);

    InitialGuess.Nn = Nn;
    InitialGuess.T_total = T_total;
    InitialGuess.dt = dt;
    if (!Limits.HasVehicle) {
        for (int k = 1; k <= Nn; ++k) {
            Sink.put(Sink.ctx, k, linspaceAt(0, T_total, Nn, k - 1), Vector15d::Zero(), Vector4d::Zero());
        }
        return PolyQPStatus::Ok;
    }

    // STREAM: X_nodes(:, k) and U_nodes(:, k) per node; U_ctrl = U_nodes(:, 1:Nc), and the last node carries U_ctrl(:, end)
    Vector4d qPrev = Vector4d::Zero();
    FlatPointT pt;
    Vector4d U_ctrl = Vector4d::Zero();
    for (int k = 1; k <= Nn; ++k) {
        const double t_k = linspaceAt(0, T_total, Nn, k - 1);
        POLYQP_TRY(flatnessInvertPoint(QPTraj, t_k, Plan.HDGRef, Limits, k > 1, qPrev, pt));
        if (k <= Nc) {
            U_ctrl = allocateControlsPoint(pt.F_body, pt.Tau_body, pt.CGz, Limits);
        }
        Sink.put(Sink.ctx, k, t_k, pt.X, U_ctrl);
    }
    return PolyQPStatus::Ok;
}

}
