// Transpiled from TOAD/Guidance/PolynomialQP/sampleOutputs.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

PolyQPStatus sampleOutputs(const QPTrajT &QPTraj, const PlanT &Plan, const LimitsT &Limits, InitialGuessT &InitialGuess) {
    // Sample Trajectory on Uniform DtSample Grid
    const Real T_total = QPTraj.T_total;
    const Real dt_target = Limits.DtSample;
    // TRANSPILE-FLAG F11: the 3999 cap silently truncates the output grid for long trajectories. Transpiled as-is per MATLAB.
    const int Nc = static_cast<int>(std::min(Real(3999), std::max(Real(1), std::ceil(T_total / dt_target))));
    POLYQP_CAPACITY(Nc <= POLYQP_MAX_GRID);
    const int Nn = Nc + 1;
    static RowGrid t_vec;
    linspace(0, T_total, Nn, t_vec);
    const Real dt = T_total / Real(Nc);

    if (!Limits.HasVehicle) {
        InitialGuess.Time = t_vec;
        InitialGuess.X = Mat15Grid::Zero(15, Nn);
        InitialGuess.U = Mat4Grid::Zero(4, Nc);
        InitialGuess.T_total = T_total;
        InitialGuess.dt_row = dt * RowGrid::Ones(Nc);
        return PolyQPStatus::Ok;
    }

    static Mat15Grid X_nodes;
    static WrenchT Wrench;
    static Mat3Grid AccelCgIgnored;
    POLYQP_TRY(flatnessInvert(QPTraj, t_vec, Plan.HDGRef, Limits, X_nodes, Wrench, AccelCgIgnored));
    static Mat4Grid U_nodes;
    POLYQP_TRY(allocateControls(Wrench, Limits, U_nodes));
    static Mat4Grid U_ctrl;
    U_ctrl = U_nodes(Eigen::all, R1(1, Nc));

    InitialGuess.Time = t_vec;
    InitialGuess.X = X_nodes;
    InitialGuess.U = U_ctrl;
    InitialGuess.T_total = T_total;
    InitialGuess.dt_row = dt * RowGrid::Ones(Nc);
    return PolyQPStatus::Ok;
}

}
