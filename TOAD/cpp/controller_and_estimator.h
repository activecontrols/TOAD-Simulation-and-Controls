#pragma once

#include "astra_structs.h"
#include "matlab_funcs.h"
#include "AbortCode/AbortCode.h"

namespace ControllerAndEstimator {
    extern Matrix15_15 Flight_P;
    extern Vector19 x_est;

    void compute_abort(float t_now, AbortTrajectory &traj);
    void init_controller_and_estimator_constants();
    Controller_Output get_controller_output(Controller_Input ci, const Vector15 &X_trg, const Vector4 &U_ff, float ideal_dT, float loop_dT, Controller_Internals *cs);
}; // namespace ControllerAndEstimator
