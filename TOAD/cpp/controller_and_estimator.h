#pragma once

#include "astra_structs.h"
#include "matlab_funcs.h"

namespace ControllerAndEstimator {
    extern Matrix15_15 Flight_P;
    extern Vector19 x_est;

    void init_controller_and_estimator_constants();
    Controller_Output get_controller_output(Controller_Input ci, const Vector15 &X_trg, const Vector4 &U_ff, float ideal_dT, float loop_dT, Controller_Internals *cs);

    // Read-only (abort harness, Abort::plan): the 15-state X_ctrl passed to ASTRAv2_Controller on the last
    // get_controller_output call, [q(4); pos(3); vel(3); bias-corrected body rates(3); 0; 0]. Zero before the first call.
    const Vector15 &get_controller_x_est();
    // Read-only (abort harness): the LESO_Position disturbance acceleration a_dist of the last call (also in co.dist_accel_*).
    const Vector3 &get_position_disturbance();
}; // namespace ControllerAndEstimator
