// MEX gateway: lets MATLAB call ControllerAndEstimator::compute_abort unchanged.
//   [t, X, U] = compute_abort_mex(t_now, X_est(15x1), last_thrust)
//   t: Nx1   X: Nx15 (q[1:4] r[5:7] v[8:10] w[11:13] 0 0)   U: Nx4 (gimbal_th, gimbal_ph, thrust, roll)
// X_est layout = ControllerAndEstimator::last_X_ctrl (x_est(0:9), omega_body, 0, 0).
#include "mex.h"
#include "AbortCode.h"

namespace ControllerAndEstimator {
extern Vector15 last_X_ctrl;   // inputs consumed by compute_abort()
extern float last_thrust;
void init_controller_and_estimator_constants();
void compute_abort(float t_now, AbortTrajectory &traj);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 3 || mxGetNumberOfElements(prhs[1]) != 15)
        mexErrMsgIdAndTxt("compute_abort_mex:args", "Usage: [t,X,U] = compute_abort_mex(t_now, X_est(15), last_thrust)");

    static bool inited = false;          // sets constantsASTRA once (mass, inertia, limits...)
    if (!inited) { ControllerAndEstimator::init_controller_and_estimator_constants(); inited = true; }

    const double *x = mxGetPr(prhs[1]);
    for (int i = 0; i < 15; ++i) ControllerAndEstimator::last_X_ctrl(i) = static_cast<float>(x[i]);
    ControllerAndEstimator::last_thrust = static_cast<float>(mxGetScalar(prhs[2]));  // after init (init zeroes it)

    static AbortTrajectory traj;
    ControllerAndEstimator::compute_abort(static_cast<float>(mxGetScalar(prhs[0])), traj);

    const int n = traj.num_points;
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n, 15, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n, 4, mxREAL);
    double *t = mxGetPr(plhs[0]), *X = mxGetPr(plhs[1]), *U = mxGetPr(plhs[2]);
    for (int k = 0; k < n; ++k) {
        t[k] = traj.time[k];
        for (int i = 0; i < 15; ++i) X[i * n + k] = traj.states[k][i];   // column-major
        for (int i = 0; i < 4;  ++i) U[i * n + k] = traj.inputs[k][i];
    }
}
