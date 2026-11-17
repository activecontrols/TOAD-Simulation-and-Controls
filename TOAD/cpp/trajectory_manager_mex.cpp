#include "mex.h"
#include <cstring>
#include "TrajectoryManager.h"

// Commands:
//   trajectory_manager_mex('init')
//   trajectory_manager_mex('reset')
//   [X_trg, U_ff, K_trans, K_rot, LA, LT] = trajectory_manager_mex('step', t, X_est, m_dry, g)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 1) {
        mexErrMsgIdAndTxt("TrajectoryManager:Usage",
                          "Usage:\n"
                          "  trajectory_manager_mex('init')\n"
                          "  trajectory_manager_mex('reset')\n"
                          "  [X, U, K_trans, K_rot, LA, LT] = trajectory_manager_mex('step', t, X_est, m_dry, g)\n");
    }

    char cmd[32];
    if (mxGetString(prhs[0], cmd, sizeof(cmd)) != 0) {
        mexErrMsgIdAndTxt("TrajectoryManager:Command", "First argument must be a command string ('init', 'reset', 'step').");
    }

    if (strcmp(cmd, "init") == 0) {
        TrajectoryManager::init();
        return;
    }

    if (strcmp(cmd, "reset") == 0) {
        TrajectoryManager::reset();
        return;
    }

    if (strcmp(cmd, "step") == 0) {
        if (nrhs < 2) {
            mexErrMsgIdAndTxt("TrajectoryManager:StepArgs", "step command requires at least t: step(t, [X_est], [m_dry], [g])");
        }

        float t = (float)mxGetScalar(prhs[1]);

        Vector15 X_est = Vector15::Zero();
        X_est(0) = 1.0f; // upright default
        if (nrhs >= 3 && !mxIsEmpty(prhs[2])) {
            double *x_pr = mxGetPr(prhs[2]);
            int n_elems = (int)mxGetNumberOfElements(prhs[2]);
            int copy_n = (n_elems < 15) ? n_elems : 15;
            for (int i = 0; i < copy_n; i++) {
                X_est(i) = (float)x_pr[i];
            }
        }

        float m_dry = 1.1690f;
        if (nrhs >= 4 && !mxIsEmpty(prhs[3])) {
            m_dry = (float)mxGetScalar(prhs[3]);
        }

        float g = 9.8015f;
        if (nrhs >= 5 && !mxIsEmpty(prhs[4])) {
            g = (float)mxGetScalar(prhs[4]);
        }

        TrajectoryManager::TrajectoryOutput out = TrajectoryManager::step(t, X_est, m_dry, g);

        // X_trg (15x1)
        if (nlhs >= 1) {
            plhs[0] = mxCreateDoubleMatrix(15, 1, mxREAL);
            double *pX = mxGetPr(plhs[0]);
            for (int i = 0; i < 15; i++) pX[i] = (double)out.X_trg(i);
        }

        // U_ff (4x1)
        if (nlhs >= 2) {
            plhs[1] = mxCreateDoubleMatrix(4, 1, mxREAL);
            double *pU = mxGetPr(plhs[1]);
            for (int i = 0; i < 4; i++) pU[i] = (double)out.U_ff(i);
        }

        // K_trans (3x6)
        if (nlhs >= 3) {
            plhs[2] = mxCreateDoubleMatrix(3, 6, mxREAL);
            double *pKT = mxGetPr(plhs[2]);
            for (int col = 0; col < 6; col++) {
                for (int row = 0; row < 3; row++) {
                    pKT[row + col * 3] = (double)out.K_trans(row, col);
                }
            }
        }

        // K_rot (3x6)
        if (nlhs >= 4) {
            plhs[3] = mxCreateDoubleMatrix(3, 6, mxREAL);
            double *pKR = mxGetPr(plhs[3]);
            for (int col = 0; col < 6; col++) {
                for (int row = 0; row < 3; row++) {
                    pKR[row + col * 3] = (double)out.K_rot(row, col);
                }
            }
        }

        // LA (6x3)
        if (nlhs >= 5) {
            plhs[4] = mxCreateDoubleMatrix(6, 3, mxREAL);
            double *pLA = mxGetPr(plhs[4]);
            for (int col = 0; col < 3; col++) {
                for (int row = 0; row < 6; row++) {
                    pLA[row + col * 6] = (double)out.LA(row, col);
                }
            }
        }

        // LT (9x6)
        if (nlhs >= 6) {
            plhs[5] = mxCreateDoubleMatrix(9, 6, mxREAL);
            double *pLT = mxGetPr(plhs[5]);
            for (int col = 0; col < 6; col++) {
                for (int row = 0; row < 9; row++) {
                    pLT[row + col * 9] = (double)out.LT(row, col);
                }
            }
        }
        return;
    }

    mexErrMsgIdAndTxt("TrajectoryManager:UnknownCmd", "Unknown command: %s", cmd);
}
