// Closed-Form MVP Abort Trajectory Generator (ASTRA only)
#include "AbortCode.h"
#include <cmath>

// Closed-form degree-9 Hermite block solution for zero terminal derivatives (v1 = a1 = j1 = s1 = 0)
// Derived from assembleQP.m: c_bot = A11^-1 * ([p1; 0; 0; 0; 0] - A10 * c_top)
static Matrix3_10 solveRestSegment(const Matrix3_5 &bc0, const Vector3 &p1, float T) {
    static const Matrix5_5 M_hermite = (Matrix5_5() <<
         126.0f,  -70.0f,  -35.0f, -15.0f,  -5.0f,
        -420.0f,  224.0f,  105.0f,  40.0f,  10.0f,
         540.0f, -280.0f, -126.0f, -45.0f, -10.0f,
        -315.0f,  160.0f,   70.0f,  24.0f,   5.0f,
          70.0f,  -35.0f,  -15.0f,  -5.0f,  -1.0f).finished();

    float T2 = T * T;
    float T3 = T2 * T;
    float T4 = T3 * T;

    Matrix3_10 coef = Matrix3_10::Zero();
    for (int d = 0; d < 3; ++d) {
        float c0 = bc0(d, 0);
        float c1 = T * bc0(d, 1);
        float c2 = 0.5f * T2 * bc0(d, 2);
        float c3 = (T3 / 6.0f) * bc0(d, 3);
        float c4 = (T4 / 24.0f) * bc0(d, 4);

        Vector5 rhs;
        rhs << (p1(d) - c0), c1, c2, c3, c4;
        Vector5 c_bot = M_hermite * rhs;

        coef(d, 0) = c0;
        coef(d, 1) = c1;
        coef(d, 2) = c2;
        coef(d, 3) = c3;
        coef(d, 4) = c4;
        coef(d, 5) = c_bot(0);
        coef(d, 6) = c_bot(1);
        coef(d, 7) = c_bot(2);
        coef(d, 8) = c_bot(3);
        coef(d, 9) = c_bot(4);
    }
    return coef;
}

// Evaluate derivatives thru R=4 in normalized time tau in [0, 1] via Horner's method (Outlined in RICO)
static void evalSegmentDerivs(const Matrix3_10 &coef, float T, float tau, Vector3 &pos, Vector3 &vel, Vector3 &acc, Vector3 &jerk, Vector3 &snap) {
    static constexpr float deriv_scale[5][10] = {
        {1.0f, 1.0f,  1.0f,   1.0f,   1.0f,   1.0f,   1.0f,   1.0f,    1.0f,    1.0f},
        {0.0f, 1.0f,  2.0f,   3.0f,   4.0f,   5.0f,   6.0f,   7.0f,    8.0f,    9.0f},
        {0.0f, 0.0f,  2.0f,   6.0f,  12.0f,  20.0f,  30.0f,  42.0f,   56.0f,   72.0f},
        {0.0f, 0.0f,  0.0f,   6.0f,  24.0f,  60.0f, 120.0f, 210.0f,  336.0f,  504.0f},
        {0.0f, 0.0f,  0.0f,   0.0f,  24.0f, 120.0f, 360.0f, 840.0f, 1680.0f, 3024.0f}
    };

    float tau_c = fmaxf(0.0f, fminf(1.0f, tau));
    float inv_T = 1.0f / fmaxf(1e-6f, T);

    Matrix3_5 out = Matrix3_5::Zero();
    float inv_T_pow = 1.0f;
    for (int k = 0; k < 5; ++k) {
        for (int d = 0; d < 3; ++d) {
            float val = 0.0f;
            for (int i = 9; i >= k; --i) {
                val = val * tau_c + deriv_scale[k][i] * coef(d, i);
            }
            out(d, k) = val * inv_T_pow;
        }
        inv_T_pow *= inv_T;
    }

    pos  = out.col(0);
    vel  = out.col(1);
    acc  = out.col(2);
    jerk = out.col(3);
    snap = out.col(4);
}

// Time scale s where upward acceleration causes thrust to reach actuator limit
static float solveThrustStretch(const Vector3 &a, const Vector3 &G_vec, float e_allow, float u_thi) {
    float g = G_vec(2);
    float a2 = a.dot(a);
    float aG = a.dot(G_vec);
    float c = g * g - e_allow * e_allow;
    if (c >= 0.0f) {
        return sqrtf(u_thi);
    }
    if (a2 < 1e-8f) {
        return 0.0f;
    }
    float disc = fmaxf(0.0f, aG * aG - a2 * c);
    float w = (-aG + sqrtf(disc)) / a2;
    if (w > 0.0f) {
        return 1.0f / sqrtf(w);
    }
    return sqrtf(u_thi);
}

// Time scale s where downward acceleration causes thrust to reach actuator limit
static float solveThrustFloorStretch(const Vector3 &a, const Vector3 &G_vec, float e_floor, float u_tlo) {
    float g = G_vec(2);
    float a2 = a.dot(a);
    float aG = a.dot(G_vec);
    float c = g * g - e_floor * e_floor;
    if (c <= 0.0f) {
        // Hover thrust m*g is already below T_lo: time dilation cannot raise thrust
        return 1.0f;
    }
    if (a2 < 1e-8f || aG >= 0.0f) {
        return 0.0f;
    }
    float disc = aG * aG - a2 * c;
    if (disc < 0.0f) {
        if (u_tlo > 1.0f) {
            return 1.0f;
        }
        return 0.0f;
    }
    float w = (-aG - sqrtf(disc)) / a2;
    if (w > 0.0f) {
        return 1.0f / sqrtf(w);
    }
    return 0.0f;
}

// Minimum-crackle free 3D arrest target (c9 = 0), XY clamped within ABORT_R_MAX and Z floored at ground
static Vector3 computeArrestTarget(const Matrix3_5 &bc0_P, const Vector3 &r0, float T0, float rho) {
    float T0_2 = T0 * T0;
    float T0_3 = T0_2 * T0;

    Vector3 P_free = bc0_P.col(0)
                   + 0.5f * bc0_P.col(1) * T0
                   + (3.0f / 28.0f) * bc0_P.col(2) * T0_2
                   + (1.0f / 84.0f) * bc0_P.col(3) * T0_3;

    Vector2 offset_cg = P_free.head<2>() - r0.head<2>();
    float dist = offset_cg.norm();
    if (dist > ABORT_R_MAX) {
        offset_cg *= (ABORT_R_MAX / dist);
    }
    float z_arrest = fmaxf(0.0f, P_free(2) - rho);
    return Vector3(r0(0) + offset_cg(0), r0(1) + offset_cg(1), z_arrest + rho);
}

// Eval of Segment 1 velocity, actuator, and tilt limits, returns dilation factor for time 
static float evaluateSegment1Dilation(const Matrix3_10 &coef1, float T1, float rho, const Vector3 &v_cap, float tilt_cap, const Vector4 &q_init, const constantsASTRA_t &constants) {
    Vector3 G_vec(0.0f, 0.0f, constants.g);
    float T_lo = (ABORT_MIN_THROTTLE + ABORT_THRUST_MARGIN) * constants.MaxThrust;
    float T_hi = (1.0f - ABORT_THRUST_MARGIN) * constants.MaxThrust;
    float GimbalLimit = (1.0f - ABORT_GIMBAL_MARGIN) * ABORT_MAX_GIMBAL_RAD;
    float RollLimit = (1.0f - ABORT_THRUST_MARGIN) * ABORT_MAX_ROLL;

    float max_s_req = 0.0f;
    Vector4 q_prev = q_init;

    for (int k = 0; k < ABORT_EVAL_GRID; ++k) {
        float tau = static_cast<float>(k) / static_cast<float>(ABORT_EVAL_GRID - 1);

        Vector3 P_pos, P_vel, P_acc, P_jerk, P_snap;
        evalSegmentDerivs(coef1, T1, tau, P_pos, P_vel, P_acc, P_jerk, P_snap);

        Vector15 X_k;
        AbortWrenchPoint Wrench_k;
        Vector3 ak;
        flatnessInvertPoint(P_pos, P_vel, P_acc, P_jerk, P_snap, rho, constants, &q_prev, X_k, Wrench_k, ak);
        Vector4 uk = allocateControlsPoint(Wrench_k, constants.rTB);
        q_prev = X_k.segment<4>(0);

        // Velocity caps
        float u_vx = fabsf(X_k(7)) / v_cap(0);
        float u_vy = fabsf(X_k(8)) / v_cap(1);
        float u_vz = fabsf(X_k(9)) / v_cap(2);
        float u_vel = fmaxf(u_vx, fmaxf(u_vy, u_vz));

        // Actuator caps 
        float T_k = uk(2);
        float u_thi = T_k / T_hi;
        float u_tlo = T_lo / fmaxf(1e-3f, T_k);
        float s_thi = solveThrustStretch(ak, G_vec, T_hi / constants.m, u_thi);
        float s_tlo = solveThrustFloorStretch(ak, G_vec, T_lo / constants.m, u_tlo);

        float u_gim = fmaxf(fabsf(uk(0)), fabsf(uk(1))) / GimbalLimit;
        float u_roll = fabsf(uk(3)) / RollLimit;

        // Tilt cap
        float b3_z = 1.0f - 2.0f * (X_k(1) * X_k(1) + X_k(2) * X_k(2));
        float tilt_k = acosf(fmaxf(-1.0f, fminf(1.0f, b3_z)));
        float u_tilt = tilt_k / tilt_cap;

        float s_gim = powf(u_gim, 0.25f);
        float s_roll = powf(u_roll, 0.25f);
        float s_tilt = sqrtf(u_tilt);

        float s_req = fmaxf(u_vel, fmaxf(s_thi, s_tlo));
        s_req = fmaxf(s_req, fmaxf(s_gim, fmaxf(s_roll, s_tilt)));
        max_s_req = fmaxf(max_s_req, s_req);
    }

    return max_s_req;
}

void computeAbort(float t_now, const Vector15 &X_est, float last_thrust, const constantsASTRA_t &constants, AbortTrajectory &traj) {
    // Unpack initial state and clamp thrust into limits
    Vector4 q_est = X_est.segment<4>(0).normalized();
    Vector3 r0 = X_est.segment<3>(4);
    Vector3 v0 = X_est.segment<3>(7);
    Vector3 w0 = X_est.segment<3>(10);

    float T_lo = (ABORT_MIN_THROTTLE + ABORT_THRUST_MARGIN) * constants.MaxThrust;
    float T_hi = (1.0f - ABORT_THRUST_MARGIN) * constants.MaxThrust;
    float T0 = fmaxf(T_lo, fminf(T_hi, last_thrust));

    // Build  boundary conditions from DCM, thrust, and gyro rates
    float rho = constants.J(0) / (constants.m * constants.rTB);
    Vector3 G_vec(0.0f, 0.0f, constants.g);
    Matrix3_3 C_IB = quatRot(q_est).transpose();
    Vector3 b3 = C_IB.col(2);
    Vector3 b3_dot = (C_IB * w0).cross(b3);

    Matrix3_5 bc0_P = Matrix3_5::Zero();
    bc0_P.col(0) = r0 + rho * b3;
    bc0_P.col(1) = v0 + rho * b3_dot;
    bc0_P.col(2) = (T0 / constants.m) * b3 - G_vec;
    bc0_P.col(3) = (T0 / constants.m) * b3_dot;
    bc0_P.col(4) = Vector3::Zero();

    // 3-segment abort: Arrest hover -> rest-to-rest to 1 m -> rest-to-rest to ground
    constexpr float c_v = 315.0f / 128.0f;
    constexpr float c_a = 9.3720f;

    // Time initialization (arrest only, decoupled from altitude drop)
    float a_up_lim = fmaxf(0.30f, 0.80f * (T_hi / constants.m - constants.g));
    float a_dn_lim = fmaxf(0.30f, 0.80f * (constants.g - T_lo / constants.m));
    float a_z_init = fminf(ABORT_DEF_ACCEL_Z, fminf(a_up_lim, a_dn_lim));
    Vector3 a_init(ABORT_DEF_ACCEL_XY, ABORT_DEF_ACCEL_XY, a_z_init);

    float T_dv = fmaxf(
        c_v * (fabsf(v0(0)) / a_init(0)),
        fmaxf(c_v * (fabsf(v0(1)) / a_init(1)), c_v * (fabsf(v0(2)) / a_init(2)))
    );
    float T_arr = fmaxf(0.4f, fminf(ABORT_MAX_T_ARREST, 0.5f * T_dv));

    // Solve Segment 0 (Arrest) with time dilation
    Vector3 v_cap(
        fmaxf(ABORT_MAX_VEL_XY, fabsf(v0(0))),
        fmaxf(ABORT_MAX_VEL_XY, fabsf(v0(1))),
        fmaxf(ABORT_MAX_VEL_Z, fabsf(v0(2)))
    );
    float tilt0 = acosf(fmaxf(-1.0f, fminf(1.0f, b3(2))));
    float tilt_cap = fmaxf(ABORT_MAX_TILT_RAD, tilt0);

    Vector3 P0 = Vector3::Zero();
    Matrix3_10 coef0 = Matrix3_10::Zero();

    float best_s = 1e9f;
    float best_T_arr = T_arr;
    int stall_count = 0;

    // Time dilation loop
    for (int iter = 0; iter <= ABORT_MAX_DILATE_ITERS; ++iter) {
        P0 = computeArrestTarget(bc0_P, r0, T_arr, rho);
        coef0 = solveRestSegment(bc0_P, P0, T_arr);

        float s_req = evaluateSegment1Dilation(coef0, T_arr, rho, v_cap, tilt_cap, q_est, constants);
        if (iter == 0 || s_req < best_s - 1e-4f) {
            best_s = s_req;
            best_T_arr = T_arr;
            stall_count = 0;
        } else {
            ++stall_count;
        }

        if (s_req <= 1.001f || T_arr >= ABORT_MAX_T_ARREST || iter == ABORT_MAX_DILATE_ITERS || stall_count >= 2) {
            break;
        }

        float s_apply = fminf(2.0f, 1.02f * s_req);
        float next_T_arr = fminf(ABORT_MAX_T_ARREST, T_arr * s_apply);
        if (next_T_arr <= T_arr + 1e-4f) {
            break;
        }
        T_arr = next_T_arr;
    }

    // Restore best iterate regardless of how dilation exited
    if (best_T_arr != T_arr) {
        T_arr = best_T_arr;
        P0 = computeArrestTarget(bc0_P, r0, T_arr, rho);
        coef0 = solveRestSegment(bc0_P, P0, T_arr);
    }

    // Assemble 1-3 segments: Arrest -> 1 m hover -> 0 m ground
    Matrix3_10 seg_coef[3];
    float seg_T[3];
    int num_segs = 0;

    seg_coef[num_segs] = coef0;
    seg_T[num_segs++]  = T_arr;

    Vector3 P_cur = P0;
    float z_cur = P_cur(2) - rho;

    if (z_cur > ABORT_HOVER_ALT + 1e-3f) {
        float dz = z_cur - ABORT_HOVER_ALT;
        float T1 = fmaxf(ABORT_MIN_SEG_TIME, fmaxf(c_v * (dz / ABORT_MAX_VEL_Z), sqrtf(c_a * (dz / a_z_init))));
        Matrix3_5 bc = Matrix3_5::Zero();
        bc.col(0) = P_cur;
        Vector3 P1(P_cur(0), P_cur(1), ABORT_HOVER_ALT + rho);
        seg_coef[num_segs] = solveRestSegment(bc, P1, T1);
        seg_T[num_segs++]  = T1;
        P_cur = P1;
        z_cur = ABORT_HOVER_ALT;
    }

    if (z_cur > 1e-3f) {
        float T2 = fmaxf(ABORT_MIN_SEG_TIME, fmaxf(c_v * (z_cur / ABORT_LAND_VEL_Z), sqrtf(c_a * (z_cur / a_z_init))));
        Matrix3_5 bc = Matrix3_5::Zero();
        bc.col(0) = P_cur;
        Vector3 P2(P_cur(0), P_cur(1), rho);
        seg_coef[num_segs] = solveRestSegment(bc, P2, T2);
        seg_T[num_segs++]  = T2;
    }

    // Sample trajectory on grid 
    float T_total = 0.0f;
    for (int s = 0; s < num_segs; ++s) {
        T_total += seg_T[s];
    }
    int Nc = static_cast<int>(ceilf(T_total / ABORT_DT_SAMPLE));
    if (Nc < 1) {
        Nc = 1;
    }
    if (Nc > ABORT_MAX_POINTS - 1) {
        Nc = ABORT_MAX_POINTS - 1;
    }
    int Nn = Nc + 1;
    traj.num_points = Nn;

    Vector4 q_prev = q_est;
    for (int k = 0; k < Nn; ++k) {
        float t_rel = (k == Nn - 1) ? T_total : (T_total * static_cast<float>(k) / static_cast<float>(Nc));
        traj.time[k] = t_now + t_rel;

        Vector3 P_pos, P_vel, P_acc, P_jerk, P_snap;
        if (k == Nn - 1) {
            P_pos = Vector3(P0(0), P0(1), rho);
            P_vel = Vector3::Zero();
            P_acc = Vector3::Zero();
            P_jerk = Vector3::Zero();
            P_snap = Vector3::Zero();
        } else {
            float t_seg = t_rel;
            int s = 0;
            while (s < num_segs - 1 && t_seg > seg_T[s]) {
                t_seg -= seg_T[s];
                ++s;
            }
            float tau = t_seg / seg_T[s];
            evalSegmentDerivs(seg_coef[s], seg_T[s], tau, P_pos, P_vel, P_acc, P_jerk, P_snap);
        }

        Vector15 X_k;
        AbortWrenchPoint Wrench_k;
        Vector3 AccelCg_k;
        flatnessInvertPoint(P_pos, P_vel, P_acc, P_jerk, P_snap, rho, constants, &q_prev, X_k, Wrench_k, AccelCg_k);
        Vector4 U_k = allocateControlsPoint(Wrench_k, constants.rTB);
        q_prev = X_k.segment<4>(0);

        for (int i = 0; i < 15; ++i) {
            traj.states[k][i] = X_k(i);
        }
        for (int i = 0; i < 4; ++i) {
            traj.inputs[k][i] = U_k(i);
        }
    }
}
