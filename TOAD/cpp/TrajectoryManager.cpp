#include "TrajectoryManager.h"
#include <cmath>

namespace TrajectoryManager {

static const AstraGuidance::TrajectoryDataset *current_dataset = nullptr;
static int cached_idx = 0;

void init() {
    current_dataset = &AstraGuidance::getDefaultTrajectory();
    cached_idx = 0;
}

void reset() {
    cached_idx = 0;
}

void set_trajectory(AstraGuidance::TrajectoryID id) {
    current_dataset = &AstraGuidance::getTrajectory(id);
    cached_idx = 0;
}

void update_controller_constants(constantsASTRA_t &constants, const TrajectoryOutput &ref) {
    constants.K_trans  = ref.K_trans;
    constants.K_rot    = ref.K_rot;
    constants.L_Att    = ref.LA;
    constants.L_Thrust = ref.LT;
}

TrajectoryOutput step(float t, const Vector15 &X_est, float m_dry, float g) {
    if (current_dataset == nullptr) {
        init();
    }

    const AstraGuidance::TrajectoryDataset &ds = *current_dataset;
    const int N = ds.meta.num_points;
    const float *T = ds.time;

    TrajectoryOutput out;

    // Helper lambdas to unpack points using Eigen::Map 
    auto getState = [&](int idx) -> Vector15 {
        return Eigen::Map<const Vector15>(ds.states + idx * 15);
    };

    auto getInput = [&](int idx) -> Vector4 {
        return Eigen::Map<const Vector4>(ds.inputs + idx * 4);
    };

    auto getKTrans = [&](int idx) -> Matrix3_6 {
        return Eigen::Map<const Matrix3_6>(ds.k_trans + idx * 18);
    };

    auto getKRot = [&](int idx) -> Matrix3_6 {
        return Eigen::Map<const Matrix3_6>(ds.k_rot + idx * 18);
    };

    auto getLA = [&](int idx) -> Matrix6_3 {
        return Eigen::Map<const Matrix6_3>(ds.la + idx * 18);
    };

    auto getLT = [&](int idx) -> Matrix9_6 {
        return Eigen::Map<const Matrix9_6>(ds.lt + idx * 54);
    };

    // Boundary Condition Checks
    if (t <= T[0]) {
        out.X_trg   = getState(0);
        out.U_ff    = getInput(0);
        out.K_trans = getKTrans(0);
        out.K_rot   = getKRot(0);
        out.LA      = getLA(0);
        out.LT      = getLT(0);
        return out;
    }

    if (t >= T[N - 1]) {
        out.X_trg = getState(N - 1);
        
        // Command upright attitude at terminal hover
        out.X_trg(0) = 1.0f;
        out.X_trg(1) = 0.0f;
        out.X_trg(2) = 0.0f;
        out.X_trg(3) = 0.0f;

        // Zero out linear velocities and angular rates
        out.X_trg.segment<6>(7).setZero();

        // Calculate total mass (m_dry + consumed propellant masses)
        float m_prop = X_est(13) + X_est(14);
        if (m_prop < 0.0f) {
            m_prop = 0.0f;
        }
        float total_mass = (m_dry > 0.0f ? m_dry : ds.meta.m_dry) + m_prop;
        float g_acc = (g > 0.0f ? g : ds.meta.g);

        // Hover thrust command 
        out.U_ff << 0.0f, 0.0f, total_mass * g_acc, 0.0f;

        // Hold terminal gains
        out.K_trans = getKTrans(N - 1);
        out.K_rot   = getKRot(N - 1);
        out.LA      = getLA(N - 1);
        out.LT      = getLT(N - 1);
        return out;
    }

    // Interval Search
    if (cached_idx < 0 || cached_idx >= N - 1) {
        cached_idx = 0;
    }

    // Check if current cached index still contains t
    if (t < T[cached_idx] || t >= T[cached_idx + 1]) {
        if (t >= T[cached_idx + 1] && (cached_idx + 2 < N) && t < T[cached_idx + 2]) {
            cached_idx++;
        } else { //gotta love binary search
            int low = 0;
            int high = N - 1;
            while (low <= high) {
                int mid = low + (high - low) / 2;
                if (T[mid] <= t) {
                    cached_idx = mid;
                    low = mid + 1;
                } else {
                    high = mid - 1;
                }
            }
            if (cached_idx >= N - 1) {
                cached_idx = N - 2;
            }
        }
    }

    const int i_low  = cached_idx;
    const int i_high = cached_idx + 1;

    // Interpolation Weight
    const float t_low  = T[i_low];
    const float t_high = T[i_high];
    const float dt     = t_high - t_low;
    const float ratio  = (dt > 1e-7f) ? ((t - t_low) / dt) : 0.0f;

    // State & Control Interpolation
    Vector15 X_low  = getState(i_low);
    Vector15 X_high = getState(i_high);
    out.X_trg = X_low + ratio * (X_high - X_low);

    // Quaternion Normalization
    float q_norm = out.X_trg.segment<4>(0).norm();
    if (q_norm > 1e-6f) {
        out.X_trg.segment<4>(0) /= q_norm;
    }

    Vector4 U_low  = getInput(i_low);
    Vector4 U_high = getInput(i_high);
    out.U_ff = U_low + ratio * (U_high - U_low);

    // Gain Matrix Interpolation
    Matrix3_6 KT_low  = getKTrans(i_low);
    Matrix3_6 KT_high = getKTrans(i_high);
    out.K_trans = KT_low + ratio * (KT_high - KT_low);

    Matrix3_6 KR_low  = getKRot(i_low);
    Matrix3_6 KR_high = getKRot(i_high);
    out.K_rot = KR_low + ratio * (KR_high - KR_low);

    Matrix6_3 LA_low  = getLA(i_low);
    Matrix6_3 LA_high = getLA(i_high);
    out.LA = LA_low + ratio * (LA_high - LA_low);

    Matrix9_6 LT_low  = getLT(i_low);
    Matrix9_6 LT_high = getLT(i_high);
    out.LT = LT_low + ratio * (LT_high - LT_low);

    return out;
}

} // namespace TrajectoryManager
