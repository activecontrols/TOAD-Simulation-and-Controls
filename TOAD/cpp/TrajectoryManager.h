#pragma once

#include "matlab_funcs.h"
#include "ASTRA_Trajectory_Data.h"

namespace TrajectoryManager {

    struct TrajectoryOutput {
        Vector15 X_trg;     // Target 15-state vector [q(4), pos(3), vel(3), omega(3), mass(2)]
        Vector4  U_ff;      // Feedforward 4-control vector [pitch, yaw, thrust, roll]
        Matrix3_6 K_trans;  // Translational LQR gain matrix (3x6)
        Matrix3_6 K_rot;    // Rotational LQR gain matrix (3x6)
        Matrix6_3 LA;       // Attitude LESO observer gain matrix (6x3)
        Matrix9_6 LT;       // Thrust LESO observer gain matrix (9x6)
    };

    /**
     * @brief Initialize the Trajectory Manager with the default trajectory.
     */
    void init();

    /**
     * @brief Reset runtime state (e.g., cached search index) for a fresh flight.
     */
    void reset();

    /**
     * @brief Select which compiled trajectory to follow.
     * @param id Trajectory ID enum defined in ASTRA_Trajectory_Data.h
     */
    void set_trajectory(AstraGuidance::TrajectoryID id);

    /**
     * @brief Execute one interpolation step at flight time t.
     * 
     * Matches TrajectoryManager.m 1:1.
     * Performs boundary checking, O(1) interval localization, linear interpolation,
     * quaternion normalization, and gain interpolation.
     * 
     * @param t Current flight time / MET in seconds
     * @param X_est Current estimated state vector (at least 15 elements)
     * @param m_dry Vehicle dry mass in kg
     * @param g Gravitational acceleration in m/s^2
     * @return TrajectoryOutput containing interpolated targets and TV-LQI/LESO gains
     */
    TrajectoryOutput step(float t, const Vector15 &X_est, float m_dry, float g);

    /**
     * @brief Helper to directly update constantsASTRA_t with the interpolated gains.
     * @param constants Reference to flight constants struct
     * @param ref The TrajectoryOutput generated from step()
     */
    void update_controller_constants(constantsASTRA_t &constants, const TrajectoryOutput &ref);

} // namespace TrajectoryManager
