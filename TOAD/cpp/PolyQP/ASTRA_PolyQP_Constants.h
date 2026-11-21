// ASTRA (vehicle 0) constants for the PolyQP abort solver.
// Copied from TOAD/LoadTOADParams.m with Vehicle = 0; the solver is ASTRA-only, so these replace the runtime
// PolyQPConstants struct (constants6DoF / constantsTOAD in the MATLAB). ASTRA carries no propellant, so the
// propellant and engine fields of LoadTOADParams (OxMass, FuMass, OxHeight, FuHeight, OxRadius, FuRadius, Ox_Z,
// Fu_Z, MaxMdot, OF) are not needed and are not carried.
//
// Compared with the flight controller (see INTEGRATION.md, "Constants"):
//   g:         9.80145 here (LoadTOADParams.m)  vs  9.8015 in ControllerAndEstimator::init_controller_and_estimator_constants
//              (constantsASTRA.g), TrajectoryManager.cpp / ASTRA_Trajectory_Data.h metadata and verify_trajectory_manager.m.
//              Not resolved here: DECISION for the lead.
//   m_dry:     1.1690 everywhere (float 1.1690f in the controller).
//   rTB:       0.26 everywhere (the controller marks it "TODO: Check if true").
//   J:         diag(0.067, 0.067, 0.02) everywhere.
//   MaxThrust: 16 everywhere.
#pragma once

namespace ASTRA {

constexpr double m_dry = 1.1690;    // kg
constexpr double g = 9.80145;       // m/s^2
constexpr double rTB = 0.26;        // m, thrust point to CG
constexpr double J_xx = 0.067;      // kg m^2
constexpr double J_yy = 0.067;      // kg m^2
constexpr double J_zz = 0.02;       // kg m^2
constexpr double MaxThrust = 16;    // N

} // namespace ASTRA
