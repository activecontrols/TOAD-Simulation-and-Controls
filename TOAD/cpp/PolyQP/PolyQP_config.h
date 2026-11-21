// PolyQP transpile support header (capacities) for TOAD/Guidance/PolynomialQP @ 5ac6bd23
// STATIC: every PolyQP array is fully fixed-size (no Eigen::Dynamic). The capacities below are sized for the ASTRA abort
// (three waypoints: ABT START, ABT GND, ABT END); a larger build (the verification MEX polyqp_mex_wide) overrides them with -D.
#pragma once

// Polynomial order: hardcoded r = 5 (minimum-crackle, degree 9, nc = 2r = 10). Every other order was removed.
#define POLYQP_MAX_R 5
#define POLYQP_MAX_NC 10
constexpr int POLYQP_R = POLYQP_MAX_R;
constexpr int POLYQP_NC = POLYQP_MAX_NC;

// STATIC: knot and segment capacities. resolveWaypoints turns the abort's three waypoints into K = 3 knots and M = 2 segments,
// always: HoldTime = 0 adds no hold segment, the free lateral landing position makes dp non-finite in MATLAB (masked here), so
// neither StraightLine nor the cruise split can fire (F06), and the first waypoint is pass-and-go. The capacities keep one spare
// segment; anything larger returns CapacityExceeded (resolveWaypoints), and the abort falls back.
#ifndef POLYQP_MAX_WP
#define POLYQP_MAX_WP 3
#endif

#ifndef POLYQP_MAX_SEG
#define POLYQP_MAX_SEG 3
#endif

#ifndef POLYQP_MAX_KNOTS
#define POLYQP_MAX_KNOTS (POLYQP_MAX_SEG + 1)
#endif

// STATIC: free unknowns per axis in solveLeg (the size of R_PP). For the abort leg: knot 1 pins position and velocity, so
// acceleration, jerk and snap are free (3); knots 2 and 3 are rest knots with only the lateral position free (1 each). x, y: 5
// free, z: 3 free. FIX-F07's reparametrisation p3 = p2 + Delta keeps the count (Tz is square). R_PP is solved as a
// POLYQP_MAX_FREE x POLYQP_MAX_FREE system padded with an identity block (solveLeg); more free unknowns return CapacityExceeded.
#ifndef POLYQP_MAX_FREE
#define POLYQP_MAX_FREE 5
#endif

// STATIC: output grid. sampleOutputs writes Nn = Nc + 1 nodes at DtSample = 0.1 s, Nc = min(3999, ceil(T_total / DtSample)).
// POLYQP_MAX_GRID caps Nc (sized from ABORT_MAX_START_ALT_M below); POLYQP_MAX_NODES = POLYQP_MAX_GRID + 1 is the node count
// the caller's trajectory arrays must hold (Abort::MAX_POINTS). The evaluation grid of evaluateConstraints (dt = 0.05 s, up
// to 3999 points) is streamed and needs no storage.
#ifndef POLYQP_MAX_GRID
#define POLYQP_MAX_GRID 600
#endif

// TRANSPILE-FLAG F24: derived capacities. sampleOutputs evaluates Nn = Nc + 1 <= POLYQP_MAX_GRID + 1 nodes; solveLeg needs nv = r * K <= POLYQP_MAX_R * POLYQP_MAX_KNOTS. Transpiled as-is per MATLAB.
#define POLYQP_MAX_NODES (POLYQP_MAX_GRID + 1)
#define POLYQP_MAX_NV (POLYQP_MAX_R * POLYQP_MAX_KNOTS)

static_assert(POLYQP_MAX_KNOTS == POLYQP_MAX_SEG + 1, "MATLAB: MAX_KNOTS = MAX_SEG + 1");
static_assert(POLYQP_NC == 2 * POLYQP_R, "nc = 2 * r");
static_assert(POLYQP_MAX_WP >= 3 && POLYQP_MAX_SEG >= 2 && POLYQP_MAX_KNOTS >= 3, "the abort needs 3 waypoints, 2 segments, 3 knots");
static_assert(POLYQP_MAX_FREE >= 5 && POLYQP_MAX_FREE <= POLYQP_MAX_NV, "the abort leg has 5 free unknowns per lateral axis");

// ------------------------------------------------------------------------------------------------------------------------------
// Grid capacity from the abort envelope (STATIC). ABORT_MAX_START_ALT_M is the verified envelope; the lead sets the real value.
//
// Worst-case initial time guess (resolveWaypoints -> initSingleSegTime, r = 5: c_v = 315/128, c_a = 9.3720):
//   Segment 1, ABT START (z0) -> ABT GND (z = ABORT_HOVER_ALT): dp = (free, free, ABORT_HOVER_ALT - z0), so only z counts.
//     v_max = [4; 4; 5] (ABT GND MaxVel default); a_max = [2.5; 2.5; a_z] with a_z = min(3.5, a_up_lim, a_dn_lim),
//     a_up_lim = max(0.3, 0.8 (T_hi / m - g)) = 2.01 m/s^2 (T_hi = 0.9 MaxThrust), a_dn_lim = max(0.3, 0.8 (g - T_lo / m)) = 4.01
//     (T_lo = 0.35 MaxThrust), so a_z = 2.01 m/s^2.
//     T_v = c_v |dz| / 5, T_a = sqrt(c_a |dz| / a_z), T_dv = c_v max_i |v0_i| / a_max_i over the start-velocity envelope
//     (|vx|, |vy| <= 5, |vz| <= 6 m/s, the verification box). T1 = max(MinSegTime, T_v, T_a, T_dv).
//   Segment 2, ABT GND -> ABT END (z = 0): rest to rest, MaxVel z = 1: T2 = max(MinSegTime, c_v ABORT_HOVER_ALT / 1, ...) = 4.92 s.
//   At 40 m: T_v = 18.70 s, T_a = 13.30 s, T_dv = 7.34 s -> T1 = 18.70 s, T_init = 23.62 s (MATLAB resolveWaypoints: 18.703 + 4.922).
// Dilation: the outer loop multiplies a violating leg's times by s_apply = min(2, 1.02 s_leg) per iteration, for up to 40
// iterations, so it has no hard bound; the allowance is measured. With the Phase 5 start boundary condition (acceleration and
// jerk pinned, FIX-F02), 1500 envelope starts (1000 upright, 700 of them at 30-40 m; 500 tilted <= 30 deg, |w| <= 1 rad/s, at
// 30-40 m) gave max T_total / T_init = 2.25 and max T_total = 48.28 s (z0 = 39.4 m). (Before Phase 5 the measured ratio was
// 1.40 and max T_total 33.08 s; the earlier allowance of 2 x T_init = 47.25 s was exceeded by a 39.8 m verification start
// needing 48.4 s.) The allowance is 2.5 x T_init = 59.05 s (11 % over the measured ratio), so Nc_max = 591 <= POLYQP_MAX_GRID
// = 600 and the caller's arrays hold Abort::MAX_POINTS = 601 rows. A plan that still needs more output nodes than
// POLYQP_MAX_NODES returns CapacityExceeded and the abort flies the fallback; the fallback itself coarsens its grid to
// fit (F34), so it cannot fail.
// ------------------------------------------------------------------------------------------------------------------------------
constexpr double ABORT_MAX_START_ALT_M = 40.0;
constexpr double ABORT_ENVELOPE_MAX_VXY = 5.0;          // m/s, verified start-velocity envelope (lateral)
constexpr double ABORT_ENVELOPE_MAX_VZ = 6.0;           // m/s, verified start-velocity envelope (vertical)
constexpr double POLYQP_ABORT_HOVER_ALT_M = 2.0;         // ABT GND altitude (DECISION-D5); AbortGuidance.cpp asserts it matches
constexpr double POLYQP_DILATION_ALLOWANCE = 2.5;        // measured max T_total / T_init = 2.25 (Phase 5 BC), plus margin

namespace PolyQPCapacity {
constexpr double csqrt(double x) {
    double y = x > 1 ? x : 1.0;
    for (int i = 0; i < 60; ++i) y = 0.5 * (y + x / y);
    return y;
}
constexpr double cmax(double a, double b) { return a > b ? a : b; }
constexpr double cmin(double a, double b) { return a < b ? a : b; }
constexpr double c_v = 315.0 / 128.0;
constexpr double c_a = 9.3720;
constexpr double MinSegTime = 0.2;
constexpr double DtSample = 0.1;
constexpr double m = 1.1690, g = 9.80145, MaxThrust = 16; // ASTRA_PolyQP_Constants.h (static_assert in buildVehicleLimits.cpp)
constexpr double a_z = cmin(3.5, cmin(cmax(0.30, 0.80 * (0.9 * MaxThrust / m - g)), cmax(0.30, 0.80 * (g - 0.35 * MaxThrust / m))));
constexpr double dz1 = ABORT_MAX_START_ALT_M - POLYQP_ABORT_HOVER_ALT_M;
constexpr double T1 = cmax(cmax(MinSegTime, c_v * dz1 / 5.0),
                           cmax(csqrt(c_a * dz1 / a_z), c_v * cmax(ABORT_ENVELOPE_MAX_VXY / 2.5, ABORT_ENVELOPE_MAX_VZ / a_z)));
constexpr double T2 = cmax(cmax(MinSegTime, c_v * POLYQP_ABORT_HOVER_ALT_M / 1.0), csqrt(c_a * POLYQP_ABORT_HOVER_ALT_M / a_z));
constexpr double T_init = T1 + T2;
constexpr double T_total_max = POLYQP_DILATION_ALLOWANCE * T_init;
constexpr int Nc_max = static_cast<int>(T_total_max / DtSample) + 1; // >= ceil(T_total_max / DtSample)
} // namespace PolyQPCapacity

static_assert(PolyQPCapacity::Nc_max <= POLYQP_MAX_GRID, "POLYQP_MAX_GRID must cover the worst-case abort from ABORT_MAX_START_ALT_M");
static_assert(POLYQP_MAX_GRID >= 41, "MATLAB: evaluateConstraints uses at least 41 grid points");
