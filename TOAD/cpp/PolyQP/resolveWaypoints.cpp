// Transpiled from TOAD/Guidance/PolynomialQP/resolveWaypoints.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

using SlabT = MatX<3, POLYQP_MAX_R>;
using SlabMaskT = Eigen::Matrix<bool, 3, Eigen::Dynamic, Eigen::ColMajor, 3, POLYQP_MAX_R>; // NO-NAN: "fixed" mask parallel to a SlabT

static PolyQPStatus wpToFix(const TOADWaypointT &wp, int r, SlabT &slab, SlabMaskT &slabFixed, bool &isRest);
static void makeCruiseKnot(const Vec3 &pos, const Vec3 &v_cruise, int r, SlabT &slab, SlabMaskT &slabFixed);
static void pinTransverse(SlabT &slab, SlabMaskT &slabFixed, const Vec3 &dp, int r);
static Real initSingleSegTime(const Vec3 &dp, const Vector3b &dpFixed, const Vec3 &vA, const Vector3b &vAFixed, const Vec3 &vB, const Vector3b &vBFixed,
                              const Vec3 &v_max, const Vec3 &a_max, int r, Real MinSegTime);

PolyQPStatus resolveWaypoints(const TOADWaypointArray &Waypoints, const LimitsT &Limits, PlanT &Plan) {
    // Validate Waypoints, Expand Holds & Cruise Profiles, and Partition Legs
    const int MAX_WP = POLYQP_MAX_WP;
    const int MAX_SEG = POLYQP_MAX_SEG;
    const int MAX_KNOTS = MAX_SEG + 1;
    const int MAX_R = POLYQP_MAX_R;
    const int r = POLYQP_R; // ASTRA-ONLY: the polynomial order is hardcoded (r = 5), so the MATLAB InvalidOrder check (r in 3..6) cannot fail and was removed

    const int nWP = Waypoints.numel();
    if (nWP < 2 || nWP > MAX_WP) {
        POLYQP_ERROR(InvalidWaypointCount, "Waypoint count must be between 2 and %d (got %d).", MAX_WP, nWP);
    }

    // Validate heading reference consistency and rest knots
    Vec3 hdg0 = Waypoints(1).HDGRef;
    if (mnorm(hdg0) < Real(1e-9)) {
        POLYQP_ERROR(InvalidHDGRef, "First waypoint HDGRef must be non-zero.");
    }
    hdg0 = hdg0 / mnorm(hdg0);

    for (int k = 1; k <= nWP; ++k) {
        const TOADWaypointT &wp = Waypoints(k);
        if (mnorm(wp.HDGRef / std::max(Real(1e-9), mnorm(wp.HDGRef)) - hdg0) > Real(1e-6)) {
            POLYQP_ERROR(InconsistentHDGRef, "Waypoint %d HDGRef must match first waypoint HDGRef.", k);
        }
        if ((wp.MaxVel.array() <= 0).any() || (wp.MaxAccel.array() <= 0).any()) {
            POLYQP_ERROR(InvalidLimits, "Waypoint %d MaxVel and MaxAccel must be strictly positive.", k);
        }
    }

    // Peak velocity factor c_v = (2r-1)! / ((r-1)!^2 * 4^(r-1)) for degree 2r-1 rest-to-rest profile
    const Real c_v = Real(315) / Real(128); // ASTRA-ONLY: the r = 5 branch; the r = 3, 4, 6 branches were removed

    // Fixed-capacity knot and segment buffers
    static Fix3 FixCg;
    FixCg.setConstant(3, MAX_R, MAX_KNOTS, Real(0));
    static Bool3 FixCgFixed; // NO-NAN: MATLAB FixCg = nan(...): every entry starts free
    FixCgFixed.setConstant(3, MAX_R, MAX_KNOTS, false);
    BoolKnot KnotIsRest = BoolKnot::Constant(MAX_KNOTS, false);
    RowSeg SegT = RowSeg::Zero(MAX_SEG);
    BoolSeg IsHold = BoolSeg::Constant(MAX_SEG, false);
    BoolSeg IsCruise = BoolSeg::Constant(MAX_SEG, false);
    Mat3Seg SegMaxVel = Mat3Seg::Constant(3, MAX_SEG, Inf());
    Mat3Seg SegMaxAccel = Mat3Seg::Constant(3, MAX_SEG, Inf());
    Mat3Knot B3Knot = Mat3Knot::Zero(3, MAX_KNOTS); // FIX-F03: per-knot B3BC (NaN: b3 = e3)
    Bool3Knot B3KnotFixed = Bool3Knot::Constant(3, MAX_KNOTS, false); // NO-NAN: B3Knot = NaN (not given)
    Mat3Knot B3DotKnot = Mat3Knot::Zero(3, MAX_KNOTS); // FIX-F03: per-knot B3DotBC (NaN: not applied)
    Bool3Knot B3DotKnotFixed = Bool3Knot::Constant(3, MAX_KNOTS, false); // NO-NAN: B3DotKnot = NaN (not given)

    // Helper to populate a knot's BCs from a waypoint
    int nKnots = 1;
    int nSegs = 0;
    SlabT slab;
    SlabMaskT slabFixed;
    bool isRest;
    POLYQP_TRY(wpToFix(Waypoints(1), r, slab, slabFixed, isRest));
    setSlab(FixCg, 1, r, slab);
    setSlab(FixCgFixed, 1, r, slabFixed);
    V1(KnotIsRest, 1) = isRest;
    Col1(B3Knot, 1) = Waypoints(1).B3BC; // FIX-F03: knot attitude from the waypoint
    Col1(B3KnotFixed, 1) = Waypoints(1).B3Fixed.array();
    Col1(B3DotKnot, 1) = Waypoints(1).B3DotBC; // FIX-F03: knot attitude from the waypoint
    Col1(B3DotKnotFixed, 1) = Waypoints(1).B3DotFixed.array();

    // If first waypoint has a positive HoldTime, insert initial hold segment
    if (!Waypoints(1).IsPassAndGo && Waypoints(1).HoldTime > 0) {
        nSegs = nSegs + 1;
        nKnots = nKnots + 1;
        POLYQP_CAPACITY(nSegs <= MAX_SEG && nKnots <= MAX_KNOTS);
        getSlab(FixCg, nKnots - 1, r, slab);
        setSlab(FixCg, nKnots, r, slab);
        getSlab(FixCgFixed, nKnots - 1, r, slabFixed);
        setSlab(FixCgFixed, nKnots, r, slabFixed);
        Col1(B3Knot, nKnots) = Col1(B3Knot, nKnots - 1); // FIX-F03: a hold keeps the knot attitude
        Col1(B3KnotFixed, nKnots) = Col1(B3KnotFixed, nKnots - 1);
        Col1(B3DotKnot, nKnots) = Col1(B3DotKnot, nKnots - 1); // FIX-F03: a hold keeps the knot attitude
        Col1(B3DotKnotFixed, nKnots) = Col1(B3DotKnotFixed, nKnots - 1);
        V1(KnotIsRest, nKnots) = true;
        V1(SegT, nSegs) = Waypoints(1).HoldTime;
        V1(IsHold, nSegs) = true;
        Col1(SegMaxVel, nSegs) = Waypoints(1).MaxVel;
        Col1(SegMaxAccel, nSegs) = Waypoints(1).MaxAccel.cwiseMin(Limits.DefaultMaxAccel);
    }

    // Expand each waypoint transition m = 1 .. nWP-1
    for (int m = 1; m <= nWP - 1; ++m) {
        const TOADWaypointT &wpB = Waypoints(m + 1);
        const Vec3 v_max = wpB.MaxVel;
        const Vec3 a_max = wpB.MaxAccel.cwiseMin(Limits.DefaultMaxAccel);

        const Vec3 pA = getFiber3(FixCg, 1, nKnots);
        const Vec3 pB = wpB.PosBC;
        const Vec3 dp = pB - pA;
        const Vector3b dpFixed = (getFiber3(FixCgFixed, 1, nKnots).array() && wpB.PosFixed.array()).matrix(); // NO-NAN: dp is NaN where either end is free

        // Check if StraightLine is requested on segment m -> m+1
        // TRANSPILE-FLAG F06: all(~isnan(dp)) is false with NaN lateral positions, so StraightLine on abort waypoint 3 is inert. Transpiled as-is per MATLAB.
        if (wpB.StraightLine && dpFixed.all() && mnorm(dp) > Real(1e-6)) {
            getSlab(FixCg, nKnots, r, slab);
            getSlab(FixCgFixed, nKnots, r, slabFixed);
            pinTransverse(slab, slabFixed, dp, r);
            setSlab(FixCg, nKnots, r, slab);
            setSlab(FixCgFixed, nKnots, r, slabFixed);
        }

        // Check if segment is long enough to saturate MaxVel and split into [Ramp-Up, Cruise, Ramp-Down]
        bool didSplit = false;
        // TRANSPILE-FLAG F06: the cruise-split guard all(~isnan(dp)) is likewise false with NaN lateral positions, so the abort never splits. Transpiled as-is per MATLAB.
        if (dpFixed.all() && mnorm(dp) > Real(1e-3)) {
            const Real T_ideal = maxv((mabs(dp).array() / v_max.array()).matrix());
            if (T_ideal > 2 * Limits.MinSegTime) {
                const Vec3 v_cruise = dp / T_ideal;
                Vec3 vA = getFiber3(FixCg, 2, nKnots);
                const Vector3b vAFixed = getFiber3(FixCgFixed, 2, nKnots);
                for (int i = 1; i <= 3; ++i) {
                    if (!V1(vAFixed, i)) V1(vA, i) = 0; // NO-NAN: isnan(vA) is !fixed
                }
                Vec3 vB = wpB.VelBC;
                for (int i = 1; i <= 3; ++i) {
                    if (!V1(wpB.VelFixed, i)) V1(vB, i) = 0; // NO-NAN: isnan(vB) is !fixed
                }

                // Bound vertical ramp acceleration by vehicle thrust band [T_lo, T_hi]
                Vec3 a_ramp = a_max;
                if (Limits.HasVehicle) {
                    const Real a_up_lim = std::max(Real(0.30), Real(0.80) * (Limits.T_hi / Limits.m_init - Limits.g));
                    const Real a_dn_lim = std::max(Real(0.30), Real(0.80) * (Limits.g - Limits.T_lo / Limits.m_dry));
                    V1(a_ramp, 3) = minv({V1(a_ramp, 3), a_up_lim, a_dn_lim});
                }

                // Check whether endpoint velocities are compatible with a zero-accel Cruise segment:
                // Only split into [Ramp-Up, Cruise, Ramp-Down] when endpoints are either rest knots
                // (HoldTime >= 0) or have an explicitly specified VelBC aligned with dp.
                // Free pass-through knots (HoldTime = NaN, VelBC = NaN) stay in a coupled multi-segment leg
                // so their velocity & curvature are solved smoothly by the QP across the corner.
                const bool hasExplicitVelA = V1(KnotIsRest, nKnots) || vAFixed.all();
                const bool hasExplicitVelB = !wpB.IsPassAndGo || wpB.VelFixed.all();
                const Vec3 u_dp = dp / mnorm(dp);
                const Real vA_trans = mnorm(vA - mdot(vA, u_dp) * u_dp);
                const Real vB_trans = mnorm(vB - mdot(vB, u_dp) * u_dp);

                if (hasExplicitVelA && hasExplicitVelB && vA_trans < Real(1e-3) && vB_trans < Real(1e-3)) {
                    const Real T_up = std::max(Limits.MinSegTime, c_v * maxv((mabs(v_cruise - vA).array() / a_ramp.array()).matrix()));
                    const Real T_dn = std::max(Limits.MinSegTime, c_v * maxv((mabs(vB - v_cruise).array() / a_ramp.array()).matrix()));

                    const Vec3 dp_up = Real(0.5) * (vA + v_cruise) * T_up;
                    const Vec3 dp_dn = Real(0.5) * (v_cruise + vB) * T_dn;
                    const Vec3 dp_rem = dp - dp_up - dp_dn;

                    // Require cruise segment to cover at least 25% of segment distance and >= 2*MinSegTime
                    const Real proj_rem = mdot(dp_rem, dp) / mdot(dp, dp);
                    const Real T_cruise = proj_rem * T_ideal;
                    if (proj_rem >= Real(0.25) && T_cruise >= 2 * Limits.MinSegTime && mnorm(dp_rem - proj_rem * dp) < Real(1e-3) * mnorm(dp)) {
                        didSplit = true;
                        const Vec3 p_in = pA + dp_up;
                        const Vec3 p_out = pB - dp_dn;

                        // 1. Ramp-Up segment (nKnots -> nKnots+1)
                        nSegs = nSegs + 1;
                        nKnots = nKnots + 1;
                        POLYQP_CAPACITY(nSegs <= MAX_SEG && nKnots <= MAX_KNOTS);
                        makeCruiseKnot(p_in, v_cruise, r, slab, slabFixed);
                        setSlab(FixCg, nKnots, r, slab);
                        setSlab(FixCgFixed, nKnots, r, slabFixed);
                        if (wpB.StraightLine) {
                            getSlab(FixCg, nKnots - 1, r, slab);
                            getSlab(FixCgFixed, nKnots - 1, r, slabFixed);
                            pinTransverse(slab, slabFixed, dp, r);
                            setSlab(FixCg, nKnots - 1, r, slab);
                            setSlab(FixCgFixed, nKnots - 1, r, slabFixed);
                        }
                        V1(KnotIsRest, nKnots) = true;
                        V1(SegT, nSegs) = T_up;
                        Col1(SegMaxVel, nSegs) = v_max;
                        Col1(SegMaxAccel, nSegs) = a_max;

                        // 2. Constant-Velocity Cruise segment (nKnots -> nKnots+1)
                        nSegs = nSegs + 1;
                        nKnots = nKnots + 1;
                        POLYQP_CAPACITY(nSegs <= MAX_SEG && nKnots <= MAX_KNOTS);
                        makeCruiseKnot(p_out, v_cruise, r, slab, slabFixed);
                        setSlab(FixCg, nKnots, r, slab);
                        setSlab(FixCgFixed, nKnots, r, slabFixed);
                        V1(KnotIsRest, nKnots) = true;
                        V1(SegT, nSegs) = T_cruise;
                        V1(IsCruise, nSegs) = true;
                        Col1(SegMaxVel, nSegs) = v_max;
                        Col1(SegMaxAccel, nSegs) = a_max;

                        // 3. Ramp-Down segment (nKnots -> nKnots+1 at wpB)
                        nSegs = nSegs + 1;
                        nKnots = nKnots + 1;
                        POLYQP_CAPACITY(nSegs <= MAX_SEG && nKnots <= MAX_KNOTS);
                        POLYQP_TRY(wpToFix(wpB, r, slab, slabFixed, isRest));
                        setSlab(FixCg, nKnots, r, slab);
                        setSlab(FixCgFixed, nKnots, r, slabFixed);
                        V1(KnotIsRest, nKnots) = isRest;
                        Col1(B3Knot, nKnots) = wpB.B3BC; // FIX-F03: knot attitude from the waypoint
                        Col1(B3KnotFixed, nKnots) = wpB.B3Fixed.array();
                        Col1(B3DotKnot, nKnots) = wpB.B3DotBC; // FIX-F03: knot attitude from the waypoint
                        Col1(B3DotKnotFixed, nKnots) = wpB.B3DotFixed.array();
                        if (wpB.StraightLine) {
                            pinTransverse(slab, slabFixed, dp, r);
                            setSlab(FixCg, nKnots, r, slab);
                            setSlab(FixCgFixed, nKnots, r, slabFixed);
                        }
                        V1(SegT, nSegs) = T_dn;
                        Col1(SegMaxVel, nSegs) = v_max;
                        Col1(SegMaxAccel, nSegs) = a_max;
                    }
                }
            }
        }

        // Standard single polynomial flight segment if not split
        if (!didSplit) {
            const Vec3 vA_prev = getFiber3(FixCg, 2, nKnots);
            const Vector3b vA_prevFixed = getFiber3(FixCgFixed, 2, nKnots);
            nSegs = nSegs + 1;
            nKnots = nKnots + 1;
            POLYQP_CAPACITY(nSegs <= MAX_SEG && nKnots <= MAX_KNOTS);
            POLYQP_TRY(wpToFix(wpB, r, slab, slabFixed, isRest));
            setSlab(FixCg, nKnots, r, slab);
            setSlab(FixCgFixed, nKnots, r, slabFixed);
            V1(KnotIsRest, nKnots) = isRest;
            Col1(B3Knot, nKnots) = wpB.B3BC; // FIX-F03: knot attitude from the waypoint
            Col1(B3KnotFixed, nKnots) = wpB.B3Fixed.array();
            Col1(B3DotKnot, nKnots) = wpB.B3DotBC; // FIX-F03: knot attitude from the waypoint
            Col1(B3DotKnotFixed, nKnots) = wpB.B3DotFixed.array();
            if (wpB.StraightLine && dpFixed.all() && mnorm(dp) > Real(1e-6)) {
                pinTransverse(slab, slabFixed, dp, r);
                setSlab(FixCg, nKnots, r, slab);
                setSlab(FixCgFixed, nKnots, r, slabFixed);
            }
            Vec3 a_seg_init = a_max;
            if (Limits.HasVehicle) {
                const Real a_up_lim = std::max(Real(0.30), Real(0.80) * (Limits.T_hi / Limits.m_init - Limits.g));
                const Real a_dn_lim = std::max(Real(0.30), Real(0.80) * (Limits.g - Limits.T_lo / Limits.m_dry));
                V1(a_seg_init, 3) = minv({V1(a_seg_init, 3), a_up_lim, a_dn_lim});
            }
            V1(SegT, nSegs) = initSingleSegTime(dp, dpFixed, vA_prev, vA_prevFixed, getFiber3(FixCg, 2, nKnots), getFiber3(FixCgFixed, 2, nKnots), v_max, a_seg_init, r,
                                                Limits.MinSegTime);
            Col1(SegMaxVel, nSegs) = v_max;
            Col1(SegMaxAccel, nSegs) = a_max;
        }

        // Append Hold segment if wpB has HoldTime > 0
        if (!wpB.IsPassAndGo && wpB.HoldTime > 0) {
            nSegs = nSegs + 1;
            nKnots = nKnots + 1;
            POLYQP_CAPACITY(nSegs <= MAX_SEG && nKnots <= MAX_KNOTS);
            getSlab(FixCg, nKnots - 1, r, slab);
            setSlab(FixCg, nKnots, r, slab);
            getSlab(FixCgFixed, nKnots - 1, r, slabFixed);
            setSlab(FixCgFixed, nKnots, r, slabFixed);
            Col1(B3Knot, nKnots) = Col1(B3Knot, nKnots - 1); // FIX-F03: a hold keeps the knot attitude
            Col1(B3KnotFixed, nKnots) = Col1(B3KnotFixed, nKnots - 1);
            Col1(B3DotKnot, nKnots) = Col1(B3DotKnot, nKnots - 1); // FIX-F03: a hold keeps the knot attitude
            Col1(B3DotKnotFixed, nKnots) = Col1(B3DotKnotFixed, nKnots - 1);
            V1(KnotIsRest, nKnots) = true;
            V1(SegT, nSegs) = wpB.HoldTime;
            V1(IsHold, nSegs) = true;
            Col1(SegMaxVel, nSegs) = v_max;
            Col1(SegMaxAccel, nSegs) = a_max;
        }
    }

    // Partition non-hold, non-cruise flight segments into decoupled legs
    // A knot decouples the QP if all 3xR derivatives are fixed
    BoolKnot KnotFullyPinned = BoolKnot::Constant(nKnots, false);
    for (int k = 1; k <= nKnots; ++k) {
        SlabMaskT blockFixed;
        getSlab(FixCgFixed, k, r, blockFixed);
        V1(KnotFullyPinned, k) = blockFixed.all(); // NO-NAN: all(~isnan(block))
    }

    IntSeg LegStartSeg = IntSeg::Zero(MAX_SEG);
    IntSeg LegEndSeg = IntSeg::Zero(MAX_SEG);
    int nLegs = 0;
    int m = 1;
    while (m <= nSegs) {
        if (V1(IsHold, m) || V1(IsCruise, m)) {
            m = m + 1;
            continue;
        }
        const int m_start = m;
        while (m < nSegs && !V1(IsHold, m + 1) && !V1(IsCruise, m + 1) && !V1(KnotFullyPinned, m + 1)) {
            m = m + 1;
        }
        nLegs = nLegs + 1;
        V1(LegStartSeg, nLegs) = m_start;
        V1(LegEndSeg, nLegs) = m;
        m = m + 1;
    }

    Plan.R = r;
    Plan.KnotCount = nKnots;
    Plan.SegCount = nSegs;
    Plan.FixCg = FixCg;
    Plan.FixCgFixed = FixCgFixed;
    Plan.KnotIsRest = KnotIsRest;
    Plan.SegT = SegT;
    Plan.IsHold = IsHold;
    Plan.IsCruise = IsCruise;
    Plan.SegMaxVel = SegMaxVel;
    Plan.SegMaxAccel = SegMaxAccel;
    Plan.LegCount = nLegs;
    Plan.LegStartSeg = LegStartSeg;
    Plan.LegEndSeg = LegEndSeg;
    Plan.HDGRef = hdg0;
    Plan.B3BC = B3Knot; // FIX-F03: per-knot attitude for the P-space BC
    Plan.B3Fixed = B3KnotFixed;
    Plan.B3DotBC = B3DotKnot; // FIX-F03: per-knot attitude rate for the P-space BC
    Plan.B3DotFixed = B3DotKnotFixed;
    return PolyQPStatus::Ok;
}

static PolyQPStatus wpToFix(const TOADWaypointT &wp, int r, SlabT &slab, SlabMaskT &slabFixed, bool &isRest) {
    // Convert TOADWaypoint to 3 x r boundary condition slab
    slab.setZero(3, r);
    slabFixed.setConstant(3, r, false); // NO-NAN: nan(3, r)
    Col1(slab, 1) = wp.PosBC;
    Col1(slabFixed, 1) = wp.PosFixed;
    Col1(slab, 2) = wp.VelBC;
    Col1(slabFixed, 2) = wp.VelFixed;
    Col1(slab, 3) = wp.AccelBC;
    Col1(slabFixed, 3) = wp.AccelFixed;
    Col1(slab, 4) = wp.JerkBC; // ASTRA-ONLY: r = 5 >= 4 always
    Col1(slabFixed, 4) = wp.JerkFixed;
    isRest = !wp.IsPassAndGo;
    if (isRest) {
        for (int k = 2; k <= r; ++k) {
            for (int d = 1; d <= 3; ++d) {
                if (!M1(slabFixed, d, k)) { // NO-NAN: isnan(slab(d, k))
                    M1(slab, d, k) = 0;
                    M1(slabFixed, d, k) = true;
                }
            }
        }
    }
    return PolyQPStatus::Ok;
}

static void makeCruiseKnot(const Vec3 &pos, const Vec3 &v_cruise, int r, SlabT &slab, SlabMaskT &slabFixed) {
    // Build fully-pinned constant-velocity transition knot
    slab.setZero(3, r);
    slabFixed.setConstant(3, r, true);
    Col1(slab, 1) = pos;
    Col1(slab, 2) = v_cruise;
}

static void pinTransverse(SlabT &slab, SlabMaskT &slabFixed, const Vec3 &dp, int r) {
    // Pin derivatives >= 1 to zero on axes transverse to straight-line displacement
    for (int d = 1; d <= 3; ++d) {
        if (std::fabs(V1(dp, d)) < Real(1e-6)) {
            for (int k = 2; k <= r; ++k) {
                if (!M1(slabFixed, d, k)) { // NO-NAN: isnan(slab(d, k))
                    M1(slab, d, k) = 0;
                    M1(slabFixed, d, k) = true;
                }
            }
        }
    }
}

static Real initSingleSegTime(const Vec3 &dp, const Vector3b &dpFixed, const Vec3 &vA, const Vector3b &vAFixed, const Vec3 &vB, const Vector3b &vBFixed,
                              const Vec3 &v_max, const Vec3 &a_max, int r, Real MinSegTime) {
    // Kinematic duration initializer for a single polynomial segment
    (void)r;
    const Real c_v = Real(315) / Real(128); // ASTRA-ONLY: the r = 5 branch; the r = 3, 4, 6 branches were removed
    const Real c_a = Real(9.3720);

    // NO-NAN: dp_abs(isnan(dp_abs)) = 0 is "use 0 where dp is not fixed"
    Vec3 dp_abs = mabs(dp);
    for (int i = 1; i <= 3; ++i) {
        if (!V1(dpFixed, i)) V1(dp_abs, i) = 0;
    }
    // NO-NAN: max(abs(vA), abs(vB)) ignores a NaN operand and dv_abs(isnan(dv_abs)) = 0: the larger of the fixed ends, 0 if neither is fixed
    Vec3 dv_abs = Vec3::Zero();
    for (int i = 1; i <= 3; ++i) {
        if (V1(vAFixed, i) && V1(vBFixed, i)) {
            V1(dv_abs, i) = std::max(std::fabs(V1(vA, i)), std::fabs(V1(vB, i)));
        } else if (V1(vAFixed, i)) {
            V1(dv_abs, i) = std::fabs(V1(vA, i));
        } else if (V1(vBFixed, i)) {
            V1(dv_abs, i) = std::fabs(V1(vB, i));
        }
    }

    const Real T_v = maxv((c_v * (dp_abs.array() / v_max.array())).matrix());
    const Real T_a = maxv(msqrt((c_a * (dp_abs.array() / a_max.array())).matrix()));
    const Real T_dv = maxv((c_v * (dv_abs.array() / a_max.array())).matrix());
    const Real T = maxv({MinSegTime, T_v, T_a, T_dv});
    return T;
}

}
