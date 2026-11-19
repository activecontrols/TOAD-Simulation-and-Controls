// Transpiled from TOAD/Guidance/PolynomialQP/resolveWaypoints.m @ 5ac6bd23
#include "PolyQP.h"
#include "matlab_compat.h"

namespace PolyQP {

using SlabT = MatX<3, POLYQP_MAX_R>;

static PolyQPStatus wpToFix(const TOADWaypointT &wp, int r, SlabT &slab, bool &isRest);
static SlabT makeCruiseKnot(const Vec3 &pos, const Vec3 &v_cruise, int r);
static SlabT pinTransverse(SlabT slab, const Vec3 &dp, int r);
static Real initSingleSegTime(const Vec3 &dp, const Vec3 &vA, const Vec3 &vB, const Vec3 &v_max, const Vec3 &a_max, int r, Real MinSegTime);

PolyQPStatus resolveWaypoints(const TOADWaypointArray &Waypoints, int r, const LimitsT &Limits, PlanT &Plan) {
    // Validate Waypoints, Expand Holds & Cruise Profiles, and Partition Legs
    const int MAX_WP = POLYQP_MAX_WP;
    const int MAX_SEG = POLYQP_MAX_SEG;
    const int MAX_KNOTS = MAX_SEG + 1;
    const int MAX_R = POLYQP_MAX_R;

    const int nWP = Waypoints.numel();
    if (nWP < 2 || nWP > MAX_WP) {
        POLYQP_ERROR(InvalidWaypointCount, "Waypoint count must be between 2 and %d (got %d).", MAX_WP, nWP);
    }
    if (r < 3 || r > MAX_R) {
        POLYQP_ERROR(InvalidOrder, "Order r must be in 3..%d (got %d).", MAX_R, r);
    }

    // Validate heading reference consistency and rest knots
    Vec3 hdg0 = Waypoints(1).HDGRef;
    if (mnorm(hdg0) < Real(1e-9)) {
        POLYQP_ERROR(InvalidHDGRef, "First waypoint HDGRef must be non-zero.");
    }
    hdg0 = hdg0 / mnorm(hdg0);

    for (int k = 1; k <= nWP; ++k) {
        const TOADWaypointT &wp = Waypoints(k);
        if (mnorm(wp.HDGRef / mmax(Real(1e-9), mnorm(wp.HDGRef)) - hdg0) > Real(1e-6)) {
            POLYQP_ERROR(InconsistentHDGRef, "Waypoint %d HDGRef must match first waypoint HDGRef.", k);
        }
        if ((wp.MaxVel.array() <= 0).any() || (wp.MaxAccel.array() <= 0).any()) {
            POLYQP_ERROR(InvalidLimits, "Waypoint %d MaxVel and MaxAccel must be strictly positive.", k);
        }
    }

    // Peak velocity factor c_v = (2r-1)! / ((r-1)!^2 * 4^(r-1)) for degree 2r-1 rest-to-rest profile
    Real c_v;
    if (r == 3) {
        c_v = Real(15) / Real(8);
    } else if (r == 4) {
        c_v = Real(35) / Real(16);
    } else if (r == 5) {
        c_v = Real(315) / Real(128);
    } else {
        c_v = Real(693) / Real(256);
    }

    // Fixed-capacity knot and segment buffers
    static Fix3 FixCg;
    FixCg.setConstant(3, MAX_R, MAX_KNOTS, NaN());
    BoolKnot KnotIsRest = BoolKnot::Constant(MAX_KNOTS, false);
    RowSeg SegT = RowSeg::Zero(MAX_SEG);
    BoolSeg IsHold = BoolSeg::Constant(MAX_SEG, false);
    BoolSeg IsCruise = BoolSeg::Constant(MAX_SEG, false);
    Mat3Seg SegMaxVel = Mat3Seg::Constant(3, MAX_SEG, Inf());
    Mat3Seg SegMaxAccel = Mat3Seg::Constant(3, MAX_SEG, Inf());
    Mat3Knot B3Knot = Mat3Knot::Constant(3, MAX_KNOTS, NaN()); // FIX-F03: per-knot B3BC (NaN: b3 = e3)
    Mat3Knot B3DotKnot = Mat3Knot::Constant(3, MAX_KNOTS, NaN()); // FIX-F03: per-knot B3DotBC (NaN: not applied)

    // Helper to populate a knot's BCs from a waypoint
    int nKnots = 1;
    int nSegs = 0;
    SlabT slab;
    bool isRest;
    POLYQP_TRY(wpToFix(Waypoints(1), r, slab, isRest));
    setSlab(FixCg, 1, r, slab);
    V1(KnotIsRest, 1) = isRest;
    Col1(B3Knot, 1) = Waypoints(1).B3BC; // FIX-F03: knot attitude from the waypoint
    Col1(B3DotKnot, 1) = Waypoints(1).B3DotBC; // FIX-F03: knot attitude from the waypoint

    // If first waypoint has a positive HoldTime, insert initial hold segment
    if (!Waypoints(1).IsPassAndGo && Waypoints(1).HoldTime > 0) {
        nSegs = nSegs + 1;
        nKnots = nKnots + 1;
        POLYQP_CAPACITY(nSegs <= MAX_SEG && nKnots <= MAX_KNOTS);
        getSlab(FixCg, nKnots - 1, r, slab);
        setSlab(FixCg, nKnots, r, slab);
        Col1(B3Knot, nKnots) = Col1(B3Knot, nKnots - 1); // FIX-F03: a hold keeps the knot attitude
        Col1(B3DotKnot, nKnots) = Col1(B3DotKnot, nKnots - 1); // FIX-F03: a hold keeps the knot attitude
        V1(KnotIsRest, nKnots) = true;
        V1(SegT, nSegs) = Waypoints(1).HoldTime;
        V1(IsHold, nSegs) = true;
        Col1(SegMaxVel, nSegs) = Waypoints(1).MaxVel;
        Col1(SegMaxAccel, nSegs) = mmin(Waypoints(1).MaxAccel, Limits.DefaultMaxAccel);
    }

    // Expand each waypoint transition m = 1 .. nWP-1
    for (int m = 1; m <= nWP - 1; ++m) {
        const TOADWaypointT &wpB = Waypoints(m + 1);
        const Vec3 v_max = wpB.MaxVel;
        const Vec3 a_max = mmin(wpB.MaxAccel, Limits.DefaultMaxAccel);

        const Vec3 pA = getFiber3(FixCg, 1, nKnots);
        const Vec3 pB = wpB.PosBC;
        const Vec3 dp = pB - pA;

        // Check if StraightLine is requested on segment m -> m+1
        // TRANSPILE-FLAG F06: all(~isnan(dp)) is false with NaN lateral positions, so StraightLine on abort waypoint 3 is inert. Transpiled as-is per MATLAB.
        if (wpB.StraightLine && (!dp.array().isNaN()).all() && mnorm(dp) > Real(1e-6)) {
            getSlab(FixCg, nKnots, r, slab);
            setSlab(FixCg, nKnots, r, pinTransverse(slab, dp, r));
        }

        // Check if segment is long enough to saturate MaxVel and split into [Ramp-Up, Cruise, Ramp-Down]
        bool didSplit = false;
        // TRANSPILE-FLAG F06: the cruise-split guard all(~isnan(dp)) is likewise false with NaN lateral positions, so the abort never splits. Transpiled as-is per MATLAB.
        if ((!dp.array().isNaN()).all() && mnorm(dp) > Real(1e-3)) {
            const Real T_ideal = maxv((mabs(dp).array() / v_max.array()).matrix());
            if (T_ideal > 2 * Limits.MinSegTime) {
                const Vec3 v_cruise = dp / T_ideal;
                Vec3 vA = getFiber3(FixCg, 2, nKnots);
                for (int i = 1; i <= 3; ++i) {
                    if (std::isnan(V1(vA, i))) V1(vA, i) = 0;
                }
                Vec3 vB = wpB.VelBC;
                for (int i = 1; i <= 3; ++i) {
                    if (std::isnan(V1(vB, i))) V1(vB, i) = 0;
                }

                // Bound vertical ramp acceleration by vehicle thrust band [T_lo, T_hi]
                Vec3 a_ramp = a_max;
                if (Limits.HasVehicle) {
                    const Real a_up_lim = mmax(Real(0.30), Real(0.80) * (Limits.T_hi / Limits.m_init - Limits.g));
                    const Real a_dn_lim = mmax(Real(0.30), Real(0.80) * (Limits.g - Limits.T_lo / Limits.m_dry));
                    V1(a_ramp, 3) = minv({V1(a_ramp, 3), a_up_lim, a_dn_lim});
                }

                // Check whether endpoint velocities are compatible with a zero-accel Cruise segment:
                // Only split into [Ramp-Up, Cruise, Ramp-Down] when endpoints are either rest knots
                // (HoldTime >= 0) or have an explicitly specified VelBC aligned with dp.
                // Free pass-through knots (HoldTime = NaN, VelBC = NaN) stay in a coupled multi-segment leg
                // so their velocity & curvature are solved smoothly by the QP across the corner.
                const bool hasExplicitVelA = V1(KnotIsRest, nKnots) || (!getFiber3(FixCg, 2, nKnots).array().isNaN()).all();
                const bool hasExplicitVelB = !wpB.IsPassAndGo || (!wpB.VelBC.array().isNaN()).all();
                const Vec3 u_dp = dp / mnorm(dp);
                const Real vA_trans = mnorm(vA - mdot(vA, u_dp) * u_dp);
                const Real vB_trans = mnorm(vB - mdot(vB, u_dp) * u_dp);

                if (hasExplicitVelA && hasExplicitVelB && vA_trans < Real(1e-3) && vB_trans < Real(1e-3)) {
                    const Real T_up = mmax(Limits.MinSegTime, c_v * maxv((mabs(v_cruise - vA).array() / a_ramp.array()).matrix()));
                    const Real T_dn = mmax(Limits.MinSegTime, c_v * maxv((mabs(vB - v_cruise).array() / a_ramp.array()).matrix()));

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
                        setSlab(FixCg, nKnots, r, makeCruiseKnot(p_in, v_cruise, r));
                        if (wpB.StraightLine) {
                            getSlab(FixCg, nKnots - 1, r, slab);
                            setSlab(FixCg, nKnots - 1, r, pinTransverse(slab, dp, r));
                        }
                        V1(KnotIsRest, nKnots) = true;
                        V1(SegT, nSegs) = T_up;
                        Col1(SegMaxVel, nSegs) = v_max;
                        Col1(SegMaxAccel, nSegs) = a_max;

                        // 2. Constant-Velocity Cruise segment (nKnots -> nKnots+1)
                        nSegs = nSegs + 1;
                        nKnots = nKnots + 1;
                        POLYQP_CAPACITY(nSegs <= MAX_SEG && nKnots <= MAX_KNOTS);
                        setSlab(FixCg, nKnots, r, makeCruiseKnot(p_out, v_cruise, r));
                        V1(KnotIsRest, nKnots) = true;
                        V1(SegT, nSegs) = T_cruise;
                        V1(IsCruise, nSegs) = true;
                        Col1(SegMaxVel, nSegs) = v_max;
                        Col1(SegMaxAccel, nSegs) = a_max;

                        // 3. Ramp-Down segment (nKnots -> nKnots+1 at wpB)
                        nSegs = nSegs + 1;
                        nKnots = nKnots + 1;
                        POLYQP_CAPACITY(nSegs <= MAX_SEG && nKnots <= MAX_KNOTS);
                        POLYQP_TRY(wpToFix(wpB, r, slab, isRest));
                        setSlab(FixCg, nKnots, r, slab);
                        V1(KnotIsRest, nKnots) = isRest;
                        Col1(B3Knot, nKnots) = wpB.B3BC; // FIX-F03: knot attitude from the waypoint
                        Col1(B3DotKnot, nKnots) = wpB.B3DotBC; // FIX-F03: knot attitude from the waypoint
                        if (wpB.StraightLine) {
                            getSlab(FixCg, nKnots, r, slab);
                            setSlab(FixCg, nKnots, r, pinTransverse(slab, dp, r));
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
            nSegs = nSegs + 1;
            nKnots = nKnots + 1;
            POLYQP_CAPACITY(nSegs <= MAX_SEG && nKnots <= MAX_KNOTS);
            POLYQP_TRY(wpToFix(wpB, r, slab, isRest));
            setSlab(FixCg, nKnots, r, slab);
            V1(KnotIsRest, nKnots) = isRest;
            Col1(B3Knot, nKnots) = wpB.B3BC; // FIX-F03: knot attitude from the waypoint
            Col1(B3DotKnot, nKnots) = wpB.B3DotBC; // FIX-F03: knot attitude from the waypoint
            if (wpB.StraightLine && (!dp.array().isNaN()).all() && mnorm(dp) > Real(1e-6)) {
                getSlab(FixCg, nKnots, r, slab);
                setSlab(FixCg, nKnots, r, pinTransverse(slab, dp, r));
            }
            Vec3 a_seg_init = a_max;
            if (Limits.HasVehicle) {
                const Real a_up_lim = mmax(Real(0.30), Real(0.80) * (Limits.T_hi / Limits.m_init - Limits.g));
                const Real a_dn_lim = mmax(Real(0.30), Real(0.80) * (Limits.g - Limits.T_lo / Limits.m_dry));
                V1(a_seg_init, 3) = minv({V1(a_seg_init, 3), a_up_lim, a_dn_lim});
            }
            V1(SegT, nSegs) = initSingleSegTime(dp, vA_prev, getFiber3(FixCg, 2, nKnots), v_max, a_seg_init, r, Limits.MinSegTime);
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
            Col1(B3Knot, nKnots) = Col1(B3Knot, nKnots - 1); // FIX-F03: a hold keeps the knot attitude
            Col1(B3DotKnot, nKnots) = Col1(B3DotKnot, nKnots - 1); // FIX-F03: a hold keeps the knot attitude
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
        SlabT block;
        getSlab(FixCg, k, r, block);
        V1(KnotFullyPinned, k) = (!block.array().isNaN()).all();
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
    Plan.B3DotBC = B3DotKnot; // FIX-F03: per-knot attitude rate for the P-space BC
    return PolyQPStatus::Ok;
}

static PolyQPStatus wpToFix(const TOADWaypointT &wp, int r, SlabT &slab, bool &isRest) {
    // Convert TOADWaypoint to 3 x r boundary condition slab
    slab.setConstant(3, r, NaN());
    Col1(slab, 1) = wp.PosBC;
    Col1(slab, 2) = wp.VelBC;
    Col1(slab, 3) = wp.AccelBC;
    if (r >= 4) {
        Col1(slab, 4) = wp.JerkBC;
    }
    isRest = !wp.IsPassAndGo;
    if (isRest) {
        for (int k = 2; k <= r; ++k) {
            for (int d = 1; d <= 3; ++d) {
                if (std::isnan(M1(slab, d, k))) {
                    M1(slab, d, k) = 0;
                }
            }
        }
    }
    return PolyQPStatus::Ok;
}

static SlabT makeCruiseKnot(const Vec3 &pos, const Vec3 &v_cruise, int r) {
    // Build fully-pinned constant-velocity transition knot
    SlabT slab = SlabT::Zero(3, r);
    Col1(slab, 1) = pos;
    Col1(slab, 2) = v_cruise;
    return slab;
}

static SlabT pinTransverse(SlabT slab, const Vec3 &dp, int r) {
    // Pin derivatives >= 1 to zero on axes transverse to straight-line displacement
    for (int d = 1; d <= 3; ++d) {
        if (std::fabs(V1(dp, d)) < Real(1e-6)) {
            for (int k = 2; k <= r; ++k) {
                if (std::isnan(M1(slab, d, k))) {
                    M1(slab, d, k) = 0;
                }
            }
        }
    }
    return slab;
}

static Real initSingleSegTime(const Vec3 &dp, const Vec3 &vA, const Vec3 &vB, const Vec3 &v_max, const Vec3 &a_max, int r, Real MinSegTime) {
    // Kinematic duration initializer for a single polynomial segment
    Real c_v;
    Real c_a;
    if (r == 3) {
        c_v = Real(15) / Real(8);
        c_a = Real(5.7735);
    } else if (r == 4) {
        c_v = Real(35) / Real(16);
        c_a = Real(7.5131);
    } else if (r == 5) {
        c_v = Real(315) / Real(128);
        c_a = Real(9.3720);
    } else {
        c_v = Real(693) / Real(256);
        c_a = Real(11.2);
    }

    Vec3 dp_abs = mabs(dp);
    for (int i = 1; i <= 3; ++i) {
        if (std::isnan(V1(dp_abs, i))) V1(dp_abs, i) = 0;
    }
    Vec3 dv_abs = mmax(mabs(vA), mabs(vB));
    for (int i = 1; i <= 3; ++i) {
        if (std::isnan(V1(dv_abs, i))) V1(dv_abs, i) = 0;
    }

    const Real T_v = maxv((c_v * (dp_abs.array() / v_max.array())).matrix());
    const Real T_a = maxv(msqrt((c_a * (dp_abs.array() / a_max.array())).matrix()));
    const Real T_dv = maxv((c_v * (dv_abs.array() / a_max.array())).matrix());
    const Real T = maxv({MinSegTime, T_v, T_a, T_dv});
    return T;
}

}
