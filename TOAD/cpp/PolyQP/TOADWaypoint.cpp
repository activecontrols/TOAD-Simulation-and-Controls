// Transpiled from TOAD/Guidance/PolynomialQP/TOADWaypoint.m @ 5ac6bd23
#include "PolyQP.h"
#include "matlab_compat.h"

namespace PolyQP {

TOADWaypointT TOADWaypoint(const Vec3 &pos, const TOADWaypointOptions &options) {
    // TOADWaypoint Builder
    // FIX-F03: thrust axis b3 at this knot for the P-space BC (NaN: b3 = e3)
    // FIX-F03: db3/dt at this knot for the P-space velocity BC (NaN: not applied)
    POLYQP_ASSERT(options.Label != nullptr);

    // Fixed-length 9-character label for codegen
    const int FIXED_LEN = 9;
    const int currentLen = static_cast<int>(std::strlen(options.Label));
    char finalLabel[10];
    if (currentLen < FIXED_LEN) {
        std::memcpy(finalLabel, options.Label, static_cast<size_t>(currentLen));
        for (int q = currentLen + 1; q <= FIXED_LEN; ++q) finalLabel[q - 1] = ' ';
    } else if (currentLen > FIXED_LEN) {
        std::memcpy(finalLabel, options.Label, static_cast<size_t>(FIXED_LEN));
    } else {
        std::memcpy(finalLabel, options.Label, static_cast<size_t>(FIXED_LEN));
    }
    finalLabel[FIXED_LEN] = '\0';

    // Rest knot check
    const bool isPassAndGo = std::isnan(options.HoldTime);

    // Resolve rest knot derivative defaults
    Vec3 velBC = options.VelBC;
    Vec3 accelBC = options.AccelBC;
    Vec3 jerkBC = options.JerkBC;
    if (!isPassAndGo) {
        for (int i = 1; i <= 3; ++i) {
            if (std::isnan(V1(velBC, i))) V1(velBC, i) = 0;
            if (std::isnan(V1(accelBC, i))) V1(accelBC, i) = 0;
            if (std::isnan(V1(jerkBC, i))) V1(jerkBC, i) = 0;
        }
    }

    // Construct fixed-layout waypoint struct
    TOADWaypointT wp;
    wp.PosBC = pos;
    wp.VelBC = velBC;
    wp.AccelBC = accelBC;
    wp.JerkBC = jerkBC;
    wp.HoldTime = options.HoldTime;
    wp.MaxVel = options.MaxVel;
    wp.MaxAccel = options.MaxAccel;
    wp.StraightLine = options.StraightLine;
    wp.HDGRef = options.HDGRef;
    std::memcpy(wp.Label, finalLabel, sizeof(wp.Label));
    wp.IsPassAndGo = isPassAndGo;
    wp.B3BC = options.B3BC; // FIX-F03: optional knot attitude for the P-space BC
    wp.B3DotBC = options.B3DotBC; // FIX-F03: optional knot attitude rate for the P-space BC
    return wp;
}

}
