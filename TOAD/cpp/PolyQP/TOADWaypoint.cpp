// Transpiled from TOAD/Guidance/PolynomialQP/TOADWaypoint.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

TOADWaypointT TOADWaypoint(const Vec3 &pos, const TOADWaypointOptions &options, const Vector3b &posFixed) {
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
    const bool isPassAndGo = !options.HoldTimeGiven; // NO-NAN: isnan(HoldTime)

    // Resolve rest knot derivative defaults
    Vec3 velBC = options.VelBC;
    Vector3b velFixed = options.VelFixed;
    Vec3 accelBC = options.AccelBC;
    Vector3b accelFixed = options.AccelFixed;
    Vec3 jerkBC = options.JerkBC;
    Vector3b jerkFixed = options.JerkFixed;
    if (!isPassAndGo) {
        // NO-NAN: isnan(x) is !fixed; a defaulted entry becomes a fixed 0
        for (int i = 1; i <= 3; ++i) {
            if (!V1(velFixed, i)) {
                V1(velBC, i) = 0;
                V1(velFixed, i) = true;
            }
            if (!V1(accelFixed, i)) {
                V1(accelBC, i) = 0;
                V1(accelFixed, i) = true;
            }
            if (!V1(jerkFixed, i)) {
                V1(jerkBC, i) = 0;
                V1(jerkFixed, i) = true;
            }
        }
    }

    // Construct fixed-layout waypoint struct
    TOADWaypointT wp;
    wp.PosBC = pos;
    wp.PosFixed = posFixed;
    for (int i = 1; i <= 3; ++i) {
        if (!V1(posFixed, i)) V1(wp.PosBC, i) = 0; // NO-NAN: a free entry holds 0
    }
    wp.VelBC = velBC;
    wp.VelFixed = velFixed;
    wp.AccelBC = accelBC;
    wp.AccelFixed = accelFixed;
    wp.JerkBC = jerkBC;
    wp.JerkFixed = jerkFixed;
    wp.HoldTime = options.HoldTimeGiven ? options.HoldTime : Real(0); // NO-NAN: pass-and-go (MATLAB NaN) holds 0
    wp.MaxVel = options.MaxVel;
    wp.MaxAccel = options.MaxAccel;
    wp.StraightLine = options.StraightLine;
    wp.HDGRef = options.HDGRef;
    std::memcpy(wp.Label, finalLabel, sizeof(wp.Label));
    wp.IsPassAndGo = isPassAndGo;
    wp.B3BC = options.B3BC; // FIX-F03: optional knot attitude for the P-space BC
    wp.B3Fixed = options.B3Fixed;
    wp.B3DotBC = options.B3DotBC; // FIX-F03: optional knot attitude rate for the P-space BC
    wp.B3DotFixed = options.B3DotFixed;
    return wp;
}

}
