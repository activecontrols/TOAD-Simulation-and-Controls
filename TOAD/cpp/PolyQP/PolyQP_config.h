// PolyQP transpile support header (capacities) for TOAD/Guidance/PolynomialQP @ 5ac6bd23
#pragma once

#ifndef POLYQP_MAX_WP
#define POLYQP_MAX_WP 16
#endif

#ifndef POLYQP_MAX_SEG
#define POLYQP_MAX_SEG 40
#endif

#ifndef POLYQP_MAX_KNOTS
#define POLYQP_MAX_KNOTS (POLYQP_MAX_SEG + 1)
#endif

#ifndef POLYQP_MAX_R
#define POLYQP_MAX_R 6
#endif

#ifndef POLYQP_MAX_NC
#define POLYQP_MAX_NC 12
#endif

#ifndef POLYQP_MAX_GRID
#define POLYQP_MAX_GRID 3999
#endif

// TRANSPILE-FLAG F24: derived capacities. sampleOutputs evaluates Nn = Nc + 1 <= POLYQP_MAX_GRID + 1 nodes; solveLeg needs nv = r * K <= POLYQP_MAX_R * POLYQP_MAX_KNOTS. Transpiled as-is per MATLAB.
#define POLYQP_MAX_NODES (POLYQP_MAX_GRID + 1)
#define POLYQP_MAX_NV (POLYQP_MAX_R * POLYQP_MAX_KNOTS)

static_assert(POLYQP_MAX_KNOTS == POLYQP_MAX_SEG + 1, "MATLAB: MAX_KNOTS = MAX_SEG + 1");
static_assert(POLYQP_MAX_R >= 3 && POLYQP_MAX_R <= 6, "MATLAB: FixOut = nan(3, 6, MAX_KNOTS) bounds r <= 6");
static_assert(POLYQP_MAX_NC >= 2 * POLYQP_MAX_R, "nc = 2 * r must fit POLYQP_MAX_NC");
static_assert(POLYQP_MAX_GRID >= 41, "MATLAB: evaluateConstraints uses at least 41 grid points");

// TRANSPILE-FLAG F25: fixed-max Eigen buffers (e.g. Rmat, nv x nv) exceed Eigen's default 128 KB static-size assert; the limit is disabled and every such buffer is static (see TRANSPILE_FLAGS.md, non-reentrancy). Transpiled as-is per MATLAB.
#if defined(EIGEN_WORLD_VERSION) && EIGEN_STACK_ALLOCATION_LIMIT
#error "Include PolyQP headers before any Eigen header (EIGEN_STACK_ALLOCATION_LIMIT must be 0 for PolyQP)"
#endif
#ifndef EIGEN_STACK_ALLOCATION_LIMIT
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#endif
