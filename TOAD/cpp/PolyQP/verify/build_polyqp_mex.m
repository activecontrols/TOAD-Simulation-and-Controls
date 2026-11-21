function build_polyqp_mex(variant)
% Verification only (never part of the flight build).
% build_polyqp_mex         -> polyqp_mex      : the flight capacities of PolyQP_config.h (abort-sized), used for every abort case
% build_polyqp_mex('wide') -> polyqp_mex_wide : the same sources with larger capacities (-D), for the missions, fuzzed legs and
%                                               grids that exceed the abort envelope by design
% (The solver is double only; the former float build, -DPOLYQP_REAL=float, was removed with the helper move.)
if nargin < 1
    variant = 'flight';
end

here      = fileparts(mfilename('fullpath'));
polyqpDir = fileparts(here);
cppDir    = fileparts(polyqpDir);
eigenDir  = fullfile(cppDir, 'Extras', 'eigen-3.4.0', 'eigen-3.4.0');

assert(isfile(fullfile(eigenDir, 'Eigen', 'Dense')), 'Eigen/Dense not found under %s.', eigenDir);
assert(isfile(fullfile(cppDir, 'matlab_funcs.h')), 'matlab_funcs.h not found under %s.', cppDir);

flightSrcs = {'PolyQP.cpp', 'resolveWaypoints.cpp', 'assembleQP.cpp', 'solveLeg.cpp', 'allocateTimes.cpp', ...
    'evaluateConstraints.cpp', 'flatnessInvert.cpp', 'allocateControls.cpp', 'sampleOutputs.cpp', 'PolyEval.cpp', ...
    'buildVehicleLimits.cpp', 'TOADWaypoint.cpp', 'AbortGuidance.cpp', 'abortFallback.cpp'};
% The double-precision MATLAB helpers (quatRot, DCM_Quat_Conversion, ComputeJtot, ...) live in cpp/matlab_helpers.cpp / matlab_funcs.h.
srcs = [{fullfile(here, 'polyqp_mex.cpp')}, fullfile(polyqpDir, flightSrcs), {fullfile(cppDir, 'matlab_helpers.cpp')}];

defs = {'-DPOLYQP_DEBUG_BOUNDS', '-DPOLYQP_ERROR_MESSAGES', '-D_USE_MATH_DEFINES'};
name = 'polyqp_mex';
if strcmpi(variant, 'wide')
    % 16 segments (17 knots, nv = 85, every unknown may be free) and MATLAB's 3999-sample grid. The 15 x 4000 grid
    % buffers exceed Eigen's 128 KB static-size check, so it is disabled in every translation unit of this build.
    defs = [defs, {'-DPOLYQP_MAX_WP=16', '-DPOLYQP_MAX_SEG=16', '-DPOLYQP_MAX_FREE=85', '-DPOLYQP_MAX_GRID=3999', '-DEIGEN_STACK_ALLOCATION_LIMIT=0'}];
    name = 'polyqp_mex_wide';
elseif ~strcmpi(variant, 'flight')
    error('build_polyqp_mex:variant', 'variant must be ''flight'' or ''wide''.');
end

if ispc
    flags = {'COMPFLAGS=$COMPFLAGS /std:c++17 /fp:precise /bigobj'};
else
    flags = {'CXXFLAGS=$CXXFLAGS -std=c++17 -ffp-contract=off'};
end

mex('-outdir', here, '-output', name, ...
    ['-I' polyqpDir], ['-I' cppDir], ['-I' fullfile(cppDir, 'Extras')], ['-I' eigenDir], ...
    defs{:}, flags{:}, srcs{:});

fprintf('Built %s.\n', name);
end
