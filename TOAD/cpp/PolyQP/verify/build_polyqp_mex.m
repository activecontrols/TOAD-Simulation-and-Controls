function build_polyqp_mex()
% Verification only (never part of the flight build).
% build_polyqp_mex -> polyqp_mex (the solver is double only; the former float build, -DPOLYQP_REAL=float, was removed
% when the helpers moved to cpp/matlab_funcs.h, which are double by construction).

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
