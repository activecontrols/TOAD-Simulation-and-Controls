function build_polyqp_mex(precision)
% Verification only (never part of the flight build).
% build_polyqp_mex          -> polyqp_mex       (Real = double)
% build_polyqp_mex('float') -> polyqp_mex_float (Real = float, -DPOLYQP_REAL=float)
if nargin < 1
    precision = 'double';
end

here      = fileparts(mfilename('fullpath'));
polyqpDir = fileparts(here);
cppDir    = fileparts(polyqpDir);
eigenDir  = fullfile(cppDir, 'Extras', 'eigen-3.4.0', 'eigen-3.4.0');

assert(isfile(fullfile(eigenDir, 'Eigen', 'Dense')), 'Eigen/Dense not found under %s.', eigenDir);
assert(isfile(fullfile(cppDir, 'matlab_funcs.h')), 'matlab_funcs.h not found under %s.', cppDir);

flightSrcs = {'PolyQP.cpp', 'resolveWaypoints.cpp', 'assembleQP.cpp', 'solveLeg.cpp', 'allocateTimes.cpp', ...
    'evaluateConstraints.cpp', 'flatnessInvert.cpp', 'allocateControls.cpp', 'sampleOutputs.cpp', 'PolyEval.cpp', ...
    'buildVehicleLimits.cpp', 'TOADWaypoint.cpp', 'ComputeJtot.cpp', 'DCM_Quat_Conversion.cpp', 'AbortGuidance.cpp', 'abortFallback.cpp', 'quatRot.cpp'};
srcs = [{fullfile(here, 'polyqp_mex.cpp')}, fullfile(polyqpDir, flightSrcs)];

defs = {'-DPOLYQP_DEBUG_BOUNDS', '-DPOLYQP_ERROR_MESSAGES', '-D_USE_MATH_DEFINES'};
name = 'polyqp_mex';
if strcmpi(precision, 'float')
    defs{end + 1} = '-DPOLYQP_REAL=float';
    name = 'polyqp_mex_float';
elseif ~strcmpi(precision, 'double')
    error('build_polyqp_mex:precision', 'precision must be ''double'' or ''float''.');
end

if ispc
    flags = {'COMPFLAGS=$COMPFLAGS /std:c++17 /fp:precise /bigobj'};
else
    flags = {'CXXFLAGS=$CXXFLAGS -std=c++17 -ffp-contract=off'};
end

mex('-outdir', here, '-output', name, ...
    ['-I' polyqpDir], ['-I' cppDir], ['-I' fullfile(cppDir, 'Extras')], ['-I' eigenDir], ...
    defs{:}, flags{:}, srcs{:});

fprintf('Built %s (%s precision).\n', name, precision);
end
