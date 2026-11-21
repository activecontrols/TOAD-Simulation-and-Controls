% abort_code_build.m  (run from inside the cpp folder)

eigenDir = 'Extras/eigen-3.4.0';

% Fail early if a path is wrong
assert(isfile('astra_mex.cpp'),                  'astra_mex.cpp not found. cd into the cpp folder first.');
assert(isfile('Extras/ArduinoEigenDense.h'),     'Extras/ArduinoEigenDense.h not found.');
assert(isfile([eigenDir '/Eigen/Dense']),        'Eigen/Dense not found under %s. Check the doubled folder.', eigenDir);

mex('AbortCode/computeAbort.cpp', ...
    'AbortCode/flatnessInvert.cpp', ...
    'AbortCode/allocateControls.cpp');

disp('Build finished.');
