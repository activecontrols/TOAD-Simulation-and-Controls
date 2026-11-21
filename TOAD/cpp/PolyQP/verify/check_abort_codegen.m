function ok = check_abort_codegen()
% Verification only (never part of the flight build).
% The abort block (Guidance/AbortGuidance.m) must stay code-generatable so it can run in a Simulink MATLAB Function block.
% It was moved out of TrajectoryManager.m (commit fd08005 replaced that file) and, since the ASTRA integration Phase 6, has
% no trigger: it plans once from X_est. This runs codegen (lib, code only) on AbortGuidance with the ASTRA constants
% (LoadTOADParams(0)), a 15-state estimate, the last command u_cmd and the disturbance a_dist. Requires MATLAB Coder.
here    = fileparts(mfilename('fullpath'));
toadDir = fileparts(fileparts(fileparts(here)));
addpath(toadDir);
addpath(fullfile(toadDir, 'Guidance'));
addpath(fullfile(toadDir, 'Guidance', 'PolynomialQP'));
addpath(fullfile(toadDir, 'Helper'));

tdir = fullfile(tempdir, 'polyqp_abort_codegen');
if isfolder(tdir)
    rmdir(tdir, 's');
end
mkdir(tdir);

c = LoadTOADParams(0);
cfg = coder.config('lib');
cfg.GenCodeOnly = true;
try
    codegen('AbortGuidance', '-args', {zeros(15, 1), c, zeros(4, 1), zeros(3, 1)}, '-config', cfg, '-d', fullfile(tdir, 'codegen'));
    ok = true;
    fprintf('check_abort_codegen: AbortGuidance generates code: OK\n');
catch ME
    ok = false;
    fprintf('check_abort_codegen: FAILED: %s\n', ME.message);
end
end
