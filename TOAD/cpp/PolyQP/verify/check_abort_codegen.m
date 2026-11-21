function ok = check_abort_codegen()
% Verification only (never part of the flight build).
% TrajectoryManager.m runs in a Simulink MATLAB Function block (Guidance/TOAD_TrajectoryGen_3.slx), so the abort
% block must stay code-generatable. The shipped file has ABORT = 0, which lets code generation drop the whole abort
% branch, so this copies TrajectoryManager.m to a temporary folder with ABORT = 1 and runs codegen (lib, code only)
% against LoadTOADParams(1) plus a small synthetic constantsTOAD.Traj. Requires MATLAB Coder.
here    = fileparts(mfilename('fullpath'));
toadDir = fileparts(fileparts(fileparts(here)));
addpath(toadDir);
addpath(fullfile(toadDir, 'Guidance', 'PolynomialQP'));
addpath(fullfile(toadDir, 'Helper'));

tdir = fullfile(tempdir, 'polyqp_abort_codegen');
if isfolder(tdir)
    rmdir(tdir, 's');
end
mkdir(tdir);
src = fileread(fullfile(toadDir, 'Guidance', 'TrajectoryManager.m'));
assert(contains(src, 'ABORT = 0;'), 'ABORT = 0 not found in TrajectoryManager.m');
src = strrep(src, 'ABORT = 0;', 'ABORT = 1;');
fid = fopen(fullfile(tdir, 'TrajectoryManager.m'), 'w');
fwrite(fid, src);
fclose(fid);
addpath(tdir);
cleanP = onCleanup(@() rmpath(tdir));

c = LoadTOADParams(1);
c.Traj.Time = (0:0.1:1)';
c.Traj.States = zeros(11, 15);
c.Traj.Inputs = zeros(11, 4);
c.Traj.KTGain = zeros(2, 3, 6);
c.Traj.KRGain = zeros(2, 3, 6);
c.Traj.LAGain = zeros(2, 6, 3);
c.Traj.LTGain = zeros(2, 9, 6);
cfg = coder.config('lib');
cfg.GenCodeOnly = true;
try
    codegen('TrajectoryManager', '-args', {0, zeros(15, 1), c}, '-config', cfg, '-d', fullfile(tdir, 'codegen'));
    ok = true;
    fprintf('check_abort_codegen: TrajectoryManager (ABORT = 1) generates code: OK\n');
catch ME
    ok = false;
    fprintf('check_abort_codegen: FAILED: %s\n', ME.message);
end
end
