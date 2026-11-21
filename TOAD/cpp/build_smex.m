function build_smex(clean)
% build_smex: High-speed incremental compiler for ASTRAv2_SMEX
% Only recompiles .cpp files that have changed since their .obj was created.
% Auto-discovers and compiles all sources/headers in AbortCode/.

if nargin < 1, clean = false; end

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end

% Locate cpp directory (works both from sandbox/deployed/cpp/ and cpp/)
cppDir = '';
cand = thisDir;
for depth = 0:4
    if isfile(fullfile(cand, 'ASTRAv2_SMEX.cpp'))
        cppDir = cand;
        break;
    elseif isfile(fullfile(cand, 'cpp', 'ASTRAv2_SMEX.cpp'))
        cppDir = fullfile(cand, 'cpp');
        break;
    end
    parent = fileparts(cand);
    if strcmp(parent, cand)
        break;
    end
    cand = parent;
end
assert(~isempty(cppDir), 'Could not locate cpp/ASTRAv2_SMEX.cpp from %s.', thisDir);

extrasDir = fullfile(cppDir, 'Extras');
abortDir  = fullfile(cppDir, 'AbortCode');

eigenDir = fullfile(extrasDir, 'eigen-3.4.0');
if ~isfile(fullfile(eigenDir, 'Eigen', 'Dense'))
    eigenDir = fullfile(extrasDir, 'eigen-3.4.0', 'eigen-3.4.0');
end

assert(isfile(fullfile(extrasDir, 'ArduinoEigenDense.h')), 'Extras/ArduinoEigenDense.h not found.');
assert(isfile(fullfile(eigenDir, 'Eigen', 'Dense')), 'Eigen/Dense not found under %s.', eigenDir);
assert(isfolder(abortDir), 'AbortCode folder not found under %s.', cppDir);

objDir = fullfile(cppDir, 'build_objs');
if ~isfolder(objDir)
    mkdir(objDir);
end

if ispc
    objExt = '.obj';
else
    objExt = '.o';
end

if clean
    fprintf('Cleaning object files in %s...\n', objDir);
    staleObjs = dir(fullfile(objDir, ['*' objExt]));
    if ~isempty(staleObjs)
        delete(fullfile(objDir, ['*' objExt]));
    end
end

sources = {
    'matlab_helpers.cpp'
    'DigitalNF.cpp'
    'Q_gen.cpp'
    'GroundEstimator.cpp'
    'FlightEstimator2.cpp'
    'LESO_Position.cpp'
    'LESO_Attitude.cpp'
    'ASTRAv2_Controller.cpp'
    'TrajectoryManager.cpp'
    'controller_and_estimator.cpp'
    'ASTRAv2_SMEX.cpp'
};

abortCpp = dir(fullfile(abortDir, '*.cpp'));
assert(~isempty(abortCpp), 'No .cpp files found in %s.', abortDir);
abortSources = cell(numel(abortCpp), 1);
for k = 1:numel(abortCpp)
    abortSources{k} = fullfile('AbortCode', abortCpp(k).name);
end
sources = [sources; abortSources];

% Check header dependencies timestamps
headers = {
    'matlab_funcs.h'
    'astra_structs.h'
    'controller_and_estimator.h'
    'TrajectoryManager.h'
    'ASTRA_Trajectory_Data.h'
};

abortHdr = dir(fullfile(abortDir, '*.h'));
abortHeaders = cell(numel(abortHdr), 1);
for k = 1:numel(abortHdr)
    abortHeaders{k} = fullfile('AbortCode', abortHdr(k).name);
end
headers = [headers; abortHeaders];

maxHeaderTime = 0;
for h = 1:numel(headers)
    hf = fullfile(cppDir, headers{h});
    if isfile(hf)
        d = dir(hf);
        if d.datenum > maxHeaderTime
            maxHeaderTime = d.datenum;
        end
    end
end

objFiles = cell(size(sources));
recompiled = 0;
maxObjTime = 0;
tTotal = tic;

for k = 1:numel(sources)
    srcFile = fullfile(cppDir, sources{k});
    assert(isfile(srcFile), 'Source file not found: %s', srcFile);
    [~, baseName] = fileparts(sources{k});
    objFile = fullfile(objDir, [baseName objExt]);
    objFiles{k} = objFile;

    srcInfo = dir(srcFile);
    objInfo = dir(objFile);

    needsCompile = isempty(objInfo) || ...
                   (srcInfo.datenum > objInfo.datenum) || ...
                   (maxHeaderTime > objInfo.datenum);

    if needsCompile
        fprintf('[%2d/%2d] Compiling %s... ', k, numel(sources), sources{k});
        t0 = tic;
        mex('-c', '-D_USE_MATH_DEFINES', ...
            ['-I' extrasDir], ...
            ['-I' eigenDir], ...
            ['-I' abortDir], ...
            ['-I' cppDir], ...
            srcFile, ...
            '-outdir', objDir);
        fprintf('(%.2f s)\n', toc(t0));
        recompiled = recompiled + 1;
        objInfo = dir(objFile);
    else
        fprintf('[%2d/%2d] Cached:    %s\n', k, numel(sources), sources{k});
    end

    if ~isempty(objInfo) && objInfo.datenum > maxObjTime
        maxObjTime = objInfo.datenum;
    end
end

% Link step
mexTarget = fullfile(cppDir, ['ASTRAv2_SMEX.', mexext]);
targetInfo = dir(mexTarget);
needsLink = (recompiled > 0) || isempty(targetInfo) || (maxObjTime > targetInfo.datenum);

if needsLink
    fprintf('Linking ASTRAv2_SMEX... ');
    if mislocked('ASTRAv2_SMEX')
        munlock('ASTRAv2_SMEX');
    end
    clear('ASTRAv2_SMEX'); % Unload MEX binary before overwriting
    t0 = tic;
    mex(objFiles{:}, '-output', 'ASTRAv2_SMEX', '-outdir', cppDir);
    fprintf('(%.2f s)\n', toc(t0));
else
    fprintf('ASTRAv2_SMEX is already up to date. Nothing to link.\n');
end

fprintf('Build complete in %.2f s (Recompiled %d of %d units).\n', toc(tTotal), recompiled, numel(sources));
end
