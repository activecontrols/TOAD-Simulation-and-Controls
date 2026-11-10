function build_smex(clean)
% build_smex: High-speed incremental compiler for ASTRAv2_SMEX
% Only recompiles .cpp files that have changed since their .obj was created.
% Cuts build time from >5 minutes down to ~5-15 seconds for single-file edits!

if nargin < 1, clean = false; end

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end

% Locate cpp directory
if exist(fullfile(thisDir, 'ASTRAv2_SMEX.cpp'), 'file')
    cppDir = thisDir;
else
    cppDir = fullfile(fileparts(fileparts(thisDir)), 'cpp');
end

eigenDir  = fullfile(cppDir, 'Extras', 'eigen-3.4.0');
extrasDir = fullfile(cppDir, 'Extras');
objDir    = fullfile(cppDir, 'build_objs');
if ~exist(objDir, 'dir'), mkdir(objDir); end

if clean
    fprintf('Cleaning object files in %s...\n', objDir);
    delete(fullfile(objDir, '*.obj'));
end

sources = {
    'matlab_helpers.cpp'
    'DigitalNF.cpp'
    'Q_gen.cpp'
    'GroundEstimator.cpp'
    'FlightEstimator.cpp'
    'FlightEstimator2.cpp'
    'LESO_Position.cpp'
    'LESO_Attitude.cpp'
    'ASTRAv2_Controller.cpp'
    'controller_and_estimator.cpp'
    'ASTRAv2_SMEX.cpp'
};

% Check header dependencies timestamps
headers = {'matlab_funcs.h', 'astra_structs.h', 'controller_and_estimator.h'};
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
tTotal = tic;

for k = 1:numel(sources)
    srcFile = fullfile(cppDir, sources{k});
    [~, baseName, ~] = fileparts(sources{k});
    objFile = fullfile(objDir, [baseName '.obj']);
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
            ['-I' cppDir], ...
            srcFile, ...
            '-outdir', objDir);
        fprintf('(%.2f s)\n', toc(t0));
        recompiled = recompiled + 1;
    else
        fprintf('[%2d/%2d] Cached:    %s\n', k, numel(sources), sources{k});
    end
end

% Link step
mexTarget = fullfile(cppDir, 'ASTRAv2_SMEX.cpp');
targetInfo = dir(mexTarget);
needsLink = (recompiled > 0) || isempty(targetInfo);

if needsLink
    fprintf('Linking ASTRAv2_SMEX... ');
    t0 = tic;
    mex(objFiles{:}, '-output', 'ASTRAv2_SMEX', '-outdir', cppDir);
    fprintf('(%.2f s)\n', toc(t0));
else
    fprintf('ASTRAv2_SMEX is already up to date. Nothing to link.\n');
end

fprintf('Build complete in %.2f s (Recompiled %d of %d units).\n', toc(tTotal), recompiled, numel(sources));
end
