function LiveSim_Runner(varargin)
% LIVESIM_RUNNER Real-Time Hardware-in-the-Loop Simulation Runner for TOAD
% Interfaces the MATLAB Lumped Parameter model with the Ground Control UI via UDP.
%
% Features:
%   - Automated C++ MEX Dependency Tracking: Checks if core source files have been
%     modified since last compilation. Recompiles automatically if changes are detected;
%     otherwise reuses last available compiled binary.
%   - Telemetry Stream (Sim -> UI): UDP 127.0.0.1:9000 (100-byte EC_FMT packet at 20 Hz)
%   - Valve Command (UI -> Sim):   UDP 127.0.0.1:9001 (12-byte packet via NIO channel)
%   - Real-Time Pacing: 100 Hz (dt = 0.010 s) synchronized with wall-clock time.
%
% Optional Name-Value Arguments:
%   'ForceCompile' - Logical (true to force C++ MEX recompile, default false)
%   'MaxSteps'     - Numeric (maximum steps to execute, default Inf for continuous run)
%   'Quiet'        - Logical (suppress routine step printouts, default false)

    % Parse optional input arguments
    p = inputParser;
    addParameter(p, 'ForceCompile', false, @islogical);
    addParameter(p, 'MaxSteps', Inf, @isnumeric);
    addParameter(p, 'Quiet', false, @islogical);
    parse(p, varargin{:});
    opts = p.Results;

    fprintf('=====================================================\n');
    fprintf('   TOAD Lumped Parameter Real-Time Simulation Runner \n');
    fprintf('=====================================================\n');

    % Setup paths robustly
    baseDir = fileparts(mfilename('fullpath'));
    
    coreCandidates = {
        fullfile(baseDir, '..', 'Core'), ...
        fullfile(baseDir, 'Core'), ...
        fullfile(pwd, 'sandbox', 'experiments', 'Core'), ...
        fullfile(pwd, 'TOAD_LumpedParam', 'Core')
    };
    coreDir = '';
    for c = 1:numel(coreCandidates)
        if exist(coreCandidates{c}, 'dir')
            coreDir = coreCandidates{c};
            break;
        end
    end

    configCandidates = {
        fullfile(baseDir, '..', 'Config'), ...
        fullfile(baseDir, 'Config'), ...
        fullfile(pwd, 'sandbox', 'experiments', 'Config'), ...
        fullfile(pwd, 'TOAD_LumpedParam', 'Config')
    };
    configDir = '';
    for c = 1:numel(configCandidates)
        if exist(configCandidates{c}, 'dir')
            configDir = configCandidates{c};
            break;
        end
    end

    thermalDir = fullfile(baseDir, '..', 'Thermal');
    propsDir   = fullfile(baseDir, '..', 'Fluid Properties');

    % Also find root TOAD_LumpedParam for FluidProperties and tables
    repoCandidates = {
        fullfile(baseDir, '..'), ...
        fullfile(baseDir, '..', '..', '..', 'TOAD_LumpedParam'), ...
        fullfile(pwd, 'TOAD_LumpedParam')
    };
    for r = 1:numel(repoCandidates)
        rp = repoCandidates{r};
        if exist(fullfile(rp, 'Core'), 'dir'), addpath(fullfile(rp, 'Core')); end
        if exist(fullfile(rp, 'Fluid Properties'), 'dir'), addpath(fullfile(rp, 'Fluid Properties')); end
    end

    if ~isempty(coreDir), addpath(coreDir, '-begin'); end
    if ~isempty(configDir), addpath(configDir, '-begin'); end
    if exist(thermalDir, 'dir'), addpath(thermalDir); end
    if exist(propsDir, 'dir'), addpath(propsDir); end
    addpath(baseDir, '-begin');

    %% --- 1. Check & Conditional C++ MEX Compilation ---
    if ~isempty(coreDir)
        checkAndCompileMEX(coreDir, opts.ForceCompile);
    end

    %% --- 2. Initialize System Topology ---
    fprintf('Initializing Flight-Ready System Topology...\n');
    if exist('BuildTOADSystem', 'file')
        [System, State] = BuildTOADSystem();
        useUniversal = true;
    else
        [System, State] = BuildTOADSystem_Merged();
        useUniversal = false;
    end

    dt = 0.010; % 100 Hz simulation rate
    psi2Pa = 6894.757;

    %% --- 3. Setup Zero-Toolbox Java UDP Sockets ---
    fprintf('Configuring UDP Sockets (Telemetry: 9000, Commands: 9001)...\n');
    sendSocket = [];
    cmdChannel = [];
    cmdBuf     = [];

    try
        sendSocket = java.net.DatagramSocket();
        destAddr   = java.net.InetAddress.getByName('127.0.0.1');
        destPort   = 9000;
    catch ME
        warning('LiveSim_Runner:SocketInit', 'Could not open telemetry socket: %s', ME.message);
    end

    try
        cmdChannel = java.nio.channels.DatagramChannel.open();
        cmdChannel.configureBlocking(false);
        cmdChannel.socket().bind(java.net.InetSocketAddress(9001));
        cmdBuf = java.nio.ByteBuffer.allocate(64);
        fprintf('Command receiver listening on UDP port 9001 (non-blocking NIO channel).\n');
    catch ME
        warning('LiveSim_Runner:CmdBind', 'Could not bind UDP port 9001 (%s). Running in broadcast-only mode.', ME.message);
    end

    % Cleanup sockets reliably on exit / Ctrl+C
    cleanupObj = onCleanup(@() cleanupSockets(sendSocket, cmdChannel));

    % Initial Link States
    LinkStates = struct();
    linkNames = fieldnames(System.Links);
    for k = 1:length(linkNames)
        fn = linkNames{k};
        LinkStates.(fn) = System.Links.(fn).State;
    end

    if isinf(opts.MaxSteps)
        fprintf('Ready! Starting 100 Hz Real-Time Execution Loop. Press Ctrl+C to stop.\n');
    else
        fprintf('Ready! Executing test batch (%d steps)...\n', opts.MaxSteps);
    end

    tStart = tic;
    simTime = 0.0;
    stepCount = 0;

    try
        while stepCount < opts.MaxSteps
            stepCount = stepCount + 1;

            %% A. Non-blocking Receive: Drain UI Valve Commands
            if ~isempty(cmdChannel)
                while true
                    cmdBuf.clear();
                    sender = cmdChannel.receive(cmdBuf);
                    if isempty(sender)
                        break; % No more packets queued
                    end
                    cmdBuf.flip();
                    if cmdBuf.remaining() >= 12
                        rawBytes = zeros(1, 12, 'int8');
                        for b = 1:12
                            rawBytes(b) = cmdBuf.get();
                        end
                        cmdMask = typecast(rawBytes(1:4), 'uint32');
                        oxThrt  = typecast(rawBytes(5:8), 'single');
                        fuThrt  = typecast(rawBytes(9:12), 'single');

                        LinkStates = applyValveMask(LinkStates, cmdMask, oxThrt, fuThrt);
                    end
                end
            end

            %% B. Step Physics Simulation
            if useUniversal
                [State, ~, System] = StepSimulation(State, dt, System, LinkStates);
            else
                [State, ~, System] = StepSimulation_Live(State, dt, System, LinkStates);
            end
            simTime = simTime + dt;

            %% C. Pack 100-byte EC_FMT Telemetry Packet
            if ~isempty(sendSocket) && mod(stepCount, 5) == 0
                % Assemble telemetry at 20 Hz
                pktBytes = assembleTelemetryPacket(State, System, stepCount, psi2Pa);
                sendPkt = java.net.DatagramPacket(typecast(pktBytes, 'int8'), 100, destAddr, destPort);
                sendSocket.send(sendPkt);
            end

            %% D. Real-Time Pacing (Wall-clock synchronization)
            elapsed = toc(tStart);
            if elapsed > simTime + 0.050
                simTime = elapsed; % Runaway catch-up protection
            end

            sleepSec = simTime - elapsed;
            if sleepSec > 0.001
                java.lang.Thread.sleep(int64(sleepSec * 1000));
            else
                java.lang.Thread.yield();
            end

            if ~opts.Quiet && mod(stepCount, 500) == 0 % Periodic 5-second status update
                pt_n2  = State.Nodes.TK_N2.P / psi2Pa;
                pt_ox  = State.Nodes.TK_O2_01.P / psi2Pa;
                pt_fu  = State.Nodes.TK_FU_01.P / psi2Pa;
                pt_pc  = State.Nodes.SKIPPER.P / psi2Pa;
                fprintf('t = %.1f s | COPV: %.0f psi | LOX: %.1f psi | FU: %.1f psi | Pc: %.1f psi\n', ...
                    simTime, pt_n2, pt_ox, pt_fu, pt_pc);
            end
        end
    catch ME
        if ~strcmp(ME.identifier, 'MATLAB:interruption')
            fprintf('Loop ended: %s\n', ME.message);
        else
            fprintf('\nSimulation stopped by user.\n');
        end
    end
end

%% =========================================================================
%% --- HELPER: Check and Conditional MEX Compilation ---
%% =========================================================================
function checkAndCompileMEX(coreDir, force)
    if nargin < 2, force = false; end

    mexBinaryName = ['CalculateLinkFlow_mex.' mexext];
    mexFile = fullfile(coreDir, mexBinaryName);
    cppFile = fullfile(coreDir, 'CalculateLinkFlow_mex.cpp');

    if ~exist(cppFile, 'file')
        % C++ source not present in coreDir
        return;
    end

    needsCompile = force;
    triggerReason = '';

    if force
        triggerReason = 'User forced recompile via ForceCompile=true';
    elseif ~exist(mexFile, 'file')
        needsCompile = true;
        triggerReason = sprintf('Compiled binary "%s" not found', mexBinaryName);
    else
        mexInfo = dir(mexFile);
        mexDate = mexInfo.datenum;

        % Check core source files that dictate the kernel implementation
        srcFiles = {
            cppFile, ...
            fullfile(coreDir, 'CalculateLinkFlow.m')
        };

        for k = 1:numel(srcFiles)
            srcPath = srcFiles{k};
            if exist(srcPath, 'file')
                srcInfo = dir(srcPath);
                if srcInfo.datenum > mexDate
                    needsCompile = true;
                    triggerReason = sprintf('Source file "%s" is newer than binary (%s vs %s)', ...
                        srcInfo.name, srcInfo.date, mexInfo.date);
                    break;
                end
            end
        end
    end

    if needsCompile
        fprintf('\n-----------------------------------------------------\n');
        fprintf('  [C++ MEX Kernel Update Check: RECOMPILE TRIGGERED]\n');
        fprintf('  Reason: %s\n', triggerReason);
        fprintf('  Compiling %s with MSVC 2022...\n', mexBinaryName);
        tStartComp = tic;
        try
            mex('-O', cppFile, '-outdir', coreDir);
            tComp = toc(tStartComp);
            fprintf('  Compilation successful in %.2f s -> %s\n', tComp, mexFile);
            % Clear in-memory symbol handles so MATLAB reloads the newly compiled binary
            clear CalculateLinkFlow_mex;
            clear CalculateLinkFlow;
        catch ME
            warning('LiveSim_Runner:MEXCompilationFailed', ...
                'MEX compilation failed: %s\nFalling back to existing binary or interpreted MATLAB solver.', ME.message);
        end
        fprintf('-----------------------------------------------------\n\n');
    else
        fprintf('C++ MEX kernel is up to date (%s). Using last available compile.\n', mexBinaryName);
    end
end

%% =========================================================================
%% --- HELPER: Assemble 100-Byte EC_FMT Telemetry Packet ---
%% =========================================================================
function pktBytes = assembleTelemetryPacket(State, System, stepCount, psi2Pa)
    % Header (4 bytes)
    hdr = uint8([mod(stepCount, 256), 0, 0, 0]);

    % Pressures (12 floats, 48 bytes)
    pt_vals = single(zeros(1, 12));
    pt_vals(1)  = single(State.Nodes.TK_N2.P / psi2Pa);
    pt_vals(2)  = single(State.Nodes.TK_O2_01.P / psi2Pa);
    pt_vals(3)  = single(State.Nodes.TK_FU_01.P / psi2Pa);
    pt_vals(4)  = single(State.Nodes.Purge_Manifold.P / psi2Pa);
    pt_vals(5)  = single(State.Nodes.OX_Manifold.P / psi2Pa);
    pt_vals(6)  = single(State.Nodes.FU_Manifold.P / psi2Pa);
    pt_vals(7)  = single(State.Nodes.DART_Chamber.P / psi2Pa);
    if isfield(State.Nodes, 'TK_N2_BULK')
        pt_vals(8) = single(State.Nodes.TK_N2_BULK.P / psi2Pa);
    end
    pt_vals(9)  = single(State.Nodes.DART_Chamber.P / psi2Pa);
    pt_vals(10) = single(State.Nodes.SKIPPER.P / psi2Pa);

    % Temperatures (6 floats, 24 bytes)
    tc_vals = single(zeros(1, 6));
    tc_vals(1) = single(State.Nodes.TK_N2.T);
    if isfield(State.Nodes.TK_O2_01, 'Liquid')
        tc_vals(2) = single(State.Nodes.TK_O2_01.Liquid.T);
    else
        tc_vals(2) = single(State.Nodes.TK_O2_01.T);
    end
    tc_vals(3) = single(State.Nodes.OX_Manifold.T);
    tc_vals(4) = single(State.Nodes.FU_Manifold.T);
    if isfield(State, 'Thermal') && isfield(State.Thermal, 'Regen')
        tc_vals(5) = single(State.Thermal.Regen.T_wall);
    elseif isfield(State.Nodes, 'SKIPPER')
        tc_vals(5) = single(State.Nodes.SKIPPER.T);
    end

    % Valve State Mask (1 uint32, 4 bytes)
    vMask = buildValveMask(System);

    % Throttle Angles (2 floats, 8 bytes)
    ox_ang = 0.0; fu_ang = 0.0;
    if isfield(System.Links, 'BV_02_04'), ox_ang = System.Links.BV_02_04.State;
    elseif isfield(System.Links, 'BV_O2_04'), ox_ang = System.Links.BV_O2_04.State;
    end
    if isfield(System.Links, 'BV_FU_04'), fu_ang = System.Links.BV_FU_04.State; end
    vAngles = single([ox_ang, fu_ang]);

    % Fill Levels (3 floats, 12 bytes)
    v_lox_frac = 0.0; v_fu_frac = 0.0;
    if isfield(State.Nodes.TK_O2_01, 'Liquid') && isfield(State.Nodes.TK_O2_01, 'V') && State.Nodes.TK_O2_01.V > 0
        v_lox_frac = State.Nodes.TK_O2_01.Liquid.V / State.Nodes.TK_O2_01.V;
    end
    if isfield(State.Nodes.TK_FU_01, 'Liquid') && isfield(State.Nodes.TK_FU_01, 'V') && State.Nodes.TK_FU_01.V > 0
        v_fu_frac = State.Nodes.TK_FU_01.Liquid.V / State.Nodes.TK_FU_01.V;
    end

    fill_levels = single([ ...
        State.Nodes.TK_N2.P / (4500 * psi2Pa), ...
        v_lox_frac, ...
        v_fu_frac]);

    % Assemble 100-byte packet
    pktBytes = [ ...
        hdr, ...
        typecast(pt_vals, 'uint8'), ...
        typecast(tc_vals, 'uint8'), ...
        typecast(vMask, 'uint8'), ...
        typecast(vAngles, 'uint8'), ...
        typecast(fill_levels, 'uint8') ];
end

%% =========================================================================
%% --- HELPER: Valve Mask Decoding & Encoding ---
%% =========================================================================
function LinkStates = applyValveMask(LinkStates, mask, oxThrt, fuThrt)
    bitMap = {
        0,  'SV_N2_01';  1,  'SV_N2_02';  2,  'SV_N2_03';  3,  'SV_N2_04';
        4,  'SV_N2_05';  5,  'SV_N2_06';  6,  'SV_N2_07';  7,  'SV_DART_OX';
        8,  'SV_DART_FU';9,  'BV_N2_01';  10, 'BV_N2_02';  11, 'BV_O2_01';
        12, 'BV_O2_02';  13, 'BV_02_03';  14, 'BV_FU_01';  15, 'BV_FU_03';
        16, 'BV_N2_FILL';
    };

    for k = 1:size(bitMap, 1)
        b = bitMap{k, 1};
        fn = bitMap{k, 2};
        isOpen = bitand(mask, bitshift(uint32(1), b)) > 0;
        val = 0.0;
        if isOpen, val = 1.0; end
        
        LinkStates.(fn) = val;
        if strcmp(fn, 'BV_02_03')
            LinkStates.BV_O2_03 = val;
        end
    end

    LinkStates.BV_02_04 = max(0.0, min(1.0, double(oxThrt)));
    LinkStates.BV_O2_04 = LinkStates.BV_02_04;
    LinkStates.BV_FU_04 = max(0.0, min(1.0, double(fuThrt)));
end

function vMask = buildValveMask(System)
    bitMap = {
        0,  'SV_N2_01';  1,  'SV_N2_02';  2,  'SV_N2_03';  3,  'SV_N2_04';
        4,  'SV_N2_05';  5,  'SV_N2_06';  6,  'SV_N2_07';  7,  'SV_DART_OX';
        8,  'SV_DART_FU';9,  'BV_N2_01';  10, 'BV_N2_02';  11, 'BV_O2_01';
        12, 'BV_O2_02';  13, 'BV_02_03';  14, 'BV_FU_01';  15, 'BV_FU_03';
        16, 'BV_N2_FILL';
    };

    vMask = uint32(0);
    for k = 1:size(bitMap, 1)
        b = bitMap{k, 1};
        fn = bitMap{k, 2};
        st = 0.0;
        if isfield(System.Links, fn)
            st = System.Links.(fn).State;
        elseif strcmp(fn, 'BV_02_03') && isfield(System.Links, 'BV_O2_03')
            st = System.Links.BV_O2_03.State;
        elseif isfield(System, 'Link') && isfield(System.Link, 'State') && isfield(System.Link.State, fn)
            st = System.Link.State.(fn);
        end
        if st > 0.5
            vMask = bitor(vMask, bitshift(uint32(1), b));
        end
    end
end

function cleanupSockets(s1, s2)
    try if ~isempty(s1), s1.close(); end; catch; end
    try if ~isempty(s2), s2.close(); end; catch; end
    fprintf('Sockets closed cleanly.\n');
end
