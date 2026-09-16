function LiveSim_Runner(varargin)
% LIVESIM_RUNNER Real-Time Hardware-in-the-Loop Simulation Runner for TOAD
% Interfaces the MATLAB Lumped Parameter model with the Ground Control UI via UDP.
%
% Ports:
%   - Telemetry Stream (Sim -> UI): UDP 127.0.0.1:9000 (100-byte EC_FMT packet)
%   - Valve Command (UI -> Sim):   UDP 127.0.0.1:9001 (12-byte command packet)
%
% Real-Time Pacing: 100 Hz (dt = 0.010 s) synchronized with wall-clock time.

    fprintf('=====================================================\n');
    fprintf('   TOAD Lumped Parameter Real-Time Simulation Runner\n');
    fprintf('=====================================================\n');

    % Setup paths
    baseDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(baseDir, '..', 'Core'));
    addpath(fullfile(baseDir, '..', 'Config'));
    addpath(fullfile(baseDir, '..', 'Thermal'));
    addpath(fullfile(baseDir, '..', 'Fluid Properties'));
    addpath(baseDir);

    % Build Flight-Ready Unified System
    fprintf('Initializing Flight-Ready System Topology...\n');
    [System, State] = BuildTOADSystem_Merged();
    dt = 0.010; % 100 Hz
    psi2Pa = 6894.757;

    % Setup Zero-Toolbox Java UDP Sockets
    fprintf('Configuring UDP Sockets (Telemetry: 9000, Commands: 9001)...\n');
    sendSocket = java.net.DatagramSocket();
    destAddr   = java.net.InetAddress.getByName('127.0.0.1');
    destPort   = 9000;

    cmdChannel = [];
    cmdBuf = [];
    try
        cmdChannel = java.nio.channels.DatagramChannel.open();
        cmdChannel.configureBlocking(false);
        cmdChannel.socket().bind(java.net.InetSocketAddress(9001));
        cmdBuf = java.nio.ByteBuffer.allocate(64);
        fprintf('Command receiver listening on UDP port 9001 (non-blocking NIO channel).\n');
    catch ME
        warning('Could not bind UDP port 9001 (%s). Running in broadcast-only mode.', ME.message);
    end

    % Cleanup on exit
    cleanupObj = onCleanup(@() cleanupSockets(sendSocket, cmdChannel));

    % Initial Link States
    LinkStates = struct();
    linkNames = fieldnames(System.Links);
    for k = 1:length(linkNames)
        fn = linkNames{k};
        LinkStates.(fn) = System.Links.(fn).State;
    end

    fprintf('Ready! Starting 100 Hz Real-Time Execution Loop. Press Ctrl+C to stop.\n');
    tStart = tic;
    simTime = 0.0;
    stepCount = 0;

    try
        while true
            stepCount = stepCount + 1;

            %% 1. Non-blocking Receive: Drain UI Valve Commands
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

            %% 2. Step Physics Simulation
            [State, FlowRates, System] = StepSimulation_Live(State, dt, System, LinkStates);
            simTime = simTime + dt;

            %% 3. Pack 100-byte EC_FMT Telemetry Packet
            % Offset 0x00: uint8 CRC (1 byte) + 3 bytes padding
            hdr = uint8([mod(stepCount, 256), 0, 0, 0]);

            % Offset 0x04: 12 floats Pressures in psia
            pt_vals = single(zeros(1, 12));
            pt_vals(1) = single(State.Nodes.TK_N2.P / psi2Pa);             % PT-N2-01 (COPV)
            pt_vals(2) = single(State.Nodes.TK_O2_01.P / psi2Pa);          % PT-O2-01 (LOX Tank)
            pt_vals(3) = single(State.Nodes.TK_FU_01.P / psi2Pa);          % PT-FU-01 (Fuel Tank)
            pt_vals(4) = single(State.Nodes.Purge_Manifold.P / psi2Pa);    % PT-N2-02 (Purge Manifold)
            pt_vals(5) = single(State.Nodes.OX_Manifold.P / psi2Pa);       % PT-O2-02 (LOX Inj Manifold)
            pt_vals(6) = single(State.Nodes.FU_Manifold.P / psi2Pa);       % PT-FU-02 (Fuel Inj Manifold)
            pt_vals(7) = single(State.Nodes.DART_Chamber.P / psi2Pa);      % PT-FU-04 (DART Igniter Pc)
            if isfield(State.Nodes, 'TK_N2_BULK')
                pt_vals(8) = single(State.Nodes.TK_N2_BULK.P / psi2Pa);     % PT-N2-BULK
            end
            pt_vals(9) = single(State.Nodes.DART_Chamber.P / psi2Pa);      % PT-DART-CHAM (Mirror)
            pt_vals(10) = single(State.Nodes.SKIPPER.P / psi2Pa);          % PT-FU-03 (SKIPPER Chamber Pc)

            % Offset 0x34: 6 floats Temperatures in Kelvin
            tc_vals = single(zeros(1, 6));
            tc_vals(1) = single(State.Nodes.TK_N2.T);                      % TC-N2-01
            tc_vals(2) = single(State.Nodes.TK_O2_01.Liquid.T);            % TC-O2-01
            tc_vals(3) = single(State.Nodes.OX_Manifold.T);                % TC-O2-02
            tc_vals(4) = single(State.Nodes.FU_Manifold.T);                % TC-FU-01
            if isfield(State.Thermal, 'Regen')
                tc_vals(5) = single(State.Thermal.Regen.T_wall);           % TC-REGEN-WALL
            end

            % Offset 0x4C: uint32 valve_state_mask
            vMask = buildValveMask(System);

            % Offset 0x50: 2 floats valve angles (throttle 0.0 to 1.0)
            ox_ang = 0.0; fu_ang = 0.0;
            if isfield(System.Links, 'BV_02_04'), ox_ang = System.Links.BV_02_04.State;
            elseif isfield(System.Links, 'BV_O2_04'), ox_ang = System.Links.BV_O2_04.State;
            end
            if isfield(System.Links, 'BV_FU_04'), fu_ang = System.Links.BV_FU_04.State; end
            vAngles = single([ox_ang, fu_ang]);

            % Offset 0x58: 3 floats fill levels (0.0 to 1.0)
            fill_levels = single([ ...
                State.Nodes.TK_N2.P / (4500 * psi2Pa), ...
                State.Nodes.TK_O2_01.Liquid.V / State.Nodes.TK_O2_01.V, ...
                State.Nodes.TK_FU_01.Liquid.V / State.Nodes.TK_FU_01.V]);

            % Send UDP datagram to UI at 20 Hz (every 5th 100 Hz physics step)
            if mod(stepCount, 5) == 0
                % Assemble 100-byte buffer
                pktBytes = [ ...
                    hdr, ...
                    typecast(pt_vals, 'uint8'), ...
                    typecast(tc_vals, 'uint8'), ...
                    typecast(vMask, 'uint8'), ...
                    typecast(vAngles, 'uint8'), ...
                    typecast(fill_levels, 'uint8') ];

                % Must use typecast(..., 'int8') to preserve bit patterns!
                % (int8() in MATLAB saturates any byte >= 128 to 127, corrupting floats)
                sendPkt = java.net.DatagramPacket(typecast(pktBytes, 'int8'), 100, destAddr, destPort);
                sendSocket.send(sendPkt);
            end

            %% 4. Real-Time Pacing
            elapsed = toc(tStart);
            % Catch-up runaway protection: clamp simTime if lagging > 50 ms
            if elapsed > simTime + 0.050
                simTime = elapsed;
            end

            sleepSec = simTime - elapsed;
            if sleepSec > 0.001
                java.lang.Thread.sleep(int64(sleepSec * 1000));
            else
                % Yield thread briefly to prevent socket starvation and UI freeze
                java.lang.Thread.yield();
            end

            if mod(stepCount, 500) == 0 % Every 5 seconds
                fprintf('t = %.1f s | COPV: %.0f psi | LOX: %.1f psi | FU: %.1f psi | Pc: %.1f psi\n', ...
                    simTime, pt_vals(1), pt_vals(2), pt_vals(3), pt_vals(7));
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

function LinkStates = applyValveMask(LinkStates, mask, oxThrt, fuThrt)
% Maps 32-bit valve mask to individual link states
    bitMap = {
        0,  'SV_N2_01';
        1,  'SV_N2_02';
        2,  'SV_N2_03';
        3,  'SV_N2_04';
        4,  'SV_N2_05';
        5,  'SV_N2_06';
        6,  'SV_N2_07';
        7,  'SV_DART_OX';
        8,  'SV_DART_FU';
        9,  'BV_N2_01';
        10, 'BV_N2_02';
        11, 'BV_O2_01';
        12, 'BV_O2_02';
        13, 'BV_02_03';
        14, 'BV_FU_01';
        15, 'BV_FU_03';
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
% Constructs 32-bit valve status mask from System.Links
    bitMap = {
        0,  'SV_N2_01';
        1,  'SV_N2_02';
        2,  'SV_N2_03';
        3,  'SV_N2_04';
        4,  'SV_N2_05';
        5,  'SV_N2_06';
        6,  'SV_N2_07';
        7,  'SV_DART_OX';
        8,  'SV_DART_FU';
        9,  'BV_N2_01';
        10, 'BV_N2_02';
        11, 'BV_O2_01';
        12, 'BV_O2_02';
        13, 'BV_02_03';
        14, 'BV_FU_01';
        15, 'BV_FU_03';
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
