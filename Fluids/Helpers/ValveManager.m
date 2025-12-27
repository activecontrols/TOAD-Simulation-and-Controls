classdef ValveManager < handle
    % ValveManager handles valve sequencing, timing, and basic dynamics

    properties
        ID_Map
        ValveParams
        NumValves

        CurrentU
        TargetU
        NextEventIdx

        SequenceTime
        SequenceID
        SequenceCv
    end
    methods
        function obj = ValveManager(System)
            % Constructor
            obj.NumValves = length(System.Links.Algebraic);
            obj.ID_Map = containers.Map('KeyType', 'char', 'ValueType', 'double');
            obj.ValveParams = struct('Tau', {});
            obj.CurrentU = zeros(obj.NumValves, 1);
            obj.TargetU = zeros(obj.NumValves, 1);

            for i = 1:obj.NumValves
                L = System.Links.Algebraic(i);

                % Define time constant
                obj.ValveParams(i).Tau = 0.015;
                obj.ID_Map(L.Name) = i;
            end

            obj.NextEventIdx = 1;
        end
        function SetTau(obj, ValveName, TauValue)
            if isKey(obj.ID_Map, ValveName)
                idx = obj.ID_Map(ValveName);
                obj.ValveParams(idx).Tau = TauValue;
            end
        end
        function LoadSequence(obj, SequenceData)
            N = size(SequenceData, 1);
            Times = zeros(N, 1); IDs = zeros(N, 1); Cvs = zeros(N, 1);

            for i = 1:N
                Times(i) = SequenceData{i, 1};
                Name = SequenceData{i, 2};
                Value = SequenceData{i, 3};

                if isKey(obj.ID_Map, Name)
                    IDs(i) = obj.ID_Map(Name);
                    Cvs(i) = Value;
                end
            end

            [obj.SequenceTime, sortIdx] = sort(Times);
            obj.SequenceID = IDs(sortIdx);
            obj.SequenceCv = Cvs(sortIdx);
            obj.NextEventIdx = 1;
        end
        function SetCommand(obj, ValveName, CvValue)
            if isKey(obj.ID_Map, ValveName)
                idx = obj.ID_Map(ValveName);
                obj.TargetU(idx) = CvValue;
            end
        end
        function U = Step(obj, t, dt)
            % Process events
            while obj.NextEventIdx <= length(obj.SequenceTime) && t >= obj.SequenceTime(obj.NextEventIdx)
                idx = obj.SequenceID(obj.NextEventIdx);
                obj.TargetU(idx) = obj.SequenceCv(obj.NextEventIdx);
                obj.NextEventIdx = obj.NextEventIdx + 1;
            end

            % Process dynamics
            Error = obj.TargetU - obj.CurrentU;
            SnapTolerance = 1e-4;
            ActiveMask = abs(Error) > SnapTolerance;
            if any(ActiveMask)
                Idxs = find(ActiveMask);
                for k = Idxs'
                    Tau = obj.ValveParams(k).Tau;
                    Alpha = 1 - exp(-dt / Tau);
                    obj.CurrentU(k) = obj.CurrentU(k) + Alpha * Error(k);
                end
            end

            % Snap to target if below tolerance
            SnapMask = abs(obj.TargetU - obj.CurrentU) < SnapTolerance;
            obj.CurrentU(SnapMask) = obj.TargetU(SnapMask);
            U = obj.CurrentU;
        end
        function PlotSequence(obj, t_end)
            % PLOTSEQUENCE Visualizes the valve schedule in Dark Mode
            % Input: t_end (Duration in seconds)
            
            if nargin < 2, t_end = 20; end 
            
            % --- Dark Mode Settings ---
            DarkBg = [0.15 0.15 0.15];
            AxColor = [0.9 0.9 0.9];
            GridAlpha = 0.2;
            
            % Create figure
            figure('Name', 'Valve Sequence Preview', 'Color', DarkBg, 'NumberTitle', 'off');
            hold on; grid on;
            xlabel('Time [s]'); ylabel('Target Cv');
            title('Valve Auto-Sequence Preview');
            
            % Generate high-res time vector
            t_plot = linspace(0, t_end, 2000);
            Y_plot = zeros(obj.NumValves, length(t_plot));
            
            % Reconstruct state
            TempTargets = zeros(obj.NumValves, 1);
            EventPtr = 1;
            
            for k = 1:length(t_plot)
                t = t_plot(k);
                
                % Process events
                while EventPtr <= length(obj.SequenceTime) && ...
                    t >= obj.SequenceTime(EventPtr)
                    id = obj.SequenceID(EventPtr);
                    val = obj.SequenceCv(EventPtr);
                    TempTargets(id) = val;
                    EventPtr = EventPtr + 1;
                end
                Y_plot(:, k) = TempTargets;
            end
            
            % --- DYNAMIC COLOR GENERATION ---
            % 1. Identify which valves are actually used
            ActiveIndices = find(max(Y_plot, [], 2) > 0);
            NumActive = length(ActiveIndices);
            
            % 2. Generate exactly 'NumActive' distinct colors
            % 'hsv' gives the widest separation of hues.
            raw_colors = hsv(NumActive);
            
            % 3. Boost brightness/saturation for Dark Mode visibility
            % (HSV is usually bright, but let's ensure they pop)
            % No extra math needed for HSV, it defaults to max brightness.
            
            LegendEntries = {};
            
            for k = 1:NumActive
                idx = ActiveIndices(k);
                
                % Grab unique color for this valve
                C = raw_colors(k, :);
                
                % Plot
                plot(t_plot, Y_plot(idx, :), 'LineWidth', 2, 'Color', C);
                
                % Reverse Lookup Name
                Keys = keys(obj.ID_Map);
                Vals = values(obj.ID_Map);
                Name = Keys{[Vals{:}] == idx};
                LegendEntries{end+1} = Name;
            end
            
            % Apply Dark Mode Styling
            set(gca, 'Color', DarkBg, ...
                     'XColor', AxColor, ...
                     'YColor', AxColor, ...
                     'GridColor', 'w', ...
                     'GridAlpha', GridAlpha);
            
            if ~isempty(LegendEntries)
                Lgnd = legend(LegendEntries, 'Location', 'best');
                set(Lgnd, 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
            end
            
            ylim([0, max(max(Y_plot)) * 1.2]); 
        end
    end
end