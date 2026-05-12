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
                else
                    error('ValveManager: Valve "%s" in AutoSequence does not match any System Links!', Name);
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
            % PLOTSEQUENCE Visualizes the valve schedule as a digital timing diagram
            % Input: t_end (Duration in seconds)
            
            if nargin < 2, t_end = 20; end 
            
            % --- Dark Mode Settings ---
            DarkBg = [0.15 0.15 0.15];
            AxColor = [0.9 0.9 0.9];
            GridAlpha = 0.2;
            
            % Create figure
            figure('Name', 'Valve Sequence Timing Diagram', 'Color', DarkBg, 'NumberTitle', 'off');
            hold on; grid on;
            xlabel('Time [s]'); 
            title('Valve Auto-Sequence Timing Diagram', 'Color', 'w');
            
            % Generate high-res time vector
            t_plot = linspace(0, t_end, 5000);
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
            
            % 1. Identify which valves are actually used
            ActiveIndices = find(max(Y_plot, [], 2) > 0);
            NumActive = length(ActiveIndices);
            
            if NumActive == 0
                disp('No active valves in sequence to plot.');
                return;
            end
            
            % 2. Generate exactly 'NumActive' distinct colors
            raw_colors = hsv(NumActive);
            
            % Setup arrays for Y-axis labels
            YTickVals = zeros(NumActive, 1);
            YTickNames = cell(NumActive, 1);
            
            Keys = keys(obj.ID_Map);
            Vals = values(obj.ID_Map);
            
            % 3. Plot each valve in its own horizontal channel
            for k = 1:NumActive
                idx = ActiveIndices(k);
                
                % Convert Cv target to binary state (1 if > 0, 0 otherwise)
                BinaryState = double(Y_plot(idx, :) > 1e-5);
                
                % Calculate vertical offset to stack channels
                % We stack top-to-bottom so the first valve is at the top
                BaseY = (NumActive - k) * 1.5; 
                Height = 1.0;
                
                ChannelY = BaseY + (BinaryState * Height);
                
                % Grab unique color for this valve
                C = raw_colors(k, :);
                
                % Plot the square wave
                plot(t_plot, ChannelY, 'LineWidth', 2, 'Color', C);
                
                % Add a subtle dashed baseline for this channel
                yline(BaseY, 'Color', [0.6 0.6 0.6], 'LineWidth', 1, 'LineStyle', '-.');
                
                % Reverse Lookup Name
                Name = Keys{[Vals{:}] == idx};
                
                % Store tick info (Center the text slightly above the baseline)
                YTickVals(k) = BaseY + (Height / 2);
                YTickNames{k} = Name;
            end
            
            % Apply Dark Mode Styling and custom Y-Axis
            set(gca, 'Color', DarkBg, ...
                     'XColor', AxColor, ...
                     'YColor', AxColor, ...
                     'GridColor', 'w', ...
                     'GridAlpha', GridAlpha);
                 
            % Apply the names to the Y-axis
            [YTickVals, sortIdx] = sort(YTickVals);
            YTickNames = YTickNames(sortIdx);
            yticks(YTickVals);
            yticklabels(YTickNames);
            
            % Remove y-axis grid lines so they don't clutter the channels
            ax = gca;
            ax.YGrid = 'off';
            
            % Set limits to tightly bound the channels
            ylim([-0.5, NumActive * 1.5]); 
            xlim([0, t_end]);
        end
    end
end