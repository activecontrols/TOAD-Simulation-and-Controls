classdef TrajectoryProgressionLogger < handle
    %% TRAJECTORYPROGRESSIONLOGGER  Progressive trajectory build and molding recorder.
    % Standalone companion logger for TrajectoryOptimizerHybrid.
    % Records Stage-1 QP splines, differential flatness guesses, and Stage-2 NLP iterates,
    % preparing ghost-linked frame datasets for AnimateTrajectoryProgression.
    %
    % Authors: PSP Active Controls & Antigravity

    properties
        opt                             % Reference to TrajectoryOptimizerHybrid
        Progression struct = struct()   % Progression container
        N double = 60                   % Control intervals
        Sx double                       % State scaling
        Su double                       % Control scaling
    end

    methods
        function obj = TrajectoryProgressionLogger(opt)
            if nargin < 1 || isempty(opt)
                return;
            end
            obj.opt = opt;
            obj.N = opt.N;
            obj.Sx = opt.Sx;
            obj.Su = opt.Su;
            obj.initContainer();
        end

        function initContainer(obj)
            %% INITCONTAINER  Initialize progression struct and mission metadata.
            opt = obj.opt;
            obj.Progression = struct();
            obj.Progression.meta = struct( ...
                'Vehicle', opt.Vehicle, ...
                'VehicleName', opt.getVehicleName(), ...
                'Maneuver', opt.Maneuver, ...
                'Version', opt.Version, ...
                'r0', opt.r0, ...
                'r_f', opt.r_f, ...
                'N', opt.N, ...
                'GlideslopeAngle', opt.GlideslopeAngle, ...
                'TakeoffClearanceAlt', opt.getTakeoffClearanceAlt(), ...
                'LandingFlareAlt', opt.getLandingFlareAlt(), ...
                'FunnelCurvature', opt.FunnelCurvature, ...
                'MaxThrust', opt.constants.MaxThrust, ...
                'm_dry', opt.constants.m_dry, ...
                'g', opt.constants.g, ...
                'MaxGimbalAngle', opt.MaxGimbalAngle);
            obj.Progression.frames = {};
        end

        function recordQP(obj, stageName, QPTraj, Fix)
            %% RECORDQP  Record polynomial QP intermediate spline geometry.
            if isempty(QPTraj), return; end
            Nn = obj.N + 1;
            Ttot = sum(QPTraj.T);
            tv = linspace(0, Ttot, Nn);
            rs = QPTraj.Eval(tv, 0);
            vs = QPTraj.Eval(tv, 1);
            as = QPTraj.Eval(tv, 2);

            g_vec = [0; 0; obj.opt.constants.g];
            m_dry = obj.opt.constants.m_dry;
            qs = zeros(4, Nn);
            for k = 1:Nn
                F_est = m_dry * (as(:,k) + g_vec);
                qs(:,k) = vectorToQuat(F_est);
            end

            Us = zeros(4, obj.N);
            Us(3, :) = m_dry * obj.opt.constants.g;

            keyPos = [];
            if nargin >= 4 && ~isempty(Fix)
                keyPos = squeeze(Fix(:, 1, :));
            end

            frame = obj.makeFrame(stageName, "Stage 1 (QP)", -1, tv, rs, vs, qs, ...
                                  zeros(3, Nn), Us(1:2, :), Us(3, :), Us(4, :), Ttot, keyPos);
            obj.appendFrame(frame);
        end

        function recordGuess(obj, stageName, guess)
            %% RECORDGUESS  Record Differential Flatness initial guess.
            if isempty(guess) || ~isfield(guess, 'X'), return; end
            Xg = guess.X; Ug = guess.U;
            frame = obj.makeFrame(stageName, "Stage 1 (Flatness Inversion)", 0, ...
                                  guess.Time, Xg(5:7, :), Xg(8:10, :), Xg(1:4, :), ...
                                  Xg(11:13, :), Ug(1:2, :), Ug(3, :), Ug(4, :), guess.T_total, []);
            obj.appendFrame(frame);
        end

        function recordCollocationIter(obj, opti, Xhat, Uhat, T_total, iter)
            %% RECORDCOLLOCATIONITER  CasADi opti.callback hook capturing NLP iterates.
            if iter <= 4
                record = true;
            elseif iter <= 16 && mod(iter, 2) == 0
                record = true;
            elseif iter <= 40 && mod(iter, 4) == 0
                record = true;
            elseif mod(iter, 10) == 0
                record = true;
            else
                record = false;
            end
            if ~record, return; end

            try
                X_val = obj.Sx .* opti.debug.value(Xhat);
                U_val = obj.Su .* opti.debug.value(Uhat);
                T_val = opti.debug.value(T_total);
                if iter == 0
                    stg = "Stage 2: Collocation Iter 0 (Initial Guess Seed)";
                else
                    stg = sprintf("Stage 2: Collocation Iter %d", iter);
                end
                obj.recordCollocationFrame(stg, X_val, U_val, T_val, iter);
            catch
            end
        end

        function recordCollocationFrame(obj, stageName, X_val, U_val, T_val, iter)
            %% RECORDCOLLOCATIONFRAME  Package and append a direct collocation iterate.
            N_nodes = obj.N;
            dt = T_val / N_nodes;
            tv = [0, cumsum(repmat(dt, 1, N_nodes))];

            frame = obj.makeFrame(stageName, "Stage 2 (NLP Collocation)", iter, tv, ...
                                  X_val(5:7, :), X_val(8:10, :), X_val(1:4, :), ...
                                  X_val(11:13, :), U_val(1:2, :), U_val(3, :), U_val(4, :), T_val, []);
            obj.appendFrame(frame);
        end

        function animOut = animate(obj, varargin)
            %% ANIMATE  Forward progression data to AnimateTrajectoryProgression.
            if isempty(obj.Progression) || ~isfield(obj.Progression, 'frames') || isempty(obj.Progression.frames)
                error('TrajectoryProgressionLogger:NoProgression', ...
                    'No progression data recorded. Run solve() with RecordProgression = true.');
            end
            animOut = AnimateTrajectoryProgression(obj.Progression, varargin{:});
        end
    end

    methods (Access = private)
        function frame = makeFrame(~, stageName, phase, iter, tv, rs, vs, qs, ws, gimb, thr, roll, Ttot, keyPos)
            %% MAKEFRAME  Construct standardized progression frame struct.
            frame = struct();
            frame.StageName = string(stageName);
            frame.Phase = string(phase);
            frame.Iter = iter;
            frame.Time = tv;
            frame.Pos = rs;
            frame.Vel = vs;
            frame.Quat = qs;
            frame.AngRate = ws;
            frame.Gimbal = gimb;
            frame.Thrust = thr;
            frame.RollCmd = roll;
            frame.T_total = Ttot;
            dr = diff(rs, 1, 2);
            frame.L_path = sum(sqrt(sum(dr.^2, 1)));
            frame.MaxGimbalDeg = rad2deg(max(sqrt(gimb(1,:).^2 + gimb(2,:).^2)));
            frame.MaxThrust = max(thr);
            frame.Keyframes = keyPos;
        end

        function appendFrame(obj, frame)
            %% APPENDFRAME  Link previous ghost trajectory and append to container.
            if isempty(obj.Progression) || ~isfield(obj.Progression, 'frames')
                obj.initContainer();
            end
            nF = numel(obj.Progression.frames);
            if nF > 0
                prev = obj.Progression.frames{nF};
                frame.GhostPos = prev.Pos;
                frame.GhostQuat = prev.Quat;
                frame.GhostGimbal = prev.Gimbal;
                frame.GhostThrust = prev.Thrust;
                frame.GhostTime = prev.Time;
            else
                frame.GhostPos = frame.Pos;
                frame.GhostQuat = frame.Quat;
                frame.GhostGimbal = frame.Gimbal;
                frame.GhostThrust = frame.Thrust;
                frame.GhostTime = frame.Time;
            end
            frame.FrameIndex = nF + 1;
            obj.Progression.frames{nF + 1} = frame;
        end
    end
end
