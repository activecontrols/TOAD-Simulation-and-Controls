function varargout = FluidsSolver2(X, System, Mode)
% FLUIDSSOLVER Unified Physics Engine
%   XDOT = FluidsSolver(X, System, 'Numerical') -> Returns vector XDOT
%   FluidsSolver([], System, 'Analytic')        -> Generates Sym_Jac.m & Sym_ODE.m

    if nargin < 3
        Mode = 'Numerical';
    end

    switch Mode
        case 'Numerical'
            %% NUMERICAL MODE (Standard)
            % Just call the physics engine with doubles
            XDOT = PhysicsEngine(X, System, [], false);
            varargout{1} = XDOT;
            
        case 'Analytic'
            %% Analytic Mode (Generates Jacobian FCN for Implicit Solver)
            fprintf('Initializing Symbolic Generation...\n');
            
            % Setup Symbolic Vector
            PLength = size(System.Nodes, 2);
            MDLength = size(System.Links.Dynamic, 2);
            MALength = size(System.Links.Algebraic, 2);
            NumSpecies = length(System.Constants.RhoArray);
            
            TotalStates = PLength + MDLength + PLength * NumSpecies;
            X_sym = sym('X', [TotalStates, 1]); 
            U_sym = sym('U', [MALength, 1]); %[Cv1; Cv2; Cv2; Cv3...];
            assume(X_sym, 'real');
            assume(U_sym, 'real');
            assume(X_sym(1:PLength), 'positive');
            
            % Run Physics Engine Symbolically
            % We pass 'true' to enable smooth differentiable approximations (tanh)
            fprintf('Computing Symbolic Physics...\n');
            FluidDynamics = PhysicsEngine(X_sym, System, U_sym, true);
            
            % Generate Optimized MATLAB Files
            fprintf('Exporting "FluidDynamics.m and FluidJacobian.m"...\n');
            matlabFunction(FluidDynamics, 'File', './Fluids/Helpers/FluidDynamics', 'Vars', {X_sym, U_sym});
            fprintf('Success! Files generated.\n');
            
            varargout{1} = true;
    end
end
