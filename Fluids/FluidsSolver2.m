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
            
            % Run Physics Engine Symbolically
            % We pass 'true' to enable smooth differentiable approximations (tanh)
            fprintf('Computing Symbolic Physics...\n');
            FluidDynamics = PhysicsEngine(X_sym, System, U_sym, true);
            
            % Compute Jacobian
            fprintf('Computing Jacobian (This may take a moment)...\n');
            FluidJacobian = jacobian(FluidDynamics, X_sym);
            
            % Generate Optimized MATLAB Files
            fprintf('Exporting "FluidDynamics.m" and "FluidJacobian.m"...\n');
            matlabFunction(FluidDynamics, 'File', './Helpers/FluidDynamics', 'Vars', {X_sym, U_sym});
            matlabFunction(FluidJacobian,  'File', './Helpers/FluidJacobian', 'Vars', {X_sym, U_sym});
            
            fprintf('Success! Files generated.\n');
            varargout{1} = true;
    end
end

%% Physics Engine
function XDOT = PhysicsEngine(X, System, U, isSym)
    %% Constants & Setup
    PLength = size(System.Nodes, 2);
    MALength = size(System.Links.Algebraic, 2);
    MDLength = size(System.Links.Dynamic, 2);
    YIdx_Start = PLength + MDLength + 1;

    % Pre-allocation (Use sym if mode is analytic, else zeros)
    if isSym
        XDOT = sym(zeros(size(X)));
        Pressure = sym(zeros(PLength, 1));
        Massflow = sym(zeros(MALength + MDLength, 1));
        RhoNode = sym(zeros(PLength, 1));
    else
        XDOT = zeros(size(X));
        Pressure = zeros(PLength, 1);
        Massflow = zeros(MALength + MDLength, 1);
        RhoNode = zeros(PLength, 1);
    end

    RhoArray = System.Constants.RhoArray;
    NumSpecies = length(RhoArray);
    LinkMap = System.LinkMap;

    %% Node Densities
    Y_All = X(YIdx_Start:end);
    Y_Matrix = reshape(Y_All, [NumSpecies, PLength])';
    
    for i = 1:PLength
        Y_Vec = Y_Matrix(i, :);
        
        % Normalize (With protection for Symbolic div-by-zero)
        % Note: For Jacobian safety, we assume Y sums to ~1 or handle it smoothly.
        % Here we strictly use the formula provided but sanitized.
        Y_Sum = sum(Y_Vec);
        if ~isSym && Y_Sum > 0
             Y_Vec = Y_Vec / Y_Sum;
        elseif isSym
             Y_Vec = Y_Vec / (Y_Sum + 1e-9); 
        end
        Y_Matrix(i, :) = Y_Vec;

        InvRho = sum(Y_Vec ./ RhoArray);
        RhoNode(i) = 1 / InvRho;
    end

    %% Dynamic Links (Pipes)
    for i = 1:MDLength
        L = System.Links.Dynamic(i);
        UpIDX = L.Up; DwIDX = L.Down;
        DeltaP = X(UpIDX) - X(DwIDX);
        FlowState = X(PLength + i);
        
        % Density Selection (Upwinding)
        % isSym -> Uses tanh blend for differentiability
        % ~isSym -> Uses hard switch for speed
        Rho = GetUpwindRho(FlowState, RhoNode(UpIDX), RhoNode(DwIDX), isSym);

        % Friction Term (Signed Square)
        % P_loss = K * m * |m|
        if isSym
            % Smooth approximation of |x|*x for derivative stability near 0
            FricTerm = FlowState * sqrt(FlowState^2 + 1e-6); 
        else
            FricTerm = abs(FlowState) * FlowState;
        end
        
        XDOT(PLength + i) = L.A/L.L * DeltaP - L.Zeta/(2*Rho*L.A*L.L) * FricTerm;
        Massflow(L.ID) = FlowState;
    end

    %% Algebraic Links (Valves)
    for i = 1:MALength
        L = System.Links.Algebraic(i);
        if strcmp(L.Type, 'Throttle')
            UpIDX = L.Up; DwIDX = L.Down;
            DeltaP = X(UpIDX) - X(DwIDX);

            % Valve Control Logic
            if isempty(U)
                Cv = L.Cv;
            else
                Cv = U(i);
            end
            
            % Upwinding based on Pressure Gradient
            Rho = GetUpwindRho(DeltaP, RhoNode(UpIDX), RhoNode(DwIDX), isSym);
            if isSym
                SqrtTerm = sqrt(abs(DeltaP)*Rho + 1e-9);
                SignTerm = tanh(100 * DeltaP);
            else
                SqrtTerm = sqrt(abs(DeltaP * Rho));
                SignTerm = sign(DeltaP);
            end
            
            Massflow(L.ID) = 2.402e-5 * Cv * SqrtTerm * SignTerm;
        end
    end

    %% 4. Pressure & Species Transport
    XDOT_Species = XDOT(YIdx_Start:end);
    XDOT_Species_Mat = zeros(PLength, NumSpecies);
    if isSym, XDOT_Species_Mat = sym(XDOT_Species_Mat); end
    for i = 1:PLength
        Node = System.Nodes(i);
        if ~Node.Fixed
            Mdot_IN = Massflow(Node.LinksIN);
            Mdot_OUT = Massflow(Node.LinksOUT);
            
            % Pressure ODE
            % dP = C^2/V * (Min - Mout)
            XDOT(i) = System.Constants.C^2 / Node.V * (sum(Mdot_IN) - sum(Mdot_OUT));
            Pressure(Node.ID) = X(i);

            % Species Transport
            % dY = (FluxIn - Y*TotalIn) / Mass
            Mass = Node.V * RhoNode(i);
            FluxIN = zeros(1, NumSpecies); 
            if isSym, FluxIN = sym(FluxIN); end
            TotalMDOT_IN = 0;

            % Process Inflows (Neighbors feeding THIS node)
            % Logic: If LinkIN has Flow > 0, it contributes.
            for k = 1:length(Node.LinksIN)
                lIdx = Node.LinksIN(k);
                flow = Massflow(lIdx);
                neighbor = LinkMap(lIdx, 1);
                
                factor = GetPosFactor(flow, isSym); % 1 if flow>0, 0 else
                FluxIN = FluxIN + (flow * factor) * Y_Matrix(neighbor, :);
                TotalMDOT_IN = TotalMDOT_IN + (flow * factor);
            end

            % Process Backflows (Neighbors feeding THIS node via OUT link)
            % Logic: If LinkOUT has Flow < 0, it contributes (reverse flow).
            for k = 1:length(Node.LinksOUT)
                lIdx = Node.LinksOUT(k);
                flow = Massflow(lIdx);
                neighbor = LinkMap(lIdx, 2); % Downstream is neighbor
                
                % We want the magnitude of negative flow: -flow
                % factor needs to be 1 if flow < 0
                factor = GetPosFactor(-flow, isSym); 
                FluxIN = FluxIN + (-flow * factor) * Y_Matrix(neighbor, :);
                TotalMDOT_IN = TotalMDOT_IN + (-flow * factor);
            end

            XDOT_Species_Mat(i, :) = (FluxIN - TotalMDOT_IN * Y_Matrix(i, :)) / Mass;
        else
            XDOT(i) = 0;
            XDOT_Species_Mat(i, :) = 0;
        end
    end
    
    XDOT(YIdx_Start:end) = reshape(XDOT_Species_Mat', [], 1);
end
function Rho = GetUpwindRho(Driver, RhoUp, RhoDw, isSym)
    if isSym
        k = 50; 
        Blend = 0.5 + 0.5 * tanh(k * Driver);
        Rho = RhoUp * Blend + RhoDw * (1 - Blend);
    else
        if Driver >= 0
            Rho = RhoUp;
        else
            Rho = RhoDw;
        end
    end
end
function factor = GetPosFactor(Val, isSym)
    if isSym
        k = 50;
        factor = 0.5 + 0.5 * tanh(k * Val);
    else
        factor = double(Val > 0);
    end
end