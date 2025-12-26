function XDOT = PhysicsEngine(X, System, U)
    %% Constants & Setup
    PLength = size(System.Nodes, 2);
    MALength = size(System.Links.Algebraic, 2);
    MDLength = size(System.Links.Dynamic, 2);
    YIdx_Start = PLength + MDLength + 1;

    % Constants (bulk modulus assumed constant between OX and IPA)
    Beta_Liq = 1.1e9;
    RhoArray = System.Constants.RhoArray;
    NumSpecies = length(RhoArray);
    LinkMap = System.LinkMap;
    Rho_Default_Liq = mean(RhoArray(1:2));

    % Pre-allocation
    if ~isreal(X)
        XDOT = complex(zeros(size(X)));
        Massflow = complex(zeros(MALength + MDLength, 1));
    else
        XDOT = zeros(size(X));
        Massflow = zeros(MALength + MDLength, 1);
    end

    %% Node Densities (Vectorized)
    % Initialize Arrays
    Nodes      = System.Nodes;
    R_Vec      = [Nodes.R]';   
    T_Vec      = [Nodes.Temp]';    
    V_Vec      = [Nodes.V]';      
    IsComb_Vec = [Nodes.IsCombustor]';

    % Process Pressures
    P_Raw    = X(1:PLength);
    Pressure = sqrt(P_Raw.^2 + 1e-6);

    % Process Mass Fractions
    Y_All = X(YIdx_Start:end);
    Y_Matrix = reshape(Y_All, [NumSpecies, PLength])';

    % Normalize
    Y_Sum    = sum(Y_Matrix, 2);
    Y_Matrix = Y_Matrix ./ (Y_Sum + 1e-9);

    % Gas Density via. Ideal Gas Law
    Rho_Gas_Vec = Pressure ./ (R_Vec .* T_Vec);

    % Extract Columns
    Y_OX = Y_Matrix(:, 1); 
    Y_FU = Y_Matrix(:, 2); 
    Y_N2 = Y_Matrix(:, 3);

    % Mixture Density
    InvRhoLiq = (Y_OX ./ RhoArray(1)) + (Y_FU ./ RhoArray(2));
    Rho_Raw   = (Y_OX + Y_FU) ./ (InvRhoLiq + 1e-12);

    % Smooth selection depending on gas composition
    IsPureGas = 0.5 + 0.5 * tanh(1e6 * (1e-9 - InvRhoLiq));
    Rho_Liq_Vec = (1 - IsPureGas) .* Rho_Raw + IsPureGas .* Rho_Default_Liq;

    % Proportion of gas (Alpha)
    Term_LiqVol = Rho_Gas_Vec ./ Rho_Liq_Vec;
    Alpha_Vec = Y_N2 ./ (Y_N2 + Term_LiqVol .* (1 - Y_N2) + 1e-9);

    % Override to pure gas if combustor
    Alpha_Vec(IsComb_Vec) = 1.0;

    % Final Node Properties
    RhoNode    = Alpha_Vec .* Rho_Gas_Vec + (1 - Alpha_Vec) .* Rho_Liq_Vec;
    Vol_Ullage = V_Vec .* Alpha_Vec;

    %% Dynamic Links (Vectorized)
    DynLinks = System.Links.Dynamic;

    if ~isempty(DynLinks)
        % Initialize Arrays
        Up_Vec   = [DynLinks.Up]';
        Down_Vec = [DynLinks.Down]';
        A_Vec    = [DynLinks.A]';
        L_Vec    = [DynLinks.L]';
        Zeta_Vec = [DynLinks.Zeta]';
        ID_Vec   = [DynLinks.ID]';

        % Pick out the Massflow states from the state vector X
        LinkIdx_Start = PLength + 1;
        LinkIdx_End   = PLength + MDLength;
        FlowState_Vec = X(LinkIdx_Start:LinkIdx_End);

        % Flow direction vector
        Dir_Vec = 0.5 + 0.5 * tanh(100 * FlowState_Vec);

        % Upstream node proprieties
        P_Up_Vec = Pressure(Up_Vec) .* Dir_Vec + Pressure(Down_Vec) .* (1 - Dir_Vec);
        T_Up_Vec = T_Vec(Up_Vec)    .* Dir_Vec + T_Vec(Down_Vec)    .* (1 - Dir_Vec);
        R_Up_Vec = R_Vec(Up_Vec)    .* Dir_Vec + R_Vec(Down_Vec)    .* (1 - Dir_Vec);

        % Gas Density per Link
        Rho_Gas_Vec = P_Up_Vec ./ (R_Up_Vec .* T_Up_Vec);

        % Smooth Selection for link composition
        Y_Up_Raw = Y_Matrix(Up_Vec, :)   .* Dir_Vec +...
                   Y_Matrix(Down_Vec, :) .* (1 - Dir_Vec);

        % Apply stratification and normalize
        LiqSum_Vec = sum(Y_Up_Raw(:, 1:2), 2);
        isLiquid      = 0.5 + 0.5 * tanh(100 * LiqSum_Vec);
        Y_Strat    = Y_Up_Raw;
        Y_Strat(:, 3) = Y_Strat(:, 3) .* (1 - isLiquid);

        % Mixture density
        % We build a density matrix and replace the third column as
        % nessesary from the above process.
        Rho_Mat_Links = repmat(RhoArray, MDLength, 1);
        Rho_Mat_Links(:, 3) = Rho_Gas_Vec;
        Rho_Mix_Vec = 1 ./ sum(Y_Strat ./ (Rho_Mat_Links + 1e-9), 2);

        % Momentum Equation
        DeltaP_Vec = X(Up_Vec) - X(Down_Vec);
        FricTerm_Vec = FlowState_Vec .* sqrt(FlowState_Vec.^2 + 1e-6);

        % ODE
        XDOT(LinkIdx_Start:LinkIdx_End) = (A_Vec ./ L_Vec) .* ...
            (DeltaP_Vec - (Zeta_Vec ./ (2 * Rho_Mix_Vec .* A_Vec.^2)) .* FricTerm_Vec);

        % Store Massflows for global use
        Massflow(ID_Vec) = FlowState_Vec;
    end

    %% Algebraic Links (Vectorized)
    AlgLinks = System.Links.Algebraic;

    if ~isempty(AlgLinks)
        % Initialize Vectors
        Up_Vec   = [AlgLinks.Up]';
        Down_Vec = [AlgLinks.Down]';
        ID_Vec   = [AlgLinks.ID]';
        Cv_Vec   = [AlgLinks.Cv]';
        A_Vec    = [AlgLinks.A]';

        % Logical Mask for Valves
        Type_Cell  = {AlgLinks.Type};
        isThrottle = strcmp(Type_Cell, 'Throttle')';

        % Control Input handling
        if ~isempty(U)
            Cv_Vec(isThrottle) = U(isThrottle); 
        end

        % Both Valves and Orifice Plates follow the same form.
        Coeff_Vec = isThrottle * 2.402e-5 + (~isThrottle) .* A_Vec;
        Sqrt_Mult = isThrottle * 1.0      + (~isThrottle) * 2.0;

        % Flow Direction & Upstream Properties
        DeltaP_Vec = X(Up_Vec) - X(Down_Vec);
        Dir_Vec    = 0.5 + 0.5 * tanh(1e3 * DeltaP_Vec);

        % Blending Upstream Properties
        P_Up_Vec = Pressure(Up_Vec) .* Dir_Vec + Pressure(Down_Vec) .* (1 - Dir_Vec);
        T_Up_Vec = T_Vec(Up_Vec)    .* Dir_Vec + T_Vec(Down_Vec)    .* (1 - Dir_Vec);
        R_Up_Vec = R_Vec(Up_Vec)    .* Dir_Vec + R_Vec(Down_Vec)    .* (1 - Dir_Vec);
        Rho_Gas_Flow = P_Up_Vec ./ (R_Up_Vec .* T_Up_Vec);
        Y_Up_Raw = Y_Matrix(Up_Vec, :)   .* Dir_Vec + ...
                   Y_Matrix(Down_Vec, :) .* (1 - Dir_Vec);

        % Stratification
        LiqSum_Vec = sum(Y_Up_Raw(:, 1:2), 2);
        isLiquid   = 0.5 + 0.5 * tanh(10 * (LiqSum_Vec - 0.1)); 
        Y_Flow     = Y_Up_Raw;
        Y_Flow(:, 3) = Y_Flow(:, 3) .* (1 - isLiquid);
        Y_Flow     = Y_Flow ./ (sum(Y_Flow, 2) + 1e-9);

        % Identify combustor sources
        IsSrcComb = (IsComb_Vec(Up_Vec) .* Dir_Vec) + ...
                    (IsComb_Vec(Down_Vec) .* (1 - Dir_Vec));

        % Mixture Density
        Rho_Mat_Links = repmat(RhoArray, MALength, 1);
        Rho_Mat_Links(:, 3) = Rho_Gas_Flow;
        Rho_Mix_Vec = 1 ./ sum(Y_Flow ./ (Rho_Mat_Links + 1e-9), 2);

        % Final Upstream Density
        Rho_Up_Final = (IsSrcComb .* Rho_Gas_Flow) + ((1 - IsSrcComb) .* Rho_Mix_Vec);

        % Choked flow handling for combustors and choked flow check
        IsLiqFlow = isLiquid .* (1 - IsSrcComb);
        NodesGamma = [System.Nodes.Gamma]';
        G_Node     = NodesGamma(Up_Vec) .* Dir_Vec + NodesGamma(Down_Vec) .* (1 - Dir_Vec);
        G_Up       = IsLiqFlow * 100 + (1 - IsLiqFlow) .* G_Node;

        % Critical Pressure Ratio
        Rc = (2 ./ (G_Up + 1)) .^ (G_Up ./ (G_Up - 1));
        DP_Choked = P_Up_Vec .* (1 - Rc);
        DP_Raw    = sqrt(DeltaP_Vec.^2 + 1e-9);

        % Smooth Min 
        DP_Eff = 0.5 * (DP_Raw + DP_Choked - sqrt((DP_Raw - DP_Choked).^2 + 1e-4));

        % Final Massflow Calculation
        SignTerm = tanh(1e3 * DeltaP_Vec);
        SqrtTerm = sqrt(Sqrt_Mult .* Rho_Up_Final .* DP_Eff + 1e-9);
        Massflow(ID_Vec) = Coeff_Vec .* Cv_Vec .* SqrtTerm .* SignTerm;
    end

    %% Conservation of Mass and Energy (Vectorized)
    % Gather link data
    Link_Up   = LinkMap(:, 1);
    Link_Down = LinkMap(:, 2);
    Link_Flow = Massflow;

    % Determine flow direction (1.0 if Up -> Down, 0.0 if Down -> Up)
    Link_Dir = 0.5 + 0.5 * tanh(100 * Link_Flow);

    % Get source compositions and perform stratification on upstream and
    % downstream nodes
    % Upstream candidates
        Y_Up_Raw   = Y_Matrix(Link_Up, :);
        LiqSum_Up  = sum(Y_Up_Raw(:, 1:2), 2);
        isLiquid_Up       = 0.5 + 0.5 * tanh(100 * LiqSum_Up);
        Y_Up_Strat = Y_Up_Raw;
        Y_Up_Strat(:, 3) = Y_Up_Strat(:, 3) .* (1 - isLiquid_Up);
        Y_Up_Strat = Y_Up_Strat ./ (sum(Y_Up_Strat, 2) + 1e-9);
    % Downstream candidates
        Y_Down_Raw   = Y_Matrix(Link_Down, :);
        LiqSum_Down  = sum(Y_Down_Raw(:, 1:2), 2);
        isLiquid_Down       = 0.5 + 0.5 * tanh(100 * LiqSum_Down);
        Y_Down_Strat = Y_Down_Raw;
        Y_Down_Strat(:, 3) = Y_Down_Strat(:, 3) .* (1 - isLiquid_Down);
        Y_Down_Strat = Y_Down_Strat ./ (sum(Y_Down_Strat, 2) + 1e-9);

    % Select composition based on flow direction
    Y_Trans = Y_Up_Strat .* Link_Dir + Y_Down_Strat .* (1 - Link_Dir);

    % Calculate Fluxes crossing each link
    Flux_Total = Link_Flow;
    Flux_Gas   = Link_Flow .* Y_Trans(:, 3);
    Flux_Liq   = Link_Flow .* (1 - Y_Trans(:, 3));
    Flux_Spec  = Link_Flow .* Y_Trans;

    % Accumulate flows to nodes
    % Net massflows
    Net_Total = accumarray(Link_Down, Flux_Total, [PLength, 1]) - ...
                accumarray(Link_Up,   Flux_Total, [PLength, 1]);

    Net_Gas   = accumarray(Link_Down, Flux_Gas, [PLength, 1]) - ...
                accumarray(Link_Up,   Flux_Gas, [PLength, 1]);
                
    Net_Liq   = accumarray(Link_Down, Flux_Liq, [PLength, 1]) - ...
                accumarray(Link_Up,   Flux_Liq, [PLength, 1]);

    % Net species flux
    Net_Spec = zeros(PLength, 3);
    Net_Spec(:,1) = accumarray(Link_Down, Flux_Spec(:,1), [PLength, 1]) - ...
                    accumarray(Link_Up, Flux_Spec(:,1), [PLength, 1]);

    Net_Spec(:,2) = accumarray(Link_Down, Flux_Spec(:,2), [PLength, 1]) - ...
                    accumarray(Link_Up, Flux_Spec(:,2), [PLength, 1]);

    Net_Spec(:,3) = accumarray(Link_Down, Flux_Spec(:,3), [PLength, 1]) - ...
                    accumarray(Link_Up, Flux_Spec(:,3), [PLength, 1]);

    % Calculate Pressure and Mass Fraction Derivatives
    Gas_For_P = Net_Gas;
    Liq_For_P = Net_Liq;

    % Handle combustor logic (all incoming mass == gas)
    Gas_For_P(IsComb_Vec) = Net_Gas(IsComb_Vec) + Net_Liq(IsComb_Vec);
    Liq_For_P(IsComb_Vec) = 0;

    % Volume Safety (combustors are pure gas volume)
    V_Safe = Vol_Ullage + 1e-6;
    V_Safe(IsComb_Vec) = V_Vec(IsComb_Vec);

    % Pressure derivatives
    dP_GasTerm = Gas_For_P .* R_Vec .* T_Vec;
    dP_LiqTerm = Pressure .* (Liq_For_P ./ max(Rho_Liq_Vec, 100));
    dP_GasMode   = (dP_GasTerm + dP_LiqTerm) ./ V_Safe;
    dP_FluidMode = (Beta_Liq ./ (V_Vec .* max(Rho_Liq_Vec, 100))) .* Net_Total;

    % Blend Modes
    BlendFactor = 0.5 + 0.5 * tanh(1e3 * (Alpha_Vec - 0.002));
    dP_Final    = BlendFactor .* dP_GasMode + (1 - BlendFactor) .* dP_FluidMode;

    % Species derivatives
    Mass_Node = V_Vec .* RhoNode;
    dY_Num    = Net_Spec - Net_Total .* Y_Matrix;
    dY_dt     = dY_Num ./ (Mass_Node + 1e-9);

    % Fixed (boundary) node constraints
    Fixed_Vec = [System.Nodes.Fixed]';
    dP_Final(Fixed_Vec == 1) = 0;
    dY_dt(Fixed_Vec == 1, :) = 0;

    % Reshape matrix and store derivatives
    XDOT(1:PLength) = dP_Final;
    XDOT(YIdx_Start:end) = reshape(dY_dt', [], 1);
end








    