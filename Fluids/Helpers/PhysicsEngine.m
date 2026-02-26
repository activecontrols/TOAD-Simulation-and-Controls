function [XDOT, Massflow] = PhysicsEngine(X, System, U)
    %% Constants & Setup
    PLength = size(System.Nodes, 2);
    MALength = size(System.Links.Algebraic, 2);
    MDLength = size(System.Links.Dynamic, 2);
    YIdx_Start = PLength + MDLength + 1;
    Combustor = isfield(System, 'Combustion') && isfield(System.Combustion, 'BurnRates');

    % Constants (bulk modulus assumed constant between OX and IPA)
    Beta_Liq = 1.1e9;
    RhoArray = System.Constants.RhoArray;
    NumSpecies = length(RhoArray);
    LinkMap = System.LinkMap;
    Rho_Default_Liq = (RhoArray(1) + RhoArray(2)) / 2;

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
    Type_Cell = {Nodes.Type}';
    IsTank_Vec = strcmp(Type_Cell, 'Tank');

    % Process Pressures
    P_Raw    = X(1:PLength);
    Pressure = sqrt(P_Raw.^2 + 1e-6);

    % Process Mass Fractions
    Y_All = X(YIdx_Start:end);
    Y_Matrix = reshape(Y_All, [NumSpecies, PLength])';

    % Normalize
    Y_Sum    = sum(Y_Matrix, 2);
    Y_Matrix = Y_Matrix ./ (Y_Sum + 1e-9);

    %% Combustor Temp Logic
    CStarEff = 0.83;
    CEA_OF = System.Combustion.CEA_OF;
    CEA_Temp = System.Combustion.CEA_Temp * CStarEff^2;

    % Current OF Ratio
    Current_OF = Y_Matrix(:, 1) ./ (Y_Matrix(:, 2) + 1e-9);
    if Combustor
        B_OX = System.Combustion.BurnRates(1, 1);
        B_FU = System.Combustion.BurnRates(1, 2);
        Current_OF(IsComb_Vec) = B_OX / (B_FU + 1e-6);
    end
    IgniterMask = 1;
    IgniterIdx = System.Combustion.IgniterID;
    if ~isempty(IgniterIdx)
        IgniterMask = U(IgniterIdx);
    end

    OF_Clamped = max(min(Current_OF, 3.3), 0.8);
    T_Dynamic = 293 * ones(size(T_Vec));
    if any(IsComb_Vec)
        OF_CombRaw = Current_OF(IsComb_Vec);
        OF_Comb = OF_Clamped(IsComb_Vec);

        T_Interp =  interp1(CEA_OF, CEA_Temp, OF_Comb, 'linear');

        MaskLean = 0.5 + 0.5 * tanh(3 * (OF_CombRaw - 0.3));
        MaskRich = 0.5 + 0.5 * tanh(3 * (4 - OF_CombRaw));
        BurnFactor = MaskLean .* MaskRich .*IgniterMask;
        T_Dynamic(IsComb_Vec) = 293 + (T_Interp - 293) .* BurnFactor;
    end
    T_Vec(IsComb_Vec) = T_Dynamic(IsComb_Vec);

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
    if ~Combustor
        Alpha_Vec(IsComb_Vec) = 1.0;
    end

    % Final Node Properties
    RhoNode    = Alpha_Vec .* Rho_Gas_Vec + (1 - Alpha_Vec) .* Rho_Liq_Vec;

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
        Dir_Vec = 0.5 + 0.5 * tanh(1e5 * FlowState_Vec);

        % Upstream node proprieties
        P_Up_Vec = Pressure(Up_Vec) .* Dir_Vec + Pressure(Down_Vec) .* (1 - Dir_Vec);
        T_Up_Vec = T_Vec(Up_Vec)    .* Dir_Vec + T_Vec(Down_Vec)    .* (1 - Dir_Vec);
        R_Up_Vec = R_Vec(Up_Vec)    .* Dir_Vec + R_Vec(Down_Vec)    .* (1 - Dir_Vec);

        % Gas Density per Link
        Rho_Gas_Link = P_Up_Vec ./ (R_Up_Vec .* T_Up_Vec);

        % Smooth Selection for link composition
        Y_Up_Raw = Y_Matrix(Up_Vec, :)   .* Dir_Vec +...
                   Y_Matrix(Down_Vec, :) .* (1 - Dir_Vec);

        % Apply stratification and normalize
        LiqSum_Vec = sum(Y_Up_Raw(:, 1:2), 2);
        isLiquid   = 0.5 + 0.5 * tanh(1e2 * (LiqSum_Vec - 0.05));
        Y_Strat    = Y_Up_Raw;
        Y_Strat(:, 3) = Y_Strat(:, 3) .* (1 - isLiquid);
        Y_Strat = Y_Strat ./ (sum(Y_Strat, 2) + 1e-9);

        % Mixture density
        % We build a density matrix and replace the third column as
        % nessesary from the above process.
        Rho_Mat_Links = repmat(RhoArray, MDLength, 1);
        Rho_Mat_Links(:, 3) = Rho_Gas_Link;
        Rho_Mix_Vec = 1 ./ sum(Y_Strat ./ (Rho_Mat_Links + 1e-9), 2);

        % Momentum Equation
        DeltaP_Vec = X(Up_Vec) - X(Down_Vec);
        FricTerm_Vec = FlowState_Vec .* sqrt(FlowState_Vec.^2 + 1e-6);

        % ODE
        Rho_Fric = max(Rho_Mix_Vec, 1.0);
        XDOT(LinkIdx_Start:LinkIdx_End) = (A_Vec ./ L_Vec) .* ...
            (DeltaP_Vec - (Zeta_Vec ./ (2 * Rho_Fric .* A_Vec.^2)) .* FricTerm_Vec);

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
        isCheck = strcmp(Type_Cell, 'Check')';
        isSignal = strcmp(Type_Cell, 'Signal')';

        % Use Cv Math if valve is Throttle or Check
        isActuated = (isCheck | isThrottle) | isSignal;

        % Control Input handling
        if ~isempty(U)
            Cv_Vec(isActuated) = U(isActuated); 
        end
        
        % Both Valves and Orifice Plates follow the same form.
        Coeff_Vec = isActuated * 2.402e-5 + (~isActuated) .* A_Vec;
        Sqrt_Mult = isActuated * 1.0      + (~isActuated) * 2.0;

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
        isLiquid   = 0.5 + 0.5 * tanh(1e2 * (LiqSum_Vec - 0.05)); 
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

        % Robust Massflow Calculation
        Eps = 0.5;
        DensityTerm = sqrt(Sqrt_Mult .* Rho_Up_Final + 1e-9);
        DPSoft = DP_Eff ./ ((DP_Eff.^2 + Eps).^0.25);
        SignTerm = tanh(1e3 * DeltaP_Vec);
        RawMassflow = Coeff_Vec .* Cv_Vec .* DensityTerm .* DPSoft .* SignTerm;
        
        % Check valve logic: Block massflows if backflow
        isBackflow = DeltaP_Vec < 0;
        BlockMask = isCheck & isBackflow;
        RawMassflow(BlockMask) = 0.0;

        Massflow(ID_Vec) = RawMassflow;
    end

    %% Conservation of Mass and Energy (Vectorized)
    % Gather link data
    Link_Up   = LinkMap(:, 1);
    Link_Down = LinkMap(:, 2);
    Link_Flow = Massflow;

    % Determine flow direction (1.0 if Up -> Down, 0.0 if Down -> Up)
    Link_Dir = 0.5 + 0.5 * tanh(1e7 * Link_Flow);

    % Get source compositions and perform stratification on upstream and
    % downstream nodes
    % Upstream candidates
        Y_Up_Raw   = Y_Matrix(Link_Up, :);
        LiqSum_Up  = sum(Y_Up_Raw(:, 1:2), 2);
        isLiquid_Up = 0.5 + 0.5 * tanh(1e2 * (LiqSum_Up - 0.05));
        isLiquid_Up = isLiquid_Up .* IsTank_Vec(Link_Up); 
        Y_Up_Strat = Y_Up_Raw;
        Y_Up_Strat(:, 3) = Y_Up_Strat(:, 3) .* (1 - isLiquid_Up);
        Y_Up_Strat = Y_Up_Strat ./ (sum(Y_Up_Strat, 2) + 1e-9);
    % Downstream candidates
        Y_Down_Raw   = Y_Matrix(Link_Down, :);
        LiqSum_Down  = sum(Y_Down_Raw(:, 1:2), 2);
        isLiquid_Down = 0.5 + 0.5 * tanh(1e2 * (LiqSum_Down - 0.05));
        isLiquid_Down = isLiquid_Down .* IsTank_Vec(Link_Down); 
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
    TotalBurn = zeros(PLength, 1);
    MassNode = V_Vec .* RhoNode;

    if Combustor
        % % Extract local mass fractions for the combustor node(s)
        Y_OX_Local = Y_Matrix(IsComb_Vec, 1);
        Y_FU_Local = Y_Matrix(IsComb_Vec, 2);

        % Calculate burn demand and limiter
        BurnDemand = B_OX + B_FU;
        ChamberLiquid = (Y_OX_Local + Y_FU_Local) .* MassNode(IsComb_Vec);
        MaxBurn = ChamberLiquid / 1e-4;
        Burn_Scale = min(1.0, MaxBurn / (BurnDemand + 1e-8));
        B_OX = B_OX * Burn_Scale * IgniterMask;
        B_FU = B_FU * Burn_Scale * IgniterMask;
        TotalBurn(IsComb_Vec) = B_OX + B_FU;

        % Combustor consumes liquid accumulated inside the chamber. Net_Liq
        % vector contains the injection massflow added to the node. We now
        % substract the burn rate (consumed liquid inside the chamber turnt
        % to gas)
        Net_Spec(IsComb_Vec, 1) = Net_Spec(IsComb_Vec, 1) - B_OX;
        Net_Spec(IsComb_Vec, 2) = Net_Spec(IsComb_Vec, 2) - B_FU;
        Net_Spec(IsComb_Vec, 3) = Net_Spec(IsComb_Vec, 3) + TotalBurn(IsComb_Vec);
    end
    
    % Effective Bulk Modulus
    % Isothermal Bulk Modulus for Ideal Gas is P
    Beta_Gas = Pressure; 
    Alpha_Vec = max(Alpha_Vec, 1e-4); 
    InvBeta  = (Alpha_Vec ./ Beta_Gas) + ((1 - Alpha_Vec) ./ Beta_Liq);
    Beta_Eff = 1 ./ (InvBeta + 1e-12);

    % Pressure ODE
    dVol_Gas = Net_Gas ./ Rho_Gas_Vec;
    dVol_Liq = Net_Liq ./ Rho_Liq_Vec;
    dVol_Tot = dVol_Liq + dVol_Gas;
    dP_Final = (Beta_Eff ./ V_Vec) .* dVol_Tot;

    % Pressure derivatives for combustors
    if any(IsComb_Vec)
        % Dynamic gas volume
        V_Gas_Actual = V_Vec(IsComb_Vec) .* Alpha_Vec(IsComb_Vec);
        V_Gas_Safe   = max(V_Gas_Actual, 0.001 * V_Vec(IsComb_Vec));

        % Gas generation in the chamber (burnrate - nozzle outflow)
        dM_Gas = TotalBurn(IsComb_Vec) + Net_Gas(IsComb_Vec);
        dP_GasTerm = (R_Vec(IsComb_Vec) .* T_Vec(IsComb_Vec) ./ V_Gas_Safe) .* dM_Gas;

        % Liquid piston (injector inflow - burnrate)
        dM_Liq = Net_Liq(IsComb_Vec) - TotalBurn(IsComb_Vec);
        dVolLiq = dM_Liq ./ max(Rho_Liq_Vec(IsComb_Vec), 1);
        dP_Piston = (Pressure(IsComb_Vec) ./ V_Gas_Safe) .* dVolLiq;

        % Overwrite combustor derivatives
        dP_Final(IsComb_Vec) = dP_GasTerm + dP_Piston;
    end

    % Species derivatives
    Mass_Node = V_Vec .* RhoNode;
    dY_Num    = Net_Spec - Net_Total .* Y_Matrix;
    dY_dt     = dY_Num ./ (Mass_Node + 1e-5);

    % Fixed (boundary) node constraints
    Fixed_Vec = [System.Nodes.Fixed]';
    dP_Final(Fixed_Vec == 1) = 0;
    dY_dt(Fixed_Vec == 1, :) = 0;

    % Reshape matrix and store derivatives
    XDOT(1:PLength) = dP_Final;
    XDOT(YIdx_Start:end) = reshape(dY_dt', [], 1);
end








    