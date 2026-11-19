function [XDOT, Massflow] = PhysicsEngine(X, System, U)
    %% Constants & Setup
    PLength = size(System.Nodes, 2);
    MALength = size(System.Links.Algebraic, 2);
    MDLength = size(System.Links.Dynamic, 2);
    YIdx_Start = PLength + MDLength + 1;

    % Constants 
    Beta_Liq = 2e7;
    RhoArray = System.Constants.RhoArray;
    NumSpecies = length(RhoArray);
    LinkMap = System.LinkMap;
    Rho_Default_Liq = (RhoArray(1) + RhoArray(2)) / 2;

    % Pre-allocation
    XDOT = zeros(size(X));
    Massflow = zeros(MALength + MDLength, 1);

    %% Node Densities (Vectorized)
    Nodes      = System.Nodes;
    R_Vec      = [Nodes.R]';   
    T_Vec      = [Nodes.Temp]';    
    V_Vec      = [Nodes.V]';      
    IsComb_Vec = [Nodes.IsCombustor]';
    Type_Cell  = {Nodes.Type}';
    IsTank_Vec = strcmp(Type_Cell, 'Tank');
    Combustor  = any(IsComb_Vec);

    % Process Pressures
    P_Raw    = X(1:PLength);
    Pressure = sqrt(P_Raw.^2 + 1e-6);

    % Process Mass Fractions
    Y_All = X(YIdx_Start:end);
    Y_Matrix = reshape(Y_All, [NumSpecies, PLength])';
    Y_Sum    = sum(Y_Matrix, 2);
    Y_Matrix = Y_Matrix ./ (Y_Sum + 1e-9);

    %% Combustor Temp Logic (Instantaneous)
    CStarEff = 0.83;
    CEA_OF = System.Combustion.CEA_OF;
    CEA_Temp = System.Combustion.CEA_Temp * CStarEff^2;

    % Current OF Ratio
    Current_OF = Y_Matrix(:, 1) ./ (Y_Matrix(:, 2) + 1e-9);
    
    if Combustor
        % Break the algebraic loop by estimating instantaneous OF ratio from current pressures
        % UPDATED IDS: OX Manifold = 11, FU Manifold = 12, Chamber = 13
        dP_OX = max(Pressure(11) - Pressure(13), 0);
        dP_FU = max(Pressure(12) - Pressure(13), 0);
        
        AlgLinks = System.Links.Algebraic;
        L_OX = AlgLinks(find([AlgLinks.ID] == 12, 1)); % Inj OX ID = 12
        L_FU = AlgLinks(find([AlgLinks.ID] == 13, 1)); % Inj FU ID = 13
        
        Est_OX = L_OX.Cv * L_OX.A * sqrt(2 * 1141 * dP_OX);
        Est_FU = L_FU.Cv * L_FU.A * sqrt(2 * 786 * dP_FU);
        
        PilotMass = 1e-3; % 1 gram/s stability baseline to prevent wild OF swings on startup
        Current_OF(IsComb_Vec) = (Est_OX + PilotMass) ./ (Est_FU + PilotMass);
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
        BurnFactor = MaskLean .* MaskRich .* IgniterMask;
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
    IsPureGas = 0.5 + 0.5 * tanh(1e3 * (1e-9 - InvRhoLiq));
    Rho_Liq_Vec = (1 - IsPureGas) .* Rho_Raw + IsPureGas .* Rho_Default_Liq;

    % Proportion of gas (Alpha)
    Term_LiqVol = Rho_Gas_Vec ./ Rho_Liq_Vec;
    Alpha_Vec = Y_N2 ./ (Y_N2 + Term_LiqVol .* (1 - Y_N2) + 1e-9);

    if ~Combustor
        Alpha_Vec(IsComb_Vec) = 1.0;
    end

    % Final Node Properties
    RhoNode = Alpha_Vec .* Rho_Gas_Vec + (1 - Alpha_Vec) .* Rho_Liq_Vec;

    %% Dynamic Links (Vectorized)
    DynLinks = System.Links.Dynamic;

    if ~isempty(DynLinks)
        Up_Vec   = [DynLinks.Up]';
        Down_Vec = [DynLinks.Down]';
        A_Vec    = [DynLinks.A]';
        L_Vec    = [DynLinks.L]';
        Zeta_Vec = [DynLinks.Zeta]';
        ID_Vec   = [DynLinks.ID]';

        LinkIdx_Start = PLength + 1;
        LinkIdx_End   = PLength + MDLength;
        FlowState_Vec = X(LinkIdx_Start:LinkIdx_End);
        Dir_Vec = 0.5 + 0.5 * tanh(1e3 * FlowState_Vec);

        P_Up_Vec = Pressure(Up_Vec) .* Dir_Vec + Pressure(Down_Vec) .* (1 - Dir_Vec);
        T_Up_Vec = T_Vec(Up_Vec)    .* Dir_Vec + T_Vec(Down_Vec)    .* (1 - Dir_Vec);
        R_Up_Vec = R_Vec(Up_Vec)    .* Dir_Vec + R_Vec(Down_Vec)    .* (1 - Dir_Vec);

        Rho_Gas_Link = P_Up_Vec ./ (R_Up_Vec .* T_Up_Vec);
        Y_Up_Raw = Y_Matrix(Up_Vec, :) .* Dir_Vec + Y_Matrix(Down_Vec, :) .* (1 - Dir_Vec);

        LiqSum_Vec = sum(Y_Up_Raw(:, 1:2), 2);
        isLiquid   = 0.5 + 0.5 * tanh(1e2 * (LiqSum_Vec - 0.05));
        Y_Strat    = Y_Up_Raw;
        Y_Strat(:, 3) = Y_Strat(:, 3) .* (1 - isLiquid);
        Y_Strat = Y_Strat ./ (sum(Y_Strat, 2) + 1e-9);

        Rho_Mat_Links = repmat(RhoArray, MDLength, 1);
        Rho_Mat_Links(:, 3) = Rho_Gas_Link;
        Rho_Mix_Vec = 1 ./ sum(Y_Strat ./ (Rho_Mat_Links + 1e-9), 2);

        DeltaP_Vec = X(Up_Vec) - X(Down_Vec);
        Laminar_Damping = 0.1; 
        FricTerm_Vec = FlowState_Vec .* sqrt(FlowState_Vec.^2 + 1e-6) + Laminar_Damping * FlowState_Vec;

        Rho_Fric = max(Rho_Mix_Vec, 1.0);
        XDOT(LinkIdx_Start:LinkIdx_End) = (A_Vec ./ L_Vec) .* ...
            (DeltaP_Vec - (Zeta_Vec ./ (2 * Rho_Fric .* A_Vec.^2)) .* FricTerm_Vec);

        Massflow(ID_Vec) = FlowState_Vec;
    end

    %% Algebraic Links (Vectorized)
    AlgLinks = System.Links.Algebraic;

    if ~isempty(AlgLinks)
        Up_Vec   = [AlgLinks.Up]';
        Down_Vec = [AlgLinks.Down]';
        ID_Vec   = [AlgLinks.ID]';
        Cv_Vec   = [AlgLinks.Cv]';
        A_Vec    = [AlgLinks.A]';

        Type_Cell  = {AlgLinks.Type};
        isThrottle = strcmp(Type_Cell, 'Throttle')';
        isCheck = strcmp(Type_Cell, 'Check')';
        isSignal = strcmp(Type_Cell, 'Signal')';

        isActuated = (isCheck | isThrottle) | isSignal;
        if ~isempty(U)
            Cv_Vec(isActuated) = U(isActuated); 
        end
        
        Coeff_Vec = isActuated * 2.402e-5 + (~isActuated) .* A_Vec;
        Sqrt_Mult = isActuated * 1.0      + (~isActuated) * 2.0;

        DeltaP_Vec = X(Up_Vec) - X(Down_Vec);
        Dir_Vec    = 0.5 + 0.5 * tanh(1e3 * DeltaP_Vec);

        P_Up_Vec = Pressure(Up_Vec) .* Dir_Vec + Pressure(Down_Vec) .* (1 - Dir_Vec);
        T_Up_Vec = T_Vec(Up_Vec)    .* Dir_Vec + T_Vec(Down_Vec)    .* (1 - Dir_Vec);
        R_Up_Vec = R_Vec(Up_Vec)    .* Dir_Vec + R_Vec(Down_Vec)    .* (1 - Dir_Vec);
        Rho_Gas_Flow = P_Up_Vec ./ (R_Up_Vec .* T_Up_Vec);
        Y_Up_Raw = Y_Matrix(Up_Vec, :) .* Dir_Vec + Y_Matrix(Down_Vec, :) .* (1 - Dir_Vec);

        LiqSum_Vec = sum(Y_Up_Raw(:, 1:2), 2);
        isLiquid   = 0.5 + 0.5 * tanh(1e2 * (LiqSum_Vec - 0.05)); 
        Y_Flow     = Y_Up_Raw;
        Y_Flow(:, 3) = Y_Flow(:, 3) .* (1 - isLiquid);
        Y_Flow     = Y_Flow ./ (sum(Y_Flow, 2) + 1e-9);

        IsSrcComb = (IsComb_Vec(Up_Vec) .* Dir_Vec) + (IsComb_Vec(Down_Vec) .* (1 - Dir_Vec));

        Rho_Mat_Links = repmat(RhoArray, MALength, 1);
        Rho_Mat_Links(:, 3) = Rho_Gas_Flow;
        Rho_Mix_Vec = 1 ./ sum(Y_Flow ./ (Rho_Mat_Links + 1e-9), 2);
        Rho_Up_Final = (IsSrcComb .* Rho_Gas_Flow) + ((1 - IsSrcComb) .* Rho_Mix_Vec);

        IsLiqFlow = isLiquid .* (1 - IsSrcComb);
        NodesGamma = [System.Nodes.Gamma]';
        G_Node     = NodesGamma(Up_Vec) .* Dir_Vec + NodesGamma(Down_Vec) .* (1 - Dir_Vec);
        G_Up       = IsLiqFlow * 100 + (1 - IsLiqFlow) .* G_Node;

        Rc = (2 ./ (G_Up + 1)) .^ (G_Up ./ (G_Up - 1));
        DP_Choked = P_Up_Vec .* (1 - Rc);
        DP_Raw    = sqrt(DeltaP_Vec.^2 + 1e-9);

        DP_Eff = 0.5 * (DP_Raw + DP_Choked - sqrt((DP_Raw - DP_Choked).^2 + 1e-4));

        Eps = 0.5;
        DensityTerm = sqrt(Sqrt_Mult .* Rho_Up_Final + 1e-9);
        DPSoft = DP_Eff ./ ((DP_Eff.^2 + Eps).^0.25);
        SignTerm = DeltaP_Vec ./ sqrt(DeltaP_Vec.^2 + 2500);
        RawMassflow = Coeff_Vec .* Cv_Vec .* DensityTerm .* DPSoft .* SignTerm;
        
        Check_Open_Factor = 0.5 + 0.5 * tanh(50 * DeltaP_Vec);
        CheckMask = (~isCheck) + (isCheck .* Check_Open_Factor);
        RawMassflow = RawMassflow .* CheckMask;

        Massflow(ID_Vec) = RawMassflow;
    end

    %% Conservation of Mass and Energy (Vectorized)
    Link_Up   = LinkMap(:, 1);
    Link_Down = LinkMap(:, 2);
    Link_Flow = Massflow;
    Link_Dir = 0.5 + 0.5 * tanh(1e4 * Link_Flow);

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

    Y_Trans = Y_Up_Strat .* Link_Dir + Y_Down_Strat .* (1 - Link_Dir);

    Flux_Total = Link_Flow;
    Flux_Gas   = Link_Flow .* Y_Trans(:, 3);
    Flux_Liq   = Link_Flow .* (1 - Y_Trans(:, 3));
    Flux_Spec  = Link_Flow .* Y_Trans;

    Net_Total = accumarray(Link_Down, Flux_Total, [PLength, 1]) - accumarray(Link_Up, Flux_Total, [PLength, 1]);
    Net_Gas   = accumarray(Link_Down, Flux_Gas, [PLength, 1]) - accumarray(Link_Up, Flux_Gas, [PLength, 1]);
    Net_Liq   = accumarray(Link_Down, Flux_Liq, [PLength, 1]) - accumarray(Link_Up, Flux_Liq, [PLength, 1]);

    Net_Spec = zeros(PLength, 3);
    Net_Spec(:,1) = accumarray(Link_Down, Flux_Spec(:,1), [PLength, 1]) - accumarray(Link_Up, Flux_Spec(:,1), [PLength, 1]);
    Net_Spec(:,2) = accumarray(Link_Down, Flux_Spec(:,2), [PLength, 1]) - accumarray(Link_Up, Flux_Spec(:,2), [PLength, 1]);
    Net_Spec(:,3) = accumarray(Link_Down, Flux_Spec(:,3), [PLength, 1]) - accumarray(Link_Up, Flux_Spec(:,3), [PLength, 1]);

    TotalBurn = zeros(PLength, 1);
    MassNode = V_Vec .* RhoNode;

    if Combustor
        % Extract INSTANTANEOUS mass flows across the injectors natively
        B_OX = max(Massflow(12), 0);
        B_FU = max(Massflow(13), 0);

        Y_OX_Local = Y_Matrix(IsComb_Vec, 1);
        Y_FU_Local = Y_Matrix(IsComb_Vec, 2);

        BurnDemand = B_OX + B_FU;
        ChamberLiquid = (Y_OX_Local + Y_FU_Local) .* MassNode(IsComb_Vec);
        MaxBurn = ChamberLiquid / 1e-4;
        Burn_Scale = min(1.0, MaxBurn / (BurnDemand + 1e-8));
        B_OX = B_OX * Burn_Scale * IgniterMask;
        B_FU = B_FU * Burn_Scale * IgniterMask;
        TotalBurn(IsComb_Vec) = B_OX + B_FU;

        Net_Spec(IsComb_Vec, 1) = Net_Spec(IsComb_Vec, 1) - B_OX;
        Net_Spec(IsComb_Vec, 2) = Net_Spec(IsComb_Vec, 2) - B_FU;
        Net_Spec(IsComb_Vec, 3) = Net_Spec(IsComb_Vec, 3) + TotalBurn(IsComb_Vec);
    end
    
    Beta_Gas = Pressure; 
    Alpha_Vec = max(Alpha_Vec, 1e-4); 
    InvBeta  = (Alpha_Vec ./ Beta_Gas) + ((1 - Alpha_Vec) ./ Beta_Liq);
    Beta_Eff = 1 ./ (InvBeta + 1e-12);

    dVol_Gas = Net_Gas ./ Rho_Gas_Vec;
    dVol_Liq = Net_Liq ./ Rho_Liq_Vec;
    dVol_Tot = dVol_Liq + dVol_Gas;
    dP_Final = (Beta_Eff ./ V_Vec) .* dVol_Tot;

    if any(IsComb_Vec)
        V_Gas_Actual = V_Vec(IsComb_Vec) .* Alpha_Vec(IsComb_Vec);
        V_Gas_Safe   = max(V_Gas_Actual, 0.001 * V_Vec(IsComb_Vec));

        dM_Gas = TotalBurn(IsComb_Vec) + Net_Gas(IsComb_Vec);
        dP_GasTerm = (R_Vec(IsComb_Vec) .* T_Vec(IsComb_Vec) ./ V_Gas_Safe) .* dM_Gas;

        dM_Liq = Net_Liq(IsComb_Vec) - TotalBurn(IsComb_Vec);
        dVolLiq = dM_Liq ./ max(Rho_Liq_Vec(IsComb_Vec), 1);
        dP_Piston = (Pressure(IsComb_Vec) ./ V_Gas_Safe) .* dVolLiq;

        dP_Final(IsComb_Vec) = dP_GasTerm + dP_Piston;
    end

    Mass_Node = V_Vec .* RhoNode;
    dY_Num    = Net_Spec - Net_Total .* Y_Matrix;
    dY_dt     = dY_Num ./ (Mass_Node + 1e-3);

    Fixed_Vec = [System.Nodes.Fixed]';
    dP_Final(Fixed_Vec == 1) = 0;
    dY_dt(Fixed_Vec == 1, :) = 0;
    % Smoothly dampen negative pressure changes
    P_Raw = X(1:PLength);
    Floor_Mask = 0.5 + 0.5 * tanh((P_Raw - 6895) / 1000); 
    
    % Allow pressure to rise freely, but clamp pressure drops near zero
    dP_Final = dP_Final .* Floor_Mask + max(dP_Final, 0) .* (1 - Floor_Mask);

    XDOT(1:PLength) = dP_Final;
    XDOT(YIdx_Start:end) = reshape(dY_dt', [], 1);
end