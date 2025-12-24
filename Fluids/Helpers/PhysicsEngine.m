function XDOT = PhysicsEngine(X, System, U, isSmooth)
    %% Constants & Setup
    PLength = size(System.Nodes, 2);
    MALength = size(System.Links.Algebraic, 2);
    MDLength = size(System.Links.Dynamic, 2);
    YIdx_Start = PLength + MDLength + 1;

    % Auto-detect Symbolic Input
    isSymbolicType = isa(X, 'sym');

    % Constants (approximate)
    Beta_Liq = 1.1e9; % Pa (Bulk Modulus of Liquid ~ Fuel)

    % Pre-allocation
    if isSymbolicType
        XDOT = sym(zeros(size(X)));
        Pressure = sym(zeros(PLength, 1));
        Massflow = sym(zeros(MALength + MDLength, 1));
        RhoNode = sym(zeros(PLength, 1));
        Vol_Ullage = sym(zeros(PLength, 1));
        Alpha_Vec = sym(zeros(PLength, 1));
        Rho_Liq_Vec = sym(zeros(PLength, 1));
    elseif ~isreal(X)
        XDOT = complex(zeros(size(X)));
        Pressure = complex(zeros(PLength, 1));
        Massflow = complex(zeros(MALength + MDLength, 1));
        RhoNode = complex(zeros(PLength, 1));
        Vol_Ullage = complex(zeros(PLength, 1));
        Alpha_Vec = complex(zeros(PLength, 1));
        Rho_Liq_Vec = complex(zeros(PLength, 1));
    else
        XDOT = zeros(size(X));
        Pressure = zeros(PLength, 1);
        Massflow = zeros(MALength + MDLength, 1);
        RhoNode = zeros(PLength, 1);
        Vol_Ullage = zeros(PLength, 1);
        Alpha_Vec = zeros(PLength, 1);
        Rho_Liq_Vec = zeros(PLength, 1);
    end

    RhoArray = System.Constants.RhoArray;
    NumSpecies = length(RhoArray);
    LinkMap = System.LinkMap;
    Rho_Default_Liq = mean(RhoArray(1:2));

    %% Node Densities
    Y_All = X(YIdx_Start:end);
    Y_Matrix = reshape(Y_All, [NumSpecies, PLength])';

    for i = 1:PLength
        % Node Parameters
        Node = System.Nodes(i);
        Temp = Node.Temp;
        R_Gas = Node.R;

        % Robust Pressure
        P = sqrt(X(i)^2 + 1e-6);
        Pressure(i) = P; 

        % Normalize Composition
        Y_Vec = Y_Matrix(i, :);
        Y_Sum = sum(Y_Vec);
        Y_Vec = Y_Vec / (Y_Sum + 1e-9);
        Y_Matrix(i, :) = Y_Vec;

        % Calculate Phase Split (Optimized for Symbolic Engine)
        Rho_Gas = P / (R_Gas * Temp);
        
        % Mass Fractions
        Y_OX = Y_Vec(1);    Y_FU = Y_Vec(2);    Y_N2 = Y_Vec(3);
        
        % Term representing Liquid Volume per unit Mass
        InvRhoLiq = (Y_OX / RhoArray(1)) + (Y_FU / RhoArray(2));
        
        % Calculate "Raw" density (Unstable near zero)
        Rho_Raw = (Y_OX + Y_FU) / (InvRhoLiq + 1e-12);
         
        % Weighting factor: 1 if gas, 0 if liquid
        % Use high-gain tanh (1e12) to switch cleanly when InvRhoLiq < 1e-9
        IsPureGas = 0.5 + 0.5 * tanh(1e6 * (1e-9 - InvRhoLiq)); 
         
        % Blend: If Gas, use Default. If Liquid, use Calculated.
        Rho_Liq_Local = (1 - IsPureGas) * Rho_Raw + IsPureGas * Rho_Default_Liq;

        % Calculate Volume Ratio term: (Rho_Gas / Rho_Liq)
        Term_LiqVol = Rho_Gas / Rho_Liq_Local;
        
        % Calculate Void Fraction (Alpha) from Mass Fraction
        Alpha = (Y_N2) / (Y_N2 + Term_LiqVol * (1 - Y_N2) + 1e-9);
        if Node.IsCombustor
            Alpha = 1;
        end
        Alpha_Vec(i) = Alpha;

        % Store
        Rho_Liq_Vec(i) = Rho_Liq_Local;
        RhoNode(i) = Alpha * Rho_Gas + (1 - Alpha) * Rho_Liq_Local;
        Vol_Ullage(i) = System.Nodes(i).V * Alpha;
    end

    %% Dynamic Links (Pipes)
    for i = 1:MDLength
        L = System.Links.Dynamic(i);
        UpIDX = L.Up; DwIDX = L.Down;
        DeltaP = X(UpIDX) - X(DwIDX);
        FlowState = X(PLength + i);
        
        % 1. Determine Flow Direction (Symbolic Safe)
        % Returns ~1 if Flow > 0, ~0 if Flow < 0, smooth blend if symbolic
        IsForward = GetPosFactor(FlowState, isSmooth);

        % 2. Get "Upstream" Properties via Blending
        Y_UpStr = Y_Matrix(UpIDX, :) * IsForward + Y_Matrix(DwIDX, :) * (1 - IsForward);
        P_UpStr = X(UpIDX) * IsForward + X(DwIDX) * (1 - IsForward);

        % 3. Apply Stratification Logic
        % This applies the "Liquid Check" to our effectively selected upstream composition
        Y_Stratified = NodeFluidFlow(Y_UpStr, isSmooth);

        % Calculate Density of this Stratified Mixture
        % Fetch Upstream Thermodynamic Proprieties
        T_Up = System.Nodes(UpIDX).Temp * IsForward + System.Nodes(DwIDX).Temp * (1 - IsForward);
        R_Up = System.Nodes(UpIDX).R    * IsForward + System.Nodes(DwIDX).R    * (1 - IsForward);
        RhoGas = P_UpStr / (R_Up * T_Up);
        
        % Build Local Density Array
        LinkRhoArray = RhoArray; 
        if isSymbolicType, LinkRhoArray = sym(LinkRhoArray); end
        LinkRhoArray(3) = RhoGas;
        
        % Calculate Mixture Density
        Rho = 1 / sum(Y_Stratified ./ (LinkRhoArray + 1e-9));
        FricTerm = FlowState * sqrt(FlowState^2 + 1e-6);

        % Massflow ODE
        XDOT(PLength + i) = L.A/L.L * (DeltaP - L.Zeta/(2*Rho*L.A*L.A) * FricTerm);
        Massflow(L.ID) = FlowState;
    end

    %% Algebraic Links (Valves)
    for i = 1:MALength
        L = System.Links.Algebraic(i);
        if strcmp(L.Type, 'Throttle')
            UpIDX = L.Up; DwIDX = L.Down;
            DeltaP = X(UpIDX) - X(DwIDX);
            
            % Determine Direction & Upstream Properties
            % (Strict upwinding is required for choked flow logic)
            % Smooth switch to select P_up and P_down
            Dir = 0.5 + 0.5 * tanh(1e3 * DeltaP); 
            P_up   = X(UpIDX) * Dir + X(DwIDX) * (1 - Dir);
            Rho_up = RhoNode(UpIDX) * Dir + RhoNode(DwIDX) * (1 - Dir);
            G_Up = System.Nodes(UpIDX).Gamma * Dir + System.Nodes(DwIDX).Gamma * (1 - Dir);

            % Calculate Effective Pressure Drop (Choked Flow Check)
            % Critical Pressure Ratio (approx 0.5 for N2)
            Rc = (2 / (G_Up + 1)) ^ (G_Up / (G_Up - 1)); 
            DP_Choked = P_up * (1 - Rc); % Max physical DeltaP before sonic choking
            
            % Smooth Min: min(a,b) approx 0.5*(a+b - sqrt((a-b)^2))
            DP_Raw = sqrt(DeltaP^2 + 1e-9);
            DP_Eff = 0.5 * (DP_Raw + DP_Choked - sqrt((DP_Raw - DP_Choked)^2 + 1e-4));
            
            % Standard valve equation terms
            SqrtTerm = sqrt(DP_Eff * Rho_up + 1e-9);
            SignTerm = tanh(1e3 * DeltaP);

            % Get Cv (Input U or fixed)
            Cv = L.Cv;
            if ~isempty(U), Cv = U(i); end

            % Compute Massflow
            Massflow(L.ID) = 2.402e-5 * Cv * SqrtTerm * SignTerm;
        end
    end

    %% Conservation of Mass and Energy
    XDOT_Species_Mat = zeros(PLength, NumSpecies);
    if isSymbolicType
        XDOT_Species_Mat = sym(zeros(PLength, NumSpecies));
    end

    for i = 1:PLength
        Node = System.Nodes(i);
        if ~Node.Fixed
            % Initialize Accumulators
            mdot_Gas_Net = 0; mdot_Liq_Net = 0; m_Net_Total = 0;
            Species_Term = zeros(1, NumSpecies);
            Mass = Node.V * RhoNode(i);
            
            if isSymbolicType
                mdot_Gas_Net = sym(0); mdot_Liq_Net = sym(0);
                m_Net_Total = sym(0); Species_Term = sym(zeros(1, NumSpecies));
            end

            % Helper to process flow (saves code duplication)
            Links = [Node.LinksIN, Node.LinksOUT];
            Dirs  = [ones(1, length(Node.LinksIN)), -1*ones(1, length(Node.LinksOUT))];
            
            for k = 1:length(Links)
                lIdx = Links(k);
                flow_raw = Massflow(lIdx);
                flow_signed = flow_raw * Dirs(k); % Positive = Entering Node
                
                % Determine Source Node for Composition
                if Dirs(k) == 1, neighbor = LinkMap(lIdx, 1);
                else,            neighbor = LinkMap(lIdx, 2); end
                
                % Check Flow Direction (Is it actually entering?)
                isInflow = GetPosFactor(flow_signed, isSmooth);
                
                % Get Composition
                Y_Neigh = NodeFluidFlow(Y_Matrix(neighbor, :), isSmooth);
                Y_Self  = NodeFluidFlow(Y_Matrix(i, :), isSmooth);
                Y_Flux  = Y_Neigh * isInflow + Y_Self * (1 - isInflow);
                
                % Split Flow by Phase (N2 is idx 3)
                GasFrac = Y_Flux(3);
                LiqFrac = 1 - GasFrac;
                mdot_Gas_Net = mdot_Gas_Net + flow_signed * GasFrac;
                mdot_Liq_Net = mdot_Liq_Net + flow_signed * LiqFrac;
                m_Net_Total  = m_Net_Total + flow_signed;
                Species_Term = Species_Term + flow_signed * (Y_Flux - Y_Matrix(i, :));
            end
            
            %% ODE's
            Rho_Liq_Local = Rho_Liq_Vec(i);
            Temp = Node.Temp;   R = Node.R;

            % Calculate Gas Volume
            if Node.IsCombustor
                mdot_Gas_Net = mdot_Gas_Net + mdot_Liq_Net;
                mdot_Liq_Net = 0; 
                V_Safe = Node.V;
            else
                V_Safe = Vol_Ullage(i);
                V_Safe = V_Safe + 1e-6;
            end
            
            dP_GasTerm = (mdot_Gas_Net * R * Temp);
            dP_LiqTerm = Pressure(i) * (mdot_Liq_Net / max(Rho_Liq_Local, 100));

            dP_GasMode = (dP_GasTerm + dP_LiqTerm) / V_Safe;
            dP_HydraulicMode = (Beta_Liq / (Node.V * max(Rho_Liq_Local, 100))) * m_Net_Total;
            
            if Node.IsCombustor
                XDOT(i) = dP_GasMode;
            else
                % Smooth Blend
                % Blend Factor: 0 = Hydraulic, 1 = Gas
                k = 1e3;
                BlendFactor = 0.5 + 0.5 * tanh(k * (Alpha_Vec(i) - 0.002));
                
                % Interpolate
                XDOT(i) = BlendFactor * dP_GasMode + (1 - BlendFactor) * dP_HydraulicMode;
            end
            
            % Species Transport ODE
            XDOT_Species_Mat(i, :) = Species_Term / (Mass + 1e-9);
        else
            XDOT(i) = 0;
            XDOT_Species_Mat(i, :) = 0;
        end
    end
    
    XDOT(YIdx_Start:end) = reshape(XDOT_Species_Mat', [], 1);
end
function factor = GetPosFactor(Val, isSym)
    if isSym
        k = 100;
        factor = 0.5 + 0.5 * tanh(k * Val);
    else
        factor = double(Val > 0);
    end
end
function YVector = NodeFluidFlow(Y, isSym)
    if isSym
        k = 100;
        F = 0.5 + 0.5 * tanh(k * sum(Y(1:2)));
        YVector = Y;
        YVector(3) = YVector(3) * (1 - F);
        YVector = YVector / sum(YVector);
    else
        if sum(Y(1:2)) > 0.05
            YVector = Y;
            YVector(3) = 0;
            YVector = YVector / sum(YVector);
        else 
            YVector = Y;
        end
    end
end