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

    %% Node Densities
    Y_All = X(YIdx_Start:end);
    Y_Matrix = reshape(Y_All, [NumSpecies, PLength])';

    for i = 1:PLength
        % Node Parameters
        Node = System.Nodes(i);
        
        % Robust Pressure
        P = sqrt(X(i)^2 + 1e-6);
        Pressure(i) = P; 

        % Normalize Composition
        Y_Vec = Y_Matrix(i, :);
        Y_Sum = sum(Y_Vec);
        Y_Vec = Y_Vec / (Y_Sum + 1e-9);
        Y_Matrix(i, :) = Y_Vec;

        % Calculate Gas Density
        Rho_Gas = P / (Node.R * Node.Temp);
        
        % Liquid Density
        Y_OX = Y_Vec(1); Y_FU = Y_Vec(2); Y_N2 = Y_Vec(3);
        InvRhoLiq = (Y_OX / RhoArray(1)) + (Y_FU / RhoArray(2));
        Rho_Raw = (Y_OX + Y_FU) / (InvRhoLiq + 1e-12);
         
        % Blend: If Gas, use Default. If Liquid, use Calculated.
        IsPureGas = Sigmoid(1e-9 - InvRhoLiq, 1e6);
        Rho_Liq_Local = (1 - IsPureGas) * Rho_Raw + IsPureGas * Rho_Default_Liq;

        % Gas Proportion
        Term_LiqVol = Rho_Gas / Rho_Liq_Local;
        Alpha = (Y_N2) / (Y_N2 + Term_LiqVol * (1 - Y_N2) + 1e-9);
        
        if Node.IsCombustor, Alpha = 1; end
        
        % Store
        Alpha_Vec(i) = Alpha;
        Rho_Liq_Vec(i) = Rho_Liq_Local;
        RhoNode(i) = Alpha * Rho_Gas + (1 - Alpha) * Rho_Liq_Local;
        Vol_Ullage(i) = Node.V * Alpha;
    end

    %% Dynamic Links (Pipes)
    for i = 1:MDLength
        L = System.Links.Dynamic(i);
        UpIDX = L.Up; DwIDX = L.Down;
        FlowState = X(PLength + i);
        
        % Determine Flow Direction
        Dir = Sigmoid(FlowState, 100);

        % Get Upstream Properties
        [RhoGas, Y_Stratified, ~] = GetUpstreamParams(System, Pressure, Y_Matrix, UpIDX, DwIDX, Dir);

        % Calculate Mixture Density
        LinkRhoArray = RhoArray; 
        LinkRhoArray(3) = RhoGas;
        Rho = 1 / sum(Y_Stratified ./ (LinkRhoArray + 1e-9));
        
        % Momentum ODE
        DeltaP = X(UpIDX) - X(DwIDX);
        FricTerm = FlowState * sqrt(FlowState^2 + 1e-6);
        XDOT(PLength + i) = L.A/L.L * (DeltaP - L.Zeta/(2*Rho*L.A*L.A) * FricTerm);
        Massflow(L.ID) = FlowState;
    end

    %% Algebraic Links (Valves & Orifices)
    for i = 1:MALength
        L = System.Links.Algebraic(i);
        UpIDX = L.Up; DwIDX = L.Down;
        DeltaP = X(UpIDX) - X(DwIDX);
        
        % Determine Direction
        Dir = Sigmoid(DeltaP, 1e3);

        % Get Upstream Properties (Unified Logic)
        [Rho_Gas_Flow, Y_Flow, P_Up] = GetUpstreamParams(System, Pressure, Y_Matrix, UpIDX, DwIDX, Dir);

        % Handle Combustor Sources (Inject Gas Only)
        IsSrcComb = (System.Nodes(UpIDX).IsCombustor * Dir) + ...
                    (System.Nodes(DwIDX).IsCombustor * (1 - Dir));
        if IsSrcComb > 0.5
            Rho_up = Rho_Gas_Flow;
        else
            LinkRhoArray = RhoArray;
            LinkRhoArray(3) = Rho_Gas_Flow;
            Rho_up = 1 / sum(Y_Flow ./ (LinkRhoArray + 1e-9));
        end

        % Choked Flow Check
        Liq_Frac = sum(Y_Flow(1:2));
        IsLiqFlow = Sigmoid(Liq_Frac - 0.1, 10) * (1 - IsSrcComb);
        
        G_Node = System.Nodes(UpIDX).Gamma * Dir + System.Nodes(DwIDX).Gamma * (1 - Dir);
        G_Up = IsLiqFlow * 100 + (1 - IsLiqFlow) * G_Node;

        % Critical Pressure Ratio & Effective DP
        Rc = (2 / (G_Up + 1)) ^ (G_Up / (G_Up - 1)); 
        DP_Choked = P_Up * (1 - Rc); 
        DP_Raw = sqrt(DeltaP^2 + 1e-9);
        DP_Eff = SmoothMin(DP_Raw, DP_Choked);
        
        % Compute Massflow
        SignTerm = tanh(1e3 * DeltaP);
        
        % Throttle Valve Logic
        if strcmp(L.Type, 'Throttle')
            SqrtTerm = sqrt(DP_Eff * Rho_up + 1e-9);
            Cv = L.Cv;
            if ~isempty(U), Cv = U(i); end
            Massflow(L.ID) = 2.402e-5 * Cv * SqrtTerm * SignTerm;
        
        % Orifice Plate Logic
        elseif strcmp(L.Type, 'Orifice')
            SqrtTerm = sqrt(2 * Rho_up * DP_Eff + 1e-9);
            Massflow(L.ID) = L.Cv * L.A * SqrtTerm * SignTerm;
        end
    end

    %% Conservation of Mass and Energy
    XDOT_Species_Mat = zeros(PLength, NumSpecies);
    if ~isreal(X), XDOT_Species_Mat = complex(zeros(PLength, NumSpecies)); end

    for i = 1:PLength
        Node = System.Nodes(i);
        if ~Node.Fixed
            % Initialize accumulators
            mdot_Gas_Net = 0; mdot_Liq_Net = 0; m_Net_Total = 0;
            Species_Term = zeros(1, NumSpecies);
            Mass = Node.V * RhoNode(i);
            
            % Process Links
            Links = [Node.LinksIN, Node.LinksOUT];
            Dirs  = [ones(1, length(Node.LinksIN)), -1*ones(1, length(Node.LinksOUT))];
            
            for k = 1:length(Links)
                lIdx = Links(k);
                flow_raw = Massflow(lIdx);
                flow_signed = flow_raw * Dirs(k); 
                
                % Determine Neighbor & Inflow Factor
                neighbor = LinkMap(lIdx, 1) * (Dirs(k) == 1) + LinkMap(lIdx, 2) * (Dirs(k) ~= 1);
                isInflow = Sigmoid(flow_signed, 100);
                
                % Upwind Composition 
                Y_Neigh = ApplyStratification(Y_Matrix(neighbor, :));
                Y_Self  = ApplyStratification(Y_Matrix(i, :));
                Y_Flux  = Y_Neigh * isInflow + Y_Self * (1 - isInflow);
                
                % Accumulate
                GasFrac = Y_Flux(3);
                LiqFrac = 1 - GasFrac;
                
                mdot_Gas_Net = mdot_Gas_Net + flow_signed * GasFrac;
                mdot_Liq_Net = mdot_Liq_Net + flow_signed * LiqFrac;
                m_Net_Total  = m_Net_Total + flow_signed;
                Species_Term = Species_Term + flow_signed * (Y_Flux - Y_Matrix(i, :));
            end
            
            %% Derivatives
            if Node.IsCombustor
                % If the node is a combustor, all the massflow is Gaseous.
                mdot_Gas_Net = mdot_Gas_Net + mdot_Liq_Net;
                mdot_Liq_Net = 0; 
                V_Safe = Node.V;
            else
                V_Safe = Vol_Ullage(i) + 1e-6;
            end
            
            dP_GasTerm = (mdot_Gas_Net * Node.R * Node.Temp);
            dP_LiqTerm = Pressure(i) * (mdot_Liq_Net / max(Rho_Liq_Vec(i), 100));

            dP_GasMode = (dP_GasTerm + dP_LiqTerm) / V_Safe;
            dP_HydraulicMode = (Beta_Liq / (Node.V * max(Rho_Liq_Vec(i), 100))) * m_Net_Total;
            
            if Node.IsCombustor
                % If Combustor, only Gas mode.
                XDOT(i) = dP_GasMode;
            else
                % Blend based on Void Fraction
                BlendFactor = Sigmoid(Alpha_Vec(i) - 0.002, 1e3);
                XDOT(i) = BlendFactor * dP_GasMode + (1 - BlendFactor) * dP_HydraulicMode;
            end
            
            XDOT_Species_Mat(i, :) = Species_Term / (Mass + 1e-9);
        end
    end
    
    XDOT(YIdx_Start:end) = reshape(XDOT_Species_Mat', [], 1);
end

%% Helper Functions

function [RhoGas, Y_Strat, P_Up, T_Up] = GetUpstreamParams(System, Pressure, Y_Matrix, uIdx, dIdx, Dir)
    % Blends thermodynamic properties from the Upstream and Downstream nodes
    % based on the flow direction 'Dir' (0 to 1).
    
    % Blend Thermo Props
    T_Up = System.Nodes(uIdx).Temp * Dir + System.Nodes(dIdx).Temp * (1 - Dir);
    R_Up = System.Nodes(uIdx).R    * Dir + System.Nodes(dIdx).R    * (1 - Dir);
    P_Up = Pressure(uIdx)          * Dir + Pressure(dIdx)          * (1 - Dir);

    % Gas Density
    RhoGas = P_Up / (R_Up * T_Up);

    % Composition
    Y_Raw = Y_Matrix(uIdx, :) * Dir + Y_Matrix(dIdx, :) * (1 - Dir);
    Y_Strat = ApplyStratification(Y_Raw);
end

function val = Sigmoid(x, k)
    val = 0.5 + 0.5 * tanh(k * x);
end

function Y_Strat = ApplyStratification(Y)
    % If liquid is present, supress gas from outlet flow
    LiqSum = sum(Y(1:2));
    F = Sigmoid(LiqSum, 100);
    
    Y_Strat = Y;
    Y_Strat(3) = Y_Strat(3) * (1 - F);
    Y_Strat = Y_Strat / sum(Y_Strat);
end

function val = SmoothMin(a, b)
    val = 0.5 * (a + b - sqrt((a - b)^2 + 1e-4));
end