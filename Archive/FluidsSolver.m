function XDOT = FluidsSolver(X, System)
    %% Solver Constants & Dimensions
    PLength = size(System.Nodes, 2);
    MALength = size(System.Links.Algebraic, 2);
    MDLength = size(System.Links.Dynamic, 2);

    % State Indexes
    YIdx_Start = PLength + MDLength + 1;

    % Initialize Vectors
    XDOT = zeros(size(X));
    Pressure = zeros(PLength, 1);
    Massflow = zeros(MALength + MDLength, 1);

    % Extract Densities
    RhoArray = System.Constants.RhoArray;
    NumSpecies = length(RhoArray);
    RhoNode = zeros(PLength, 1);
    LinkMap = System.LinkMap;

    %% Calculate Node Densities
    Y_All = X(YIdx_Start:end);
    Y_Matrix = reshape(Y_All, [NumSpecies, PLength])';
    for i = 1:PLength
        Y_Vec = Y_Matrix(i, :);
        Y_Vec = max(0, Y_Vec);
        if sum(Y_Vec) > 0
            Y_Vec = Y_Vec / sum(Y_Vec);
        end
        Y_Matrix(i, :) = Y_Vec;

        % Calculate mix density
        InvRho = sum(Y_Vec ./ RhoArray);
        RhoNode(i) = 1 / InvRho;
    end

    %% Solve the dynamic massflow states
    for i = 1:MDLength
        % Link proprieties
        UpIDX = System.Links.Dynamic(i).Up;
        DwIDX = System.Links.Dynamic(i).Down;
        DeltaP = X(UpIDX) - X(DwIDX);
        A = System.Links.Dynamic(i).A;
        L = System.Links.Dynamic(i).L;
        Zeta = System.Links.Dynamic(i).Zeta;
        ID = System.Links.Dynamic(i).ID;

        % Pick Link Density depending on flow direction
        if X(PLength + i) >= 0
            Rho = RhoNode(UpIDX);
        else
            Rho = RhoNode(DwIDX);
        end

        % ODE
        XDOT(PLength + i) = A / L * DeltaP - Zeta / (2 * Rho * A * L) * abs(X(PLength + i)) * X(PLength + i);
        Massflow(ID) = X(PLength + i);
    end

    %% Solve the algebraic massflow states
    for i = 1:MALength
        % Throttle Valve Massflow Equation
        if strcmp(System.Links.Algebraic(i).Type, 'Throttle')
            % Link Proprieties
            UpIDX = System.Links.Algebraic(i).Up;
            DwIDX = System.Links.Algebraic(i).Down;
            ID = System.Links.Algebraic(i).ID;
            Cv = System.Links.Algebraic(i).Cv;
            DeltaP = X(UpIDX) - X(DwIDX);

            % Pick Link Density depending on flow direction
            if DeltaP >= 0
                Rho = RhoNode(UpIDX);
            else
                Rho = RhoNode(DwIDX);
            end

            % Equation
            Massflow(ID) = 0.0076 * Cv * sqrt(abs(DeltaP * Rho)) * sign(DeltaP)...
            * System.Links.Algebraic(i).Active;
        end
    end
    
    %% Solve the pressure states and the Prop Transport
    XDOT_Species = zeros(PLength, NumSpecies);
    for i = 1:PLength
        % Node Proprieties
        LinksIN = System.Nodes(i).LinksIN;
        LinksOUT = System.Nodes(i).LinksOUT;
        Mdot_IN = Massflow(LinksIN);
        Mdot_OUT = Massflow(LinksOUT);
        C = System.Constants.C;
        V = System.Nodes(i).V;

        % ODE Pressure
        if ~System.Nodes(i).Fixed
            XDOT(i) = C^2 / V * (sum(Mdot_IN) - sum(Mdot_OUT));
        else
            XDOT(i) = 0;
        end
        ID = System.Nodes(i).ID;
        Pressure(ID) = X(i);

        % ODE Transport 
        if ~System.Nodes(i).Fixed
            FluxIN = zeros(1, NumSpecies);
            Mass = System.Nodes(i).V * RhoNode(i);
            TotalMDOT_IN = 0;
    
            % Inflows
            Mask = Mdot_IN > 0;
            if any(Mask)
                Flows = Mdot_IN(Mask);
                IDs = LinksIN(Mask);
                Neighbors = LinkMap(IDs, 1);
                TotalMDOT_IN = TotalMDOT_IN + sum(Flows);
                FluxIN = FluxIN + Flows' * Y_Matrix(Neighbors, :);
            end

            % Backflows
            Mask = Mdot_OUT < 0;
            if any(Mask)
                Flows = abs(Mdot_OUT(Mask));
                IDs = LinksOUT(Mask);
                Neighbors = LinkMap(IDs, 2);
                TotalMDOT_IN = TotalMDOT_IN + sum(Flows);
                FluxIN = FluxIN + Flows' * Y_Matrix(Neighbors, :);
            end

            % Transport Derivative
            XDOT_Species(i, :) = (FluxIN - TotalMDOT_IN * Y_Matrix(i, :)) / Mass;
        else
            XDOT_Species(i, :) = 0;
        end
    end

    %% Package State Vector
    XDOT(YIdx_Start:end) = reshape(XDOT_Species', [], 1);
end