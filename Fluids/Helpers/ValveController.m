function [Cv_OX, Cv_FU] = ValveController(Pc, OF)
    % This function returns the optimal flow coefficient for the main
    % propellant valves given a chamber pressure target and a mass flow
    % ratio.

    PSI_2_Pa = 6895;

    % Tank set pressures (assumed constant for now)
    P_Tank = 550 * PSI_2_Pa;

    % Fluid proprieties
    Rho_OX = 1141;
    Rho_FU = 786;

    % Nozzle and Combustor Proprieties (IRL you replace Cstar with actual
    % measured)
    A_th = 0.0011;
    T_comb = 2666 * 0.83^2;
    R_gas = 442;
    Gamma = 1.207;
    Rc = (2 / (Gamma + 1)) ^ (Gamma / (Gamma - 1));

    % This adjusts the theoretical C* to match the simulation's
    % specific orifice physics. IRL you replace this C_Star value with the
    % real measured Cstar.
    Cd_Nozzle = 0.79; 
    C_Star = sqrt(R_gas * T_comb) / (Cd_Nozzle * sqrt(2 * (1 - Rc)));

    % Calculate Required total mdot
    Pc = Pc * PSI_2_Pa;
    MD_Total = Pc * A_th / C_Star;

    % Mass balance
    MD_OX = MD_Total * (OF / (OF + 1));
    MD_FU = MD_Total * (1 / (OF + 1));

    % Injector Geometry (IRL, replace CdA values with measured)
    InjOX_A  = 3.146e-5; InjOX_Cd = 0.45;
    InjFU_A  = 2.6e-5;   InjFU_Cd = 0.77;
    InjOX_CdA = InjOX_A * InjOX_Cd; InjFU_CdA = InjFU_A * InjFU_Cd;

    % Desired injector pressure drops and manifold pressures
    DP_InjOX = MD_OX^2 / (2 * Rho_OX * InjOX_CdA^2);
    DP_InjFU = MD_FU^2 / (2 * Rho_FU * InjFU_CdA^2);
    P_ManOX = Pc + DP_InjOX;
    P_ManFU = Pc + DP_InjFU;
    
    %% Valve Cv Targets
    % Conversion factor
    Conv = 2.402e-5;
    DP_ValveOX = P_Tank - P_ManOX;
    DP_ValveFU = P_Tank - P_ManFU;

    % Cv Commands
    Cv_OX = MD_OX / (Conv * sqrt(Rho_OX * DP_ValveOX));
    Cv_FU = MD_FU / (Conv * sqrt(Rho_FU * DP_ValveFU));
end