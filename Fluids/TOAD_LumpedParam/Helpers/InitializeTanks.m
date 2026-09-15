function [Y0_LOX, Y0_IPA] = InitializeTanks(TankVolOX, TankVolIPA, StartP)
    %% Calculate Initial Tank Conditions
    % Fluid Properties [kg/m^3]
    rho_OX = 1141; % LOX
    rho_FU = 786;  % IPA

    % Constants
    R_N2_init = 296.8; 
    T_amb_init = 293;
    P_init_Pa = StartP * 6895; 

    % Tank Geometry
    TargetUllage = 0.10; 

    % Initial N2 Density (1 atm)
    Rho_Gas_Init = P_init_Pa / (R_N2_init * T_amb_init);

    %% LOX TANK CALCULATIONS
    Vol_Gas_OX = TankVolOX * TargetUllage;
    Vol_Liq_OX = TankVolOX * (1 - TargetUllage);
    
    Mass_Gas_OX = Vol_Gas_OX * Rho_Gas_Init;
    Mass_Liq_OX = Vol_Liq_OX * rho_OX;
    Mass_Total_OX = Mass_Gas_OX + Mass_Liq_OX;
    
    Y_OX_Liq = Mass_Liq_OX / Mass_Total_OX;
    Y_OX_N2  = Mass_Gas_OX / Mass_Total_OX;
    Y0_LOX = [Y_OX_Liq, 0, Y_OX_N2]; % [OX, FU, N2]

    %% IPA TANK CALCULATIONS
    Vol_Gas_FU = TankVolIPA * TargetUllage;
    Vol_Liq_FU = TankVolIPA * (1 - TargetUllage);
    
    Mass_Gas_FU = Vol_Gas_FU * Rho_Gas_Init;
    Mass_Liq_FU = Vol_Liq_FU * rho_FU;
    Mass_Total_FU = Mass_Gas_FU + Mass_Liq_FU;
    
    Y_FU_Liq = Mass_Liq_FU / Mass_Total_FU;
    Y_FU_N2  = Mass_Gas_FU / Mass_Total_FU;
    Y0_IPA = [0, Y_FU_Liq, Y_FU_N2]; % [OX, FU, N2]

    %% Print Results
    fprintf('Init LOX Tank: %.2f kg (Y_N2: %.6f)\n', Mass_Total_OX, Y_OX_N2);
    fprintf('Init IPA Tank: %.2f kg (Y_N2: %.6f)\n', Mass_Total_FU, Y_FU_N2);
end