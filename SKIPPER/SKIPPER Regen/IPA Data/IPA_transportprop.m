function [mu, k, rho, cp, boiling_temp, state] = IPA_properties(T,P)

    

    %% Constants     
    alpha = 5.4574*10^(-9); % Viscosity Pressure correlation coeff
    C = 2.9531289036*10^(-9); % Thermal conductivity correlation coeff

    %% Import Tables
    pressurearray = [1,1.5,2,3,4,5,7,10,15,20,30,40,50,70,100,150,200,300,400,500,700,1000,1500,2000,3000,4000,5000,7000,10000,15000,20000,30000,40000,50000]; % pressures in dataset (singlephase) [kpa]
    numVars = 35;
    % boiling point table
    opts = detectImportOptions("/bin/ipadata.xlsx");
    opts.Sheet = 'SLVM';
    SLVM = readtable("/bin/ipadata.xlsx",opts)
    boiling_point_table = table2array([SLVM(2:end,2), SLVM(2:end,1)])
    
    % transport properties (k & mew) tables
    opts = detectImportOptions("/bin/ipadata.xlsx");
    preview("/bin/ipadata.xlsx",opts)
    opts.Sheet = 'transport';
    %opts.SelectedVariableNames = [1:5]; 
    %opts.DataRange = '2:11';
    properties = readtable(pwd + "/bin/ipadata.xlsx",opts);
    visc_table = str2double(table2array([properties(2:35,2) properties(2:35,3)]));
    k_table = str2double(table2array([properties(2:end,2) properties(2:end,4)]));

    % density table
    opts = detectImportOptions("/bin/ipadata.xlsx");
    preview("/bin/ipadata.xlsx",opts);
    opts.Sheet = 'Density';
    opts = spreadsheetImportOptions('NumVariables',numVars);
    opts.DataRange = "A2:AI64"; 
    opts.Sheet = 'Density';
    density_import = readtable(pwd + "/bin/ipadata.xlsx",opts);
    density_table = str2double(table2array(density_import));

    % Cp table 
    opts.Sheet = 'Cp';
    cp_import = readtable(pwd + "/bin/ipadata.xlsx",opts);
    cp_table = str2double(table2array(cp_import));

    %% Interpolate Tables 
    boiling_temp = interp1(str2double(boiling_point_table(:,1)), str2double(boiling_point_table(:,2)), P/1000, 'linear','extrap'); % Boiling temperature [Kelvin]
    rho = interp2(pressurearray,density_table(:,1),density_table(:,2:end),P/1000,T); % Density [kg/m^3]
    cp = interp2(pressurearray,cp_table(:,1),cp_table(:,2:end),P/1000,T); %Specific heat [kJ/Kg]
    mu_data = interp1(visc_table(:,1), visc_table(:,2), T, 'cubic','extrap'); % Viscosity from dataset [Pa*s]
    k_data = interp1(k_table(:,1), k_table(:,2), T, 'linear', 'extrap');  % Thermal conductivity from dataset [W/mk]

    %% Determine State 
    if T >= boiling_temp 
        state = "boiling"
    else 
        state = "liquid"
    end

    %% Calculate Curve Fits

    % Viscosity temperature Curve fit
    if (T - 273.15) < 36.5
        mu_fit = 4.5054*exp(-.031*(T-273.15));
    else
        mu_fit = 3.3724*exp(-.023*(T-273.15));
    end
    %visc_tp =  visc_temp*exp(alpha*(P-101325));
    mu_fit = visc_tp /1000; % Viscosity from curve fits [Pa*s}

    % Viscosity Average between data and fit
    if isnan(mu_data) 
        mu_avg = mu_fit;
    else
        mu_avg = (mu_fit + mu_data)./2;
    end

    % Pressure curve fits
    mu = mu_avg*exp(alpha*(P-101325)); % Viscosity [Pa*s]
    k = k_data*(1+C*(P-101325)); % Thermal Conductivity [W/mK]
end