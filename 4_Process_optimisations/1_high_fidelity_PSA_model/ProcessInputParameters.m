function InputParams = ProcessInputParameters(process_vars, iso_params, material, N)
%InputParameters: Specify input parameters for PSA simulation
%   Number of finite volumes is always an input from the calling script/
%   function. In addition, x is used to allow varying specific parameters 
%   (times, pressures, heats of adsorptions, etc.) from the calling script/
%   function. It is not necessary for simple simulations.
    
%% State all input parameters for the simulation
    
    % Retrieve required process variables provided as inputs of the function
    P_0         = process_vars(1)       ;   % Adsorption pressure [Pa]
    ndot_0      = process_vars(2)       ;   % Inlet molar flux [mol/s/m^2]
    t_ads       = process_vars(3)       ;   % Time of adsorption step [s]
    alpha       = process_vars(4)       ;   % Light product reflux ratio [-]
    epsilon     = process_vars(5)       ;   % Void fraction
    e_p         = process_vars(6)       ;   % porosity of the pellet [-]
    r_p         = process_vars(7)       ;   % radius of pellets

    % Retrieving Params that depend on the adsorbent material
    IsothermPar        = iso_params    ; 
    
    % Operating bed parameters
    t_pres      = 20                    ;   % Maximum/time of pressurization step [s]
    t_CnCdepres = 15                    ;   % Maximum/time of depressurization step [s]
    t_CoCdepres = 70                    ;   % Maximum/time of depressurization step [s]
    t_LR        = t_ads                 ;   % Time of light reflux step [s]
    t_HR        = t_LR                  ;   % Time of heavy reflux step [s]
    tau         = 0.5                   ;   % Parameter used for determining speed of pressure change
    P_inlet     = 1.02                  ;   % Pressure of feed gas at the inlet of the adsorption step
    P_I         = 0.1e5                 ;   % Intermediate pressure [Pa] (Not used in modified skarstrom cycle)
    P_l         = 0.1e5                 ;   % Purge Pressure [Pa]
    beta        = 1                     ;   % Heavy product reflux ratio [-]

    % Flue gas parameters and constants
    R          = 8.314                  ;   % Universal gas constant [J/mol/K : Pa*m^3/mol/K]
    T_0        = 313.15                 ;   % Feed temperature of flue gas [K]
    y_0        = 0.15                   ;   % Inlet gas CO2 mole fraction[-]
    Ctot_0     = P_0/R/T_0              ;   % Inlet total concentration [mol/m^3]
    v_0        = ndot_0/Ctot_0          ;   % Inlet velocity and scaling parameter [m/s]
    mu         = 1.72e-5                ;   % Viscosity of gas [Pa*s]
    D_m        = 1.2995e-5              ;   % Molecular diffusivity [m^2/s]
    K_z        = 0.09                   ;   % Thermal conduction in gas phase [W/m/k]
    C_pg       = 30.7                   ;   % Specific heat of gas [J/mol/k]
    C_pa       = 30.7                   ;   % Specific heat of adsorbed phase [J/mol/k]
    MW_CO2     = 0.04402                ;   % Molecular weight of CO2 [kg/mol]
    MW_N2      = 0.02802                ;   % Molecular weight of N2 [kg/mol]
	%feed_gas  = 'Constant Pressure'    ;   % Whether flue gas during the feed step has a constant pressure or velocity
    feed_gas   = 'Constant Velocity'    ;   % Whether flue gas during the feed step has a constant pressure or velocity
    
    % Adsorbent parameters
    ro_crys     = material(1)           ;   % crystal density [kg/m^3]
    C_ps        = material(2)           ;   % Specific heat capacity of the adsorbent [J/kg/K]
    ro_s        = ro_crys*(1-e_p)       ;   % Density of the adsorbent [kg/m^3]
    tau_s       = 3                     ;   % Tortuosity factor [-]
    q_s         = 5.84                  ;   % Molar loading scaling factor [mol/kg]
    q_s0        = q_s*ro_s              ;   % Molar loading scaling factor [mol/m^3]
    
    % Isotherm parameters
    q_s_b      = [IsothermPar(1),  IsothermPar(7)]   ;   % Saturation loading on site b [mol/kg]
    q_s_d      = [IsothermPar(2),  IsothermPar(8)]   ;   % Saturation loading on site d [mol/kg]
    b          = [IsothermPar(3),  IsothermPar(9)]   ;   % Pre-exponential factor for site b [Pa-1]
    d          = [IsothermPar(4),  IsothermPar(10)]  ;   % Pre-exponential factor for site d [Pa-1]
    deltaU_b   = [IsothermPar(5),  IsothermPar(11)]  ;   % Heat of adsorption for site b [J/mol]
    deltaU_d   = [IsothermPar(6),  IsothermPar(12)]  ;   % Heat of adsorption for site d [J/mol]
    
    KH_b_1 =q_s_b(1)*b(1)*exp(deltaU_b(1)./(8.314*T_0));
    KH_d_1 =q_s_d(1)*d(1)*exp(deltaU_d(1)./(8.314*T_0));
    KH_b_2 =q_s_b(2)*b(2)*exp(deltaU_b(2)./(8.314*T_0));
    KH_d_2 =q_s_d(2)*d(2)*exp(deltaU_d(2)./(8.314*T_0));
    deltaU_1 = (deltaU_b(1)*KH_b_1 + deltaU_d(1)*KH_d_1)/(KH_b_1+KH_d_1); % Isosteric Heat of adsorption (component 1)
    deltaU_2 = (deltaU_b(2)*KH_b_2 + deltaU_d(2)*KH_d_2)/(KH_b_2+KH_d_2); % Isosteric Heat of adsorption (component 2)
    deltaU = [deltaU_1, deltaU_2];  
    IsothermParams = [q_s_b, q_s_d, b, d, deltaU_b, deltaU_d, IsothermPar(13)] ;

    % Column properties
    L          = 1        ;   % Length of the column [m]
    radius     = 0.1445   ;   % radius of the column [m] 

    % Linear Driving force approximation
    k_p         = (15*e_p*D_m)/(tau_s*(r_p^2));                  % pellet LDF coefficient
    q_in        = Isotherm(y_0, P_0, T_0, IsothermParams) ;
    k_CO2_LDF   = (k_p*y_0*P_0)/(q_in(1)*ro_s.*8.314.*T_0) ;     % Mass transfer coefficient for CO2 [1/s]
    k_N2_LDF    = (k_p*(1-y_0)*P_0)/(q_in(2)*ro_s.*8.314.*T_0) ; % Mass transfer coefficient for N2 [1/s]
    
%% Distribute the values to the necessary variables
    Params     = zeros(42, 1) ;
    Params(1)  = N			  ;
    Params(2)  = deltaU(1)    ;
    Params(3)  = deltaU(2)    ;
    Params(4)  = ro_s		  ;
    Params(5)  = T_0		  ;
    Params(6)  = epsilon	  ;
    Params(7)  = r_p		  ;
    Params(8)  = mu			  ;
    Params(9)  = R			  ;
    Params(10) = v_0		  ;
    Params(11) = q_s0		  ;
    Params(12) = C_pg		  ;
    Params(13) = C_pa		  ;
    Params(14) = C_ps		  ;
    Params(15) = D_m		  ;
    Params(16) = K_z		  ;
    Params(17) = P_0		  ;
    Params(18) = L			  ;
    Params(19) = MW_CO2		  ;
    Params(20) = MW_N2		  ;
    Params(21) = k_CO2_LDF	  ;
    Params(22) = k_N2_LDF	  ;
    Params(23) = y_0		  ;
    Params(24) = tau		  ;
    Params(25) = P_l		  ;
    Params(26) = P_inlet	  ;
    Params(27) = 1			  ;   % Place for y at outlet of Adsorption = y at inlet of Light Reflux: y_LP
                                  % y_LR = 1 - No initial guess for inlet CO2 mole fraction in Light Reflux step
    Params(28) = 1			  ;   % Place for T at outlet of Adsorption = T at inlet of Light Reflux: T_LP
                                  % T_LR = 1 - No initial guess for inlet temperature in Light Reflux step
    Params(29) = 1			  ;   % Place for ndot at outlet of Adsorption = ndot at inlet of Light Reflux
                                  % ndot_LR = 1 - No initial guess for inlet ndotin Light Reflux step
    Params(30) = alpha    	  ;
    Params(31) = beta         ;
    Params(32) = P_I          ;
    Params(33) = y_0          ;   % Place for y at outlet of CnC depressurization = y at inlet of of Heavy Reflux: y_HP
                                  % y_HR = y_0 - Initial guess for inlet CO2 mole fraction in Heavy Reflux step
    Params(34) = T_0          ;   % Place for T at outlet of CnC depressurization = T at inlet of of Heavy Reflux: T_HP
                                  % T_HR = T_0 - Initial guess for inlet temperature in Heavy Reflux step
    Params(35) = ndot_0*beta  ; % 300/30   % Place for ndot at outlet of CnC depressurization = ndot at inlet of Heavy Reflux
                                  % ndot_HR = ndot_0*beta - Initial guess for inlet ndotin Heavy Reflux step
    Params(36) = 0.01    	  ;   % Place for y at outlet of Adsorption = y at inlet of CnC pressurization: y_LP
                                  % y_LR = 0.01 - Initial guess for inlet CO2 mole fraction in CnC pressurization step
    Params(37) = T_0    	  ;   % Place for T at outlet of Adsorption = T at inlet of CnC pressurization: T_LP
                                  % T_LR = T_0 - Initial guess for inlet temperature in CnC pressurization step
    Params(38) = ndot_0  	  ;   % Place for ndot at outlet of Adsorption = ndot at inlet of CnC pressurization 
                                  % NOTE: not used, seems not necessary. ndot_LR = ndot_0 - Initial guess for inlet ndot
                                  % in CnC pressurization step
    Params(39) = radius       ;   % Radius of the column [m]
    Params(40) = e_p          ;   % porosity of the pellet [-]
    Params(41) = tau_s        ;   % tortuosity factor [-] 
    
    
    if strcmpi(feed_gas, 'Constant Pressure') == 1
        Params(end) = 1 ;
    elseif strcmpi(feed_gas, 'Constant Velocity') == 1
        Params(end) = 0 ;
    else
        error('Please specify whether inlet velocity or pressure is constant for the feed step')
    end
	
    Times          = [ t_pres; t_ads; t_CnCdepres; t_LR; t_CoCdepres; t_HR ] ;

%% Economic Parameters (depricated - not used in the current model).
    
    desired_flow                          = 100          ;   % Desired flow rate of flue gas per column [mol/s]
    electricity_cost                      = 0.07         ;   % Cost of electricity [$/kWh]
    hour_to_year_conversion               = 8000         ;   % total hours in a year, the remaining time is assumed to be down for maitenance  [hr/year]
    life_span_equipment                   = 20           ;   % Life Span of all equuipment besides adsorbent [years]
    life_span_adsorbent                   = 5            ;   % Life Span of adsorbent [years]
    CEPCI                                 = 536.4        ;   % CEPCI of present month year (Jan 2016).
    
    % change this according to the cycle to be simulated
    cycle_time = t_pres + t_ads + t_HR + t_CnCdepres + t_LR ;   % total time required for 1 cycle [s]
    
    EconomicParams    = zeros(6, 1)              ;
    EconomicParams(1) = desired_flow             ;
    EconomicParams(2) = electricity_cost         ;
    EconomicParams(3) = cycle_time               ;
    EconomicParams(4) = hour_to_year_conversion  ;
    EconomicParams(5) = life_span_equipment      ;
    EconomicParams(6) = life_span_adsorbent      ;
    EconomicParams(7) = CEPCI                    ;
%   
%% Combine all lists into one variable of cells that can easily be passed
    InputParams{1} = Params          ;
    InputParams{2} = IsothermParams  ;
    InputParams{3} = Times           ;
    InputParams{4} = EconomicParams  ;
%     InputParams{5} = economic_class  ;
%   
end 