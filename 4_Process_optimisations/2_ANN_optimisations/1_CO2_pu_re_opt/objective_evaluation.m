function [objectives, constraints] = objective_evaluation(x,input,type,net_Pu,net_Re)

%% Design Parameters Input:                                  Index Position
% P_H         = Adsorption Pressure (Pa)                            1
% v_feed      = velocity of the inlet (m/s)                         2
% t_ads       = Time of adsorption step (s)                         3
% alpha_LR    = Light reflux ratio (-)                              4
% epsilon     = Bed void fraction (-)                               5
% e_p         = Pellet void fraction (-)                            6
% r_p         = Radius of the pellet (m)                            7 
%% Material-specific Parameters Input:
% ro_crys     = Crystal density (kg/m3)                             8 
% C_ps        = Specific Cp of adsorbent (J/kg/K)                   9
%% DSL Model Parameters Input:
%************************* CO2 **************************
% q_s1        = Langmuir saturation @ site 1 (mmol/g)               10
% q_s2        = Langmuir saturation @ site 2 (mmol/g)               11
% b01         = Equilibrium parameter (bar^-1)                      12
% b02         = Equilibrium parameter (bar^-1)                      13
% deltaH1     = Heat of Ads @ site 1 (J/mol)                        14
% deltaH2     = Heat of Ads @ site 2 (J/mol)                        15
% ************************* N2 ***************************
% q_s1        = Langmuir saturation @ site 1 (mmol/g)               16
% q_s2        = Langmuir saturation @ site 2 (mmol/g)               17
% b02         = Equilibrium parameter (bar^-1)                      18
% b01         = Equilibrium parameter (bar^-1)                      19
% deltaH1     = Heat of Ads @ site 1/2 (bar)                        20
% deltaH2     = Heat of Ads @ site 1/2 (bar)                        21

%% Initialize objectives and constraints output
switch type
    case 'ProcessEvaluation'
        constraints = [0, 0];
    case 'EconomicEvaluation'
        constraints = [0, 0];
    otherwise
        error('Error. %s is not a recognizable type of operation.' , type) ;
end
objectives  = [0, 0];


%% Defining the Optimisable parameters in the input matrix
input(1) = x(1); % P_H
input(2) = x(2); % t_ads
input(3) = x(3); % v_feed
input(4) = x(4); % alpha_LR 
input(5) = x(5); % epsilon 
input(6) = x(6); % e_p
input(7) = x(7); % r_p

%% Surrogate model predictions
% We have to supply the predictors as transposed arrays, as the neural net 
% function expects the number of variables to be equal to the number of rows
eps = 0.0051;
pu_comp = net_Pu(input');
re_comp = net_Re(input');
purity = ((1+eps).*exp(pu_comp))./(1+exp(pu_comp));
recovery= ((1+eps).*exp(re_comp))./(1+exp(re_comp));

    
% Compile constraint and objective evaluations    
con = recovery - 0.9       ;
if con < 0
    constraints(1) = abs(con) ;
%     constraints(1) = 0;
end

switch type
    case 'ProcessEvaluation'
        objectives(1) = -purity       ;
        objectives(2) = -recovery  ;
        
        con = recovery - 1            ; % outlet must be more pure in CO2 than inlet
        if con > 0
            constraints(2) = abs(con) ;
            constraints(2) = 0        ;
        end
        
    case 'EconomicEvaluation'
        objectives(1) = -productivity       ;
        objectives(2) = energy_requirements  ;
        
        con = recovery - 0.9       ;
        if con < 0
            constraints(1) = abs(con) ;
%             constraints(1) = 0;
        end
        con = purity - 0.9            ;
        if con < 0
            constraints(2) = abs(con) ;
%             constraints(2) = 0 ;
        end
end