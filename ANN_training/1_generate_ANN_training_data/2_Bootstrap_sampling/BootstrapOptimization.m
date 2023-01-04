%% Script to run an economic optimisation, maximising the productivity and
% minimising the energy of the PSA process. This method is used to generate 
% training data for ANN model via a bootstrap optimisation

clc; clear; close all;
format long
addpath('CycleSteps\')
addpath('GA_files\')
addpath('utils')
load('DSL_samples_200k.mat'); IsothermParams = DSL_samples;
N = 30 ;
type = 'EconomicEvaluation' ;

%% *********************** FOR GENERATING SAMPLES ************************
% Design Parameters Input: [lower_bound, upper_bound]
P_H         = [1e5, 5e5];       % Adsorption Pressure (Pa)             1
v_feed      = [0.1, 2];         % velocity of the inlet (m/s)          2
t_ads       = [5,500];          % Time of adsorption step (s)          3
alpha_LR    = [0.01, 0.2];      % Light reflux ratio (-)               4

% Material Parameters Input: [lower_bound, upper_bound]
e_b         = [0.35, 0.45];     % Bed void fraction [-]                5
e_p         = [0.3, 0.7];       % Pellet void fraction [-]             6
r_p         = [0.5e-3, 2.5e-3]; % Pellet radius [m]                    7
ro_crys     = [600, 2600];      % Crystal density (kg/m3)              8 
C_ps        = [500, 1500];      % Specific Cp of adsorbent (J/kg/K)    9

% Isotherm Parameters Input
DSL         = [1, 200000];

% Lower bound of decision variables
lower_bound = [P_H(1),v_feed(1),t_ads(1),alpha_LR(1),e_b(1),e_p(1),...
    r_p(1),ro_crys(1),C_ps(1),DSL(1)];

% Upper bound of decision variables
upper_bound = [P_H(2),v_feed(2),t_ads(2),alpha_LR(2),e_b(2),e_p(2),...
    r_p(2),ro_crys(2),C_ps(2),DSL(2)];


%% Economic Optimisation
Function = @(x) PSACycleSimulation( x, IsothermParams, type, N, []) ; % Function to simulate the PSA cycle

options            = nsgaopt();                          % create default options structure
options.popsize    = 150;                                % population size
options.outputfile = 'run1.txt';
options.maxGen     = 45;                                 % max generation
options.vartype    = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ;
options.numObj     = 2 ;                                 % number of objectives
options.numVar     = 10 ;                                % number of design variables
options.numCons    = 3 ;                                 % number of constraints
options.lb         = lower_bound;                        % lower bound of x
options.ub         = upper_bound;                        %  upper bound of x
options.nameObj    = {'-productivity','energy'} ;        % the objective names are showed in GUI window.
options.objfun     = Function                   ;        % objective function handle

options.useParallel = 'yes' ;                            % parallel computation is non-essential here
options.poolsize     = 28   ;                            % number of worker processes

result = nsga2(options)     ;                            % begin the optimization!


