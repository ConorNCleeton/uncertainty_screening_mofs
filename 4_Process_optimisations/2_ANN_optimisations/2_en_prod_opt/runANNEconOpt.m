% Run a Energy-Productivity optimisation for CRAFTED MOF materials which meet
% the US DOE constraints consistently 
clear all; close all; clc

% Load in path dependencies and .mat files
addpath('GA_files'); 
addpath('sorted_material_dataframes\'); % csv files containing the material input data
load('net_Pu_logb0.mat');    % Purity ANN 
load('net_Re_logb0.mat');    % Recovery ANN
load('net_Prod_logb0.mat');  % Productivity ANN
load('net_En_logb0.mat');    % Energy ANN 
load('CL1_MOFs.mat');        % consistent DOE materials

% Load in material property dataframes for each Forcefield (FF_dfs)
myFolder = sprintf('%s/sorted_material_dataframes/',pwd);
filePattern = fullfile(myFolder, '*.csv');
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    FF_dfs{k} = readtable(fullFileName);
end
type = 'EconomicEvaluation';

% Define the ANN inputs: optimisable parameters (placeholders)
P_H         = [1];              % Adsorption Pressure (Pa)             1
v_feed      = [1];              % velocity of the inlet (m/s)          2
t_ads       = [1];              % Time of adsorption step (s)          3
alpha_LR    = [1];              % Light reflux ratio (-)               4
epsilon     = [1];              % Void fraction of bed (-)             5
e_p         = [1];              % Void fraction of the pellet (-)      6
r_p         = [1];              % radius of the pellet (m)             7


% DSL Model Parameters Input: [lower_bound, upper_bound]
%************************* CO2 **************************
% q_s1CO2     = [0, 12.5];         Langmuir saturation @ site 1 (mmol/g) 10
% q_s2CO2     = [0, 12.5];         Langmuir saturation @ site 2 (mmol/g) 11
% b01         = [1e-9, 1e-2];      Equilibrium parameter (bar^-1)        12
% b02         = [1e-9, 1e-2];      Equilibrium parameter (bar^-1)        13
% deltaH1     = [15000, 45000];    Heat of Ads @ site 1 (J/mol)          14
% deltaH2     = [15000, 45000];    Heat of Ads @ site 2 (J/mol)          15
% % ************************* N2 ***************************
% q_s1N2      = [0, 12.5];         Langmuir saturation @ site 1 (mmol/g) 16
% q_s2N2      = [0, 12.5];         Langmuir saturation @ site 2 (mmol/g) 17
% b01          = [1e-9, 1e-2];     Equilibrium parameter (bar^-1)        18
% b02          = [1e-9, 1e-2];     Equilibrium parameter (bar^-1)        19
% deltaH1     = [5000, 25000];     Heat of Ads @ site 1/2 (J/mol)        20
% deltaH2     = [5000, 25000];     Heat of Ads @ site 1/2 (J/mol)        21
 

for material = 2%:length(alwaysDOE)
    iso_param_df = [];
    material_name = alwaysDOE{material,1};
    disp(material_name)
    objectivesEcon = cell(length(FF_dfs),5);
    variablesEcon  = cell(length(FF_dfs),5);
    for FF = 1:length(FF_dfs)
        data = FF_dfs{FF};
        material_data = table2array(data(strcmp(data.material,material_name),2:15));
        FF_name       = data(strcmp(data.material,material_name),:).Forcefield{1};
        disp(FF_name)
        iso_param_df = [iso_param_df;material_data(3:end)];
        for iter = 1
            rng(iter*10); % different population each iteration
            
            % Define the ANN inputs: material-specific properties
            ro_crys     = max(min(material_data(1),2600),600);      % Crystal density (kg/m3)       8
            C_ps        = max(min(material_data(2),1500),500);      % Specific Cp of adsorbent (J/kg/K)    9
            IsothermParams     = material_data(3:end);

            % ********** for net_logb0 ANN models ************************
            IsothermParams(3:4)  = log(IsothermParams(3:4));
            IsothermParams(9:10) = log(IsothermParams(9:10));
            % ************************************************************

            % Gather the inputs
            inputs =[P_H,v_feed,t_ads,alpha_LR,epsilon,e_p,r_p,...
                ro_crys,C_ps,IsothermParams];
            
            if ismember(0,inputs) || ismember(NaN,inputs)
                objectivesEcon{FF,iter} = zeros(2,70);
                variablesEcon{FF,iter} = zeros(7,70);
            else
                
                %% Perform Process Optimisation
                Function = @(x) objective_evaluation(x,inputs,type,net_Pu,net_Re,net_Prod,net_En) ; % Function to surrogate model
                
                options         = nsgaopt() ;                            % create default options structure
                options.popsize = 70        ;                            % population size
                options.maxGen  = 250        ;                            % max generation
                
                options.vartype    = [1, 1, 1, 1, 1, 1, 1]         ;
                options.outputfile = sprintf('%s_output.txt',material_name) ;
                options.numObj  = 2 ;                                    % number of objectives
                options.numVar  = 7 ;                                    % number of design variables
                options.numCons = 2 ;                                    % number of constraints
                
                % x = [P_H, v_feed, t_feed, alpha_LR, epsilon, e_p, r_p];
                options.lb      = [1e5, 0.1, 5, 0.01, 0.35, 0.3, 0.5e-3];     % lower bound of x
                options.ub      = [5e5, 2, 500, 0.2, 0.45, 0.7, 2.5e-3];      % upper bound of x
                options.nameObj = {'-purity','-recovery'} ;               % the objective names are showed in GUI window.
                options.objfun  = Function               ;                % objective function handle
%                 options.useParallel = 'yes' ;                             % parallel computation is non-essential here
%                 options.poolsize     = 4   ;                              % number of worker processes
                result = nsga2(options)     ;                             % begin the optimization!
                
                [obj, vars] = sortt(loadpopfile(sprintf('%s_output.txt',material_name)));
                objectivesEcon{FF,iter} = obj;
                variablesEcon{FF,iter} = vars;
                
                delete(sprintf('%s_output.txt',material_name));
                % delete .mat file and then resave with every iteration. Just in case
                % the for loop terminates with some error and the information from
                % previous iterations are not lost.
                if ~exist('opt_results', 'dir')
                    mkdir('opt_results')
                end
                savetofile(objectivesEcon,sprintf('%s/opt_results/%s_objectivesEcon.mat',pwd,material_name))
                savetofile(variablesEcon,sprintf('%s/opt_results/%s_variablesEcon.mat',pwd,material_name))
            end
        end
    end
    visualise_results(objectivesEcon,variablesEcon,material_name,iso_param_df,[1,2,4,5,6,8]);
end

function savetofile(data,fullfilename)
    save(fullfilename,'data');
end

