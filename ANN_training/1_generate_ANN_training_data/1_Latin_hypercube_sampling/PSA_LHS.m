%% Script: generates training data for ANN model using high fidelity PSA 
% process model. Latin hypercube sampling is used to sample the input 
% phase space. 

clc; clear; close all;
parpool('local',28);
format long
addpath('CycleSteps')
addpath('utils')

N = 30 ;
type = 'ProcessEvaluation' ;

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

% Lower bound of decision variables
lower_bound = [P_H(1),v_feed(1),t_ads(1),alpha_LR(1),e_b(1),e_p(1),...
    r_p(1),ro_crys(1),C_ps(1)];

% Upper bound of decision variables
upper_bound = [P_H(2),v_feed(2),t_ads(2),alpha_LR(2),e_b(2),e_p(2),...
    r_p(2),ro_crys(2),C_ps(2)];

% LHS sampling
x = lhsdesign_mod(100000,lower_bound,upper_bound);
%% ***********************************************************************


% Extracting the appropriate DSL model parameters
load('LHS_samples_200k.mat');
IsothermParameters = [DSL_samples,zeros(size(DSL_samples,1),1)];


%% Evaluation of smaller subsets 
run = 1;
if run == 1
    x = x(1:10000,:);
    IsothermParameters = IsothermParameters(1:10000,:);
else
    x = x((run-1)*10000 + 1:run*10000,:);
    IsothermParameters = IsothermParameters((run-1)*10000 + 1:run*10000,:);
end

% For Printing to output for machine learning training purposes
to_print = [x,IsothermParameters(:,1:12)];
[samples,vars] = size(to_print);

% Input matrix (Same order as decision variables, except v_feed converted
% to n_dot0 using P_H).
input_to_sim = [x(:,1),x(:,2).*x(:,1)/8.314/313.15,x(:,3),x(:,4),x(:,5),...
    x(:,6),x(:,7),x(:,8),x(:,9)];

% Opening a parallel pool constant object to create a temporary file on
% each worker within the parallel loop, which will then be cleaned up and
% concatenated into a single output file after the operation
c = parallel.pool.Constant(@() fopen(tempname(pwd),'wt'),@fclose);

% Each worker can operate on a different data set or different portion of 
% distributed data, and can communicate with other participating workers
% while performing the parallel computations using spmd 
spmd
    A=(fopen(c.Value));
end
parfor i = 1:samples 
    [purity,recovery,prod,energy,cycle_iter] = PSACycle(input_to_sim(i,:), IsothermParameters(i,:), [], type, N);

    % Printing the outputs first
    fprintf(c.Value, '%9.6f\t%9.6f\t%9.8f\t%9.3f\t',purity,recovery,prod,energy,cycle_iter);
    
    % Then printing the inputs
    for j = 1:11
        fprintf(c.Value, '%9.6f\t',to_print(i,j));
    end
    for j = 12:13
        fprintf(c.Value, '%9.14f\t',to_print(i,j));
    end
    for j = 14:15
        fprintf(c.Value, '%9.4f\t',to_print(i,j));
    end
    for j = 16:17
        fprintf(c.Value, '%9.6f\t', to_print(i,j));
    end
    for j = 18:19
        fprintf(c.Value, '%9.14f\t',to_print(i,j));
    end
    fprintf(c.Value, '%9.4f\t%9.4f\n',to_print(i,20),to_print(i,21));
    
end
% Clear the object c 
clear c;
%%
% Moving temporary files to usable file names
for i = 1: length(A)
    a= char(A(1,i)); 
    NAMES(i,:) = a(1, length(pwd)+2:length(a));  
    movefile(NAMES(i,:),sprintf('%d.txt',i),'f');
end

% concatenating the separate output files into a master file 
fileout=sprintf('output%d.txt',run);
fout=fopen(fileout,'w');
for cntfiles=1:length(A)
  fin=fopen(sprintf('%d.txt',cntfiles));
  while ~feof(fin)
    fprintf(fout,'%s \n',fgetl(fin));
  end
  fclose(fin);
end
fclose(fout);

% delete the unnecesary files by the following loop
fclose('all');
for i = 1:length(A)
    delete(sprintf('%d.txt',i))
end
delete(gcp('nocreate'));


