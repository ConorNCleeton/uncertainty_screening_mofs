% Script File to train an artificial neural network for
% Energy and Productivity performance indicators of a modified skarstrom
% cycle.

parpool('local',28);        
clear all; close all; clc

%% Load in the data
% supply the appropriate file path to the training data file. 
load('../1_generate_ANN_training_data/3_training_data/total_raw_data.mat');
train_data = curated_data;
out = train_data(:,3:end);

%% Extracting responses and setting up loop variables
% 1: Productivity, 2: Energy
responses_raw = out(:,1:2);                    

% Transformation 1: Transfoming the responses with a log transformation 
% (reduces skewness of the data)
responses_log = log(responses_raw)/log(10);

% Predictors
predictors_raw = out(:,3:end);

% Setting up the appropriate responses -
responses = responses_log;
predictors = predictors_raw;


%% Train / Test Split
% Use all of the data, define splits in the ANN model setup
Train = [responses,predictors];
%% ******************** BEGIN THE TRAINING ******************************
% ANN Architecture - 3 hidden Layers, 45 neurons each
hiddenLayerSize = [45, 45, 45];

%% Productivity
net_Prod = fitnet(hiddenLayerSize,'trainlm');

% Define the net parameters you want to be fine-tuned
rng(100); % Reproducible random partition
net_Prod.divideParam.trainRatio = 90/100;
net_Prod.divideParam.valRatio = 5/100;
net_Prod.divideParam.testRatio = 5/100;   % to check for overfitting
net_Prod.trainParam.showWindow = false;   % open the nntraintool window
net_Prod.trainParam.epochs = 2500;        % number of training iterations
net_Prod.trainParam.max_fail = 20;

% Train the network net using the training data. We have to supply the
% predictors as transposed arrays, as the train function expects the number
% of variables to be = to the number of rows
[net_Prod,tr_Prod] = train(net_Prod,Train(:,3:end)',Train(:,1)','useParallel','yes','showResources','yes');
fprintf('Productivity: DONE\n');
%% Energy
net_En = fitnet(hiddenLayerSize,'trainlm');

% Define the net parameters you want to be fine-tuned
rng(100); % Reproducible random partition
net_En.divideParam.trainRatio = 90/100;
net_En.divideParam.valRatio = 5/100;
net_En.divideParam.testRatio = 5/100;    % to check for overfitting
net_En.trainParam.showWindow = false;    % open the nntraintool window
net_En.trainParam.epochs = 2500;         % number of training iterations
net_En.trainParam.max_fail = 25;        % stop trainining if mse does not reduce

% Train the network net using the training data. We have to supply the
% predictors as transposed arrays, as the train function expects the number
% of variables to be = to the number of rows
[net_En, tr_En] = train(net_En,Train(:,3:end)',Train(:,2)','useParallel','yes','showResources','yes');
fprintf('Energy: DONE\n');

%% save required nets and data 
save('net_Prod.mat','net_Prod');
save('net_En.mat','net_En');
TrainData = Train(tr_En.trainInd,:); dlmwrite('EnProdtrainData.txt',TrainData);
TestData = Train(tr_En.testInd,:); dlmwrite('EnProdtestData.txt',TestData);

%% Ancillary Functions
% Create train/test split
function [inputTrain, inputTest] = train_test_split(input,split)
    rng(10);
    split = cvpartition(size(input,1),'HoldOut',split);
    idx = split.test;
    % Separate to training and test data
    inputTrain = input(~idx,:);
    inputTest  = input(idx,:);
end

% Create k-partitions
function [prtA]=randDivide(data,kfold)
    rng(100);
    n=length(data);
    np=(n-rem(n,kfold))/kfold; % r = rem(a,b) returns the remainder after division of a by b
    B=data;
    [c,idx]=sort(rand(n,1));
    C=data(idx,:);
    i=1;
    j=1;
    ptrA={};
    idxo={};
%     n-mod(n,K)
    while i<n-mod(n,kfold)
        prtA{j}=C(i:i+np-1,:);
        idxo{i}=idx(i:i+np-1,1);
        i=i+np;
        j=j+1;
    end
    prtA{j}=C(i:n,:);
end

% Randomly sample N times from the data file
function [response_out,predictor_out,test_set] = sample(response,predictor, N)
    rng(3);
    idx = randperm(length(response),N);
    response_out = response(idx,:);
    predictor_out = predictor(idx,:);
    filter = ismember(response,response_out,'rows');
    test_set = [response(~filter,:),predictor(~filter,:)];
end


