% Script File to train an artificial neural network for Purity, Recovery,
% performance indicators of a modified skarstrom cycle.

clear all; close all; clc
parpool('local',16);

%% Load in the data
load('curated_data_linearANDnonlinear.mat');
out = curated_data;

%% Extracting responses and setting up loop variables
% 1: Purity, 2: Recovery
responses_raw = out(:,1:2);                    

% Transdormation 2: This transformation is necessary for the Purity
% and Recovery predictors in order to be coupled with the optimisation
% algorithm; it is to obtain physically meaningful values of the Recovery
% when we try and maximise it (constrains the output to be max 1.001).
eps = 0.0051;
responses_process = log(responses_raw(:,1:2)./((1+eps)-responses_raw(:,1:2)));

% Predictors
predictors_raw = out(:,5:end);

% Setting up the appropriate responses -
responses = responses_process;
predictors = predictors_raw;

%% Train / Test Split
% Use all of the data, define splits in the ANN model architecture
Train = [responses,predictors];
%% ******************** BEGIN THE TRAINING ******************************
% ANN Architecture - 2 hidden Layers, 21 neurons each
hiddenLayerSize = [45, 45, 45];

%% Purity
net_Pu = fitnet(hiddenLayerSize,'trainlm');

% Define the net parameters you want to be fine-tuned
rng(100);
net_Pu.divideParam.trainRatio = 90/100;
net_Pu.divideParam.valRatio = 5/100;
net_Pu.divideParam.testRatio = 5/100;  % to check for overfitting
net_Pu.trainParam.showWindow = false;  % open the nntraintool window
net_Pu.trainParam.epochs = 2500;       % number of training iterations
net_Pu.trainParam.max_fail = 25;       % Default for LM, but 0 for br so in this we impose a stopping criteria for br

% Train the network net using the training data. We have to supply the
% predictors as transposed arrays, as the train function expects the number
% of variables to be = to the number of rows
[net_Pu, tr_Pu] = train(net_Pu,Train(:,3:end)',Train(:,1)','useParallel','yes','showResources','yes');
fprintf('Purity: DONE\n');
%% Recovery
net_Re = fitnet(hiddenLayerSize,'trainlm');

% Define the net parameters you want to be fine-tuned
rng(100);
net_Re.divideParam.trainRatio = 90/100;
net_Re.divideParam.valRatio = 5/100;
net_Re.divideParam.testRatio = 5/100;    % to check for overfitting
net_Re.trainParam.showWindow = false;    % open the nntraintool window
net_Re.trainParam.epochs = 2500;         % number of training iterations
net_Re.trainParam.max_fail = 25;

% Train the network net using the training data. We have to supply the
% predictors as transposed arrays, as the train function expects the number
% of variables to be = to the number of rows
[net_Re, tr_Re] = train(net_Re,Train(:,3:end)',Train(:,2)','useParallel','yes','showResources','yes');
fprintf('Recovery: DONE\n');

%% save required nets and data 
save('net_Pu.mat','net_Pu');
save('net_Re.mat','net_Re');
TrainData = Train(tr_Pu.trainInd,:); dlmwrite('PuRetrainData.txt',TrainData);
TestData = Train(tr_Pu.testInd,:); dlmwrite('PuRetestData.txt',TestData);


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


