% Plot performance of ANN model using training and testing data
close all; clear all; clc;

% Load in the trained ANN models
addpath('trained_ANN_models\');
load('net_En_logb0.mat'); load('net_Prod_logb0.mat');
load('net_Pu_logb0.mat'); load('net_Re_logb0.mat');

% Load in the testing data (supply the appropriate file path as needed)
test_Proc = readmatrix('../1_generate_ANN_training_data/3_training_data/testDataProc_transformed.txt'); % Test data for purity / recovery
test_Econ = readmatrix('../1_generate_ANN_training_data/3_training_data/testDataEcon_transformed.txt'); % Test data for energy / prod

% Load in the full dataset (i.e. the training and testing data)
% NB: this is the raw data before applying log-like tansformations to the
% target KPIs (see supplementary note 4 of publication).
load('../1_generate_ANN_training_data/3_training_data/total_raw_data.mat');

% we train the ANN models on the logarithm of the b0 parameters. Here, we
% are just converting the raw data to the appropriate format for the ANN
% model input
Proc_data = [curated_data(:,1:2),curated_data(:,5:end)]; Proc_data(:,14:15) = log(Proc_data(:,14:15)); Proc_data(:,20:21) = log(Proc_data(:,20:21)); 
Econ_data = [curated_data(:,3:4),curated_data(:,5:end)]; Econ_data(:,14:15) = log(Econ_data(:,14:15)); Econ_data(:,20:21) = log(Econ_data(:,20:21));
eps = 0.0051;

% Create a tiledlayout
figure('Color','w');
set(gcf,'Position',[2000 100 500 400])
t = tiledlayout('flow','TileSpacing','tight');
nexttile;
% predict the purity using the training and testing data
scatter(((1+eps).*exp(net_Pu(Proc_data(:,3:end)')))./(1+exp(net_Pu(Proc_data(:,3:end)'))),Proc_data(:,1),2,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerEdgeAlpha',0.2,'MarkerFaceAlpha',.2); hold on;
scatter(((1+eps).*exp(net_Pu(test_Proc(:,3:end)')))./(1+exp(net_Pu(test_Proc(:,3:end)'))),((1+eps).*exp(test_Proc(:,1)))./(1+exp(test_Proc(:,1))),2,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerEdgeAlpha',0.35,'MarkerFaceAlpha',.35); hold on;
plot([0,1],[0,1],'LineWidth',1,'Color','k');
lin = fitlm(((1+eps).*exp(net_Pu(test_Proc(:,3:end)')))./(1+exp(net_Pu(test_Proc(:,3:end)'))),((1+eps).*exp(test_Proc(:,1)))./(1+exp(test_Proc(:,1)))); R2_Pu = lin.Rsquared.Adjusted;
axis([0 1 0 1]);
text(0.55,0.1,sprintf('R^2_{adj} = %1.3f',R2_Pu),'Units','Normalized');
text(0.03,0.94,'Purity','Units','Normalized','FontWeight','bold')
box on; 

nexttile;
% predict the recovery using the training and testing data
scatter(((1+eps).*exp(net_Re(Proc_data(:,3:end)')))./(1+exp(net_Re(Proc_data(:,3:end)'))),Proc_data(:,2),2,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerEdgeAlpha',0.2,'MarkerFaceAlpha',.2); hold on;
scatter(((1+eps).*exp(net_Re(test_Proc(:,3:end)')))./(1+exp(net_Re(test_Proc(:,3:end)'))),((1+eps).*exp(test_Proc(:,2)))./(1+exp(test_Proc(:,2))),2,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerEdgeAlpha',0.35,'MarkerFaceAlpha',.35); hold on;
plot(Proc_data(:,2),Proc_data(:,2),'LineWidth',1,'Color','k');
lin = fitlm(((1+eps).*exp(net_Re(test_Proc(:,3:end)')))./(1+exp(net_Re(test_Proc(:,3:end)'))),((1+eps).*exp(test_Proc(:,2)))./(1+exp(test_Proc(:,2)))); R2_Re = lin.Rsquared.Adjusted;
axis([0 1 0 1]);
text(0.55,0.1,sprintf('R^2_{adj} = %1.3f',R2_Re),'Units','Normalized');
text(0.03,0.94,'Recovery','Units','Normalized','FontWeight','bold')
box on; 

nexttile;
% predict the prod using the training and testing data
scatter(Econ_data(:,1),10.^(net_Prod(Econ_data(:,3:end)')),2,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerEdgeAlpha',0.2,'MarkerFaceAlpha',.2,'HandleVisibility','off'); hold on; 
scatter(10.^(test_Econ(:,1)),10.^(net_Prod(test_Econ(:,3:end)')),2,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerEdgeAlpha',0.35,'MarkerFaceAlpha',.35,'HandleVisibility','off'); hold on; plot(Econ_data(:,1),Econ_data(:,1),'LineWidth',1,'Color','k','HandleVisibility','off')
lin = fitlm(10.^(net_Prod(test_Econ(:,3:end)')),10.^(test_Econ(:,1))); R2_Prod = lin.Rsquared.Adjusted;
scatter(Econ_data(1,1),10.^(net_Prod(Econ_data(1,3:end)')),2,'MarkerFaceColor','r','MarkerEdgeColor','r')
scatter(Econ_data(1,1),10.^(net_Prod(Econ_data(1,3:end)')),2,'MarkerFaceColor','b','MarkerEdgeColor','b')
set(gca,'Xscale','log'); set(gca,'Yscale','log');
axis([min(Econ_data(:,1)),max(Econ_data(:,1)),min(Econ_data(:,1)),max(Econ_data(:,1))]);
text(0.55,0.1,sprintf('R^2_{adj} = %1.3f',R2_Prod),'Units','Normalized');
text(0.03,0.83,{'Productivity','[{\itmol_{CO_2} kg^{-1} s^{-1}}]'},'Units','Normalized','FontWeight','bold')
legend('train data','test data','FontWeight','bold');
legend('boxoff')
legend('Position',[0.26 0.21 0.25 0.075])
box on; 

nexttile;
% predict the energy using the training and testing data
scatter(Econ_data(:,2),10.^(net_En(Econ_data(:,3:end)')),2,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerEdgeAlpha',0.2,'MarkerFaceAlpha',.2); hold on; plot(Econ_data(:,2),Econ_data(:,2),'LineWidth',1,'Color','k','HandleVisibility','off')
scatter(10.^(test_Econ(:,2)),10.^(net_En(test_Econ(:,3:end)')),2,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerEdgeAlpha',0.35,'MarkerFaceAlpha',.35); hold on; plot(Econ_data(:,2),Econ_data(:,2),'LineWidth',1,'Color','k')
lin = fitlm(10.^(net_En(test_Econ(:,3:end)')),10.^(test_Econ(:,2))); R2_En = lin.Rsquared.Adjusted;
set(gca,'Xscale','log'); set(gca,'Yscale','log');
axis([min(Econ_data(:,2)),max(Econ_data(:,2)),min(Econ_data(:,2)),max(Econ_data(:,2))]);
text(0.55,0.1,sprintf('R^2_{adj} = %1.3f',R2_En),'Units','Normalized');
text(0.03,0.83,{'Energy','[{\itkWh ton^{-1}_{CO_2}}]'},'Units','Normalized','FontWeight','bold')
box on; 

xlabel(t, 'ANN Predictions')
ylabel(t, 'PSA Model Data')
