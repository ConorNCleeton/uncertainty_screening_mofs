% Plot the data distributions used to train the final surrogate model
clear; clc; close all;

% load in the required data and partition to the appropriate variable types
load('../1_generate_ANN_training_data/3_training_data/total_raw_data.mat');
train_data = curated_data; 
design_vars = train_data(:,5:8);          % PSA cycle variables
material_vars = train_data(:,9:13);       % pellet and crystal variables 
Isotherm_Params = train_data(:,14:end);   % CO2 and N2 DSL model params

%% plot target distributions
figure('Color','w');
t = tiledlayout('flow','TileSpacing','compact');
nexttile; histogram(train_data(:,1),'NumBins',30,'Normalization','probability','HandleVisibility','off');
xlabel('Purity [-]'); set(gca,'YTick',[]); set(gca,'YTicklabel',[]); 
nexttile; histogram(train_data(:,2),'NumBins',30,'Normalization','probability','HandleVisibility','off');
xlabel('Recovery [-]'); set(gca,'YTick',[]); set(gca,'YTicklabel',[]); 
nexttile; histogram(train_data(logical(train_data(:,3)<0.01),3),'NumBins',30,'Normalization','probability','HandleVisibility','off');
xlabel('Productivity [\it{mol_{CO_2} kg^{-1} s^{-1}}]');set(gca,'YTick',[]); set(gca,'YTicklabel',[]); 
nexttile; histogram(train_data(logical(train_data(:,4)<1e4),4),'NumBins',30,'Normalization','probability','HandleVisibility','off');
xlabel('Energy [\it{kWh ton_{CO_2}^{-1}}]');set(gca,'YTick',[]); set(gca,'YTicklabel',[]); 



%% Plot the design variables
idx_titles = ["P_H [Pa]" "\nu_{Feed} [m s^{-1}]" "t _{Ads} [s]" "\alpha_{LR} [-]"];
idx_length = (size(design_vars,2)).^2;
idx_position = 1:1:idx_length;
idx_position = reshape(idx_position,[size(design_vars,2),size(design_vars,2)])';
figure('Color','w');
for i = 1:size(design_vars,2)
    for j = 1:size(design_vars,2)
        subplot(size(design_vars,2),size(design_vars,2),idx_position(i,j));
        if i == j
            histogram(design_vars(:,i),15);
        else
            scatter(design_vars(:,j),design_vars(:,i),0.1,'o');
        end
        if j == 1
            ylabel(idx_titles(i),'FontWeight','bold');
        end
        if i == 4
            xlabel(idx_titles(j),'FontWeight','bold');
        end
        box on; 
    end
    box on; 
end

%% Plot the material properties
idx_titles = ["\epsilon_{B} [-]" "\epsilon_{p} [-]" "r_p [m]" "\rho_{s} [kg m^{-3}]" "c_{p,s} [J kg^{-1} K^{-1}]"];
idx_length = (size(material_vars,2)).^2;
idx_position = 1:1:idx_length;
idx_position = reshape(idx_position,[size(material_vars,2),size(material_vars,2)])';
figure('Color','w');
for i = 1:size(material_vars,2)
    for j = 1:size(material_vars,2)
        subplot(size(material_vars,2),size(material_vars,2),idx_position(i,j));
        if i == j
            histogram(material_vars(:,i),15);
        else
            scatter(material_vars(:,j),material_vars(:,i),0.1,'o');
        end
        if j == 1
            ylabel(idx_titles(i),'FontWeight','bold');
        end
        if i == 5
            xlabel(idx_titles(j),'FontWeight','bold');
        end
        box on;
    end
    box on;
end

%% Plot the isotherms
pts = 50;
n_samples = 100;%length(Isotherm_Params);
idx =randi(length(Isotherm_Params),100000,1);
P = linspace(0,1,pts);
figure('Color','w'); hold on;
subplot(1,2,1)
hold on;
for i = idx'%length(Isotherm_Params)
    q = DSL(Isotherm_Params(i,1:6),298,P);
    p1 = plot(P,q,'-b');
    p1.Color(4) = 0.1;
end
xlabel('P [bar]'); ylabel('q^*_{CO_2} [mmol g^{-1}]');
hold off;
box on;

subplot(1,2,2)
hold on;
for i = idx'%length(Isotherm_Params)
    q = DSL(Isotherm_Params(i,7:end),298,P);
    p1 = plot(P,q,'-r');
    p1.Color(4) = 0.1;
end
xlabel('P [bar]'); ylabel('q^*_{N_2} [mmol g^{-1}]');
box on; 


%% Utility function
function q = DSL(params,T,P)

q_s1 = params(1); 
q_s2 = params(2);
b01 = params(3);
b02 = params(4);
deltaH1 = params(5);
deltaH2 = params(6);

b1 = b01*exp(deltaH1/(8.314*T));
b2 = b02*exp(deltaH2/(8.314*T));

site1 = (q_s1*b1.*P)./(1+b1.*P);
site2 = (q_s2*b2.*P)./(1+b2.*P);

q = site1+site2;

end