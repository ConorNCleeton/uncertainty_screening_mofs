function [] = visualise_results(obj_struct,vars_struct,name,iso_params)
% Function to visualise the results of the optimisation. Each figure will
% produce a purity - recovery pareto front plot with subplots for each
% different optimisation (using different random seeds for the
% initialisation to ensure no local minima are observed).
FF = {'DRE+DDEC','DRE+EQeq','DRE+Neutral','DRE+Qeq','UFF+DDEC','UFF+EQeq','UFF+Neutral','UFF+Qeq'};
obj = obj_struct;
f = figure('Color','w','visible','off');
set(groot,'defaulttextinterpreter','latex');

t = tiledlayout('flow','TileSpacing','tight');
nexttile;
c= [1 0 0;0.0235 0.8000 0.2431;0 0 0;0 0 1];
set(gcf,'Position',[2000 100 1000 750])

for k = 1:5
    obj_data = obj(:,k);
    makeplot(obj_data, c, FF)
    nexttile;
end

% Get max purity output 
hold on;
hatch=fill([0.9,1,1,0.9],[0,0,1.2,1.2],[0.5 0.5 0.5],'HandleVisibility','off');
hatch.FaceAlpha=0.3;
for k = 1:length(FF)  
    obj_data = obj(k,:);
    cat_Data = cat(1,obj{k,:});
    [val,idx] = max(cat_Data(:,1));
    [rows, cols] = size(obj{k,1});
    iter = floor(idx/rows)+1; % choose iteration with the highest purity    
    if k<length(FF)/2 + 1
        try
            scatter(obj{k,iter}(:,1),obj{k,iter}(:,2),15,'s','MarkerEdgeColor',c(k,:),'MarkerFaceColor',c(k,:),'DisplayName',FF{k})
        catch
            scatter([0 0],[0 0],15,'s','MarkerEdgeColor',c(k,:),'MarkerFaceColor',c(k,:),'DisplayName',FF{k})
        end
        
    else
        try
            scatter(obj{k,iter}(:,1),obj{k,iter}(:,2),15,'MarkerEdgeColor',c(k-(length(FF)/2),:),'MarkerFaceColor',c(k-(length(FF)/2),:),'DisplayName',FF{k})
        catch
            scatter([0 0],[0 0],15,'MarkerEdgeColor',c(k-(length(FF)/2),:),'MarkerFaceColor',c(k-(length(FF)/2),:),'DisplayName',FF{k})
        end
    end
end
ylabel('Recovery');xlabel('Purity');
title('Final Pareto Fronts')
axis([0.15 1 0.9 1.01]);
box on; 


% Plot CO2 isotherms
nexttile;
pts = 10000;
P = linspace(0,10,pts);
hold on;
for i = 1:length(FF)
    q = DSL(iso_params(i,1:6),298,P);
    if i < length(FF)/2 + 1
        plot(P,q,'--','Color',c(i,:),'LineWidth',1.5,'DisplayName',FF{i});
    else
        plot(P,q,'Color',c(i-(length(FF)/2),:),'LineWidth',1.5,'DisplayName',FF{i});
    end
end
xlabel('P [bar]'); ylabel('$q_{CO_2}^*\ [mmol\ g^{-1}]$');
% set(gca,'Xscale','log')
box on; 
xlim([0 1]);
lgd2 = legend; 
lgd2.Layout.Tile = 9;

nexttile;
P = linspace(0,1,pts);
hold on;
for i = 1:length(FF)
    q = DSL(iso_params(i,7:end),298,P);
    if i < length(FF)/2 + 1
        plot(P,q,'--','Color',c(i,:),'LineWidth',1.5);
    else
        plot(P,q,'Color',c(i-(length(FF)/2),:),'LineWidth',1.5);
    end
end
xlabel('P [bar]'); ylabel('$q_{N_2}^*\ [mmol\ g^{-1}]$');
box on; 

%% SAVE FIGURE
% exportgraphics(f,sprintf('%s/visualisations/%s.png',pwd,name),'Resolution',300)
if ~exist('visualisations', 'dir')
    mkdir('visualisations')
end
print(sprintf('%s/visualisations/%s',pwd,name),'-dpng','-r300')

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

    function []= makeplot(obj,c,FF)
    hold on;
    h=fill([0.9,1,1,0.9],[0,0,1.2,1.2],[0.5 0.5 0.5],'HandleVisibility','off');
    h.FaceAlpha=0.3;
    for j = 1:length(FF)
        if j<length(FF)/2 + 1
            try
                scatter(obj{j,1}(:,1),obj{j,1}(:,2),15,'MarkerEdgeColor',c(j,:),'MarkerFaceColor',c(j,:),'DisplayName',FF{j})
            catch
                scatter([0 0],[0 0],15,'MarkerEdgeColor',c(j,:),'MarkerFaceColor',c(j,:),'DisplayName',FF{j})
            end

        else
            try
                scatter(obj{j,1}(:,1),obj{j,1}(:,2),15,'s','MarkerEdgeColor',c(j-(length(FF)/2),:),'MarkerFaceColor',c(j-(length(FF)/2),:),'DisplayName',FF{j})
            catch
                scatter([0 0],[0 0],15,'s','MarkerEdgeColor',c(j-(length(FF)/2),:),'MarkerFaceColor',c(j-(length(FF)/2),:),'DisplayName',FF{j})
            end
        end
    end
    ylabel('Recovery');xlabel('Purity');
    axis([0.15 1 0.9 1.01]);
    box on; 
    end
end
