function [results, parameters]=pareto_progression(A, ploton,type)
%% Function file for extracting the objective function with every generation
% A = load population file from optimisation 
% -ploton: Optional flag to indicate whether the data should be plotted (0 = no plot, 1 = plot).

a=A.pops(:,:); % Extracting the population structure array

[generations, pop_size] = size(a);

objectives=length(a(1, 1).obj);
variables=length(a(1, 1).var);

results = zeros(generations, pop_size, objectives);

for i= 1:generations
    for j = 1:pop_size
        for k = 1:objectives
            results(i,j,k) = a(i,j).obj(k);
        end
    end
end
for i= 1:generations
    for j = 1:pop_size
        for k = 1:variables
            parameters(i,j,k) = a(i,j).var(k);
        end
    end
end

results=abs(results);
if ploton
    figure('Color','w');
    for i = 1:generations
        hold on;
        if i == generations
            plot(results(i,:,1),results(i,:,2),'or')
        else
            plot(results(i,:,1),results(i,:,2),'ob')
        end
    end
    switch type
        case 'Process_Optimisation'
            xlabel('Purity');
            ylabel('Recovery')
        case 'Economic_Optimisation'
            xlabel('Productivity / mol_{CO_2} kg_{adsorbent}^{-1} s^{-1}');
            ylabel('Energy / kWh ton_{CO_2}^{-1}')
    end
end

% Extracting the population distribution
% b=A.pops(end, :);
% [~, n]=size(b);
% q=length(b(1, 1).obj);
% qq=b(1, 1).var;
% [~, nn]=size(qq);
% parameters=zeros(n, nn);
% for i= 1:n
%    parameters(i, 1:nn)=b(1, i).var;
% end

end