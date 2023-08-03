% Visualise all results contained in the /opt_results/ subfolder
clear all; close all; clc;

addpath('./opt_results');

% get optimised outputs
myFolder = sprintf('%s/opt_results/',pwd);
opt_results_path = fullfile(myFolder, '*_objectivesProc.mat');
var_results_path = fullfile(myFolder, '*_variablesProc.mat');
myFolder = sprintf('%s/sorted_material_dataframes/',pwd);
filePattern = fullfile(myFolder, '*.csv');
theFiles = dir(filePattern);
optFiles = dir(opt_results_path);
varFiles = dir(var_results_path);

for k = 1:length(theFiles)
    % get iso parameters
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    FF_dfs{k} = readtable(fullFileName);
end

for k = 1:length(optFiles)
    try
        % get outputs
        baseFileName = optFiles(k).name;
        optFileName = fullfile(optFiles(k).folder, baseFileName);
        baseFileName = varFiles(k).name;
        varFileName = fullfile(varFiles(k).folder, baseFileName);
        opt_data = load(optFileName);
        opt_data = opt_data.data; 
        var_data = load(varFileName);
        var_data = var_data.data; 
        
        % get name 
        n = split(baseFileName,'_');
        name = n{1};
        
        % get iso params
        iso_params = [];  
        for FF = 1:length(FF_dfs)
            data = FF_dfs{FF};
            material_data = data(strcmp(data.material,name),:);
            iso_params = [iso_params;table2array(material_data(1,4:end-2))];
        end
        visualise_results(opt_data,var_data,name,iso_params);
    catch
        continue
    end
end