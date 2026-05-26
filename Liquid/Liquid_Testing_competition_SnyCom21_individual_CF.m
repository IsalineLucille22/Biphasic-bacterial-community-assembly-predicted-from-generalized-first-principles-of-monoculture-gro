%Compare observed individual growth in soil of Syncom21 members with observed growth in Syncom20 mixture, and with modeled SynCom20 growth without cross-feeding.

%clearvars -except z 
clear
close all

%data path for individual growth
pwd_init = pwd;
pwd_init = fileparts(pwd_init);
addpath(strcat(pwd_init, '/Data'))
addpath(pwd_init)
name_Exp = 'Zenodo_Senka_21_CF_Liquid';

%% plot individual biomass in soil %Doesn't make sense here because not the same environment


%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
Data_Evol = readtable(strcat('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data/Liquid_Data/','abund_cfu_se.xlsx'), 'Sheet', 3, 'Range','26:47', 'Format','auto');%Data in liquid from Clara's experiments

biomass = readtable('MergedData_copy.xlsx','sheet','Growth kinetic params','Range','U3:X24','VariableNamingRule','preserve');
strains = readtable('MergedData_copy.xlsx','sheet','Growth kinetic params','Range','C3:C24','VariableNamingRule','preserve');

combined = [strains, biomass]; %Individual biomass in soil
combined = sortrows(combined,'Name','ascend');
combined = [combined, table(mean(table2array(combined(:,2:5)),2))];
combined = renamevars(combined,{'Var1'},{'mean'});
S = size(strains); S = S(1);


%actual absolute abundance data, has quadriplicate values for each time point

name = string(table2array(Parameters_set(1:S,1)));

% get single cell weight list

single_cell_mass = readtable(strcat('Data/','MergedData_copy.xlsx'), 'Sheet', 'Growth kinetic params','Range','C3:D24','Format','auto');
single_cell_mass = sortrows(single_cell_mass,'Name','ascend');

%sort alphabetically to strains
[name, idx2]= sortrows(name);

Data_Evol = Data_Evol(idx2,:);

%Automatized that
max_obs_biomass = mean(table2array(Data_Evol(:,18:S)),2); %timepoint 7

max_obs_biomass = [table(name), array2table(max_obs_biomass), Data_Evol(:,18:S)];%[table(name),array2table(max_obs_biomass), Data_Evol(:,22:25)];


%make a detour to have the strain names from top to bottom in alphabetical order

inv_names = sortrows(combined,'Name','descend');
xvalues = categorical(inv_names.Name);
X = reordercats(xvalues, inv_names.Name);
yvalues = max_obs_biomass.max_obs_biomass./(single_cell_mass.assumedCellDwInFgC*1e-15);

%make sure to flip the yvalues to have the inverted order

figH = figure;

subplot('Position',[0.3 0.1 0.2 0.62]);

barh(X,flip(log10(yvalues))) %observed biomass in Syncom21
hold on
barh(X,flip(log10(combined.mean)), 'FaceAlpha', 0.5)
%Scatter for each replicate
for i = 2:5
    scatter(flip(log10(table2array(combined(:,i)))),1:S,50,'.','r')
end
%Automatized this
for i = 3:6
    scatter(flip(log10(table2array(max_obs_biomass(:,i)./(single_cell_mass.assumedCellDwInFgC*1e-15)))),1:21,50,'.','b')
end

xlim([3 9]);
grid on
iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'comparison-individual-SynCom21-biomass-in-soil.pdf');
saveas(figH,FigName,'pdf');

% ttest differences

ind_syncom_test = mattest(table2array(combined(:,2:5)),table2array(max_obs_biomass(:,3:6)./(single_cell_mass.assumedCellDwInFgC*1e-15)));
FDR = mafdr(ind_syncom_test, 'BHFDR', true); % Estimate positive false discovery rate for multiple hypothesis testing
ind_syncom_test=[combined.Name table(ind_syncom_test) table(FDR)];
writetable(ind_syncom_test, 'Data/individual-syncom-mattest.csv');

% take the difference of the means as fraction of the SynCom20 biomass

diff_vector = max_obs_biomass.max_obs_biomass-(combined.mean.*single_cell_mass.assumedCellDwInFgC*1e-15);
diff_fraction = diff_vector./max_obs_biomass.max_obs_biomass;

close all

figH = figure;
subplot('Position',[0.1 0.3 0.4 0.3]);
stem(diff_vector);
xtickangle(90)
set(gca,'xtick',1:S,'xticklabel',combined.Name)
ylabel('biomass difference');
title('Syncom - individual');
grid on
ylim([-4e-4 4e-4]);

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'syncom-individual-biomass-diff.pdf');
saveas(figH,FigName,'pdf');
%% comparison of Syncom actual growth to Syncom prediction with/without cross-feeding
% run script 'Plot_syncom_from_multi_simulation_data_no_CF'
% variable 'z' has the predicted species' biomass wiht/without cross-feeding 

load(strcat('Data/', name_Exp, 'z.mat'))
CF_biomass = z; %Reorder z according to the names
CF_biomass(CF_biomass < 0) = 0;
nb_data_set = size(CF_biomass, 2) - 1;
%the variable 'name' corresponds to the strain order on no_CF_biomass

CF_biomass = [table(name) array2table(CF_biomass)];
CF_biomass = sortrows(CF_biomass,'name','ascend');

no_CF_pop_counts=[table(name) array2table(max(z./(single_cell_mass.assumedCellDwInFgC*1e-15), 1))]; %Smaller than the dry biomass per gram per mL correct it
no_CF_pop_counts=sortrows(no_CF_pop_counts,'name','ascend');

% plot on top of each other

X = reordercats(xvalues,inv_names.Name);

close all

figH = figure;

subplot('Position',[0.3 0.1 0.2 0.62]);

barh(X,flip(log10(yvalues))) %observed biomass in Syncom20
hold on
barh(X,flip(log10(mean(table2array(no_CF_pop_counts(:, 2:(nb_data_set + 1))),2))), 'FaceAlpha', 0.5)
for i = 3:6
    scatter(flip(log10(table2array(max(max_obs_biomass(:,i)./(single_cell_mass.assumedCellDwInFgC*1e-15), 1)))),1:S,50,'.','b')
end
for i = 2:(nb_data_set + 1)
    scatter(flip(log10(table2array(no_CF_pop_counts(:,i)))),1:S,50,'.','r')
end

xlim([3 9]);
grid on
iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'comparison-SynCom20-model-CF-biomass-in-soil.pdf');
saveas(figH,FigName,'pdf');

% ttest differences

syncom_CF_test = mattest(table2array(CF_biomass(:, 2:(nb_data_set + 1))),table2array(max_obs_biomass(:,3:6)));

FDR = mafdr(syncom_CF_test, 'BHFDR', true); % Estimate positive false discovery rate for multiple hypothesis testing

syncom_CF_test=[combined.Name table(syncom_CF_test) table(FDR)];

writetable(syncom_CF_test, 'Data/syncom-model-CF-mattest.csv');

% take the difference of the means as fraction of the SynCom20 biomass
diff_vector = max_obs_biomass.max_obs_biomass - mean(table2array(CF_biomass(:,2:(nb_data_set + 1))),2);

diff_fraction = diff_vector./max_obs_biomass.max_obs_biomass;

close all

figH = figure;

subplot('Position',[0.1 0.3 0.4 0.3]);

stem(diff_vector);
xtickangle(90)
set(gca,'xtick',1:S,'xticklabel',combined.Name)
ylabel('biomass difference');
title('Syncom - model-CF');
grid on
ylim([-4e-4 4e-4]);
iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'syncom-model_CF-biomass-diff.pdf');
saveas(figH,FigName,'pdf');