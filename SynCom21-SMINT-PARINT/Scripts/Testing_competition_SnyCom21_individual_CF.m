%Compare observed individual growth in soil of Syncom21 members with observed growth in Syncom20 mixture, and with modeled SynCom20 growth without cross-feeding.
%Be carful to adapt the variable Data_Range (l.52) depending on the experiment
%name and time of the analysis (1wk, 21 days, etc)


%clearvars -except z 
clear
close all

% data path for individual growth
scriptFolder = fileparts(mfilename('fullpath'));
projectFolder = fileparts(scriptFolder);
dataFolder = strcat(fullfile(projectFolder, 'Data'), '/');
FiguresFolder = strcat(fullfile(projectFolder, 'Figures'), '/');
name_Exp = 'Zenodo_Senka_Average_Interaction';

%% plot individual biomass in soil
%Addition stat

biomass = readtable(strcat(dataFolder, 'MergedData.xlsx'),'sheet','Growth kinetic params with rep','Range','U3:X24','VariableNamingRule','preserve');
strains = readtable(strcat(dataFolder, 'MergedData.xlsx'),'sheet','Growth kinetic params with rep','Range','C3:C24','VariableNamingRule','preserve');

combined = [strains, biomass]; %Individual biomass in soil
combined = sortrows(combined,'Name','ascend');
combined = [combined, table(mean(table2array(combined(:,2:5)),2))];
combined = renamevars(combined,{'Var1'},{'mean'});
S = size(strains); S = S(1);

%actual absolute abundance data, has quadriplicate values for each time point
%Loadind data
Parameters_set = readtable(strcat(dataFolder,'MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');

%actual absolute abundance data, has quadriplicate values for each time point

%Data_Evol = readtable(strcat('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data/Liquid_Data/','abund_cfu_IG.xlsx'), 'Sheet', 2, 'Range','46:67', 'Format','auto');%1:22 without correction for 0. Data in soil extract from Phil and Clara's experiments
% Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 7, 'Range','1:22', 'Format','auto'); %Senka's data
Data_Evol = readtable(strcat(dataFolder,'SSC21_genera_relative-abundances.xlsx'), 'Sheet', 4, 'Range','1:22', 'Format','auto'); %Bruna's data
Time_step = [0 12 22 38 70 168 504];%Bruna %[0 1 3 7 10 21]*24;%Senka %
name = string(table2array(Parameters_set(1:21,1)));
[~, nb_obs] = size(Data_Evol);
nb_obs = nb_obs - 1;
nb_time_step = length(Time_step);
nb_rep = nb_obs/nb_time_step;

% get single cell weight list

single_cell_mass = readtable(strcat(dataFolder,'MergedData.xlsx'), 'Sheet', 'Growth kinetic params with rep','Range','C3:D24','Format','auto');
single_cell_mass = sortrows(single_cell_mass,'Name','ascend');

%sort alphabetically to strains
[name, idx2]= sortrows(name);
Data_Evol = Data_Evol(idx2,:);

%CHANGE SELECTED ROW ACCORDING TO TIME, 1wk vs 21 days or other
%possibilities
Data_Range = Data_Evol(:,32:36); %Bruna 21 days;%Data_Evol(:, 22:25); %Senka 21-days% mean(table2array(Data_Evol(:,14:17)),2);%timepoint 7 %mean(table2array(Data_Evol(:,27:31)),2); %timepoint 21 %
mean_obs_biomass = mean(table2array(Data_Range),2); 
mean_obs_biomass = [table(name), array2table(mean_obs_biomass), Data_Range];%[table(name),array2table(max_obs_biomass), Data_Evol(:,14:17)];%


%make a detour to have the strain names from top to bottom in alphabetical order

inv_names = sortrows(combined,'Name','descend');
xvalues = categorical(inv_names.Name);
X = reordercats(xvalues, inv_names.Name);
yvalues = mean_obs_biomass.mean_obs_biomass./(single_cell_mass.assumedCellDwInFgC*1e-15);

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
% %Automatized this
for i = 3:6
    scatter(flip(log10(table2array(mean_obs_biomass(:,i)./(single_cell_mass.assumedCellDwInFgC*1e-15)))),1:S,50,'.','b')
end
xlim([4 9]);
grid on
iFolderName = FiguresFolder;
FigName = strcat(iFolderName, name_Exp, 'comparison-individual-SynCom21-biomass-in-soil.pdf');
% saveas(figH,FigName,'pdf');
print(figH,FigName, '-dpdf', '-painters');

% ttest differences
ind_syncom_test = mattest(table2array(combined(:,2:5)),table2array(mean_obs_biomass(:,3:6)./(single_cell_mass.assumedCellDwInFgC*1e-15)));
FDR = mafdr(ind_syncom_test, 'BHFDR', true); % Estimate positive false discovery rate for multiple hypothesis testing
ind_syncom_test=[combined.Name table(ind_syncom_test) table(FDR)];
filename = fullfile(dataFolder, 'individual-syncom-mattest.csv');
writetable(ind_syncom_test, filename);

% take the difference of the means as fraction of the SynCom20 biomass

mean_diff_vector = mean_obs_biomass.mean_obs_biomass-(combined.mean.*single_cell_mass.assumedCellDwInFgC*1e-15);
diff_fraction = mean_diff_vector./mean_obs_biomass.mean_obs_biomass;

close all

figH = figure;
subplot('Position',[0.1 0.3 0.4 0.3]);
stem(mean_diff_vector);
xtickangle(90)
set(gca,'xtick',1:21,'xticklabel',combined.Name)
ylabel('biomass difference');
title('Syncom - individual');
grid on
ylim([-3e-4 3e-4]);

iFolderName = FiguresFolder;
FigName = strcat(iFolderName, name_Exp, 'syncom-individual-biomass-diff.pdf');
print(figH,FigName, '-dpdf', '-painters');
%% comparison of Syncom actual growth to Syncom prediction with/without cross-feeding
% run script 'Plot_syncom_from_multi_simulation_data_no_CF'
% variable 'z' has the predicted species' biomass wiht/without cross-feeding 

load(strcat(dataFolder, name_Exp, 'z.mat'))
CF_biomass = z; %Reorder z according to the names
CF_biomass(CF_biomass < 0) = 0;
nb_rep_sim = length(CF_biomass(1, :)) - 1;
%the variable 'name' corresponds to the strain order on no_CF_biomass

CF_biomass = [table(name) array2table(CF_biomass)];
CF_biomass = sortrows(CF_biomass,'name','ascend');

no_CF_pop_counts=[table(name) array2table(max(z./(single_cell_mass.assumedCellDwInFgC*1e-15), 1))]; %Smaller than the dry biomass per gram per mL correct it
no_CF_pop_counts=sortrows(no_CF_pop_counts,'name','ascend');
nb_data_set = size(CF_biomass, 2) - 1;

% plot on top of each other
X = reordercats(xvalues,inv_names.Name);

close all
figH = figure;
subplot('Position',[0.3 0.1 0.2 0.62]);

barh(X,flip(log10(yvalues))) %observed biomass in Syncom21
hold on
barh(X,flip(log10(mean(table2array(no_CF_pop_counts(:,2:(nb_data_set + 1))),2))), 'FaceAlpha', 0.5)
for i = 3:6
    scatter(flip(log10(table2array(max(mean_obs_biomass(:,i)./(single_cell_mass.assumedCellDwInFgC*1e-15), 1)))),1:S,50,'.','b')
end
for i = 2:(nb_data_set + 1)
    scatter(flip(log10(table2array(no_CF_pop_counts(:,i)))),1:S,50,'.','r')
end

xlim([5 9]);
grid on
iFolderName = FiguresFolder;
FigName = strcat(iFolderName, name_Exp, 'comparison-SynCom20-model-CF-biomass-in-soil.pdf');
print(figH,FigName, '-dpdf', '-painters');

%Permanova and permutation-test on Bray-Curtis distance to get on statistic 

% A, B: n x m (m=4)
m_1 = nb_rep;
m_2 = nb_rep_sim;
mat_Comm_1 = table2array(mean_obs_biomass(:, 3:6));
mat_Comm_2 = table2array(CF_biomass(:,2:(nb_data_set + 1)));
X = [mat_Comm_1.'; mat_Comm_2.'];          % (m_1 + m_2)xn   rows = replicates, cols = species
group = [ones(m_1, 1); 2*ones(m_2, 1)];  % labels: 1 = SynCom20, 2 = SnyCom21
X = log1p(X);  % log(1+x)
D = BrayCurtisDistance(X);
nperm = 9999; %Choice to get the smallest p-value = 0.0001
[F_1, R2_1, p_1] = permanova1(D, group, nperm); %R2 = SS_bet/SS_tot, fraction mult variation explained by the difference between the 2 communities
%Permutation test entire community
p_val_tot_1 = Permutation_Test_Comm(D, group, nperm);

% ttest differences

syncom_CF_test = mattest(table2array(CF_biomass(:,2:(nb_data_set + 1))),table2array(mean_obs_biomass(:, 3:6)));
FDR = mafdr(syncom_CF_test, 'BHFDR', true); % Estimate positive false discovery rate for multiple hypothesis testing
syncom_CF_test=[combined.Name table(syncom_CF_test) table(FDR)];
filename = fullfile(dataFolder, 'syncom-model-CF-mattest.csv');
writetable(ind_syncom_test, filename);


% take the difference of the means as fraction of the SynCom21 biomass
rand_res_sim = randperm(min(nb_rep_sim, nb_rep), min(nb_rep_sim, nb_rep)) + 1;
diff_vector = table2array(mean_obs_biomass(:, 3:end) - table2array(CF_biomass(:, rand_res_sim)));
mean_diff_vector = mean_obs_biomass.mean_obs_biomass - mean(table2array(CF_biomass(:,2:(nb_data_set + 1))),2);
diff_fraction = mean_diff_vector./mean_obs_biomass.mean_obs_biomass;
close all
figH = figure;
subplot('Position',[0.1 0.3 0.4 0.3]);

stem(mean_diff_vector);
xtickangle(90)
set(gca,'xtick',1:S,'xticklabel',combined.Name)
ylabel('biomass difference');
title('Syncom - model-CF');
grid on
ylim([-3e-4 2e-4]);
iFolderName = FiguresFolder;
FigName = strcat(iFolderName, name_Exp, 'syncom-model_CF-biomass-diff.pdf');
% saveas(figH,FigName,'pdf');
print(figH,FigName, '-dpdf', '-painters');

%%%New stacked 

figK = figure;
bar(mean_diff_vector)
hold on
scatter(1:S, diff_vector, 60, '.', 'r')
grid on
axis([0 22 -3e-04 8e-05])
xticks(1:length(name)); 
xticklabels(name);
ylabel('measured - simulated')

iFolderName = FiguresFolder;
FigName = strcat(iFolderName, name_Exp, 'stacked_plot_diff.pdf');
print(figK,FigName, '-dpdf', '-painters');
%% Comparison individual Observed vs simulated. We simulate the biomass after one week and not 21 days. Senka (trained) data
%Simulations soil mono-cultures

num_fig = 1;
Name_file = 'struct_tot_SynCom21_Soil_only_Lyso_Parfor';%'data_to_save_SynCom21_Soil_only_Lyso_Parfor_Predation';%

biomass = readtable('MergedData.xlsx','sheet','Growth kinetic params with rep','Range','U3:X24','VariableNamingRule','preserve');
strains = readtable('MergedData.xlsx','sheet','Growth kinetic params with rep','Range','C3:C24','VariableNamingRule','preserve');

combined = [strains, biomass]; %Individual biomass in soil
combined = sortrows(combined,'Name','ascend');
combined = [combined, table(mean(table2array(combined(:,2:5)),2))];
combined = renamevars(combined,{'Var1'},{'mean'});
S = size(strains); S = S(1);
[~,name_sorted] = sort(name);
name_sorted = flip(name_sorted);

%actual absolute abundance data, has quadriplicate values for each time point
%Loadind data
Parameters_set = readtable(strcat(dataFolder,'MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
Dry_set = readtable(strcat(dataFolder,'MergedData.xlsx'), 'Sheet', 5,'Range','29:50', 'Format','auto');

%actual absolute abundance data, has quadriplicate values for each time point
Data_Evol = readtable(strcat(dataFolder,'S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 7, 'Range','1:22', 'Format','auto'); %Senka's data
% Data_Evol = readtable(strcat('Data/','SSC21_genera_relative-abundances.xlsx'), 'Sheet', 4, 'Range','1:22', 'Format','auto'); %Bruna's data
name = string(table2array(Parameters_set(1:21,1)));
Time_step =[0 1 3 7 10 21]*24;%Senka %
[~, nb_obs] = size(Data_Evol);
nb_obs = nb_obs - 1;
nb_time_step = length(Time_step);
nb_rep = nb_obs/nb_time_step;
Data_Evol_temp = table2array(Data_Evol(:, 2:end));
mean_y_0 = mean(Data_Evol_temp(:,1:nb_rep), 2); 
Dry_set = table2array(Dry_set(:, 8));


% get single cell weight list

single_cell_mass = readtable(strcat(dataFolder,'MergedData.xlsx'), 'Sheet', 'Growth kinetic params with rep','Range','C3:D24','Format','auto');
single_cell_mass = sortrows(single_cell_mass,'Name','ascend');

%sort alphabetically to strains
[name, idx2] = sortrows(name);
Data_Evol = Data_Evol(idx2,:);
Data_Range = Data_Evol(:,14:17);%timepoint 7 %Data_Evol(:, 22:25); %Data_Evol(:,32:36); %Bruna 21 days;%Senka 21-days% mean(table2array(Data_Evol(:,27:31)),2); %timepoint 21 %
mean_obs_biomass = mean(table2array(Data_Range),2); 
mean_obs_biomass = [table(name), array2table(mean_obs_biomass), Data_Evol(:,14:17)];%[table(name), array2table(mean_obs_biomass), Data_Range];%


%make a detour to have the strain names from top to bottom in alphabetical order

inv_names = sortrows(combined,'Name','descend');
xvalues = categorical(inv_names.Name);
X = reordercats(xvalues, inv_names.Name);
yvalues = mean_obs_biomass.mean_obs_biomass./(single_cell_mass.assumedCellDwInFgC*1e-15);

%make sure to flip the yvalues to have the inverted order

figL = figure;
subplot('Position',[0.3 0.1 0.2 0.62]);
%SnyCom21 after 1week
barh(X, flip(log10(yvalues))) %observed biomass in Syncom21
hold on
Comm_1wk = table2array(mean_obs_biomass(:, 3:6)./(single_cell_mass.assumedCellDwInFgC*1e-15));
for i = 3:6
    scatter(flip(log10(table2array(mean_obs_biomass(:,i)./(single_cell_mass.assumedCellDwInFgC*1e-15)))),1:S,50,'.','b')
end
xlim([4 9]);

%Individual growth after one week
figV = figure;
subplot('Position',[0.3 0.1 0.2 0.62]);
barh(X, flip(log10(combined.mean)), 'FaceAlpha', 0.5)
hold on
Mono_Obs_1wk = table2array(combined(:,2:5));
for i = 2:5
    scatter(flip(log10(table2array(combined(:,i)))),1:S,50,'.','r')
end
xlim([4 9]);


%Load parameters
load(strcat(dataFolder, Name_file));
data_to_save = struct_tot_SynCom21_Soil_only_Lyso_Parfor;%'data_to_save_SynCom21_Soil_only_Lyso_Parfor_Predation';%

%Initialization of the colors
colors_ind = [4, 1, 10, 15, 12, 3, 8, 6, 20, 18, 19, 2, 5, 21, 13, 9, 14, 7, 16, 17, 11];
colors_init = {'#B35806', '#E08214', '#D53E4F', '#B2182B', '#B6604B', '#C51B7D', ...
               '#DE77AE', '#F1B6DA', '#FDAE61', '#FEE090', '#A6D96A', '#5AAE61', ...
               '#01665E', '#35978F', '#1B7837', '#C2A5CF', '#9970AB', '#762A83', ...
               '#80CDC1', '#C7EAE5', '#2166AC', '#4393C3', '#B35806'};%distinguishable_colors(S);
colors = {};
for i = 1:S
    colors{i} = colors_init{colors_ind(i)};
end

%Setting for Matlab ODE solver
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
nb_data_set = size(data_to_save, 2);

%Initialization of the model parameters fixed for all replicates 
t_0 = 0; %Time 0
num_fig = num_fig + 1;
end_biomass = zeros(S, nb_rep);
yield_Pred = 0.2;
Time_step = [0 1 3 7]*24;%One-week time step
tspan = [0, max(Time_step)]; %[0, 22*24]; %Time interval in hours
for j = 1:nb_rep
    ind_sim = randi(nb_data_set);
    %Initialization of the interactions 
    
    kappa_mat = data_to_save(ind_sim).kappa_mat;
    CrossFeed_Mat = data_to_save(ind_sim).CrossFeed_Mat_Temp;
    Death_Mat_Temp = data_to_save(ind_sim).Death_Mat_Temp; 
    Threshold_CF = data_to_save(ind_sim).Threshold_CF;
    Threshold_death = data_to_save(ind_sim).Threshold_death; 
    Lag_time_Cons = data_to_save(ind_sim).Lag_time_Cons; 
    Lag_time_Pred = data_to_save(ind_sim).Lag_time_Pred;
    Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
    death_rate =  data_to_save(ind_sim).death_rate;
    name = string(table2array(Parameters_set(1:21,1)));
    Pred_Mat_Lyso = data_to_save(ind_sim).Pred_Mat_Lyso;
    Threshold_Pred = data_to_save(ind_sim).Threshold_Pred;
    R = data_to_save(ind_sim).R; 
    nb_Res = length(R); %Number of resources (group of resources)
    Resource_Matrix = data_to_save(ind_sim).Resource_Matrix; %Addition of a line for Lysobacter
    
    %Modifications upper part
    for i = 1:S
        S_sim = i; %Species present into the subsystem
        nb_temp = length(S_sim); %Number of species
       
        kappa_mat_temp = kappa_mat(S_sim,:);
        CrossFeed_Mat_temp = CrossFeed_Mat(S_sim, S_sim);
        Mat_kappa_3_temp = Mat_kappa_3(S_sim, S_sim);
        Resource_Matrix_temp = Resource_Matrix(S_sim,:);
        Pred_Mat_Lyso_temp = 0*Pred_Mat_Lyso(S_sim,S_sim); %No pred in mono
        Death_Mat_Temp_temp = Death_Mat_Temp(S_sim,S_sim);
        Threshold_CF_temp = Threshold_CF(S_sim);
        Threshold_Pred_temp = Threshold_Pred(S_sim);
        Threshold_death_temp = Threshold_death(S_sim);
        Lag_time_Cons_temp = Lag_time_Cons(S_sim,:);
        Lag_time_Pred_temp = 0*Lag_time_Pred(S_sim,S_sim); %No lag time
        name_temp = name(S_sim);
        mean_y_0_temp = mean_y_0(S_sim);
        death_rate_temp = death_rate(S_sim);
        
        %Setting for Matlab ODE solver
        nb_tot_Species = S*3 + nb_Res;
        opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
        
        %Initialization of the colors
        
        colors = distinguishable_colors(60);
        
        nb_replicates = 1;
        
        %Initial concentrations using a normal distribution
        mat_y_0 = mean_y_0_temp;
        
        mat_y_0 = [mat_y_0 zeros(nb_temp,1) zeros(nb_temp,1)]; 
        
        y_0 = sum(mat_y_0(:,1:2),2);
        mat_y_0 = reshape(mat_y_0', 1, []);
        mat_y_0 = [mat_y_0 R];
        
        sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat_temp, CrossFeed_Mat_temp, Mat_kappa_3_temp, Resource_Matrix_temp,...
        Threshold_CF_temp, Threshold_death_temp, Threshold_Pred_temp, Death_Mat_Temp_temp, death_rate_temp,...
        Pred_Mat_Lyso_temp, yield_Pred, nb_temp, Lag_time_Cons_temp, Lag_time_Pred_temp, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
        z_temp = deval(sol, Time_step);
        X = z_temp(mod(1:nb_temp*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:nb_temp*3,3) == 2,:); %Complexes'biomass
        W = z_temp(mod(1:nb_temp*3,3) == 0,:); %Byproducts'biomass
        R_temp = z_temp(nb_temp*3 + 1: end,:); %Byproducts'biomass
        
        z_temp = z_temp(1:(end-nb_Res), end);
        z_temp = reshape(z_temp',3, nb_temp);
        z_temp = z_temp';
        z = sum(z_temp(:,1:2), 2);%Total biomass of all species after 1 week
        nb_cell_per_mL = X/(10*Dry_set(i)*1e-15);%log10(X/(10*Dry_set(i)*1e-15));%nb cells per gram dry soil%X/(10*Dry_set(i)*1e-15);%Log values for the plot, linear values for the ttest
        end_biomass(i, j) = nb_cell_per_mL(end);%X(end);%
    end
end
Mono_Sim_1wk = end_biomass(:, 1:nb_rep);
Mono_Sim_1wk = Mono_Sim_1wk(name_sorted,:);
S = S_sim;

figG = figure;
subplot('Position',[0.3 0.1 0.2 0.62]);
mean_end_biomass = mean(end_biomass, 2);

barh(name(name_sorted), log10(mean_end_biomass(name_sorted)))
hold on
for i = 1:nb_rep
    scatter(log10(end_biomass(name_sorted, i)), 1:S_sim, 50, '.', 'r')
end
xlim([4 9]);
grid on
iFolderName = FiguresFolder;
FigName = strcat(iFolderName, name_Exp, 'individual-SynCom21-simulated-biomass-in-soil.pdf');
saveas(figG, FigName,'pdf');
print(figG, FigName, '-dpdf', '-painters');

delta_Mono = zeros(S,1);
p_value_Mono = zeros(S,1);
delta_Mono_Co = zeros(S,1);
p_value_Mono_Co = zeros(S,1);
pvals_Mono = zeros(S,1); hvals_Mono = zeros(S,1);
pvals_Co = zeros(S,1); hvals_Co = zeros(S,1);
p_val_ttest_Mono = mattest(Mono_Obs_1wk, Mono_Sim_1wk);
corrected_FDR_Mono = mafdr(p_val_ttest_Mono, 'BHFDR', true); % Estimate positive false discovery rate for multiple hypothesis testing
p_val_ttest_Co = mattest(Mono_Obs_1wk, Comm_1wk);
corrected_FDR_Co = mafdr(p_val_ttest_Co, 'BHFDR', true); % Estimate positive false discovery rate for multiple hypothesis testing
for i = 1:S
    [pvals_Mono(i), hvals_Mono(i)] = ranksum(Mono_Obs_1wk(i,:)', Mono_Sim_1wk(i,:)');
    [pvals_Co(i), hvals_Co(i)] = ranksum(Mono_Obs_1wk(i,:)', Comm_1wk(i,:)');
    delta_Mono(i) = cliffsDelta(Mono_Obs_1wk(i,:), Mono_Sim_1wk(i,:));
    p_value_Mono(i) = Permutation_Test(Mono_Obs_1wk(i,:), Mono_Sim_1wk(i, :), delta_Mono(i));
    delta_Mono_Co(i) = cliffsDelta(Mono_Obs_1wk(i,:), Comm_1wk(i,:));
    p_value_Mono_Co(i) = Permutation_Test(Mono_Obs_1wk(i,:), Comm_1wk(i,:), delta_Mono_Co(i));
end
ind_syncom_test_Cliff_1wk_Fig2 = table( ...
    name(flip(name_sorted)), p_val_ttest_Mono, corrected_FDR_Mono, p_val_ttest_Co, corrected_FDR_Co, ...
    'VariableNames', { ...
    'Name','p_val_ttest_Mono','FDR_Mono','p_val_ttest_Co','FDR_Co'});
filename = fullfile(dataFolder, 'ttest_stat_Fig2.csv');
writetable(ind_syncom_test, filename);

%% Comparison between SynCom20-21.
%We compare the time point 21 days of the SynCom20 to time point 21 days of
%the Syncom21 (this can be modified). Can be run independently from the
%other sections

% data path for individual growth
scriptFolder = fileparts(mfilename('fullpath'));
projectFolder = fileparts(scriptFolder);
dataFolder = strcat(fullfile(projectFolder, 'Data'), '/');
FiguresFolder = strcat(fullfile(projectFolder, 'Figures'), '/');

name_Exp = 'Zenodo_Senka_Average_Interaction';
%actual absolute abundance data, has quadriplicate values for each time point
Data_Evol_21 = readtable(strcat(dataFolder,'S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 7, 'Range','1:22', 'Format','auto'); %Senka's data
Data_Evol_20 = readtable(strcat(dataFolder,'S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 4, 'Range','1:21', 'Format','auto');
Parameters_set = readtable(strcat(dataFolder,'MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
Data_Evol_20 = Data_Evol_20(:, 1:53);
newRow = Data_Evol_20(1,:);             
newRow{1,1} = cellstr('Lysobacter'); 
newRow{1,2:end} = 0;
Data_Evol_20 = [Data_Evol_20(1:5, :); newRow; Data_Evol_20(6:end, :)];
nb_rep_20 = 4;
nb_rep_21 = 4;
biomass = readtable(strcat(dataFolder, 'MergedData.xlsx'),'sheet','Growth kinetic params with rep','Range','U3:X24','VariableNamingRule','preserve');
strains = readtable(strcat(dataFolder, 'MergedData.xlsx'),'sheet','Growth kinetic params with rep','Range','C3:C24','VariableNamingRule','preserve');
combined = [strains, biomass]; %Individual biomass in soil
combined = sortrows(combined,'Name','ascend');
combined = [combined, table(mean(table2array(combined(:,2:5)),2))];
combined = renamevars(combined,{'Var1'},{'mean'});
single_cell_mass = readtable(strcat(dataFolder,'MergedData.xlsx'), 'Sheet', 'Growth kinetic params with rep','Range','C3:D24','Format','auto');
single_cell_mass = sortrows(single_cell_mass,'Name','ascend');

name = string(table2array(Parameters_set(1:21,1)));
[name, idx2]= sortrows(name);
S = length(name);

% name_wo_Lysobacter = string(table2array(Parameters_set(1:20,1)));
% [name_wo_Lysobacter, idx2_wo_Lysobacter]= sortrows(name_wo_Lysobacter);
% S_wo_Lysobacter = length(name_wo_Lysobacter);

%sort alphabetically to strains. idx2 previously computed
Data_Evol_21 = Data_Evol_21(idx2,:);
Data_Evol_20 = Data_Evol_20(idx2,:);
Data_Evol_21_wo_Lysobacter = Data_Evol_21;
% Data_Evol_21_wo_Lysobacter{11,2:end} = 0;
Data_Evol_20_wo_Lysobacter = Data_Evol_20;
% Data_Evol_21_wo_Lysobacter = Data_Evol_21_wo_Lysobacter(idx2_wo_Lysobacter,:);
% Data_Evol_20_wo_Lysobacter = Data_Evol_20_wo_Lysobacter(idx2_wo_Lysobacter,:);
Data_Evol_21_wo_Lysobacter_end = table2array(Data_Evol_21_wo_Lysobacter(:,(end - nb_rep_21 + 1):end));
Data_Evol_20_wo_Lysobacter_end = table2array(Data_Evol_20_wo_Lysobacter(:,(end - nb_rep_20 + 1):end));
Data_Evol_21_wo_Lysobacter_end = Data_Evol_21_wo_Lysobacter_end;%./sum(Data_Evol_21_wo_Lysobacter_end).*sum(Data_Evol_20_wo_Lysobacter_end); %.*21/20;%Normalization to correspond to SynCom20 biomass

mean_obs_biomass_21 = mean(table2array(Data_Evol_21(:, (end - nb_rep_21 + 1):end)), 2); 
mean_obs_biomass_21 = [table(name), array2table(mean_obs_biomass_21), Data_Evol_21(:, (end - nb_rep_21 + 1):end)];
mean_obs_biomass_20 = mean(table2array(Data_Evol_20(:,(end - nb_rep_20 + 1):end)), 2); 
mean_obs_biomass_20 = [table(name), array2table(mean_obs_biomass_20), Data_Evol_20(:, (end - nb_rep_20 + 1):end)];

%make a detour to have the strain names from top to bottom in alphabetical order
inv_names = sortrows(combined,'Name','descend');
xvalues = categorical(inv_names.Name);
X = reordercats(xvalues, inv_names.Name);
yvalues_21 = mean_obs_biomass_21.mean_obs_biomass_21./(single_cell_mass.assumedCellDwInFgC*1e-15);
yvalues_20 = mean_obs_biomass_20.mean_obs_biomass_20./(single_cell_mass.assumedCellDwInFgC*1e-15);
yvalues_21_wo_Lysobacter = mean(Data_Evol_21_wo_Lysobacter_end, 2)./(single_cell_mass.assumedCellDwInFgC*1e-15);
yvalues_20_wo_Lysobacter = mean(Data_Evol_20_wo_Lysobacter_end, 2)./(single_cell_mass.assumedCellDwInFgC*1e-15);

%make sure to flip the yvalues to have the inverted order

figH = figure;
subplot('Position',[0.3 0.1 0.2 0.62]);

barh(X,flip(log10(yvalues_21))) %observed biomass in Syncom21
hold on
barh(X,flip(log10(yvalues_20)), 'FaceAlpha', 0.5)
%Scatter for each replicate
for i = 3:(3 + nb_rep_21 - 1)
    scatter(flip(log10(table2array(mean_obs_biomass_21(:,i)./(single_cell_mass.assumedCellDwInFgC*1e-15)))),1:S, 50, '.', 'r')
end
%Automatized this
for i = 3:(3 + nb_rep_20 - 1)
    scatter(flip(log10(table2array(mean_obs_biomass_20(:,i)./(single_cell_mass.assumedCellDwInFgC*1e-15)))),1:S, 50, '.', 'b')
end

xlim([3 9]);
grid on
iFolderName = FiguresFolder;
FigName = strcat(iFolderName, name_Exp, 'comparison-SynCom20-SynCom21-biomass-in-soil.pdf');
print(figH,FigName, '-dpdf', '-painters');
mat_21 = table2array(mean_obs_biomass_21(:,3:(3 + nb_rep_21 - 1))./(single_cell_mass.assumedCellDwInFgC*1e-15));
mat_20 = table2array(mean_obs_biomass_20(:,3:(3 + nb_rep_20 - 1))./(single_cell_mass.assumedCellDwInFgC*1e-15));


figK = figure;
subplot('Position',[0.3 0.1 0.2 0.62]);

barh(X, flip(log10(yvalues_21_wo_Lysobacter))) %observed biomass in Syncom21
hold on
barh(X,flip(log10(yvalues_20_wo_Lysobacter)), 'FaceAlpha', 0.5)
%Scatter for each replicate
for i = 1:(nb_rep_21 - 1)
    scatter(flip(log10(Data_Evol_21_wo_Lysobacter_end(:,i)./(single_cell_mass.assumedCellDwInFgC*1e-15))),1:S, 50, '.', 'r')
end
%Automatized this
for i = 1:(nb_rep_20 - 1)
    scatter(flip(log10(Data_Evol_20_wo_Lysobacter_end(:,i)./(single_cell_mass.assumedCellDwInFgC*1e-15))),1:S, 50, '.', 'b')
end

xlim([3 9]);
grid on
iFolderName =  FiguresFolder;
FigName = strcat(iFolderName, name_Exp, 'comparison-SynCom20-SynCom21-biomass-in-soil.pdf');
print(figH,FigName, '-dpdf', '-painters');
mat_21_wo_Lyso = Data_Evol_21_wo_Lysobacter_end./(single_cell_mass.assumedCellDwInFgC*1e-15);
mat_20_wo_Lyso = Data_Evol_20_wo_Lysobacter_end./(single_cell_mass.assumedCellDwInFgC*1e-15);

%Normalized 
mat_21_renorm = mat_21_wo_Lyso;

% ttest differences. Can't use this test for abundances and only 4
% replicates. Assume normality and same variances what is not likely false
ind_syncom_test = mattest(mat_21, mat_20);
FDR = mafdr(ind_syncom_test, 'BHFDR', true); % Estimate positive false discovery rate for multiple hypothesis testing

% Mann-Whitney test, much appropriate for 4 replicates. No assumption about
% the distribution
pvals = zeros(S,1);
hvals = zeros(S,1);
delta = zeros(S,1);
p_value = zeros(S,1);
for i = 1:S
    [pvals(i), hvals(i)] = ranksum(mat_21(i,:)', mat_20(i,:)');
    delta(i) = cliffsDelta(mat_21(i,:), mat_20(i,:));
    p_value(i) = Permutation_Test(mat_21(i,:), mat_20(i,:), delta(i));
end
ind_syncom_test = [table(flip(X)) table(ind_syncom_test) table(FDR) table(hvals) table(pvals) table(delta) table(p_value)];
filename = fullfile(dataFolder, 'SynCom21-SynCom20-mattest.csv');
writetable(ind_syncom_test, filename);

%Permanova and permutation-test on Bray-Curtis distance to get on statistic 

% A, B: n x m (m=4)
m = nb_rep_20;
X = [mat_20.'; mat_21.'];          % (2m)xn   rows = replicates, cols = species
group = [ones(m,1); 2*ones(m,1)];  % labels: 1 = SynCom20, 2 = SnyCom21
X = log1p(X);  % log(1+x)
D = BrayCurtisDistance(X);
nperm = 9999; %Choice to get the smallest p-value = 0.0001
[F, R2, p] = permanova1(D, group, nperm); %R2 = SS_bet/SS_tot, fraction mult variation explained by the difference between the 2 communities
%Permutation test entire community
p_val_tot = Permutation_Test_Comm(D, group, nperm);


mean_obs_biomass_21_renorm = table2array(mean_obs_biomass_21(:, 3:end)); %With Lysobacter %Data_Evol_21_wo_Lysobacter_end;%Without Lysobacter% 
diff_vector = mean_obs_biomass_21_renorm - table2array(mean_obs_biomass_20(:, 3:end));
mean_diff_vector = mean(mean_obs_biomass_21_renorm, 2) - mean_obs_biomass_20.mean_obs_biomass_20;%mean_obs_biomass_21.mean_obs_biomass_21 - mean_obs_biomass_20.mean_obs_biomass_20;
diff_fraction = mean_diff_vector./mean_obs_biomass_21.mean_obs_biomass_21;
ratio = 1 - mean(mean_obs_biomass_21_renorm, 2)./mean_obs_biomass_20.mean_obs_biomass_20;

close all

figH = figure;
subplot('Position',[0.1 0.3 0.4 0.3]);
stem(diff_vector);
xtickangle(90)
set(gca,'xtick',1:S,'xticklabel',combined.Name)
ylabel('biomass difference');
title('Syncom21 - Syncom20');
grid on
ylim([-2e-4 2e-4]);

iFolderName =  FiguresFolder;
FigName = strcat(iFolderName, name_Exp, 'syncom20-21-biomass-diff.pdf');
print(figH,FigName, '-dpdf', '-painters');

%%%New stacked 

figK = figure;
bar(mean_diff_vector)
hold on
scatter(1:S, diff_vector, 60, '.', 'r')
grid on
axis([0 22 -3e-04 2.5e-04])
xticks(1:length(name)); 
xticklabels(name);
ylabel('measured - simulated')

iFolderName =  FiguresFolder;
FigName = strcat(iFolderName, name_Exp, 'stacked_plot_diff_20_21.pdf');
print(figK,FigName, '-dpdf', '-painters');