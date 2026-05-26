%Simulations with fitted interspecific interactions
%use here the data_to_save_no_CF with 4 replicate simulations
%input is S20_S21_abs_abund_cfu_Senka.xlsx, Sheet 4, with 8 replicates 
%starting abundances for Bra, Phe, Mes, Cau, Coh, Tar, Bur, Chi and Muc are corrected
%Parameters set used with assumed preys: 'struct_tot_SynCom21_Soil_only_Lyso_Parfor';
%Parameters set used for predation: 'data_to_save_SynCom21_Soil_only_Lyso_Parfor_Newv6';

clear
close all

pwd_init = pwd;
pwd_init = fileparts(pwd_init);
addpath(strcat(pwd_init, '/Data'))
addpath(pwd_init)


name_Exp = 'Zenodo_Senka_Average_Interaction';
Name_file = 'struct_tot_SynCom21_Soil_only_Lyso_Parfor';%'data_to_save_SynCom21_Soil_only_Lyso_Parfor_Predation';%
Name_file_Resources_Death = 'data_to_save_SynCom21_Reduced'; %Data fitted on moncultures

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','118:139', 'Format','auto');

%actual absolute abundance data, has quadriplicate values for each time point

% Data_Evol = readtable(strcat('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data/Liquid_Data/','abund_cfu_IG.xlsx'), 'Sheet', 2, 'Range','46:67', 'Format','auto');%1:22 without correction for 0%'46:67' with correction for disappearing species. %73:94 removing disappearing species Data in soil extract from Phil and Clara's experiments
Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 7, 'Range','1:22', 'Format','auto'); %Senka's data
% Data_Evol = readtable(strcat('Data/','SSC21_genera_relative-abundances.xlsx'), 'Sheet', 4, 'Range','1:22', 'Format','auto'); %Bruna's data
% Data_Evol = readtable(strcat('Data/','SynComNX_data_SenkaCB_IG.xlsx'), 'Sheet', 4, 'Range','30:51', 'Format','auto'); %Subcommunities. modified data 57:78. Not modified 30:51. Sheets 2-4
% Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 13, 'Range','1:22', 'Format','auto'); %Phil's weekly data
% Data_Evol = Data_Evol(:, [1:3]);
mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8));
S = height(Data_Evol);
Time_step = [0 1 3 7 10 21]*24;%Senka %[0 12 22 38 70 168 504];%Bruna %[0 12 48 96 7*24 21*24];%Clara %[0 24 72 168];%Subcommunities %[0 168];%Phil's data%Measured time step in hours
Time_step_BP = 0:10:500;
nb_res = 12;
name = string(table2array(Parameters_set(1:S,1)));

%Load parameters
load(strcat('Data/', Name_file));
data_to_save = struct_tot_SynCom21_Soil_only_Lyso_Parfor;%data_to_save_SynCom21_Soil_only_Lyso_Parfor_Predation;%
load(strcat('Data/', Name_file_Resources_Death));
data_to_save_Res_Death = data_to_save_Mono;
out_struct = Weighted_Struct(data_to_save);
nb_data_set = size(data_to_save,2);
Data_Evol_temp = table2array(Data_Evol(:, 2:end));
nb_obs = length(Data_Evol_temp(1,:));
nb_time_step = length(Time_step);
nb_rep = nb_obs/nb_time_step;
nb_rep_sim = max(nb_rep, 8); %size(data_to_save,2)%Number of replicates to be simulated

% Fitted parameters
StackPlotTot = zeros(S, length(Time_step));
StackPlot_Meas_added_errors = zeros(S, length(Time_step)); %We perturbate the observed values accordingling to the simulated variance between data
z = zeros(S, nb_rep_sim);
z_no_CF = zeros(S, nb_rep_sim);
z_obs = zeros(S, nb_rep);
z_obs_added_errors = zeros(S, nb_rep_sim); %We perturbate the observed values accordingling to the simulated variance between data
exp_biomass = zeros(S, length(Time_step), nb_rep_sim);
exp_resources = zeros(nb_res, length(Time_step), nb_rep_sim);
exp_byproducts = zeros(S, length(Time_step), nb_rep_sim);
dexp_biomass = zeros(S, length(Time_step_BP), nb_rep_sim);
dexp_resources = zeros(nb_res, length(Time_step_BP), nb_rep_sim);
dexp_byproducts = zeros(S, length(Time_step_BP), nb_rep_sim);
exp_biomass_BP = zeros(S, length(Time_step_BP), nb_rep_sim);
exp_byproducts_BP = zeros(S, length(Time_step_BP), nb_rep_sim);
exp_biomass_wo_CF = zeros(S, length(Time_step), nb_rep_sim);
Measured_Abund_mod_var = zeros(S, length(Time_step), nb_rep_sim);

%Initialization of the colors
colors_ind = [4, 1, 10, 15, 12, 3, 8, 6, 20, 18, 19, 2, 5, 21, 13, 9, 14, 7, 16, 17, 11];

[sorted_name, idx] = sortrows(name);
colors = {'#B35806', '#E08214', '#D53E4F', '#B2182B', '#D6604B', '#C51B7D', ...
               '#DE77AE', '#F1B6DA', '#FDAE61', '#FEE090', '#A6D96A', '#5AAE61', ...
               '#01665E', '#35978F', '#1B7837', '#C2A5CF', '#9970AB', '#762A83', ...
               '#80CDC1', '#C7EAE5', '#2166AC', '#4393C3', '#B35806'};



Measured_Abund = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
for i = 1:nb_rep
    Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
    z_obs(:, i) = Measured_Abund(:, end, i);
end    
StackPlot_Meas = Measured_Abund./sum(Measured_Abund);
StackPlot_Meas = mean(StackPlot_Meas, 3);
%display the biomass variations among the 4 observed replicates
num_fig = 1;
figure(num_fig)
for i = 1:nb_rep
    for j = 1:S
        plot(Time_step, Measured_Abund(j,:,i), '-*', 'Color', colors{j});%, Time_step, R_temp, 'o');
        hold on
    end
end
axis([0 max(Time_step) 0 4e-04])
num_fig = num_fig + 1;
mean_y_0 = mean(Measured_Abund(:,1,:), 3);
var_mat = var(Measured_Abund, 0, 3);
Measured_Abund_average = mean(Measured_Abund,3);

figK = figure;
b = bar(mean(Measured_Abund_average(flip(idx), 2:end), 3)', 'stacked');
grid on
hold on
for i = 1:numel(b)
    b(i).FaceColor = colors{21 - i + 1};
end
axis([0 (nb_time_step + 1) 0 6e-04])
legend(flip(sorted_name), 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked observed')
iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'Transcriptomics_BarPlot.pdf');
print(figK,FigName, '-dpdf', '-painters');
close all

%Setting for Matlab ODE solver
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-12);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
[Shann_sim, Simp_sim, Shann_obs, Simp_obs] = deal(zeros(nb_rep_sim, nb_time_step));
stacked_plot_diff_CF = zeros(S, nb_rep_sim);
for zz = 1:nb_rep_sim %nb_rep %nb_data_set

    % out_struct = Weighted_Struct(data_to_save); %Compute weigthed interaction over all the datasets

    ind_sim = randi(nb_data_set);
    kappa_mat = data_to_save(ind_sim).kappa_mat;
    CrossFeed_Mat = data_to_save(ind_sim).CrossFeed_Mat_Temp;
    Death_Mat_Temp = data_to_save(ind_sim).Death_Mat_Temp; %data_to_save_Res_Death.Death_Mat_Temp; %data_to_save(ind_sim).Death_Mat_Temp; 
    Threshold_CF = data_to_save(ind_sim).Threshold_CF;
    Threshold_death = data_to_save(ind_sim).Threshold_death; 
    Lag_time_Cons = data_to_save(ind_sim).Lag_time_Cons; %data_to_save_Res_Death.Lag_time_Cons; %data_to_save(ind_sim).Lag_time_Cons; 
    Lag_time_Pred = data_to_save(ind_sim).Lag_time_Pred;
    Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
    death_rate =  data_to_save(ind_sim).death_rate;
    Pred_Mat_Lyso = data_to_save(ind_sim).Pred_Mat_Lyso;
    Threshold_Pred = data_to_save(ind_sim).Threshold_Pred;
    R = data_to_save(ind_sim).R; %0.475 subcommunities wo Lysobacter or Mucilaginibacter. 0.6 wo Pseudomonas1. Reduce 1.25 + init Senka. 1.5 just reduce %1.6 for Phil's data
    nb_Res = length(R); %Number of resources (group of resources)
    Resource_Matrix = data_to_save(ind_sim).Resource_Matrix; %data_to_save_Res_Death.Resource_Matrix; %data_to_save(ind_sim).Resource_Matrix; %Addition of a line for Lysobacter

    % kappa_mat = out_struct.kappa_mat;
    % CrossFeed_Mat = out_struct.CrossFeed_Mat_Temp;
    % Death_Mat_Temp = out_struct.Death_Mat_Temp; %data_to_save_Res_Death.Death_Mat_Temp; %data_to_save(ind_sim).Death_Mat_Temp; 
    % Threshold_CF = out_struct.Threshold_CF;
    % Threshold_death = out_struct.Threshold_death; 
    % Lag_time_Cons = out_struct.Lag_time_Cons; %data_to_save_Res_Death.Lag_time_Cons; %data_to_save(ind_sim).Lag_time_Cons; 
    % Lag_time_Pred = out_struct.Lag_time_Pred;
    % Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
    % name = string(table2array(Parameters_set(1:21,1)));
    % Pred_Mat_Lyso = out_struct.Pred_Mat_Lyso;
    % Threshold_Pred = out_struct.Threshold_Pred;
    % R = 0.4*out_struct.R; %0.4 subcommunities wo Lysobacter or Mucilaginibacter. 0.5 wo Pseudomonas1. Reduce 1.25 + init Senka. 1.5 just reduce
    % nb_Res = length(R); %Number of resources (group of resources)
    % Resource_Matrix = data_to_save(ind_sim).Resource_Matrix; %data_to_save_Res_Death.Resource_Matrix; %data_to_save(ind_sim).Resource_Matrix; %Addition of a line for Lysobacter
    
    % Measured_Abund = table2array(Data_Evol(1:20, 2:7));
    % Number of surviving species after 8 weeks
    Threshold_Surviving = 1e-10;
    nb_Surv_Obs = sum(Measured_Abund_average > Threshold_Surviving);
    
    %Initialization of the model parameters fixed for all replicates 
    t_0 = 0; %Time 0
    tspan = [0, max(Time_step)*24]; %Time interval in hours
    yield_Pred = 0.2;% 20% of yield for predation
    
    num_fig = 3;
    %Creation of a starting point based on a weighted combination of the
    %observed initial biomass
    rand_rep = rand(1, nb_rep);%randi(5, 1, nb_rep) - 1;
    rand_rep = rand_rep./min(sum(rand_rep));
    noise = normrnd(0, 1e-8, S, 1);
    mean_y_0 = sum(rand_rep.*squeeze(Measured_Abund(:, 1, :)), 2);
    noise(mean_y_0 == 0) = 0;
    mean_y_0 = mean_y_0 + noise;
    mean_y_0(mean_y_0 < 0) = 0;
    % mean_y_0(16) = 0; %For Bruna's data, no microbacterium

    %Initial concentrations using a normal distribution
    mat_y_0 = mean_y_0;
    
    mat_y_0 = [mat_y_0 zeros(S,1) zeros(S,1)];


    y_0 = sum(mat_y_0(:,1:2),2);
    mat_y_0 = reshape(mat_y_0', 1, []);
    mat_y_0 = [mat_y_0 R];

    sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat, CrossFeed_Mat, Mat_kappa_3, Resource_Matrix,...
        Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp, death_rate,...
        Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
    dz = @(t,y) fun_CF_Death_Lyso(t, y, kappa_mat, CrossFeed_Mat, Mat_kappa_3, Resource_Matrix, ...
        Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp, death_rate, ...
        Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10);
    
    z_temp = deval(sol, Time_step);
    X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
    P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
    W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
    R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass

    z_temp_BP = deval(sol, Time_step_BP);
    X_BP = z_temp_BP(mod(1:S*3,3) == 1,:); %Species'biomass
    W_BP = z_temp_BP(mod(1:S*3,3) == 0,:); %Byproducts'biomass

    %%%No CFP
    sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat, 0.*CrossFeed_Mat, 0.*Mat_kappa_3, Resource_Matrix,...
        0.*Threshold_CF, Threshold_death, 0.*Threshold_Pred, Death_Mat_Temp, death_rate,...
        0.*Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
    z_temp_no_CF = deval(sol, Time_step);
    X_no_CF = z_temp_no_CF(mod(1:S*3,3) == 1,:); %Species'biomass
    P_no_CF = z_temp_no_CF(mod(1:S*3,3) == 2,:); %Complexes'biomass
    X_no_CF = X_no_CF + P_no_CF;

    
    dz_temp = zeros(size(z_temp_BP));   % same size
    
    for j = 1:nb_time_step
        dz_temp(:,j) = dz(Time_step(j), z_temp_BP(:,j));
        dX = dz_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        dP = dz_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
        dW = dz_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
        dR_temp = dz_temp(S*3 + 1: end,:); %Byproducts'biomass
    end


    figure(num_fig)
    for j = 1:S
        plot(Time_step, X(j,:), '-o', 'Color', colors{j});%, Time_step, R_temp, 'o');
        hold on
    end
    legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
    axis([0 max(Time_step) 0 4e-04]) %axis([0 600 0 2e-03])
    num_fig = num_fig + 1;
    figure(num_fig)
    Measured_Abund_disp = Measured_Abund_average + normrnd(zeros(S, length(Time_step)), sqrt(var_mat));
    Measured_Abund_disp(Measured_Abund_disp < 0) = 0;
    for j = 1:S
        plot(Time_step, Measured_Abund_disp(j,:), '-*', 'Color', colors{j});%, Time_step, R_temp, 'o');
        hold on
    end
    legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
    axis([0 max(Time_step) 0 4e-04])

    num_fig = num_fig + 1;
    z_temp = z_temp(1:(end-nb_Res), end);
    z_temp = reshape(z_temp', 3, S);
    z_temp = z_temp';
    z(:, zz) = X(:, end);%X(:, 4);%Total biomass of all species after 1 week or at the end of the experiement
    z_no_CF(:, zz) = X_no_CF(:, end);
    z_obs_added_errors(:,zz) = Measured_Abund_disp(:, end);
    StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
    %Diversity indices
    [Shann_sim(zz, :), Simp_sim(zz, :)] = Shannon_Simpson_Indices(S, StackPlot);
    [Shann_obs(zz, :), Simp_obs(zz, :)] = Shannon_Simpson_Indices(S, mean(StackPlot_Meas,3));
    u = X(:, end) - X(:, 1); v = Measured_Abund(:, end) - Measured_Abund(:,1);
    % angle(:, zz) = u
    StackPlotTot = StackPlotTot + X./sum(X);
    StackPlot_Meas_added_errors = StackPlot_Meas_added_errors + Measured_Abund_disp./sum(Measured_Abund_disp);

    exp_biomass(:,:, zz) = X;
    exp_resources(:,:, zz) = R_temp; 
    exp_byproducts(:,:, zz) = W;
    exp_biomass_BP(:,:, zz) = X_BP;
    exp_byproducts_BP(:,:, zz) = W_BP;
    dexp_biomass(:,:, zz) = dX;
    dexp_resources(:,:, zz) = dR_temp; 
    dexp_byproducts(:,:, zz) = dW;
    exp_biomass_wo_CF(:,:, zz) = X_no_CF;
    stacked_plot_diff_CF(:, zz) = Measured_Abund_average(:, end) - X(:, end);
    Measured_Abund_mod_var(:,:, zz) = Measured_Abund_disp;
end
mean_Shann_sim = mean(Shann_sim); mean_Shann_obs = mean(Shann_obs);
mean_Simp_sim = mean(Simp_sim); mean_Simp_obs = mean(Simp_obs);

%%%New stacked figure. Diff observed vs modelled 
mean_stacked_plot_CF = mean(stacked_plot_diff_CF, 2);

%Permanova and permutation-test on Bray-Curtis distance to get on statistic 

m_1 = nb_rep;
m_2 = nb_rep_sim;
mat_Comm_1 = z_obs;
mat_Comm_2 = z;
X = [mat_Comm_1.'; mat_Comm_2.'];        % (m_1 + m_2)xn   rows = replicates, cols = species
group = [ones(m_1, 1); 2*ones(m_2, 1)];  % labels: 1 = SynCom20, 2 = SnyCom21
X = log1p(X);  % log(1+x)
D = BrayCurtisDistance(X);
[F_CF_2, R2_CF_2, p_CF_2] = permanova1_exact(D, group);
%Permutation test entire community
nperm = 9999;
p_val_tot_CF = Permutation_Test_Comm(D, group, nperm);

m_1 = nb_rep;
m_2 = nb_rep_sim;
mat_Comm_1 = z_obs;
mat_Comm_2 = z_no_CF;
X = [mat_Comm_1.'; mat_Comm_2.'];        % (m_1 + m_2)xn   rows = replicates, cols = species
group = [ones(m_1, 1); 2*ones(m_2, 1)];  % labels: 1 = SynCom20, 2 = SnyCom21
X = log1p(X);  % log(1+x)
D = BrayCurtisDistance(X);
[F_1_no_CF, R2_no_CF, p_no_CF] = permanova1_exact(D, group); %R2 = SS_bet/SS_tot, fraction mult variation explained by the difference between the 2 communities
%Permutation test entire community
p_val_tot_no_CF = Permutation_Test_Comm(D, group, nperm);

%MDS analysis
XB = reshape(permute(exp_biomass, [2 3 1]), [], S);  % (mA*nb_time_step)xn
XA = reshape(permute(Measured_Abund, [2 3 1]), [], S);  % (mB*nb_time_step)xn
XC = reshape(permute(exp_biomass_wo_CF, [2 3 1]), [], S);  % (mB*nb_time_step)xn
X  = [XA; XB; XC];

%metadata
env  = [repmat("A", nb_rep*nb_time_step, 1); repmat("B", nb_rep_sim*nb_time_step, 1); repmat("C", nb_rep_sim*nb_time_step, 1)];
repA = repelem((1:nb_rep)', nb_time_step);                 
repB = repelem((1:nb_rep_sim)', nb_time_step) + nb_rep;     
repC = repelem((1:nb_rep_sim)', nb_time_step) + nb_rep + nb_rep_sim;  
rep  = [repA; repB; repC];
time = repmat((1:nb_time_step)', nb_rep + 2*nb_rep_sim, 1);          
% transform (recommended with Bray–Curtis + biomass)
X = log1p(X);

% Bray–Curtis distances + PCoA
D = BrayCurtisDistance(X);
[Y, eigvals] = cmdscale(D);  % Convert Bray-Curtis distance into Euclidean coordinates
env_A = env(env == 'A'); env_B = env(env == 'B'); env_C = env(env == 'C');
Y_AB = [Y(env == 'A', :); Y(env == 'B', :)]; Y_AC = [Y(env == 'A', :); Y(env == 'C', :)];
rep_AB = [rep(env == 'A'); rep(env == 'B')]; rep_AC = [rep(env == 'A'); rep(env == 'C')];
time_AB = [time(env == 'A'); time(env == 'B')]; time_AC = [time(env == 'A'); time(env == 'C')];
[PCoA_Metric_AB, angle_PCoA_AB, p_val_AB] = TrajectoryDist(Y_AB, rep_AB, time_AB, [env_A; env_B]);
[PCoA_Metric_AC, angle_PCoA_AC, p_val_AC] = TrajectoryDist(Y_AC, rep_AC, time_AC, [env_A; env_C]);

eigpos = eigvals(eigvals > 0);
explained = eigpos/sum(eigpos);
disp(explained(1:5))
disp(cumsum(explained(1:5)))

k = min(5, size(Y,2));
tbl = table(env, rep, time);
for i = 1:k
    tbl.(sprintf("PC%d", i)) = Y(:,i);
end

% Mixed model per axis (repeat for i=1...k)
% Fixed: env*time, Random: replicate intercept
lme1 = fitlme(tbl, "PC1 ~ env*time + (1|rep)");
disp(anova(lme1));

figG = figure;
gscatter(Y(:,1), Y(:,2), env)
xlabel('PCoA1')
ylabel('PCoA2')
title('PCoA (Bray-Curtis)')
axis equal
hold on

envs = unique(env);
colors_MDS = lines(length(envs));

for e = 1:length(envs)
    idx_env = env == envs(e);   
    reps = unique(rep(idx_env));   
    for r = reps'
        temp = (env == envs(e)) & (rep == r);
        
        [~, order] = sort(time(temp));
        coords = Y(temp,:);
        coords = coords(order,:);
        
        alpha_line = 0.3;
        alpha_marker = 0.8;
        
        % trajectory line
        plot(coords(:,1), coords(:,2), '--', ...
             'Color', [colors_MDS(e,:) alpha_line], ...
             'LineWidth', 1.5);
        
        % markers
        scatter(coords(:,1), coords(:,2), 40, ...
                'MarkerFaceColor', colors_MDS(e,:), ...
                'MarkerEdgeColor', colors_MDS(e,:), ...
                'MarkerFaceAlpha', alpha_marker, ...
                'MarkerEdgeAlpha', alpha_marker);
    end
end

xlabel('PCo1')
ylabel('PCo2')
title('Community trajectories in PCoA space')
axis equal
legend(envs)
iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'SimCF_noCF_MDS_Plot.pdf');
print(figG, FigName, '-dpdf', '-painters');

figK = figure;
bar(mean_stacked_plot_CF(idx, :))
hold on
scatter(1:S , stacked_plot_diff_CF(idx, :),60,'.','r')
grid on
axis([0 22 -3e-04 0.7e-04]) %1.5e-04%1e02% For prop analysis
xticks(1:length(sorted_name)); 
xticklabels(sorted_name);
ylabel('measured - simulated')

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'stacked_plot_diff_CF.pdf');
print(figK,FigName, '-dpdf', '-painters');

var_data = var(exp_biomass, 0, 3);

StackPlotTot = StackPlotTot/zz;
StackPlot_Meas_added_errors = StackPlot_Meas_added_errors/zz;

z_fin_sim = mean(z, 2); %Absolute stationary abundances
z_fin_obs_added_errors = mean(z_obs_added_errors, 2);%z_fin_obs = Measured_Abund(1:S, end); %Absolute stationary abundances
nb_Surv_Sim = sum(StackPlotTot > Threshold_Surviving);
disp(sum(z_fin_sim)/sum(Measured_Abund_average(:,end)))

h_vect = zeros(1,S);
p_vect = zeros(1,S);
for i = 1:S
    [h_vect(i),p_vect(i)] = ttest2(StackPlot_Meas_added_errors(i,:)', StackPlotTot(i,:)');
end
% name_diff = name(logical(h_vect));

figure(num_fig);
bar(StackPlotTot', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked all replicates')
num_fig = num_fig + 1;

figure(num_fig); 
bar(mean(StackPlot_Meas_added_errors, 3)', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked observed')
num_fig = num_fig + 1;

[diff_vect, I, fact] = Sort_diff(StackPlotTot(:,end), StackPlot_Meas_added_errors(:,end));
name_order = name(I);

figure(num_fig);
for j = 1:S
    scatter(StackPlot_Meas_added_errors(j, end), StackPlotTot(j, end), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors{j}, 'MarkerFaceColor', colors{S+1-j});%col(S+1-j,:))      
    hold on
end
axis([0 0.7 0 0.7]);
axis square
reflin = refline(1,0);
reflin.Color = 'r';
xlabel('Experiment'); 
ylabel('Simulation');
legend(name);
title('Scatter Tot');

num_fig = num_fig + 1;

figure(num_fig);
for j = 1:S
    scatter(z_fin_obs_added_errors(j), z_fin_sim(j), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors{j}, 'MarkerFaceColor', colors{S+1-j});%col(S+1-j,:))      
    hold on
end
axis([0 3*10^(-4) 0 4*10^(-4)]);
reflin = refline(1,0);
axis square
reflin.Color = 'r';
xlabel('Experiment'); 
ylabel('Simulation');
legend(name);
title('Scatter absolute abundance');
num_fig = num_fig + 1;

figure(num_fig);
plot(1:length(Time_step), Simp_obs, 'b--o')
hold on
plot(1:length(Time_step), Simp_sim, 'r--o')
num_fig = num_fig + 1;

iFolderName = strcat(cd, '/Figures/');
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName = num2str(get(FigHandle, 'Number'));
    FigName = strcat('Fig', FigName);
    FigName = strcat(iFolderName, FigName, name_Exp, 'SynCom_21_');
    set(0, 'CurrentFigure', FigHandle);
    saveas(FigHandle, FigName, 'pdf');
end

%Save measured and obs average 
Struct_Temp = struct();

Struct_Temp.(strcat(name_Exp, 'StackPlot_Meas')) = StackPlot_Meas;
Struct_Temp.(strcat(name_Exp, 'StackPlot_Tot'))  = StackPlotTot;

save(fullfile(cd, 'Data', [name_Exp 'Stacked_abundances.mat']), '-struct', 'Struct_Temp');
%% Computation of the rates distribution for each species

% [~, nb_iter] = size(data_to_save_SynCom21_Soil_only_Lyso_Parfor_7);
% Mat_Struct = {data_to_save_SynCom21_Soil_only_Lyso_Parfor_7.CrossFeed_Mat_Temp};
% [Sample_rate, LN_fit, nb_neg_val] = Rates_Distribution(12, Mat_Struct, num_fig, name);

%% plot figure with the observed or simulated biomass per species at the end

%cd('/Users/jvanderm/Library/CloudStorage/OneDrive-UniversitédeLausanne/SynCom model paper/Data')

close all
figK = figure;
subplot('Position',[0.3 0.2 0.62 0.2]);
bar(mean(Measured_Abund(:, end, :), 3))

hold on

for i = 1:nb_rep %number of replicates
    scatter(1:S, Measured_Abund(:, end, i),60,'.','r')
end

set(gca, 'YScale', 'log')
ylim([1e-7 3e-4])
yticks([1e-6 1e-5 1e-4]);

xticks(1:length(name)); 
xticklabels(name);
grid on
ylabel('Biomass g/mL')
title('observed');

subplot('Position',[0.3 0.5 0.62 0.2]);

bar(mean(z,2))
hold on
for i=1:nb_rep %number of replicates
 scatter(1:S, z(:,i), 60, '.', 'k')
end

set(gca, 'YScale', 'log')
ylim([1e-7 3e-4])
grid on
ylabel('Biomass g/mL')
yticks([1e-6 1e-5 1e-4]);
title('simulated');

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'expected-end-biomass-SynCom20-in-soil-data-CF.pdf');
% saveas(figK,FigName,'pdf');
print(figK,FigName, '-dpdf', '-painters');


% sort of double scatter plot for simulated and experimental data
figL = figure;
subplot('Position',[0.1 0.1 0.55 0.55]);

for i=1:nb_rep %number of replicates
    scatter(median(z,2),  Measured_Abund(:, end, i), 60,'.','k')
    hold on
end
axis([0 2e-4 0 2e-4])
axis square
hold on
for i = 1:nb_rep_sim %nb_data_set 
    y =  median(Measured_Abund(:, end, :),3);
    scatter(z(:,i), y, 60,'.','b')
    hold on
end

set(gca, 'YScale', 'log', 'XScale','log')
ylim([1e-7 3e-4])
xlim([1e-7 3e-4])
reflin = refline(1,0);
axis square
reflin.Color = 'r';
ylabel('Experiment'); 
xlabel('Simulation');
grid on
xticks([1e-6 1e-5 1e-4]);

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'expected-end-biomass-SynCom21-in-soil-scatter.pdf');
% saveas(figL,FigName,'pdf');
print(figL,FigName, '-dpdf', '-painters');

%individual boxcharts

close all

figL = figure;
subplot('Position',[0.1 0.1 0.55 0.55]);
boxchart(z')
set(gca, 'YScale', 'log')
ylim([1e-7 3e-4])
%saveas(figL,'boxchart-z.pdf','pdf');


close all

figK = figure;

tiledlayout(3,7)

ax = gobjects(S,1);
Time_step_days = Time_step;%/24;
for i = 1:S  
    ax(i) = nexttile(i);

    for k = 1:nb_rep_sim
        plot(Time_step_days, exp_biomass(i, :, k),'m--');
        hold on
    end
    hold on

    for kk = 1:nb_rep %number of replicate observations
        plot(Time_step_days, Measured_Abund(i, :, kk),'g--');
    end

    set(ax(i), 'YScale', 'log');
    ylim(ax(i), [1e-9 1e-2]);  
    ylim(ax(i), [0 2e-4]); 
    xlim([0 max(Time_step_days)]);
    % yticklabels({}); %turn off the ytick labels, but they go from 1e-8 to 1e-4
    grid on
    bact_name = name(i);
    title(bact_name);
end

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'no-log-ind-growth-plots-model-data_CF-start.pdf');
saveas(figK,FigName,'pdf');
print(figK, FigName, '-dpdf', '-painters');

%% display mean and stdev of the cross-feeding matrices over time

tmp = zeros(S, S, nb_rep_sim);

for i = 1:nb_rep_sim
    tmp (:,:,i) = data_to_save(i).CrossFeed_Mat_Temp;
end

mean_CF = mean(tmp,3);
var_CF = std(tmp,[],3);

close all

figK = figure;

subplot('Position',[0.1 0.1 0.3 0.3]);

imagesc(mean_CF');
colorbar

title('crossfeeding');

subplot('Position',[0.5 0.1 0.3 0.3]);

imagesc(var_CF');
colorbar

title('crossfeeding std');

%saveas(figK,'heatmap-sim-cross-feeding-data_to_save_w_corr.pdf','pdf');
%%Useless if considering SynCom21 because the variations only occur for
%%Lysobacter

% multi stack bar plot of relative abundances in the simulations
% compared to the mean of 8 observations

[strains, idx] = sortrows(name); %to sort alphabetically
mean_CF = mean_CF(idx,:); %sort the byproducts rows accordingly
mean_CF = mean_CF(:,idx); %sort the columns accordingly

var_CF=var_CF(idx,:);
var_CF=var_CF(:,idx);

close all

figK=figure;

subplot('Position',[0.15 0.2 0.3 0.4]);
h = heatmap(strains, strains, mean_CF, Colormap=hot(8));

h.Title = 'mean cross-feeding';

subplot('Position',[0.6 0.2 0.3 0.4]);
h = heatmap(var_CF, Colormap=pink);

h.Title = 'cross-feeding standard error';

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'heatmap-names-sim-cross-feeding-data_to_save_w_corr2.pdf');
saveas(figK,FigName,'pdf');
print(figK, FigName, '-dpdf', '-painters');

close all

figK = figure;

%tiledlayout(3,7)
tiledlayout(6,7)

for i = 1:nb_rep_sim
    nexttile(i)
    b = bar((exp_biomass(:,:,i)./sum(exp_biomass(:,:,i)))','stacked','FaceColor','flat');%b = bar((exp_biomass{i}./sum(exp_biomass{i}))','stacked','FaceColor','flat');
    for j=1:S
        b(j).FaceColor = colors{j};
    end
    title(strcat('simulation',num2str(i)));
end

StackPlot_Meas_m = mean(StackPlot_Meas, 3);

nexttile
b = bar(StackPlot_Meas_m','stacked','FaceColor','flat');
for j=1:S
    b(j).FaceColor = colors{j};
end
title('mean obs');

legend(name, 'Orientation', 'vertical', 'Location', 'none','Position',[0.75 0.1 0.2 0.4])

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'rel-abund-model-data-CF.pdf');
% saveas(figK,FigName,'pdf');
print(figK, FigName, '-dpdf', '-painters');

%%simple plot of the maximum biomass attained in the simulations versus observations

sim_fin_biomass = zeros(nb_rep_sim, 1);
for i = 1:nb_rep_sim
    sim_fin_biomass(i,:) = max(sum(exp_biomass(:,:,i)));
    % sim_fin_biomass(i,:) = max(sum(exp_biomass(:,end,i)));
end

sim_fin_biomass_wo_CF = zeros(nb_rep_sim, 1);
for i = 1:nb_rep_sim
    sim_fin_biomass_wo_CF(i,:) = max(sum(exp_biomass_wo_CF(:,:,i)));
    % sim_fin_biomass_wo_CF(i,:) = max(sum(exp_biomass_wo_CF(:,end,i)));
end

%repeat calculation of the measured abundances
obs_fin_biomass = zeros(nb_rep, 1);
for i = 1:nb_rep
    obs_fin_biomass(i,:) = max(sum(Measured_Abund(:,:,i)));
end

close all
figH = figure;

subplot('Position',[0.1 0.1 0.2 0.4])

y = [mean(sim_fin_biomass); mean(sim_fin_biomass_wo_CF); mean(obs_fin_biomass)];
disp(mean(sim_fin_biomass)/mean(obs_fin_biomass))
disp(mean(sim_fin_biomass_wo_CF)/mean(obs_fin_biomass))

X = categorical({'Sim','Sim_no_CFP','Obs'});
X = reordercats(X,{'Sim','Sim_no_CFP', 'Obs'});
b = bar(X,y,'FaceColor','flat');
set(gca,'FontSize',7);
b.CData(1,:) = [1 0 1];
b.CData(2,:) = [0 0 1];
b.CData(2,:) = [0 1 1];

hold on
for j = 1:nb_rep_sim
    scatter(X(1), sim_fin_biomass(j),50,'k','.')
end
hold on
for j = 1:nb_rep_sim
    scatter(X(2), sim_fin_biomass_wo_CF(j),50,'k','.')
end
hold on
for j = 1:nb_rep
    scatter(X(3), obs_fin_biomass(j),50,'k','.')
end
ylabel('biomass (gC DW/mL)')
ylim([0 1e-3])

delta = cliffsDelta(sim_fin_biomass, obs_fin_biomass);
p_value = Permutation_Test(sim_fin_biomass, obs_fin_biomass, delta);

%No CF
delta_wo_CF = cliffsDelta(sim_fin_biomass_wo_CF, obs_fin_biomass);
p_value_wo_CF = Permutation_Test(sim_fin_biomass, obs_fin_biomass, delta_wo_CF);
title('Max biomass','FontSize',9);

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'max-data-CF-sim-obs-biomass.pdf');
saveas(figH,FigName,'pdf');
print(figH,FigName, '-dpdf', '-painters');

[Sortedname, idx] = sort(name);
z = z(idx, :);
save(strcat('Data/',name_Exp, 'z'), 'z')

%% plot resource usage

close all

for i = 1:nb_rep_sim
p = plot(Time_step, exp_resources(:,:,i)');
	for j = 1:nb_Res
	    p(j).Color = colors{j};
	end
hold on
end


%% plot byproduct formation by strain for all replicates

close all

figK = figure;

% tiledlayout(4,5)

tiledlayout(3, 7)
% tiledlayout(6, 7)
sum_BP = squeeze(sum(exp_byproducts_BP));
sum_dBP = squeeze(sum(dexp_byproducts));
sum_X = squeeze(sum(exp_biomass_BP));
sum_dX = squeeze(sum(dexp_biomass));

ax = gobjects(S,1);
Time_step_days = Time_step/24;
exp_byproducts(exp_byproducts == 0) = 1e-15;
for i=1:S

    nexttile(i)
    ax(i) = nexttile(i);
	for k = 1:nb_rep_sim
	    % plot(Time_step, exp_byproducts{k}(i,:),'m-');
        plot(Time_step_days, exp_byproducts(i,:, k),'m-');
        set(ax(i), 'YScale', 'log');
        ylim(ax(i), [1e-07 1e-2]);  
        xlim([0 max(Time_step_days)]);
        hold on
	end
	hold on

    grid on
    bact_name=name(i);
    title(bact_name);
end

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'Byproduct_Growth2.pdf');
saveas(figK,FigName,'pdf');
print(figK,FigName, '-dpdf', '-painters');

% as heatmap with the mean of the byproducts at each time point
tmp = exp_byproducts;

mean_byproducts = mean(tmp, 3);

[strains, idx] = sortrows(name); %to sort alphabetically
mean_byproducts = mean_byproducts(idx,:); %sort the byproducts matrix accordingly

close all

figK=figure;

subplot('Position',[0.2 0.2 0.2 0.55]);
h = heatmap(Time_step/24, strains, mean_byproducts);

h.Title = 'Byproduct formation';
h.XLabel = 'Time (days)';


iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'byproducts-heatmap.pdf');
% saveas(figK,FigName,'pdf');
print(figK, FigName, '-dpdf', '-painters');

%% sum byproducts in the system plus the modeled resource utilization

figK = figure;

subplot('Position',[0.1 0.1 0.3 0.3]);

	for k = 1:nb_rep_sim
	    byproduct_sum = sum(exp_byproducts(:, :, k),1);
	    resource_sum = sum(exp_resources(:, :, k),1);
    
	    plot(Time_step/24,byproduct_sum,'k-')
    
	    hold on
	    plot(Time_step/24,resource_sum,'m-')
	end

xlabel('Time (days)');
ylabel('Concentration (mg C/ml)');
grid on


iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'sum-byproducts-resources.pdf');
% saveas(figK,FigName,'pdf');
print(figK,FigName, '-dpdf', '-painters');

%% Growth rates plots 

close all

figK = figure;

% tiledlayout(3, 7)
tiledlayout(6, 7)

for i = 1:S

    nexttile(i)

	for k = 1:nb_rep_sim
        plot(Time_step, exp_biomass(i,:,k),'m-');
	    hold on
	end
	hold on

    for kk = 1:nb_rep %number of replicate observations
        plot(Time_step, Measured_Abund(i, :, kk),'r--');
    end

    set(gca,'YScale','log');
    ylim([1e-8 5e-4]);
    yticklabels({}); %turn off the ytick labels, but they go from 1e-8 to 1e-4
    grid on
    bact_name=name(i);
    title(bact_name);
end

%%calculate mu growth rates

%%Post processing analysis with moving interval on LN-transformed area data to calculate mu-Monod

mu_sim = zeros(nb_rep_sim, S);

for i = 1:nb_rep_sim
    A = exp_biomass(:,:,i)';%exp_biomass{i}';
    
    %%extract mu and rsq with a gliding scale - here 4 points  
    t = Time_step; %time scale in hours for the length of the table 
    int = 4; %number of points for the floating interval
    
    for k = 1:S %size(A,2)
	    for z = 1:length(t) - int
	        [mu_max(z), slope(z), intercept(z)] = regression(t(z:z+int), log(A(z:z+int, k))'); 
	    end
        max_m = max(mu_max);
        mu_sim(i,k) = max_m;
    end
end 


mu_exp = zeros(nb_rep, S);
for i = 1:nb_rep 

    % A = table2array(Data_Evol(:,(obs_X)-1+i))';
    A = Measured_Abund(:, :, i)';
    
    %%extract mu and rsq with a gliding scale - here 4 points
    t = Time_step; %Time_step2; %time scale in hours for the length of the table
    
    int = 3; %number of points for the floating interval
    figure(2)
    for k = 1:S %size(A,2)
	    for z = 1:length(t)-int
            temp_data = A(z:z+int,k);
            temp_data(temp_data <= 0) = 1e-07;
	        [mu_max, slope, intercept] = regression(t(z:z+int), log(temp_data)');%regression(t(z:z+int), log(A(z:z+int,k))'); 
            plot(t(z:z+int), log(temp_data), '-o')
            hold on
	    end
        max_m = max(mu_max);
        mu_exp(i,k) = max_m;
    end
end 


% %% Growth rates computation 
% 
% addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/Shared spatial/SurfaceModels/Figures generation')
% 
% 
% %For observations
% for i = 1:nb_rep 
% 
%     A = Measured_Abund(:, :, i)';
%     [time_init, mu_log] = deal(zeros(S, 1));
%     for k = 1:S
% 
%         Evol_mass = A(:,k);
%         lag_time = 0;
% 
%         x0 = [1, max(Evol_mass), Evol_mass(1), lag_time];
%         opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-10);
%         fun_Logistic = @(x, Time_step) x(2)./(1 + (x(2) - x(3))/x(3).*exp(-x(1)*Time_step));
% 
%         [pars, ~] = fminsearch(@(x) -logLikelihood(Evol_mass, 1e02, Time_step, 0, x(1), x(4), x(2), opts_1),...
%                 x0,...%Initial guess
%                 optimset('MaxFunEvals', 20000));
% 
%         sol = ode45(@(t,x) dfun(t, x, Evol_mass(1), 0, pars(1), pars(4), pars(2)), [0 max(Time_step)], Evol_mass(1), opts_1);
%         x_est = deval(sol, Time_step);
% 
%         time_init(k) = pars(4);
%         mu_log(k) = pars(1);
% 
%     end     
%        pos.time_init{i} = time_init; %fitted lag time (h)
%        pos.mu_log{i} = mu_log; %fitted logaritmic growth rate - to compare to mu_monod; multiply by log(2)
% end
% 
% 
% %For simulations
% [mu_max_Monod_LT, LT_Monod] = deal(zeros(nb_rep_sim, S));
% for k = 1:nb_rep_sim %number of replicate simulations
% 
%     data = exp_biomass(:, :, k); 
%     for i = 1:S
% 
%         data_Evol = data(:,1:end);
%         Evol_mass = data_Evol(i, :);
%         Evol_mass_R = sum(R_temp); %mean_R_0*ones(1, length(Time_step));
%         gamma = std(data_Evol); %Standard deviation of the observations. 
%         opt_1 = odeset('RelTol',1e-9,'AbsTol',1e-10); %Setting for ODEs solver
% 
%         %Implicit Monod estimation with lag time
%         r_1 = 1 ;
%         Ks = max(Evol_mass)/2; %Half velocity constant converted in equivalent OD value 0.0581 for OD.
%         x0 = [0.05, 0.1, 0.5, Ks]; %[mu_max, rho, t_lag, Ks]
%         pars_2 = fminsearch(@(x) -logLikelihood_2(Evol_mass, gamma, Time_step, Evol_mass_R, x(1), x(2), x(3),  x(4), opt_1),...
%             x0,...%Initial guess
%             optimset('MaxFunEvals', 20000));
%         sol = ode45(@(t,x) dfun_2(t, x, pars_2(2), pars_2(1), pars_2(3),  pars_2(4)), [0 max(Time_step)], [Evol_mass(1) Evol_mass_R(1)], opt_1);
%         x_est_2 = deval(sol, Time_step);
%         x_est_2 = x_est_2(1,:);
% 
%         mu_max_Monod_LT(k,i) = pars_2(1); %to have per hour
%         LT_Monod(k,i) = pars_2(3); %to have hour    
%     end
% end
% 
% %%remove any mu-prediction above 1
% mu_max_Monod_LT(mu_max_Monod_LT > 1)=nan;
% [mu_max_Monod_exp, LT_Monod_exp] = deal(zeros(nb_rep, S));
% for k = 1:nb_rep %number of experimental repetitions
% 
%     % data=table2array(Data_Evol(:,(obs_X)-1+k))';
%     data = Measured_Abund(:, :, k)';
% 
%     for i = 1:S %the strains
% 
%         data_Evol = data(:, 1:end); %Needed?
%         Evol_mass = data_Evol(:,i);
%         Evol_mass_R = sum(R_temp);%mean_R_0*ones(1, length(Time_step2)); %We d
%         gamma = std(data_Evol,0,2); %Standard deviation of the observations. 
%         opt_1 = odeset('RelTol',1e-9,'AbsTol',1e-10); %Setting for ODEs solver
% 
%         %Implicit Monod estimation with lag time
% 
%         Ks = max(Evol_mass)/2;% %Half velocity constant converted 
% 
%         x0 = [0.05, 0.1, 1, Ks]; %[mu_max, rho, t_lag, Ks]
%         pars_2 = fminsearch(@(x) -logLikelihood_2(Evol_mass, gamma, Time_step, Evol_mass_R, x(1), x(2), x(3),  x(4), opt_1),...
%             x0,...%Initial guess
%             optimset('MaxFunEvals', 20000));
%         sol = ode45(@(t,x) dfun_2(t, x, pars_2(2), pars_2(1), pars_2(3),  pars_2(4)), [0 max(Time_step)], [Evol_mass(1) Evol_mass_R(1)], opt_1);
%         x_est_2 = deval(sol, Time_step);
%         x_est_2 = x_est_2(1,:);
% 
%         mu_max_Monod_exp(k, i) = pars_2(1); %to have per hour
%         LT_Monod_exp(k, i) = pars_2(3); %to have hour  
%     end
% end
% 
% %%remove any mu-prediction above 1
% mu_max_Monod_exp(mu_max_Monod_exp > 1) = nan;