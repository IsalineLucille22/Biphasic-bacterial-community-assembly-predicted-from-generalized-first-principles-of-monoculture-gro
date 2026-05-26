%Simulations with fitted interspecific interactions for the transferred
%data
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

name_Exp = 'Phil_SynCom21_21wk_Paper_Transfer_Predation_Sim_y0';%'Senka_Comp_Ind_SynCom21_21wk_woP_Paper';
Name_file = 'data_to_save_SynCom21_Soil_only_Lyso_Parfor_Predation';%'struct_tot_SynCom21_Soil_only_Lyso_Parfor';
Name_file_Resources_Death = 'data_to_save_SynCom21_Reduced'; %Data fitted on moncultures

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','118:139', 'Format','auto');

%Initial abundances transfers 
y_0_mean_T = load(strcat('Data/', name_Exp, 'mean_y_0_T.mat'));%load(strcat('Data/', 'Phil_SynCom21_21wk_Paper_Transfer_Predation_Set_Obs_y0mean_y_0_T.mat'));%load(strcat('Data/', name_Exp, 'mean_y_0_T.mat'));
y_0_mean_T = y_0_mean_T.mean_y_0_T;

%actual absolute abundance data, has quadriplicate values for each time point
Data_Evol_tot = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 13, 'Range','1:22', 'Format','auto'); %Phil's weekly data
Data_Evol = Data_Evol_tot(:, 1:3);
Data_Evol_tot = table2array(Data_Evol_tot(:, 2:end));
mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8));
S = height(Data_Evol);
Time_step = [0 168];%Phil's data%[0 1 3 7 10 21]*24;%Senka %[0 12 22 38 70 168 504];%Bruna %[0 24 72 168];%Subcommunities %[0 12 48 96 7*24 21*24];%Clara %Measured time step in hours
Time_step_BP = 0:10:500;
nb_res = 12;
name = string(table2array(Parameters_set(1:S,1)));

%Load parameters
load(strcat('Data/', Name_file));
data_to_save = data_to_save_SynCom21_Soil_only_Lyso_Parfor_Predation;%struct_tot_SynCom21_Soil_only_Lyso_Parfor;
load(strcat('Data/', Name_file_Resources_Death));
data_to_save_Res_Death = data_to_save_Mono;
out_struct = Weighted_Struct(data_to_save);
nb_data_set = size(data_to_save,2);
Data_Evol_temp = table2array(Data_Evol(:, 2:end));
nb_obs = length(Data_Evol_temp(1,:));
nb_time_step = length(Time_step);
nb_rep = nb_obs/nb_time_step;
nb_rep_sim = nb_rep; %size(data_to_save,2)%Number of replicates to be simulated

% Fitted parameters
StackPlotTot = zeros(S, length(Time_step));
StackPlot_Meas_added_errors = zeros(S, length(Time_step)); %We perturbate the observed values accordingling to the simulated variance between data
z = zeros(S, nb_rep_sim);
z_obs_added_errors = zeros(S, nb_rep_sim); %We perturbate the observed values accordingling to the simulated variance between data
exp_biomass = zeros(S, length(Time_step), nb_rep_sim);

Measured_Abund = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
for g = 1:nb_rep
    Measured_Abund(:,:,g) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (g - 1));
end    

mean_y_0 = mean(Measured_Abund(:,1,:), 3);
var_mat = var(Measured_Abund, 0, 3);
Measured_Abund_average = mean(Measured_Abund,3);
%Initial concentrations using a normal distribution

%Initialization of the colors
[sorted_name, idx] = sortrows(name);
colors = {'#B35806', '#E08214', '#D53E4F', '#B2182B', '#B6604B', '#C51B7D', ...
               '#DE77AE', '#F1B6DA', '#FDAE61', '#FEE090', '#A6D96A', '#5AAE61', ...
               '#91665E', '#35978F', '#1B7837', '#C2A5CF', '#9970AB', '#762A83', ...
               '#80CDC1', '#C7EAE5', '#2166AC', '#4393C3', '#B35806'};



%Setting for Matlab ODE solver
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-12);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
stacked_plot_diff_CF = zeros(S, nb_rep_sim);
mean_y_0_iter = zeros(S, nb_rep_sim);

%Initialization figures
nb_exp = 8;
figH = gobjects(nb_exp,1);

for j = 1:nb_exp
    figH(j) = figure('Name', sprintf('Exp %d', j));
    t = tiledlayout(3,7);
    for k = 1:21
        ax(k) = nexttile;
        hold(ax(k),'on');
        grid(ax(k),'on');
        set(ax(k),'YScale','log');
        % xlabel(ax(k),'Time (days)');
        % ylabel(ax(k),'Biomass');
        title(ax(k), sprintf('Panel %d', k));
    end
    title(t, sprintf('Experiment %d', j));  % global title
end
angle = zeros(S, nb_rep_sim, nb_exp);
mean_y_0_T = zeros(S, nb_exp);
for zz = 1:nb_rep_sim

    out_struct = Weighted_Struct(data_to_save); %What is that?

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
    R = 1.6*data_to_save(ind_sim).R; %0.4 subcommunities wo Lysobacter or Mucilaginibacter. 0.6 wo Pseudomonas1. Reduce 1.25 + init Senka. 1.5 just reduce %1.6 for Phil's data
    nb_Res = length(R); %Number of resources (group of resources)
    Resource_Matrix = data_to_save(ind_sim).Resource_Matrix; %data_to_save_Res_Death.Resource_Matrix; %data_to_save(ind_sim).Resource_Matrix; %Addition of a line for Lysobacter
   
   
    %Initialization of the model parameters fixed for all replicates 
    t_0 = 0; %Time 0
    n_exp = 1; %No transfer in Philip experiment 
    tspan = [0, max(Time_step)*24]; %Time interval in hours
    yield_Pred = 0.2;% 20% of yield for predation
    
    for i = 1:nb_exp
        Data_Evol_temp = Data_Evol_tot(:, 2*(i - 1) + 1:(2*i));
        Measured_Abund = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
        for g = 1:nb_rep
            Measured_Abund(:,:,g) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (g - 1));
        end    
        
        mean_y_0 = mean(Measured_Abund(:,1,:), 3);
        var_mat = var(Measured_Abund, 0, 3);
        Measured_Abund_average = mean(Measured_Abund,3);
        noise = normrnd(0, 0, S, 1);
        mean_y_0 = mean_y_0 + noise;
        mat_y_0 = mean_y_0;
        mat_y_0([6 21]) = 0*10*mat_y_0(21);

        % if i > 1
        %     mat_y_0 = y_0_mean_T(:, i); %mean_y_0_T;
        % end
        
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

        %%%No CFP
        sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat, 0.*CrossFeed_Mat, 0.*Mat_kappa_3, Resource_Matrix,...
            0.*Threshold_CF, Threshold_death, 0.*Threshold_Pred, Death_Mat_Temp, death_rate,...
            0.*Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
        z_temp_no_CF = deval(sol, Time_step);
        X_no_CF = z_temp_no_CF(mod(1:S*3,3) == 1,:); %Species'biomass
        P_no_CF = z_temp_no_CF(mod(1:S*3,3) == 2,:); %Complexes'biomass
        W_no_CF = z_temp_no_CF(mod(1:S*3,3) == 0,:); %Byproducts'biomass
        R_temp_no_CF = z_temp_no_CF(S*3 + 1: end,:); %Byproducts'biomass
        X_no_CF = X_no_CF + P_no_CF;

        z_temp_BP = deval(sol, Time_step_BP);
        X_BP = z_temp_BP(mod(1:S*3,3) == 1,:); %Species'biomass
        P_BP = z_temp_BP(mod(1:S*3,3) == 2,:); %Complexes'biomass
        W_BP = z_temp_BP(mod(1:S*3,3) == 0,:); %Byproducts'biomass
        R_temp_BP = z_temp_BP(S*3 + 1: end,:); %Byproducts'biomass
        
        dz_temp = zeros(size(z_temp_BP));   % same size
        
        for j = 1:nb_time_step
            dz_temp(:,j) = dz(Time_step(j), z_temp_BP(:,j));
            dX = dz_temp(mod(1:S*3,3) == 1,:); %Species'biomass
            dP = dz_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
            dW = dz_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
            dR_temp = dz_temp(S*3 + 1: end,:); %Byproducts'biomass
        end

        Measured_Abund_disp = Measured_Abund_average + normrnd(zeros(S, length(Time_step)), 0.1*sqrt(var_mat));
        Measured_Abund_disp(Measured_Abund_disp < 0) = 0;
        z_temp = z_temp(1:(end-nb_Res), end);
        z_temp = reshape(z_temp', 3, S);
        z_temp = z_temp';
        z(:,zz) = X(:, end);%Total biomass of all species after 1 week or at the end of the experiement
        z_obs_added_errors(:,zz) = Measured_Abund_disp(:, end);
        StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
        u = [7*ones(1, S); (log10(X(:, end)) - log10(X(:, 1)))'];
        u(isnan(u)) = 0;
        v = [7*ones(1, S); (log10(Measured_Abund(:, end)) - log10(Measured_Abund(:,1)))'];
        v(isnan(v)) = 0;
        for s = 1:S
            angle(s, zz, i) = acos(u(:, s)'*v(:, s)/(norm(u(:, s))*norm(v(:, s))));
        end
        mean_y_0_T(:, i) = mean_y_0_T(:, i) + X(:, end)/10; %In the new transferred data only 10% of each species is kept

        StackPlotTot = StackPlotTot + X./sum(X);
        StackPlot_Meas_added_errors = StackPlot_Meas_added_errors + Measured_Abund_disp./sum(Measured_Abund_disp);
    
        exp_biomass(:,:, zz) = X;
        stacked_plot_diff_CF(:, zz) = Measured_Abund_average(:, end) - X(:, end);
        
        figure(figH(i));
        Time_step_days = Time_step/24;
        

        ax = zeros(1, S);
        for ii = 1:S  
            ax(ii) = nexttile(ii);
            
            plot(Time_step_days, X(ii, :),'m--');
            hold on
            
            plot(Time_step_days, Measured_Abund(ii, :),'g--');
          
            set(ax(ii), 'YScale', 'log');
            ylim(ax(ii), [1e-9 1e-2]);  
            xlim([0 max(Time_step_days)]);
            grid on
            bact_name = name(ii);
            title(bact_name);
        end
        drawnow limitrate
    end
    %reinitialize obs and mat_y0
end

mean_y_0_T = mean_y_0_T/nb_rep_sim;
% save(strcat('Data/',name_Exp, 'mean_y_0_T'), 'mean_y_0_T')

iFolderName = fullfile(cd,'Figures2');
if ~exist(iFolderName,'dir')
    mkdir(iFolderName)
end

for j = 1:nb_exp
    FigName = fullfile(iFolderName, ...
        sprintf('%s_ind-growth-plots-model-data_CF-start_exp%d.pdf', ...
                name_Exp, j));

    saveas(figH(j), FigName);                     % quick save
    print(figH(j), FigName, '-dpdf', '-painters');% publication-quality
end