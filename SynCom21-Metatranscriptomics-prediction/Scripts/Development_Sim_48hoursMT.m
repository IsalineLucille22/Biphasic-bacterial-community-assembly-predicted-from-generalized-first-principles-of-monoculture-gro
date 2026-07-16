%Simulations with fitted interspecific interactions
%use here the data_to_save_no_CF with 4 replicate simulations
%input is S20_S21_abs_abund_cfu_Senka.xlsx, Sheet 4, with 8 replicates 
%starting abundances for Bra, Phe, Mes, Cau, Coh, Tar, Bur, Chi and Muc are corrected


clear
close all

scriptFolder = fileparts(mfilename('fullpath'));
projectFolder = fileparts(scriptFolder);
dataFolder = strcat(fullfile(projectFolder, 'Data'), '/');
FiguresFolder = strcat(fullfile(projectFolder, 'Figures'), '/');


name_Exp = 'Zenodo_Bruna_48_hours_Predation';
Name_file = 'struct_tot_SynCom21_Soil_only_Lyso_Parfor';%'data_to_save_SynCom21_Soil_only_Lyso_Parfor_Predation';
Name_file_Resources_Death = 'data_to_save_SynCom21_Reduced'; %Data fitted on moncultures

%Loadind data
Parameters_set = readtable(strcat(dataFolder,'MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat(dataFolder,'MergedData.xlsx'), 'Sheet', 9, 'Range','118:139', 'Format','auto');

%actual absolute abundance data, has quadriplicate values for each time point

Data_Evol = readtable(strcat(dataFolder,'SSC21_genera_relative-abundances.xlsx'), 'Sheet', 4, 'Range','1:22', 'Format','auto'); %Metatranscriptomics'data
mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8));
S = height(Data_Evol);
Time_step = [0 12 22 38 70 168 504];%Metatranscriptomics'data %Measured time step in hours
nb_res = 12;
name = string(table2array(Parameters_set(1:S,1)));
Time_step_sim = 0:1:48;

%Load parameters
load(strcat(dataFolder, Name_file));
data_to_save = struct_tot_SynCom21_Soil_only_Lyso_Parfor;%data_to_save_SynCom21_Soil_only_Lyso_Parfor_Predation;
out_struct = Weighted_Struct(data_to_save);


%Initialization of the colors
colors_ind = [4, 1, 10, 15, 12, 3, 8, 6, 20, 18, 19, 2, 5, 21, 13, 9, 14, 7, 16, 17, 11];
colors_init = {'#B35806', '#E08214', '#D53E4F', '#B2182B', '#B6604B', '#C51B7D', ...
               '#DE77AE', '#F1B6DA', '#FDAE61', '#FEE090', '#A6D96A', '#5AAE61', ...
               '#91665E', '#35978F', '#1B7837', '#C2A5CF', '#9970AB', '#762A83', ...
               '#80CDC1', '#C7EAE5', '#2166AC', '#4393C3', '#B35806'};%distinguishable_colors(S);
colors = {};
for i = 1:S
    colors{i} = colors_init{colors_ind(i)};
end

Data_Evol_temp = table2array(Data_Evol(:, 2:end));
std_y_0 = 0;
nb_obs = length(Data_Evol_temp(1,:));
nb_time_step = length(Time_step);
nb_rep = nb_obs/nb_time_step;
nb_sim = 2*nb_rep;
Measured_Abund = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
for i = 1:nb_rep
    Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
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

% Fitted parameters
StackPlotTot = zeros(S, length(Time_step_sim));
StackPlot_Meas_added_errors = zeros(S, length(Time_step)); %We perturbate the observed values accordingling to the simulated variance between data
z = zeros(S, nb_sim);
z_obs_added_errors = zeros(S, nb_sim); %We perturbate the observed values accordingling to the simulated variance between data
exp_biomass = zeros(S, length(Time_step_sim), nb_sim);
exp_resources = zeros(nb_res, length(Time_step_sim), nb_sim);
exp_byproducts = zeros(S, length(Time_step_sim), nb_sim);


figure(num_fig); 
bar(mean(Measured_Abund_average, 3)', 'stacked');
axis([0 (nb_time_step + 1) 0 2.5*(sum(Measured_Abund_average(:, end)))])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked observed')
num_fig = num_fig + 1;

%Setting for Matlab ODE solver
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-12);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
nb_data_set = size(data_to_save,2);
[Shann_sim, Simp_sim, Shann_obs, Simp_obs] = deal(zeros(nb_data_set, nb_time_step));
Props_3hours = zeros(S, 5);
for zz = 1:nb_sim 

    out_struct = Weighted_Struct(data_to_save); %What is that?

    ind_sim = randi(nb_data_set);%zz; 
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
    R = data_to_save(ind_sim).R; %0.4 subcommunities wo Lysobacter or Mucilaginibacter. 0.5 wo Pseudomonas1. Reduce 1.25 + init Senka. 1.5 just reduce
    nb_Res = length(R); %Number of resources (group of resources)
    Resource_Matrix = data_to_save(ind_sim).Resource_Matrix; %data_to_save_Res_Death.Resource_Matrix; %data_to_save(ind_sim).Resource_Matrix; %Addition of a line for Lysobacter
    % Resource_Matrix(21,:) = 0; %Assuming Pseudomonas2 can't consume the initial resources of soil 
    
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
    n_exp = 1; %No transfer in Philip experiment 
    tspan = [0, max(Time_step)*24]; %Time interval in hours
    yield_Pred = 0.2;% 20% of yield for predation
    
    nb_replicates = 1;
    num_fig = 3;
    rand_rep = randi(5, 1, nb_rep) - 1;
    rand_rep = rand_rep./sum(rand_rep);
    noise = normrnd(0, 1e-9, S, 1);
    mean_y_0 = sum(rand_rep.*squeeze(Measured_Abund(:, 1, :)), 2);
    noise(mean_y_0 == 0) = 0;
    mean_y_0 = mean_y_0 + noise;
    mean_y_0(mean_y_0 < 0) = 0;
    mean_y_0(16) = 0; %Metatranscriptomics'data, no Microbacterium
    for i = 1:nb_replicates
        %Initial concentrations using a normal distribution
        mat_y_0 = mean_y_0;
        
        mat_y_0 = [mat_y_0 zeros(S,1) zeros(S,1)];
    
        for k = 1:n_exp
            y_0 = sum(mat_y_0(:,1:2),2);
            mat_y_0 = reshape(mat_y_0', 1, []);
            mat_y_0 = [mat_y_0 R];
       
            sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat, CrossFeed_Mat, Mat_kappa_3, Resource_Matrix,...
                Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp, death_rate,...
                Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
            
            z_temp = deval(sol, Time_step_sim);
            X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
            P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
            W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
            R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass

            figure(num_fig)
            for j = 1:S
                plot(Time_step_sim, X(j,:), '-o', 'Color', colors{j});%, Time_step, R_temp, 'o');
                hold on
            end
            legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
            axis([0 max(Time_step_sim) 0 4e-04]) %axis([0 600 0 2e-03])
            num_fig = num_fig + 1;

            z_temp = z_temp(1:(end-nb_Res), end);
            z_temp = reshape(z_temp', 3, S);
            z_temp = z_temp';
            z(:,zz) = X(:, end);%X(:, 4);%Total biomass of all species after 1 week or at the end of the experiement
            StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
        end
        StackPlotTot = StackPlotTot + X./sum(X);
    
        exp_biomass(:,:, zz) = X;
        exp_resources(:,:, zz) = R_temp; 
        exp_byproducts(:,:, zz) = W;
        Props_3hours(:, zz) = X(:, 4)/sum(X(:, 4));
    end

end
mean_Shann_sim = mean(Shann_sim); mean_Shann_obs = mean(Shann_obs);
mean_Simp_sim = mean(Simp_sim); mean_Simp_obs = mean(Simp_obs);
 


var_data = var(exp_biomass, 0, 3);

StackPlotTot = StackPlotTot/zz;
StackPlot_Meas_added_errors = StackPlot_Meas_added_errors/zz;

StackPlotTot = StackPlotTot./nb_replicates;
z_fin_sim = mean(z, 2); %Absolute stationary abundances
z_fin_obs_added_errors = mean(z_obs_added_errors, 2);

h_vect = zeros(1,S);
p_vect = zeros(1,S);
for i = 1:S
    [h_vect(i),p_vect(i)] = ttest2(StackPlot_Meas_added_errors(i,:)', StackPlotTot(i,:)');
end

iFolderName = FiguresFolder;
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

save(fullfile(dataFolder, [name_Exp 'Stacked_abundances.mat']), '-struct', 'Struct_Temp');

%% plot figure with the observed or simulated biomass per species at the end

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
for i = 1:nb_sim %number of replicates
 scatter(1:S, z(:,i), 60, '.', 'k')
end

set(gca, 'YScale', 'log')
ylim([1e-7 3e-4])
grid on
ylabel('Biomass g/mL')
yticks([1e-6 1e-5 1e-4]);
title('simulated');

iFolderName = FiguresFolder;
FigName = strcat(iFolderName, name_Exp, 'expected-end-biomass-SynCom20-in-soil-data-CF.pdf');
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
for i = 1:nb_sim 
    y = median(Measured_Abund(:, end, :), 3);
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

iFolderName = FiguresFolder;
FigName = strcat(iFolderName, name_Exp, 'expected-end-biomass-SynCom21-in-soil-scatter.pdf');
print(figL,FigName, '-dpdf', '-painters');

%individual boxcharts
close all

figL = figure;
subplot('Position',[0.1 0.1 0.55 0.55]);
boxchart(z')
set(gca, 'YScale', 'log')
ylim([1e-7 3e-4])

close all

figK = figure;

tiledlayout(3,7)

ax = gobjects(S,1);
for i = 1:S
    ax(i) = nexttile(i);
    hold(ax(i),'on')

    Ysim = squeeze(exp_biomass(i,:,:));   

    sim_med = median(Ysim, 2, 'omitnan');
    sim_q25 = prctile(Ysim, 25, 2);
    sim_q75 = prctile(Ysim, 75, 2);

    x = Time_step_sim(:)/24;
    fill(ax(i), [x; flipud(x)], [sim_q25; flipud(sim_q75)], ...
         'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(ax(i), x, sim_med, 'm-', 'LineWidth', 1.5);

    ylim(ax(i), [0 1.5e-4]); %ylim(ax(i), [0 1.5e-4]); for Pseudomonas2 %ylim(ax(i), [0 2e-5]); for Rahnella
    xlim(ax(i), [0 max(x)]);
    grid(ax(i),'on')
    title(ax(i), name(i));
end

iFolderName = FiguresFolder;
FigName = strcat(iFolderName, name_Exp, 'ind-growth-plots-model-data_CF-start.pdf');
saveas(figK,FigName,'pdf');
print(figK, FigName, '-dpdf', '-painters');