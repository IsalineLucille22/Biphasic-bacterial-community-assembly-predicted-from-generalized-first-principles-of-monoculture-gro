%Simulations with fitted interspecific interactions
%use here the data_to_save_no_CF with 4 replicate simulations
%input is S20_S21_abs_abund_cfu_Senka.xlsx, Sheet 4, with 8 replicates 
%starting abundances for Bra, Phe, Mes, Cau, Coh, Tar, Bur, Chi and Muc are corrected


clear
close all

pwd_init = pwd;
addpath(strcat(pwd_init, '/Data'))
addpath(pwd_init)
addpath(strcat(pwd_init, '/21 SynCom Script (with Lysobacter)'))

%Save or Not
save_data = 0; %1 if save, 0 otherwise
name_Exp = 'Zenodo_SynCom20_Paper_Random_interaction';
Name_file = 'data_to_save_Inter_v2';
Name_file_Saved = 'Fixed_Res_12_Groups_v2';

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','25:46','Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','71:91', 'Format','auto');

%actual absolute abundance data, has quadriplicate values for each time point

%Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'),
%'Sheet', 4, 'Range','1:21', 'Format','auto'); %Combination data fitting
%and test
Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 3, 'Range','1:21', 'Format','auto'); %Test with the new SynCom 20
mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8));
S = height(Data_Evol);

%Load parameters
load(strcat('Data/', Name_file));
data_to_save = data_to_save_Inter_v2;
Time_step = [0 1 3 7 10 21]*24; %Time in hours; %[0 1 3 7 10 15 21 22]*24; %Time in hours
Time_step_HTP = 0:0.5:21*24; %Time in hours; %More simulated time points than in observations
nb_res = 12;
Data_Evol_temp = table2array(Data_Evol(:, 2:end));
Data_Evol_temp = Data_Evol_temp(:, ismember(mod(1:length(Data_Evol_temp(1,:)),4), [1, 2, 3]));%20% of the data to train. Hold-out method.
nb_obs = length(Data_Evol_temp(1,:));
nb_time_step = length(Time_step);
nb_rep = nb_obs/nb_time_step;
nb_rep_sim = nb_rep;

% Fitted parameters
StackPlotTot = zeros(S, length(Time_step));
StackPlot_Meas_added_errors = zeros(S, length(Time_step)); %We perturbate the observed values accordingling to the simulated variance between data
z = zeros(S, nb_rep_sim);
z_no_CF = zeros(S, nb_rep_sim);
z_obs = zeros(S, nb_rep);
% z_obs_added_errors = zeros(S, size(data_to_save,2)); %We perturbate the observed values accordingling to the simulated variance between data
nb_data_set = size(data_to_save,2); 
exp_byproducts = zeros(S, length(Time_step), nb_rep_sim);
exp_resources = zeros(nb_res, length(Time_step), nb_rep_sim);
exp_biomass = zeros(S, length(Time_step), nb_rep_sim);
exp_biomass_wo_CF = zeros(S, length(Time_step), nb_rep_sim);
exp_byproducts_HTP = zeros(S, length(Time_step_HTP), nb_rep_sim);
exp_resources_HTP = zeros(nb_res, length(Time_step_HTP), nb_rep_sim);
name = string(table2array(Parameters_set(1:S, 1)));
[sorted_name, idx_sorted] = sort(name);

%Initialization of the colors
colors_ind = [4, 1, 10, 15, 12, 8, 6, 20, 18, 19, 2, 5, 21, 13, 9, 14, 7, 16, 17, 11];
colors_init = {'#B35806', '#E08214', '#D53E4F', '#B2182B', '#B6604B', '#C51B7D', ...
               '#DE77AE', '#F1B6DA', '#FDAE61', '#FEE090', '#A6D96A', '#5AAE61', ...
               '#91665E', '#35978F', '#1B7837', '#C2A5CF', '#9970AB', '#762A83', ...
               '#80CDC1', '#C7EAE5', '#2166AC', '#4393C3', '#B35806'};
colors = {};
for i = 1:S
    colors{i} = colors_init{colors_ind(i)};
end



Measured_Abund = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
for i = 1:nb_rep
    Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
    z_obs(:, i) = Measured_Abund(:, end, i);
end
StackPlot_Meas = Measured_Abund./sum(Measured_Abund);
StackPlot_Meas = mean(StackPlot_Meas, 3);
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


% Measured_Abund = table2array(Data_Evol(1:20, 2:7));
%Number of surviving species after 8 weeks
Threshold_Surviving = 1e-10;
nb_Surv_Obs = sum(Measured_Abund_average > Threshold_Surviving);

%Initialization of the model parameters fixed for all replicates 
t_0 = 0; %Time 0
tspan = [0, max(Time_step)*24]; %Time interval in hours

%Setting for Matlab ODE solver
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-12);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
[Shann_sim, Simp_sim, Shann_obs, Simp_obs] = deal(zeros(nb_data_set, nb_time_step));

% %To generate random interactions
CrossFeed_Mat = lognrnd((log(0.8) + mu_max_dist(:,1)).*ones(S,S), mu_max_dist(:,2).*ones(S,S));%.val val% of zeros. ConsumerxProducer
CrossFeed_Mat_Temp = zeros(S,S); %Temporal predation matrix
rand_indice = rand(S,S) > 0.5; %Percentage of zeros. Larger is the value smaller is the number of interaction
CrossFeed_Mat_Temp(rand_indice) = CrossFeed_Mat(rand_indice); %Put some element to 0
Death_Mat = zeros(S,S);
CrossFeed_Mat_Temp(1:1+size(CrossFeed_Mat_Temp,1):end) = 0;
CrossFeed_Mat = CrossFeed_Mat_Temp;
Threshold_CF = max(normrnd(1.0e-03, 0*0.1e-05, S, 1),0); 


stacked_plot_diff_CF = zeros(S, nb_rep_sim);
stacked_plot_diff_no_CF = zeros(S, nb_rep_sim);
Threshold_prop = 0.05;
for zz = 1:nb_rep_sim %nb_iter

    ind_sim = randi(nb_data_set);
    kappa_mat = data_to_save(ind_sim).kappa_mat;
    Resource_Matrix = data_to_save(ind_sim).Resource_Matrix;
    Threshold_death = data_to_save(ind_sim).Threshold_death;
    Lag_time_Cons = data_to_save(ind_sim).Lag_time_Cons; 
    Lag_time_Pred = data_to_save(ind_sim).Lag_time_Pred;
    R = data_to_save(ind_sim).R;
    nb_Res = length(R); %Number of resources (group of resources)
    Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
    death_rate =  data_to_save(ind_sim).death_rate;
    Death_Mat(1:1+size(Death_Mat,1):end) = death_rate; %For Random
       
    
    num_fig = 2;
    rand_rep = randi(5, 1, nb_rep) - 1;
    if sum(rand_rep) == 0
        rand_rep(1) = 1; %To avoid having NaN value. If zero, then consider the 1st replicate as only obs.
    end
    rand_rep = rand_rep./sum(rand_rep);
    noise = normrnd(0, 1e-9, S, 1);
    mean_y_0 = sum(rand_rep.*squeeze(Measured_Abund(:, 1, :)), 2);
    noise(mean_y_0 == 0) = 0;
    mean_y_0 = mean_y_0 + noise;
    mean_y_0(mean_y_0 < 0) = 0;
    %Initial concentrations using a normal distribution
    mat_y_0 = mean_y_0;
    
    mat_y_0 = [mat_y_0 zeros(S,1) zeros(S,1)];

    y_0 = sum(mat_y_0(:,1:2),2);
    mat_y_0 = reshape(mat_y_0', 1, []);
    mat_y_0 = [mat_y_0 R];

    %With CF
    sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat, CrossFeed_Mat, Mat_kappa_3, Resource_Matrix, Threshold_CF, Threshold_death, Death_Mat, death_rate, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1); %Multiple resource groups
    z_temp = deval(sol, Time_step);
    z_temp_HTP = deval(sol, Time_step_HTP);
    X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
    P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
    W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
    R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass
    W_HTP = z_temp_HTP(mod(1:S*3,3) == 0,:); %Byproducts'biomass
    R_temp_HTP = z_temp_HTP(S*3 + 1: end,:); %Byproducts'biomass
    % X = X + P;
    %Wo CF
    sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat, 0.*CrossFeed_Mat, Mat_kappa_3, Resource_Matrix, Threshold_CF, Threshold_death, Death_Mat, death_rate, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1); %Multiple resource groups
    z_temp_no_CF = deval(sol, Time_step);
    X_no_CF = z_temp_no_CF(mod(1:S*3,3) == 1,:); %Species'biomass
    P_no_CF = z_temp_no_CF(mod(1:S*3,3) == 2,:); %Complexes'biomass
    W_no_CF = z_temp_no_CF(mod(1:S*3,3) == 0,:); %Byproducts'biomass
    R_temp_no_CF = z_temp_no_CF(S*3 + 1: end,:); %Byproducts'biomass
    X_no_CF = X_no_CF + P_no_CF;
    figure(num_fig)
    for j = 1:S
        plot(Time_step, X(j,:), '-o', 'Color', colors{j});%, Time_step, R_temp, 'o');
        hold on
    end
    %legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
    axis([0 600 0 4e-04])
    num_fig = num_fig + 1;
    figure(num_fig)
    Measured_Abund_disp = Measured_Abund_average + normrnd(zeros(S, length(Time_step)), sqrt(var_mat));
    for j = 1:S
        plot(Time_step, Measured_Abund_disp(j,:), '-*', 'Color', colors{j});%, Time_step, R_temp, 'o');
        hold on
    end
    legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
    axis([0 600 0 4e-04])

    num_fig = num_fig + 1;
    z_temp = z_temp(1:(end-nb_Res), end);
    z_temp = reshape(z_temp',3, S);
    z_temp = z_temp';
    z(:,zz) = sum(z_temp(:,1:2), 2);%Total biomass of all species after 1 week
    z_temp_no_CF = z_temp_no_CF(1:(end-nb_Res), end);
    z_temp_no_CF = reshape(z_temp_no_CF',3, S);
    z_temp_no_CF = z_temp_no_CF';
    z_no_CF(:,zz) = sum(z_temp_no_CF(:,1:2), 2);%Total biomass of all species after 1 week
    StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
    [Shann_sim, Simp_sim] = Shannon_Simpson_Indices(S,StackPlot);
    [Shann_obs, Simp_obs] = Shannon_Simpson_Indices(S, mean(StackPlot_Meas,3));
    
    exp_byproducts(:,:, zz) = W;
    exp_byproducts_HTP(:,:, zz) = W_HTP;
    exp_biomass(:,:, zz) = X;
    exp_biomass_wo_CF(:,:, zz) = X_no_CF;
    exp_resources(:,:, zz) = R_temp; 
    exp_resources_HTP(:,:, zz) = R_temp_HTP; 
    StackPlotTot = StackPlotTot + X./sum(X);
    prop_obs = Measured_Abund_average(:, end)/sum(Measured_Abund_average(:, end));
    stacked_plot_diff_CF(:, zz) = Measured_Abund_average(:, end) - X(:, end);%100*abs(Measured_Abund_average(:, end) - X(:, end))./Measured_Abund_average(:, end);
    % stacked_plot_diff_CF(prop_obs < Threshold_prop, zz) = 0; %For propos
    % computations
    stacked_plot_diff_no_CF(:, zz) = Measured_Abund_average(:, end) - X_no_CF(:, end);%100*abs(Measured_Abund_average(:, end) - X_no_CF(:, end))./Measured_Abund_average(:, end);
end
%Consider only species > 0.05 of the total final biomass. 
tot_props = sum(prop_obs(prop_obs >= Threshold_prop));
mean_stacked_plot_CF = mean(stacked_plot_diff_CF, 2);
mean_stacked_plot_no_CF = mean(stacked_plot_diff_no_CF, 2);
% perc_diff_obs = mean_stacked_plot_CF./Measured_Abund_average(:, end);
% perc_diff_obs_no_CF = mean_stacked_plot_no_CF./Measured_Abund_average(:, end);
average_CF = sum(mean(stacked_plot_diff_CF,2))/sum(prop_obs >= Threshold_prop);
average_no_CF = sum(mean(stacked_plot_diff_no_CF,2))/sum(prop_obs >= Threshold_prop);

%Permanova and permutation-test on Bray-Curtis distance to get on statistic 

m_1 = nb_rep;
m_2 = nb_rep_sim;
mat_Comm_1 = z_obs;
mat_Comm_2 = z;
X = [mat_Comm_1.'; mat_Comm_2.'];        %(m_1 + m_2)xn   rows = replicates, cols = species
group = [ones(m_1, 1); 2*ones(m_2, 1)];  %labels: 1 = obs, 2 = sim
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
xlabel('PCo1')
ylabel('PCo2')
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

xlabel('PCoA1')
ylabel('PCoA2')
title('Community trajectories in PCoA space')
axis equal
legend(envs)
iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'SimCF_noCF_MDS_Plot.pdf');
print(figG, FigName, '-dpdf', '-painters');

figK = figure;
bar(mean_stacked_plot_CF(idx_sorted, :))
hold on
scatter(1:S, stacked_plot_diff_CF(idx_sorted, :),60,'.','r')
grid on
axis([0 21 -4e-04 1.5e-04]) %1.5e-04%1e02% For prop analysis
xticks(1:length(sorted_name)); 
xticklabels(sorted_name);
ylabel('measured - simulated')

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'stacked_plot_diff_CF.pdf');
print(figK,FigName, '-dpdf', '-painters');

figK = figure;
bar(mean_stacked_plot_no_CF(idx_sorted, :))
hold on
scatter(1:S, stacked_plot_diff_no_CF(idx_sorted, :),60,'.','r')
grid on
axis([0 21 -0.8e-04 1.5e-04])
xticks(1:length(sorted_name)); 
xticklabels(sorted_name);
ylabel('measured - simulated')

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'stacked_plot_diff_no_CF.pdf');
print(figK,FigName, '-dpdf', '-painters');

close all 

StackPlotTot = StackPlotTot./nb_rep_sim;
z_fin_sim = z(1:S,end); %Absolute stationary abundances
z_fin_obs = Measured_Abund(1:S, end); %Absolute stationary abundances
nb_Surv_Sim = sum(StackPlotTot > Threshold_Surviving);
disp(sum(z_fin_sim)/sum(z_fin_obs))

h_vect = zeros(1,S);
p_vect = zeros(1,S);
for i = 1:S
    [h_vect(i), p_vect(i)] = ttest2(StackPlot_Meas(i,:)', StackPlotTot(i,:)');
end
name_diff = name(logical(h_vect));

figure(num_fig);
bar(StackPlotTot', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked all replicates')
num_fig = num_fig + 1;

figure(num_fig); 
bar(mean(StackPlot_Meas,3)', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked observed')
num_fig = num_fig + 1;

[diff_vect, I, fact] = Sort_diff(StackPlotTot(:,end), StackPlot_Meas(:,end));
name_order = name(I);
stem(diff_vect);
xtickangle(90)
set(gca,'xtick',1:21,'xticklabel',name_order)
ylabel('Sim - Obs')
num_fig = num_fig + 1;

figure(num_fig);
for j = 1:S
    scatter(StackPlot_Meas(j, end), StackPlotTot(j,end), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors{j}, 'MarkerFaceColor', colors{S+1-j});%col(S+1-j,:))      
    hold on
end
axis([0 0.5 0 0.5]);
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
    scatter(z_fin_obs(j), z_fin_sim(j), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors{j}, 'MarkerFaceColor', colors{S+1-j});%col(S+1-j,:))      
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

X = X'; %For R sript comparison

%% Computation of the rates distribution for each species

% [~, nb_iter] = size(data_to_save_Inter_v2);
% Mat_Struct = {data_to_save_Inter_v2.CrossFeed_Mat_Temp};
% [Sample_rate, LN_fit, nb_neg_val] = Rates_Distribution(12, Mat_Struct, num_fig, name);

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
ylabel('Biomass mg/L')
title('observed');

subplot('Position',[0.3 0.5 0.62 0.2]);

bar(mean(z,2))
hold on
for i=1:nb_rep %number of replicates
    scatter(1:S,z(:,i),60,'.','k')
end

set(gca, 'YScale', 'log')
ylim([1e-7 3e-4])
grid on
ylabel('Biomass mg/L')
yticks([1e-6 1e-5 1e-4]);
title('simulated');

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'expected-end-biomass-SynCom20-in-soil-data-CF.pdf');
print(figK,FigName, '-dpdf', '-painters');


% sort of double scatter plot for simulated and experimental data
figL = figure;
subplot('Position',[0.1 0.1 0.55 0.55]);

for i = 1:nb_rep %number of replicates
    scatter(median(z,2),  Measured_Abund(:, end, i), 60,'.','k')
    hold on
end
hold on
for i = 1:nb_rep_sim %strains
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
FigName = strcat(iFolderName, name_Exp, 'expected-end-biomass-SynCom20-in-soil-scatter.pdf');
print(figL,FigName, '-dpdf', '-painters');

close all

figL = figure;
subplot('Position',[0.1 0.1 0.55 0.55]);
boxchart(reshape(Measured_Abund(:, end, :), S, nb_rep)')
set(gca, 'YScale', 'log')
print(figL, strcat(iFolderName, name_Exp, 'boxchart-observed.pdf'), '-dpdf', '-painters');

close all

figK = figure;

tiledlayout(4,5)

ax = gobjects(S,1);
Time_step_days = Time_step/24;
for i = 1:S  
    ax(i) = nexttile(i);
    
    for k = 1:nb_rep_sim
        plot(Time_step_days, exp_biomass(i, :, k),'g--');
        hold on
    end
    hold on

    for k = 1:nb_rep_sim
        plot(Time_step_days, exp_biomass_wo_CF(i, :, k),'LineStyle','--','Color',[0 0 1 0.1]);
        hold on
    end
    hold on
    
    for kk = 1:nb_rep %number of replicate observations
        plot(Time_step_days, Measured_Abund(i, :, kk),'r--');
    end
    
    set(ax(i), 'YScale', 'log');
    ylim(ax(i), [1e-9 1e-3]);  
    xlim([0 max(Time_step_days)]);
    grid on
    bact_name = name(i);
    title(bact_name);
end

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'ind-growth-plots-model-data_CF-start.pdf');
print(figK, FigName, '-dpdf', '-painters');

%% display mean and stdev of the cross-feeding matrices over time
tmp = zeros(S, S, nb_rep_sim);

for i = 1:nb_data_set
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

% multi stack bar plot of relative abundances in the simulations
% compared to the mean of 8 observations

[strains, idx]=sortrows(name); %to sort alphabetically
mean_CF=mean_CF(idx,:); %sort the byproducts rows accordingly
mean_CF=mean_CF(:,idx); %sort the columns accordingly

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
print(figK, FigName, '-dpdf', '-painters');

close all

figK = figure;

tiledlayout(3, 4)

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
print(figK, FigName, '-dpdf', '-painters');

%%simple plot of the maximum/final biomass attained in the simulations versus observations

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
    % obs_fin_biomass(i,:) = max(sum(Measured_Abund(:,end,i)));
end

close all
figH = figure;

subplot('Position',[0.1 0.1 0.2 0.4])

y = [mean(sim_fin_biomass); mean(sim_fin_biomass_wo_CF); mean(obs_fin_biomass)];

X = categorical({'Sim_CF', 'Sim_wo_CF', 'Obs'});
X = reordercats(X,{'Sim_CF', 'Sim_wo_CF', 'Obs'});
b = bar(X,y,'FaceColor','flat');
set(gca,'FontSize',7);
b.CData(1,:) = [1 0 1];
b.CData(2,:) = [0 0 1];
b.CData(3,:) = [1 0 0];

hold on
for j = 1:nb_rep_sim
    scatter(X(1), sim_fin_biomass(j),50,'k','.')
end
hold on
for j = 1:nb_rep
 scatter(X(2), sim_fin_biomass_wo_CF(j),50,'k','.')
end
hold on
for j = 1:nb_rep
 scatter(X(3), obs_fin_biomass(j),50,'k','.')
end
ylabel('biomass (gC/mL)')
ylim([0 2e-3])

title('max biomass','FontSize',9);

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'max-data-CF-sim-obs-biomass.pdf');
print(figH,FigName, '-dpdf', '-painters');
disp(mean(sim_fin_biomass)/mean(obs_fin_biomass))
disp(mean(sim_fin_biomass_wo_CF)/mean(obs_fin_biomass))

[Sortedname, idx] = sort(name);
z = z(idx, :);
save(strcat('Data/',name_Exp, 'z'), 'z')

%%%%Permutation tests
delta = cliffsDelta(sim_fin_biomass, obs_fin_biomass);
p_value = Permutation_Test(sim_fin_biomass, obs_fin_biomass, delta);

%No CF
delta_wo_CF = cliffsDelta(sim_fin_biomass_wo_CF, obs_fin_biomass);
p_value_wo_CF = Permutation_Test(sim_fin_biomass, obs_fin_biomass, delta_wo_CF);
title('Max biomass','FontSize',9);

%%%%Heatmap byproducts formation
tmp = exp_byproducts;

mean_byproducts = mean(tmp, 3);
mean_byproducts = mean_byproducts(idx,:); %sort the byproducts matrix accordingly

close all

figK=figure;

subplot('Position',[0.2 0.2 0.2 0.55]);
h = heatmap(Time_step/24, strains, mean_byproducts);

h.Title = 'Byproduct formation';
h.XLabel = 'Time (days)';

%% sum byproducts in the system plus the modeled resource utilization

figK = figure;

subplot('Position',[0.1 0.1 0.3 0.3]);

for k = 1:nb_rep_sim
    % byproduct_sum_t = sum(exp_byproducts_t{k},1);
    % resource_sum_t = sum(exp_resources_t{k},1);
    % byproduct_sum = sum(exp_byproducts(:, :, k),1);
    % resource_sum = sum(exp_resources(:, :, k),1);
    byproduct_sum_HTP = sum(exp_byproducts_HTP(:, :, k),1);
    resource_sum_HTP = sum(exp_resources_HTP(:, :, k),1);


    plot(Time_step_HTP/24, byproduct_sum_HTP,'k-')

    hold on
    plot(Time_step_HTP/24, resource_sum_HTP,'m-')
end

xlabel('Time (days)');
ylabel('Concentration (mg C/ml)');
grid on


iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'sum-byproducts-resources.pdf');
print(figK,FigName, '-dpdf', '-painters');