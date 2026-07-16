function [PCoA_Metric, cos_theta, p_val]  = TrajectoryDist(Y_PCoA, rep, time, env)
unique_time = length(unique(time));
unique_env = unique(env);
ind_A = env == unique_env(1);
ind_B = env == unique_env(2);
env_A = Y_PCoA(ind_A, :); 
env_B = Y_PCoA(ind_B, :);
nb_env = length(unique(env));
rep_A = rep(ind_A); rep_B = rep(ind_B);
nb_rep_A = numel(unique(rep_A)); nb_rep_B = numel(unique(rep_B));
time_A = time(ind_A); time_B = time(ind_B);
temp_tot = 0; disp_A = 0; disp_B = 0;
cos_theta = 0;
[~, m] = size(Y_PCoA);
mu_A = zeros(unique_time, m);
mu_B = zeros(unique_time, m);
for i = 2:unique_time
    temp_A = env_A(time_A == i, :); temp_B = env_B(time_B == i, :);
    temp_tot = temp_tot + norm(mean(temp_A, 1) - mean(temp_B, 1));
    disp_A = disp_A  + sum(vecnorm(temp_A - mean(temp_A), 2, 2))/nb_rep_A;
    disp_B = disp_B  + sum(vecnorm(temp_B - mean(temp_B), 2, 2))/nb_rep_B;
    mu_A(i - 1, :) = mean(env_A(time_A == i - 1, :), 1); mu_B(i - 1, :) = mean(env_B(time_B == i - 1, :), 1);
    delta_A_temp = mean(env_A(time_A == i, :), 1) - mean(env_A(time_A == i - 1, :), 1);
    delta_B_temp = mean(env_B(time_B == i, :), 1) - mean(env_B(time_B == i - 1, :), 1);
    cos_theta = cos_theta + delta_A_temp*delta_B_temp'/(norm(delta_A_temp)*norm(delta_B_temp));
end
mu_A(i, :) = mean(env_A(time_A == i, :), 1); mu_B(i, :) = mean(env_B(time_B == i, :), 1);
% delta_A = mean(env_A(time_A == i, :), 1) - mean(env_A(time_A == 1, :), 1);
% delta_B = mean(env_B(time_B == i, :), 1) - mean(env_B(time_B == 1, :), 1);
cos_theta = cos_theta/(unique_time - 1);%delta_A*delta_B'/(norm(delta_A)*norm(delta_B));%close to 1 similar dynamics, 0 unrelated, <0, opposit direction
p_val = Permutt_cosin(unique_time, mu_A, mu_B, cos_theta);
WithIn_tot = (disp_A + disp_B)/(nb_env*(unique_time - 1));
PCoA_Metric = (1/unique_time)*temp_tot/WithIn_tot; %<1 env similar, =1 envi differs as much as rep var, >1 env has a large effect
end