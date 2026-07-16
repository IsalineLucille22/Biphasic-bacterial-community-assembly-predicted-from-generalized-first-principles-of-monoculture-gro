function p_val = Permutt_cosin(unique_time, mu_A, mu_B, obs_cosin)
%H0: trajectories are no more aligned than random vs H1: trajectories are more aligned/similar than random
nperm = 5000;
perm_vals = zeros(nperm,1);
for p = 1:nperm
    perm_idx = randperm(unique_time);
    mu_B_perm = mu_B(perm_idx, :);
    cos_vals_perm = zeros(unique_time-1,1);

    for t = 1:(unique_time - 1)
        delta_A = mu_A(t+1, :) - mu_A(t, :);
        delta_B = mu_B_perm(t+1, :) - mu_B_perm(t, :);

        if norm(delta_A) > 1e-12 && norm(delta_B) > 1e-12
            cos_vals_perm(t) = (delta_A*delta_B')/(norm(delta_A)*norm(delta_B));
        else
            cos_vals_perm(t) = NaN;
        end
    end

    perm_vals(p) = mean(cos_vals_perm, 'omitnan');
end
p_val = (1 + sum(perm_vals >= obs_cosin)) / (nperm + 1);
end