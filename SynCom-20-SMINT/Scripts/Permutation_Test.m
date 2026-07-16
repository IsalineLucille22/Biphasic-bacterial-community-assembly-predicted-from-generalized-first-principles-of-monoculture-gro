function p_value = Permutation_Test(Community_1, Community_2, delta_obs)
    combined = [Community_1; Community_2];
    len_1 = length(Community_1);
    len_2 = length(Community_2);
    n_obs = len_1 + len_2;
    nb_perm = 100000;
    delta = zeros(1,nb_perm);
    for i = 1:nb_perm
        temp_perm = combined(randperm(n_obs));
        group_A = temp_perm(1:len_1);
        group_B = temp_perm(len_1 + 1: end);
        delta(i) = cliffsDelta(group_A, group_B);
    end
    p_value = (sum(abs(delta) >= abs(delta_obs)) + 1) / (nb_perm + 1);
end