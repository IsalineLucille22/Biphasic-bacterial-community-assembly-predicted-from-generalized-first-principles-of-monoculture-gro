function p_val = Permutation_Test_Comm(Dist, group, nperm)
% observed mean between-group distance
idxA = find(group == 1);
idxB = find(group == 2);
obs = mean(Dist(idxA, idxB), 'all');
stat = zeros(nperm, 1);
N = numel(group);
for b = 1:nperm
    gp = group(randperm(N));
    iA = find(gp==1); 
    iB = find(gp==2);
    stat(b) = mean(Dist(iA, iB), 'all');
end
p_val = (1 + sum(stat >= obs))/(nperm + 1);
end