function [Fobs, R2obs, p_exact] = permanova1_exact(D, group)
group = group(:);
N = size(D,1);
ug = unique(group);

if numel(ug) ~= 2
    error('Only works for 2 groups.');
end

D2 = D.^2;
J = eye(N) - ones(N)/N;
G = -0.5*J*D2*J;
G = (G + G')/2;  %numerical symmetry
SST = trace(G);
%observed statistic
[Fobs, R2obs] = permanova_stats(G, group, SST);

%enumerate all possible splits
nA = sum(group == ug(1));
combA = nchoosek(1:N, nA);

Fperm = zeros(size(combA,1),1);

for b = 1:size(combA,1)
    gp = ug(2)*ones(N,1);
    gp(combA(b,:)) = ug(1);
    Fperm(b) = permanova_stats(G, gp, SST);
end

p_exact = (1 + sum(Fperm >= Fobs))/(length(Fperm) + 1);
end