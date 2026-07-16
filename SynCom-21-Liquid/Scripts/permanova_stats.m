function [F, R2] = permanova_stats(G, group, SST)
N = length(group);
X = [ones(N,1), double(group==2)];
H = X * ((X.'*X)\X.');

SSB = trace(H*G);
SSW = SST - SSB;

dfb = 1;
dfw = N - 2;

F = (SSB/dfb) / (SSW/dfw);
R2 = SSB / SST;
end