function [F, R2, p] = permanova1(D, group, nperm)
% Function created by Anderson 2001
% PERMANOVA one-factor
% D: N x N distance matrix
% group: N x 1 group labels (e.g., 1/2)
% nperm: number of permutations

N = size(D,1);
if size(D,2) ~= N
    error('D must be square');
end
group = group(:);
ug = unique(group);
k = numel(ug);

% Gower-centered matrix from distances
D2 = D.^2;
J = eye(N) - (1/N)*ones(N);
G = -0.5*J*D2*J;

% Total sum of squares = trace(G)
SST = trace(G);

% Compute fitted (between-group) sum of squares via hat matrix
% Design matrix with intercept + (k-1) dummies
X = ones(N,1);
for ii = 2:k
    X = [X, double(group == ug(ii))];
end
H = X*((X.'*X)\X.');  % hat matrix
SSB = trace(H * G);     % between
SSW = SST - SSB;        % within

dfb = k - 1;
dfw = N - k;

F = (SSB/dfb)/(SSW/dfw);
R2 = SSB/SST;

% Permutations: shuffle group labels
Fperm = zeros(nperm,1);
for b = 1:nperm
    gp = group(randperm(N));
    % rebuild Xp
    Xp = ones(N,1);
    for ii = 2:k
        Xp = [Xp, double(gp == ug(ii))];
    end
    Hp = Xp*((Xp.'*Xp)\Xp.');
    SSBp = trace(Hp * G);
    SSWp = SST - SSBp;
    Fperm(b) = (SSBp/dfb)/(SSWp/dfw);
end

% p-value with +1 correction
p = (1 + sum(Fperm >= F))/(nperm + 1);
end