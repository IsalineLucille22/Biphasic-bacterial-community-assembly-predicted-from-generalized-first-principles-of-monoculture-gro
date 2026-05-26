clear
close all

A = load(strcat('Data/', 'data_to_save_SynCom21_Soil_only_Lyso_Parfor_Newv5'));
B = load(strcat('Data/', 'data_to_save_SynCom21_Soil_only_Lyso_Parfor_Newv6'));

% If you don't know the variable names inside the .mat:
nameA = fieldnames(A);  s1 = A.(nameA{1});
nameB = fieldnames(B);  s2 = B.(nameB{1});

% 1) Check same set of fields
f1 = fieldnames(s1);
f2 = fieldnames(s2);
assert(isequal(sort(f1), sort(f2)), 'The two structs do NOT have the same fields.');

% 2) Force same field order (this is often the real reason concatenation "doesn't work")
s2 = orderfields(s2, s1);

% 3) Concatenate observations
data_to_save_SynCom21_Soil_only_Lyso_Parfor_Predation = [s1(:); s2(:)];   % result is 60x1 struct (use s.' if you insist on 1x60)
data_to_save_SynCom21_Soil_only_Lyso_Parfor_Predation = data_to_save_SynCom21_Soil_only_Lyso_Parfor_Predation.';              % optional: make it 1x60


save('data_to_save_SynCom21_Soil_only_Lyso_Parfor_Predation.mat','data_to_save_SynCom21_Soil_only_Lyso_Parfor_Predation')

