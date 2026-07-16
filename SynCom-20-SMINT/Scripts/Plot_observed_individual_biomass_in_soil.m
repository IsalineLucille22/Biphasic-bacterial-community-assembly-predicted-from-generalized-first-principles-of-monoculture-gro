%observed individual growth in soil of Syncom21 members. Script written by
%Jan Roelof van der Meer

clear
close all

scriptFolder = fileparts(mfilename('fullpath'));
projectFolder = fileparts(scriptFolder);
dataFolder = strcat(fullfile(projectFolder, 'Data'), '/');
FiguresFolder = strcat(fullfile(projectFolder, 'Figures'), '/');

%% plot individual biomass in soil

biomass = readtable(strcat(dataFolder, 'MergedData.xlsx'),'sheet','Growth kinetic params with rep','Range','U3:X24','VariableNamingRule','preserve');

strains = readtable(strcat(dataFolder, 'MergedData.xlsx'),'sheet','Growth kinetic params with rep','Range','C3:C24','VariableNamingRule','preserve');

combined = [strains, biomass];

combined = sortrows(combined,'Name','descend');

xvalues = categorical(combined.Name);

X = reordercats(xvalues,combined.Name);

yvalues = mean(table2array(combined(:,2:5)),2);

close all

figH = figure;

subplot('Position',[0.3 0.1 0.2 0.62]);
barh(X,log10(yvalues))
xlim([6 9]);
hold on
for i = 2:5
    scatter(log10(table2array(combined(:, i))), 1:21, 50, '.', 'k')
end
grid on

iFolderName = FiguresFolder;
FigName = strcat(iFolderName, 'observed-individual-biomass-in-soil.pdf');
print(figH, FigName, '-dpdf', '-painters');


% compare total biomass of the individual monocultures in soil and in community

%sum of all monocultures

sum_yvalues = sum(table2array(combined(:,2:5)));

%SynCom21 values

Senka = readtable(strcat(dataFolder, 'S20_S21_abs_abund_cfu_Senka.xlsx'),'sheet','absolute abundances','VariableNamingRule','preserve');

%filter for the syncom 21

SynCom21 = Senka(matches(Senka.condition,'S21'),:);

%filter for the 7d sample

SynCom21_7d = SynCom21(matches(SynCom21.("time.d"),'7d'),:);


%take sum of day7 values, multiply by 10 to have cells per g

day7_sum = 10*sum(table2array(SynCom21_7d(:,7:27)),2);

Xnames = {'community','monocultures'};
Xvalues = categorical(Xnames);


figH = figure;

subplot('Position',[0.1 0.3 0.3 0.4]);

bar(Xvalues(1),mean(day7_sum))
hold on
for i = 1:4
    scatter(1,day7_sum(i,:),200,'.','k')
end
grid on

bar(Xvalues(2),mean(sum_yvalues))
for i=1:4
    scatter(2,sum_yvalues(:,i),200,'.','k')
end

[h, p] = ttest(day7_sum,sum_yvalues');


[p1, h1]=ranksum(day7_sum,sum_yvalues');

iFolderName = FiguresFolder;
FigName = strcat(iFolderName, 'Fig2D-total-size-monocultures-syncom21-v1.pdf');
print(figH, FigName, '-dpdf', '-painters');

%% plot individual biomass in soil

biomass = readtable(strcat(dataFolder, 'MergedData.xlsx'),'sheet','Growth kinetic params with rep','Range','U3:X24','VariableNamingRule','preserve');

strains = readtable(strcat(dataFolder, 'MergedData.xlsx'),'sheet','Growth kinetic params with rep','Range','C3:C24','VariableNamingRule','preserve');

%SynCom21 values

Senka = readtable(strcat(dataFolder, 'S20_S21_abs_abund_cfu_Senka.xlsx'),'sheet','absolute abundances','VariableNamingRule','preserve');

%filter for the syncom 21

SynCom21=Senka(matches(Senka.condition,'S21'),:);

%filter for the 7d sample

SynCom21_7d=SynCom21(matches(SynCom21.("time.d"),'7d'),:);
[nb_rep, nb_info_species] = size(SynCom21_7d);


%correct names

SynCom21_7d = renamevars(SynCom21_7d,{'Pseudomonas1_R2','Pseudomonas2_R2'},{'Pseudomonas1','Pseudomonas2'});

SynCom21_biomass = [];

for i = 1:height(strains)
    tmp = SynCom21_7d(:,matches(SynCom21_7d.Properties.VariableNames,table2cell(strains(i,:))));
    
    SynCom21_biomass(i,:)=table2array(tmp);
end

%multiply by 10 to have cells per g

SynCom21_biomass = SynCom21_biomass.*10;

%combine the data

combined = [strains, biomass, array2table(SynCom21_biomass)];

combined = sortrows(combined,'Name','descend');

%sum of all monocultures

sum_yvalues = sum(table2array(combined(:,2:5)));

%take sum of SynCom21 day7 values

day7_sum = sum(table2array(combined(:,6:9)));


Xnames={'community','monocultures'};
Xvalues=categorical(Xnames);

%Color code associated to each species
syncom21_colors = {'#B35806', '#E08214', '#D53E4F', '#B2182B', '#D6604B', '#C51B7D', ...
               '#DE77AE', '#F1B6DA', '#FDAE61', '#FEE090', '#A6D96A', '#5AAE61', ...
               '#01665E', '#35978F', '#1B7837', '#C2A5CF', '#9970AB', '#762A83', ...
               '#80CDC1', '#C7EAE5', '#2166AC'};

syncom21_colors = flip(syncom21_colors);

close all

figH = figure;

ax1 = subplot('Position',[0.1 0.3 0.3 0.4], 'Parent', figH);
cats = Xvalues;                
xNum = 1:numel(cats);           

SynCom_val = [combined.SynCom21_biomass1, combined.SynCom21_biomass2, ...
          combined.SynCom21_biomass3, combined.SynCom21_biomass4];

Sum_Isolation = [combined.Rep1, combined.Rep2, combined.Rep3, combined.Rep4];

SynMean = mean(SynCom_val,2);
RepMean = mean(Sum_Isolation,2);

bar(xNum(1), SynMean, 'stacked'); 
hold on
bar(xNum(2), RepMean, 'stacked');
grid on

markerSize = 80;
jitter = 0.08; 

for i = 1:nb_rep
    y = sum(SynCom_val(:, i));
    x = xNum(1) + (rand(1) - 0.5)*jitter; 
    scatter(x, y, markerSize, 'k', 'filled')
end

for i = 1:nb_rep
    y = sum(Sum_Isolation(:, i));
    x = xNum(2) + (rand(1) - 0.5)*jitter;
    scatter(x, y, markerSize, 'k', 'filled')
end

colororder(syncom21_colors)
%legend

ax2 = subplot('Position',[0.6 0.3 0.2 0.4], 'Parent', figH);
hold(ax2,'on'); grid(ax2,'on')

bar(Xvalues(1), mean(day7_sum))
hold on
for i=1:nb_rep
    scatter(1, day7_sum(:,i),200,'.','k')
end
grid on

bar(Xvalues(2), mean(sum_yvalues))
for i=1:nb_rep
    scatter(2, sum_yvalues(:,i),200,'.','k')
end

delta = cliffsDelta(day7_sum, sum_yvalues);
p_value = Permutation_Test(day7_sum, sum_yvalues, delta);

iFolderName = FiguresFolder;
FigName = strcat(iFolderName, 'Fig2D-total-size-monocultures-syncom21-v2.pdf');
print(figH, FigName, '-dpdf', '-painters');