close all
clear 

load(strcat('Data/', 'Subcommunity_wo_LysoStacked_abundances'));
load(strcat('Data/', 'Subcommunity_wo_MuciStacked_abundances'));
load(strcat('Data/', 'Subcommunity_wo_Pseudo1and2Stacked_abundances'));
load(strcat('Data/', 'SynCom21_SenkaStacked_abundances'));
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
[S, m] = size(Parameters_set);

obs_1234 = [Subcommunity_wo_Pseudo1and2StackPlot_Meas(:, end), Subcommunity_wo_LysoStackPlot_Meas(:,end),...
    Subcommunity_wo_MuciStackPlot_Meas(:, end), SynCom21_SenkaStackPlot_Meas(:, end)];
sim_1234 = [Subcommunity_wo_Pseudo1and2StackPlot_Tot(:, end), Subcommunity_wo_LysoStackPlot_Tot(:,end),...
    Subcommunity_wo_MuciStackPlot_Tot(:, end), SynCom21_SenkaStackPlot_Tot(:, end)];

name = string(table2array(Parameters_set(1:S,1)));

figK = figure;
subplot(1, 2, 1)
bar(obs_1234', 'stacked')
axis([0, 5, 0, 1])
xticklabels({'C1 Obs', 'C2 Obs', 'C3 Obs','C4 Obs'})
subplot(1, 2, 2)
bar(sim_1234', 'stacked')
xticklabels({'C1 Sim', 'C2 Sim', 'C3 Sim','C4 Sim'})
axis([0, 5, 0, 1])
%legend(name, 'Orientation', 'vertical', 'Location', 'southeast')

iFolderName = strcat(cd, '/Figures2/');
FigName = strcat(iFolderName, 'CompFinalSub.pdf');
saveas(figK,FigName,'pdf');
print(figK, FigName, '-dpdf', '-painters');