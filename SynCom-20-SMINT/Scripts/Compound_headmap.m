%Script to create the resource heatmap and log-normal fitted distribution +
%histograms made by Jan

cd('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/SynCom model paper/Data/Data IG')


strains=readtable('MergedData.xlsx','sheet','Heatmap resources Monod','Range','A74:AI95','VariableNamingRule','preserve');

strains=sortrows(strains,'Names');
strains(:,36)=strains(:,3);
strains=renamevars(strains,'Var36','succinate');
strains(:,3)=[];

compounds=strains.Properties.VariableNames;

group_sums=table2array(strains(:,2:end));

% plot


xvalues=compounds(2:end);
yvalues=strains.Names;

close all

figH=figure;

h=heatmap(xvalues, yvalues, group_sums)
h.MissingDataColor=[0.5,0.5,0.5];

%h.ColorScaling = 'scaledcolumns';

saveas(figH,'compound-growthrates-heatmap.pdf','pdf');


% individual bar plots

figH=figure;

t=tiledlayout(5,5)

for i=1:length(yvalues)

nexttile

toplot=sort(group_sums(i,:),'descend');

histogram(toplot,'binWidth',0.1,'normalization','probability')
ylim([0 1]);
xlim([0 1]);
grid on
title(yvalues{i});

end
title(t,'growth rate distributions');

saveas(figH,'growth-rate-distributions.pdf','pdf');

%plot with the log normal distributions to scale

group_sums(group_sums==0)=nan;

close all

figH=figure;

t=tiledlayout(5,5)

for i=1:length(yvalues)

nexttile

toplot=sort(group_sums(i,:),'descend');

h=histogram(toplot,'binWidth',0.05,'normalization','probability');
hold on
pd=fitdist(toplot','LogNormal');
x_values = 0:0.05:1;
y = pdf(pd,x_values);
plot(x_values,(y./max(y)).*max(h.Values),'LineWidth',2);

ylim([0 0.6]);
yticks=[0.2,0.4];
xlim([0 1]);
grid minor
title(yvalues{i});

end
title(t,'growth rate distributions');

saveas(figH,'growth-rate-distributions-lognormal.pdf','pdf');
