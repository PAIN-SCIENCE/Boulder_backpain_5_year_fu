%% set up
clear all

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu';
figdir = fullfile(basedir,'figs');
fiveyr = readtable(fullfile(basedir,'data',"5YearsFollowUp deidentified reIDed.csv"));

colors = [0 0 1; .5 0 .5; .5 .5 .5];


%% violin plots of group diffs at 5 years 
create_figure('davio')


opts = {'box',3, 'colors', colors, 'violin', 'half', ...
    'scatter',1,'jitter', 2, 'jitterspacing', 1, 'scattercolors','same',...
    'scattersize',22,'bins',12,'xtlabels',{'PRT' 'Placebo' 'Usual Care'}};

h = daviolinplot(fiveyr.pain_avg, fiveyr.group, opts{:});

ylabel('Pain Intensity');
ylim([-1 8]);
xlim([0 4]);
set(gca,'FontSize', 20, 'ytick',0:8);

set(gcf,'Position', [  440   287   689   455]);
print(gcf, fullfile(figdir, 'pain_violins.pdf'), '-dpdf');

%% save the colors
cols{1} = h.ds(1).FaceColor;
cols{2} = h.ds(2).FaceColor;
cols{3} = h.ds(3).FaceColor;

%% compute % PF / NPF at 5 years, across 3 groups
clc
fiveyr.pfnpf = fiveyr.pain_avg <= 1;
[percimp, chi, p] = crosstab(fiveyr.pfnpf, fiveyr.group);
n = histcounts(fiveyr.group);
pfnpf = percimp(2,:) ./ n

% chi square test of PF NPF is significant
chi, p

%% compute % PF / NPF at 5 years, PRT vs one group
clc
fiveyr_temp = fiveyr(fiveyr.group==1 | fiveyr.group==2, :); % change the 3 to a 2 to compare PRT vs which group
fiveyr_temp.pfnpf = fiveyr_temp.pain_avg <= 1;
[percimp, chi, p] = crosstab(fiveyr_temp.pfnpf, fiveyr_temp.group);
n = histcounts(fiveyr_temp.group);
if n(2)==0, n=[n(1) n(3)]; end % hack to make it work for groups 1 and 3
pfnpf = percimp(2,:) ./ n

% chi square test of PF NPF is significant
chi, p

%% bar graph % PF / NPF at 5 years

create_figure('bar2')
gray = [.6 .6 .6];
b = bar(pfnpf * 100, 'FaceColor', 'flat');

xtips1 = b.XEndPoints;
ytips1 = b.YEndPoints;    
labels1 = strcat(string(round(b.YData)), '%');
text(xtips1+.08,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'FontSize', 18)

ylabel({'% with pain <= 1'})
set(gca, 'YTick', 0:25:75)
ylim([0 75])

set(gca, 'FontSize', 18, 'xtick', 1:3, 'XTickLabel', {'PRT' 'Placebo' 'Usual Care'})
set(gcf,'Position', [ 1224        1186         300         304])


% Apply the custom colors to each bar
for i = 1:3
    b.CData(i, :) = cols{i};
end
grid on

a = gca; a.XTickLabelRotation = 30
saveas(gcf, fullfile(figdir, 'percent_PFNPF_3groups.pdf'))