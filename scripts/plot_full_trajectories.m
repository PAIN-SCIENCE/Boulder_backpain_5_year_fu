clear all

% load in data through 1 yr
datadir = '/Users/yoni/Repositories/OLP4CBP/data';
load(fullfile(datadir, 'final_12mo_outcomes_wide.mat'));

% load in data at 5 yr
basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu';
fiveyr = readtable(fullfile(basedir,'data',"five_yr_long_final.csv"));
fiveyr = fiveyr(fiveyr.time==1,{'id', 'pain_avg'});

colors = [0 0 1; .5 0 .5; .5 .5 .5];
figdir = fullfile(basedir,'figs');

%% drop 6 Ps w/ missing baseline
nobaseline = setdiff(fiveyr.id, outcomes_wide.id);
fiveyr = fiveyr(~ismember(fiveyr.id, nobaseline),:);
fiveyr.Properties.VariableNames{2} = 'pain_avg_5yr'; 

%% merge
fullfive = outerjoin(outcomes_wide, fiveyr, 'Keys', 'id');

wh1 = fullfive.group==1;
wh2 = fullfive.group==2;
wh3 = fullfive.group==3;

%% plot

create_figure('main outcome by group')

timepoints = strcat('pain_avg_', {'baseline' 't2_arm_1' 'x1_month_follow_up_arm_1' 'x2_month_follow_up_arm_1' 'x3_month_follow_up_arm_1' 'x6_month_follow_up_arm_1' 'x12_month_follow_up_arm_1', '5yr'});

xvals = [1:5 7 10 14]; %[1:(length(timepoints)-1) (length(timepoints)+3)];

xlab = strrep(timepoints, 'pain_avg', '');
xlab = strrep(xlab, '_arm_1', '');
xlab = strrep(xlab, '_follow_up', '');
xlab = strrep(xlab, '_', ' ');
xlab = strrep(xlab, 'x', '');
xlab = strrep(xlab, 'baseline', 'Baseline');
xlab = strrep(xlab, 't2', 'Post-tx');

set(gca, 'XTick', xvals, 'XTickLabel', xlab, 'FontSize', 22)
ylim([0 5]) 
ylabel('Pain intensity (0 - 10)'), xlim([xvals(1) xvals(end)])
ax = gca;
ax.XTickLabelRotation = 40;

h1 = lineplot_columns(table2array(fullfive(wh1, timepoints)), 'shade', 'x', xvals, 'color', colors(1,:));
h2 = lineplot_columns(table2array(fullfive(wh2, timepoints)), 'shade', 'x', xvals, 'color', colors(2,:));
h3 = lineplot_columns(table2array(fullfive(wh3, timepoints)), 'shade', 'x', xvals, 'color', colors(3,:));

legend([h1.line_han h2.line_han h3.line_han],{'PRT' 'Placebo' 'Usual Care'}, 'Location', 'Northeast')
set(gcf, 'Position', [     399   164   563   391])

print(gcf, fullfile(figdir, 'pain_trajectories.pdf'), '-dpdf');