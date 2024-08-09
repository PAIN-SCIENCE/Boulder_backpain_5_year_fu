%% set up
clear all

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu';
figdir = fullfile(basedir,'figs');
fiveyr = readtable(fullfile(basedir,'data',"5YearsFollowUp deidentified reIDed.csv"));

%% violin plots of group diffs at 5 years 
create_figure('davio')

opts = { 'box',3,'boxwidth',0.8,'boxcolors','k',...
    'scatter',2,'jitter', 2,'scattercolors','same',...
    'scattersize',14,'bins',12,'xtlabels',{'PRT' 'Placebo' 'Usual care'}};

daviolinplot(fiveyr.pain_worst, fiveyr.group, opts{:});

set(gcf,'Position', [  440   287   418   411]);
ylabel('Pain Intensity at 5 years')
ylim([-1 8])
set(gca,'FontSize', 20, 'ytick',0:8)

save(fullfile(figdir, 'pain_violins.pdf'))