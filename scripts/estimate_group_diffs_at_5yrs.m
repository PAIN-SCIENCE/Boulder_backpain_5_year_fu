%% set up
clear all

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu/data';
fiveyr = readtable(fullfile(basedir,"5YearsFollowUp deidentified reIDed.csv"));

%% model outcomes -- need to think through!

% model_olp4cbp_outcomes_at5yr('pain_avg', fiveyr)

%% plot



violin_plot_groups(fiveyr.pain_avg, fiveyr.group, 'Avg Pain')


function violin_plot_groups(outcome, grp, label)
    wh1 = grp==1; wh2 = grp==2; wh3 = grp==3;
    create_figure('change scores')
    violinplot({outcome(wh1), outcome(wh2), outcome(wh3)}, 'mc', 'k') 
    set(gca, 'XTickLabel', {'PRT', 'PLA', 'UC'}, 'FontSize', 24), ylabel(label)
    legend off
    hline(0, 'k--')
%    saveas(gcf, fullfile(figdir, 'change_scores_pain.svg'))
end