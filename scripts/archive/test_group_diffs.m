%% set up
clear all

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu';
fiveyr = readtable(fullfile(basedir,'data',"5YearsFollowUp deidentified reIDed scale scores only.csv"));

wh1 = fiveyr.group==1;
wh2 = fiveyr.group==2;
wh3 = fiveyr.group==3;

%% test group diffs on pain intensity at 5 yr

outcome = 'pain_avg';
mdl = test_group_diffs_at_5yr(outcome, fiveyr)

%% Outcomes Table: for each outcome, show M (SD) for each group and test of group diffs
clc
fprintf('outcome\tPRT\tPLA\tUC\tPRTvsCom\tp\tPLAvsUC\tp\n')
for i = [3 21:35]
    outcome_name = fiveyr.Properties.VariableNames{i};
    outcome = fiveyr{:,outcome_name};
    var1 = outcome(wh1); var2 = outcome(wh2); var3 = outcome(wh3);

    % if we have baseline data
    if i < 30
        mdl = test_group_diffs_at_5yr(outcome_name, fiveyr);
    
        fprintf(['%s\t%3.2f (%3.2f)\t%3.2f (%3.2f)\t%3.2f (%3.2f)\t' ...
            '%3.2f\t%0.3f\t%3.2f\t%0.3f\n'], ...
        outcome_name, nanmean(var1), nanstd(var1), nanmean(var2), nanstd(var2), nanmean(var3), nanstd(var3), ...
        mdl.Coefficients.Estimate(3), mdl.Coefficients.pValue(3), mdl.Coefficients.Estimate(4), mdl.Coefficients.pValue(4));
        %gOLPvsTAU.hedgesg, asts3, gOLPvsTAU.hedgesgCi, bayesF, bayesInt)
    
    else % no baseline data
        mdl = test_group_diffs_at_5yr_nobl(outcome_name, fiveyr);
    
        fprintf(['%s\t%3.2f (%3.2f)\t%3.2f (%3.2f)\t%3.2f (%3.2f)\t' ...
            '%3.2f\t%0.3f\t%3.2f\t%0.3f\n'], ...
        outcome_name, nanmean(var1), nanstd(var1), nanmean(var2), nanstd(var2), nanmean(var3), nanstd(var3), ...
        mdl.Coefficients.Estimate(2), mdl.Coefficients.pValue(2), mdl.Coefficients.Estimate(3), mdl.Coefficients.pValue(3));
        
    end
end


%% test for group diffs at 5 years, controlling for baseline
% I am not using the MEM b/c I would need to convert from wide to long
% also a MEM seems less apt when I want to test for effects at a specific
% timepoint (5 years), controlling for baseline. I don't see why I would
% control for other timepoints (eg 1 year)

function mdl = test_group_diffs_at_5yr(outcome, fiveyr)

bl = [outcome '_baseline'];

% contrast coding
fiveyr.PRT_vs_combined = fiveyr.group == 1;
fiveyr.PLA_vs_UC = .5 * (fiveyr.group == 2) + -.5 * (fiveyr.group == 3);

mdl = fitglm(fiveyr, [outcome ' ~ PRT_vs_combined + PLA_vs_UC + ' bl]);

end


%% test for group diffs at 5 years, no baseline
% if baseline not available

function mdl = test_group_diffs_at_5yr_nobl(outcome, fiveyr)

% contrast coding
fiveyr.PRT_vs_combined = fiveyr.group == 1;
fiveyr.PLA_vs_UC = .5 * (fiveyr.group == 2) + -.5 * (fiveyr.group == 3);

mdl = fitglm(fiveyr, [outcome ' ~ PRT_vs_combined + PLA_vs_UC']);

end