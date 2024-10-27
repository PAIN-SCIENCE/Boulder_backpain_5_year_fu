%% set up
clear all

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu';
fiveyr = readtable(fullfile(basedir,'data',"five_yr_long_final.csv"));

%% test group diffs on pain intensity at 5 yr
clc
outcome = 'odi';
[mdl1, mdl2] = test_group_diffs_at_5yr(outcome, fiveyr)

%% Outcomes Table: for each outcome, show M (SD) for each group and test of group diffs
clc
fprintf('outcome\tPRT\tPLA\tUC\tPRTvsCom\tp\tPLAvsUC\tp\n')
for i = [3]% 21:35]
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


%% MEM testing for group diffs at 5 years, controlling for baseline
function [mdl1, mdl2] = test_group_diffs_at_5yr(outcome, fiveyr)

% contrast coding
fiveyr.PRT_vs_PLA = .5 * (fiveyr.group == 1) + -.5 * (fiveyr.group == 2);
fiveyr.PRT_vs_UC = .5 * (fiveyr.group == 1) + -.5 * (fiveyr.group == 3);

mdl1 = fitlme(fiveyr, [outcome ' ~ time*PRT_vs_UC + age + gender + (1|id)']);
mdl2 = fitlme(fiveyr, [outcome ' ~ time*PRT_vs_PLA + age + gender +(1|id)']);

end