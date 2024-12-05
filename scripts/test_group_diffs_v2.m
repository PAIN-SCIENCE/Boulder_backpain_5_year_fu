%% set up
clear all

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu';
fiveyr = readtable(fullfile(basedir,'data',"five_yr_long_final.csv"));

%% test group diffs on pain intensity at 5 yr
clc
outcome = 'pain_avg';
[mdl1, mdl2] = test_group_diffs_at_5yr(outcome, fiveyr)

%% Outcomes Table: for each outcome, show M (SD) for each group and test of group diffs
clc
fprintf('outcome\tPRT\tPLA\tUC\tPRTvsPLA\tp\tPLAvsUC\tp\n')
for i = 3:15
    outcome_name = fiveyr.Properties.VariableNames{i};
    
    % prep for M,SD per group
    var1 = fiveyr{fiveyr.time==1 & fiveyr.group==1,outcome_name};
    var2 = fiveyr{fiveyr.time==1 & fiveyr.group==2,outcome_name};
    var3 = fiveyr{fiveyr.time==1 & fiveyr.group==3,outcome_name};

    % if we have baseline data
    if i < 30
        [mdl1, mdl2] = test_group_diffs_at_5yr(outcome_name, fiveyr);
    
        fprintf(['%s\t%3.2f (%3.2f)\t%3.2f (%3.2f)\t%3.2f (%3.2f)\t' ...
            '%3.2f\t%0.3f\t%3.2f\t%0.3f\n'], ...
        outcome_name, nanmean(var1), nanstd(var1), nanmean(var2), nanstd(var2), nanmean(var3), nanstd(var3), ...
        mdl1.Coefficients.Estimate(6), mdl1.Coefficients.pValue(6), mdl2.Coefficients.Estimate(6), mdl2.Coefficients.pValue(6));
        %gOLPvsTAU.hedgesg, asts3, gOLPvsTAU.hedgesgCi, bayesF, bayesInt)
    
    else % no baseline data
        warning('no baseline?')
        mdl = test_group_diffs_at_5yr_nobl(outcome_name, fiveyr);
    
        fprintf(['%s\t%3.2f (%3.2f)\t%3.2f (%3.2f)\t%3.2f (%3.2f)\t' ...
            '%3.2f\t%0.3f\t%3.2f\t%0.3f\n'], ...
        outcome_name, nanmean(var1), nanstd(var1), nanmean(var2), nanstd(var2), nanmean(var3), nanstd(var3), ...
        mdl.Coefficients.Estimate(2), mdl.Coefficients.pValue(2), mdl.Coefficients.Estimate(3), mdl.Coefficients.pValue(3));
        
    end
end


%% MEM testing for group diffs at 5 years, controlling for baseline
function [mdl1, mdl2] = test_group_diffs_at_5yr(outcome, fiveyr)

% use reference (dummy) coding, as we have a meaningful reference point:
% pre-treatment for time; control group for group. Group gets automatically
% dummy coded in fitlme.

% reference code time to 0/1 for pre/post
fiveyr.time = fiveyr.time > 0;

mdl1 = fitlme(fiveyr, [outcome ' ~ time*group + age + gender +(1|id)'], 'Exclude', fiveyr.group==3);
mdl2 = fitlme(fiveyr, [outcome ' ~ time*group + age + gender + (1|id)'], 'Exclude', fiveyr.group==2);

% old; contrast coding
%fiveyr.PRT_vs_PLA = .5 * (fiveyr.group == 1) + -.5 * (fiveyr.group == 2);
%fiveyr.PRT_vs_UC = .5 * (fiveyr.group == 1) + -.5 * (fiveyr.group == 3);
% old
% fiveyr.time = (fiveyr.time / 2) + .5; % to go from -1/1 to -.5/.5
%fiveyr.PRT_vs_PLA = fiveyr.PRT_vs_PLA + .5;

end