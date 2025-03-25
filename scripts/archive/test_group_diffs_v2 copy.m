%% set up
clear all

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu';
fiveyr = readtable(fullfile(basedir,'data',"five_yr_long_final.csv"));

%% test group diffs on pain intensity at 5 yr
clc
outcome = 'pain_avg';
[mdl1, mdl2] = test_group_diffs_at_5yr(outcome, fiveyr)
[smd1, smd2] = smd_at_5yr(outcome, fiveyr)


%% Outcomes Table: for each outcome, show M (SD) for each group and test of group diffs
clc
fprintf('outcome\tPRT\tPLA\tUC\tPRTvsPLA\tp\tSMD\tPLAvsUC\tp\tSMD\n')
for i = 5%3:15
    outcome_name = fiveyr.Properties.VariableNames{i};
    
    % prep for M,SD per group
    var1 = fiveyr{fiveyr.time==1 & fiveyr.group==1,outcome_name};
    var2 = fiveyr{fiveyr.time==1 & fiveyr.group==2,outcome_name};
    var3 = fiveyr{fiveyr.time==1 & fiveyr.group==3,outcome_name};

    % if we have baseline data
    if i < 30
        [mdl1, mdl2] = test_group_diffs_at_5yr(outcome_name, fiveyr);
        [smd1, smd2] = smd_at_5yr(outcome_name, fiveyr);
        adjustedMeans = baseline_adjusted_means_at_5yr(outcome_name, fiveyr)
    
        fprintf(['%s\t%3.2f (%3.2f)\t%3.2f (%3.2f)\t%3.2f (%3.2f)\t' ...
            '%3.2f\t%0.3f\t%3.2f\t%3.2f\t%0.3f\t%3.2f\n'], ...
        outcome_name, nanmean(var1), nanstd(var1), nanmean(var2), nanstd(var2), nanmean(var3), nanstd(var3), ...
        mdl1.Coefficients.Estimate(6), mdl1.Coefficients.pValue(6), smd1, ...
        mdl2.Coefficients.Estimate(6), mdl2.Coefficients.pValue(6), smd2);
        %gOLPvsTAU.hedgesg, asts3, gOLPvsTAU.hedgesgCi, bayesF, bayesInt)

    
        % refit with baseline set to avg baseline -- ?? and say predicted
        % value per group?

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

%% estimate effect size at 5 years
function [smd1, smd2] = smd_at_5yr(outcome, fiveyr)

mytabl = fiveyr(:,{outcome 'group', 'time'});

% SD across all participants at baseline. N is balanced across groups.
pooledbaselineSD = std(mytabl{mytabl.time==-1,1},'omitnan');

g1t1 = mean(fiveyr{fiveyr.time==-1 & fiveyr.group==1,outcome},'omitnan');
g2t1 = mean(fiveyr{fiveyr.time==-1 & fiveyr.group==2,outcome},'omitnan');
g3t1 = mean(fiveyr{fiveyr.time==-1 & fiveyr.group==3,outcome},'omitnan');

g1t2 = mean(fiveyr{fiveyr.time==1 & fiveyr.group==1,outcome},'omitnan');
g2t2 = mean(fiveyr{fiveyr.time==1 & fiveyr.group==2,outcome},'omitnan');
g3t2 = mean(fiveyr{fiveyr.time==1 & fiveyr.group==3,outcome},'omitnan');

% PRT vs PLA
smd1 = ((g1t2 - g1t1) - (g2t2 - g2t1)) / pooledbaselineSD;
smd2 = ((g1t2 - g1t1) - (g3t2 - g3t1)) / pooledbaselineSD;

end

function [adjustedMeans] = baseline_adjusted_means_at_5yr(outcome_name, data)

% to estimate baseline adjusted means, first convert data to short form
% Separate baseline (-1) and post-intervention (1) data
baselineData = data(data.time == -1, {'id' outcome_name 'group', 'time'});
postData = data(data.time == 1, {'id' outcome_name 'group', 'time'});

% Merge baseline values with post-intervention data for ANCOVA
mergedData = outerjoin(postData, baselineData, 'Keys', 'id', ...
    'MergeKeys', true, 'Type', 'Left');
mergedData.Properties.VariableNames{2} = 'Post';
mergedData.Properties.VariableNames{3} = 'Group';
mergedData.Properties.VariableNames{5} = 'Baseline';

% Handle missing baseline or post-intervention values
mergedData = mergedData(~isnan(mergedData.Baseline) & ~isnan(mergedData.Post), :);

% Fit the ANCOVA model
mergedData.Group = categorical(mergedData.Group);
mdl = fitlm(mergedData, 'Post ~ Group + Baseline');
disp(mdl);

% Calculate the mean baseline value
meanBaseline = mean(mergedData.Baseline, 'omitnan') % Handle missing data safely

% Get the coefficients from the model
coefficients = mdl.Coefficients.Estimate;

% Calculate adjusted means for each group
groups = categories(mergedData.Group);
adjustedMeans = table(groups, 'VariableNames', {'Group'});

for i = 1:length(groups)
    if i == 1
        % Reference group
        adjustedMean = coefficients(1) + coefficients(4) * meanBaseline;
    else
        % Other groups (account for the group effect)
        groupEffect = coefficients(i) % Group effect starts at 2nd coefficient; for group 2, it is 2nd coefficient, etc
        adjustedMean = coefficients(1) + groupEffect + coefficients(4) * meanBaseline;
    end
    adjustedMeans.AdjustedMean(i) = adjustedMean;
end

end