% b/c there is no baseline (post-treatment) data for attribution ratings
% and rankings, the code here is different.

% note there IS baseline data on the text-based attributions -- that is not
% address in this script

% EXECUTIVE SUMMARY OF FINDINGS
% A) on VAS questions, we see clear, super strong group diffs
% B) on rank-order questions, of the 7 categories pertaining to mind/brain,
% each individual category is never ranked high for a group on avg (4.5 is
% the highest avg rank). There are no sig group diffs for PRT vs PLA. We do
% see sig group diffs for PRT vs UC on three items: as expected, items are
% ranked lower for PRT.
% C) across the above two outcomes, we DO see placebo effects on
% attributions!
% 

%% set up
clear all

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu';
fiveyrAttr = readtable(fullfile(basedir,'data',"5YearsFollowUp deidentified reIDed.csv"));

% begin and end of attributions
a=find(strcmp(fiveyrAttr.Properties.VariableNames, 'chronic_pain_attribution_scale_v2_timestamp'));
b=find(strcmp(fiveyrAttr.Properties.VariableNames, 'chronic_pain_attribution_scale_v2_complete'));

fiveyrAttr = fiveyrAttr(:, [1 2 a+1:b-1]);

wh1 = fiveyrAttr.group==1;
wh2 = fiveyrAttr.group==2;
wh3 = fiveyrAttr.group==3;

%% rename a few vars
fiveyrAttr.Properties.VariableNames(9:11) = {'STR_rating' 'MB_rating', 'STRvsMB_rating'};

%% merge in age and gender
forAG = readtable(fullfile(basedir,'data',"five_yr_long_final.csv"));
forAG = forAG(forAG.time==1, {'id', 'age', 'gender', 'tsk11', 'pain_avg', 'pcs'});

fiveyrAttr = join(fiveyrAttr, forAG);


%% group diffs in rank order values

ranked_cats_MB = {
'changes_in_my_brain_s_pain'
'stress'
'fear_of_anxiety'
'some_other_emotion_e_g_fea'
'personality_traits_e_g_per'
'relationships_with_other_p'
'childhood_experiences_that'
};

%% test for group diffs
fprintf('outcome\tPRT\tPLA\tUC\tPRTvsUC\tp\tPLAvsUC\tp\n')
outcomes = {'STR_rating' 'MB_rating', 'STRvsMB_rating'};
%outcomes = ranked_cats_MB;

clc
fprintf('outcome\tPRT\tPLA\tUC\tPRTvsPLA\tp\tPLAvsUC\tp\n')

for i=1:length(outcomes)
    outcome_name = outcomes{i};
    outcome = fiveyrAttr{:,outcome_name};
    var1 = outcome(wh1); var2 = outcome(wh2); var3 = outcome(wh3);
    
    [mdl1, mdl2] = test_group_diffs_at_5yr_nobl(outcome_name, fiveyrAttr);
    [smd1, smd2] = smd_at_5yr_nobl(outcome_name, fiveyrAttr);
    
    fprintf(['%s\t%3.2f (%3.2f)\t%3.2f (%3.2f)\t%3.2f (%3.2f)\t' ...
                '%3.2f\t%0.3f\t%3.2f\t%0.3f\t%3.2f\t%0.3f\n'], ...
            outcome_name, nanmean(var1), nanstd(var1), nanmean(var2), nanstd(var2), nanmean(var3), nanstd(var3), ...
            mdl1.Coefficients.Estimate(4), mdl1.Coefficients.pValue(4), smd1, ...
            mdl2.Coefficients.Estimate(4), mdl2.Coefficients.pValue(4), smd2);
end

%% test whether pain at 5 yrs is correlated within group with two mechs: TSK, MBattr
clc
fiveyrAttr.group = nominal(fiveyrAttr.group);
fitlme(fiveyrAttr, 'pain_avg ~ pcs + group')


%% test group diffs at 5 years; no pre-treatment values, no effect of time
% beta represents group difference magnitude, controlling for age and
% gender
function [mdl1, mdl2] = test_group_diffs_at_5yr_nobl(outcome, fiveyr)

% contrast coding
fiveyr.PRT_vs_PLA = .5 * (fiveyr.group == 1) + -.5 * (fiveyr.group == 2);
fiveyr.PRT_vs_UC = .5 * (fiveyr.group == 1) + -.5 * (fiveyr.group == 3);

fiveyr.PRT_vs_PLA = fiveyr.PRT_vs_PLA + .5;

mdl1 = fitlme(fiveyr, [outcome ' ~ PRT_vs_PLA + age + gender']);
mdl2 = fitlme(fiveyr, [outcome ' ~ PRT_vs_UC + age + gender']);

end

% compute SMD at 5 years 
function [smd1, smd2] = smd_at_5yr_nobl(outcome, fiveyr)

g1 = fiveyr{fiveyr.group==1,outcome};
g2 = fiveyr{fiveyr.group==2,outcome};
g3 = fiveyr{fiveyr.group==3,outcome};

%smd1=meanEffectSize(g1,g2,effect="glass");
%smd2=meanEffectSize(g1,g3,effect="cohen");

% for consistency w/ how i compute SMD for other 2ndary outcomes
smd1=(mean(g1,"omitnan")-mean(g2,"omitnan"))/std([g1; g2; g3],"omitnan");
smd2=(mean(g1,"omitnan")-mean(g3,"omitnan"))/std([g1; g2; g3],"omitnan");

end