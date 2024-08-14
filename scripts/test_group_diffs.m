%% set up
clear all

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu';
fiveyr = readtable(fullfile(basedir,'data',"5YearsFollowUp deidentified reIDed.csv"));

%% test group diffs on pain intensity at 5 yr

outcome = 'pain_avg';
test_group_diffs_at_5yr(outcome, fiveyr)

%% test for group diffs at 5 years, controlling for baseline
% I am not using the MEM b/c I would need to convert from wide to long
% also a MEM seems less apt when I want to test for effects at a specific
% timepoint (5 years), controlling for baseline. I don't see why I would
% control for other timepoints (eg 1 year)

function test_group_diffs_at_5yr(outcome, fiveyr)

bl = [outcome '_baseline'];

% contrast coding
fiveyr.PRT_vs_combined = fiveyr.group == 1;
fiveyr.PLA_vs_UC = .5 * (fiveyr.group == 2) + -.5 * (fiveyr.group == 3);

fitglm(fiveyr, [outcome ' ~ PRT_vs_combined + PLA_vs_UC + ' bl])

end