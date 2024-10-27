%% set up
clear all

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu';
fiveyr = readtable(fullfile(basedir,'data',"5YearsFollowUp deidentified reIDed scale scores only.csv"));

load('/Users/yoni/Repositories/OLP4CBP/data/final_12mo_outcomes_long.mat');

%% prep 5 yr data

% drop all baseline values from five yr
wh = contains(fiveyr.Properties.VariableNames, '_baseline');
fiveyr(:,wh) = [];

% for now, drop attributions data; we don't have that at baseline
wh = contains(fiveyr.Properties.VariableNames, 'attribution');
fiveyr(:,wh) = [];

% add NaN age and gender
fiveyr.age = NaN(height(fiveyr),1);
fiveyr.gender = NaN(height(fiveyr),1);

%% prep baseline data

% retain only baseline from acute data and remove columns we don't have at
% 5 yrs
outcomes_long = outcomes_long(outcomes_long.redcap_event_name == 'baseline',[1 3:17 43 45])
outcomes_long(:,'ceq') = [];

% rearrange order
columnsToMove = {'bpi_intensity', 'bpi_interference'};
columnsToKeep = setdiff(outcomes_long.Properties.VariableNames, columnsToMove, 'stable');
outcomes_long = outcomes_long(:, [columnsToKeep(1:13), columnsToMove columnsToKeep(14:end)]);
 
columnsToMove = {'promis_sleep'};
columnsToKeep = setdiff(outcomes_long.Properties.VariableNames, columnsToMove, 'stable');
outcomes_long = outcomes_long(:, [columnsToKeep(1:5) columnsToMove columnsToKeep(6:end)]);

% rescale -- no need
% outcomes_wide.pcs = outcomes_wide.pcs * 13;

% rename sopa
outcomes_long.Properties.VariableNames{13} = 'sopa_control';

clc
head(outcomes_long)
head(fiveyr)

%% anyone at 5 yr but not baseline??
setdiff(fiveyr.id, outcomes_long.id)

%% add age and gender to five yr

% Loop through subjects
for i = 1:height(fiveyr)
    fiveyr.age(i) = outcomes_long.age(outcomes_long.id == fiveyr.id(i)); 
    fiveyr.gender(i) = outcomes_long.gender(outcomes_long.id == fiveyr.id(i)); 
end


%% concat the baseline and the 5 year data

% add time variable
outcomes_long.time = zeros(height(outcomes_long), 1) - 1;
fiveyr.time = ones(height(fiveyr), 1);

% concat
fiveyr_long = [outcomes_long; fiveyr];

%% write to file
writetable(fiveyr_long, fullfile(basedir,'data',"five_yr_long_final.csv"));
