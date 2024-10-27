%% set up
clear all

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu';
fiveyr = readtable(fullfile(basedir,'data',"5YearsFollowUp deidentified reIDed scale scores only.csv"));

load('/Users/yoni/Repositories/OLP4CBP/data/final_acute_outcomes_wide.mat');
load('/Users/yoni/Repositories/OLP4CBP/data/final_12mo_outcomes_long.mat');

%% prep 5 yr data

% drop all baseline values from five yr
wh = contains(fiveyr.Properties.VariableNames, '_baseline');
fiveyr(:,wh) = [];

% for now, drop attributions data; we don't have that at baseline
wh = contains(fiveyr.Properties.VariableNames, 'attribution');
fiveyr(:,wh) = [];

fiveyr.time = zeros(height(fiveyr), 1);

% for now, drop sopa -- don't have baseline data for that in this baseline file...


%% prep baseline data

% retain only baseline from acute data
wh=contains(outcomes_wide.Properties.VariableNames, '_baseline');
wh(1:2) = 1; % id and group

outcomes_wide(:,~wh) = [];

% remove '_baseline' from names
outcomes_wide.Properties.VariableNames = strrep(outcomes_wide.Properties.VariableNames, '_baseline', '');

% remove last three: drug use
outcomes_wide(:, end-4:end) = [];

% rearrange order
columnsToMove = {'bpi_intensity', 'bpi_interference'};
outcomes_wide = outcomes_wide(:, [setdiff(outcomes_wide.Properties.VariableNames, columnsToMove, 'stable'), columnsToMove]);

columnsToMove = {'promis_sleep'};
columnsToKeep = setdiff(outcomes_wide.Properties.VariableNames, columnsToMove, 'stable');
outcomes_wide = outcomes_wide(:, [columnsToKeep(1:5) columnsToMove columnsToKeep(6:end)]);

% rescale -- no need
% outcomes_wide.pcs = outcomes_wide.pcs * 13;

outcomes_wide.time = ones(height(outcomes_wide), 1) * 5;

clc
head(outcomes_wide)
head(fiveyr)

%% concat the baseline and the 5 year data

