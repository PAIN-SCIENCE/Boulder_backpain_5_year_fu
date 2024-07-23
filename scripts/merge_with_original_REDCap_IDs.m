% This script will merge the 5 yr data with the REDCap IDs from the
% original study. Once we have those IDs, we can merge with various other
% prepared datasets from that study

clear all

%% load both datasets

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu/data';

orig = readtable(fullfile(basedir,"OriginalProject Names RecordIds IDENTIFIED.csv"));
orig = orig(:, {'record_id', 'first_name', 'last_name'});
orig.name = lower(strcat(orig.first_name, '_', orig.last_name));


fiveyr = readtable(fullfile(basedir,"5YearsFollowUp IDENTIFIED.csv"));
fiveyr.name = lower(strcat(fiveyr.first_name, '_', fiveyr.last_name));

%% find duplicate names in original data and remove

myCellArray = orig.name;

% Find duplicates
[uniqueValues, ~, idxDuplicates] = unique(myCellArray, 'stable');
duplicates = uniqueValues(hist(idxDuplicates, unique(idxDuplicates)) > 1);

% Display duplicates
disp('Duplicates:');
disp(duplicates);

% are any of the duplicates also in 5 yr data? NO.
any(ismember(fiveyr.name, duplicates))

% OK to drop all duplicates
orig = orig(~ismember(orig.name, duplicates), :);


%% Who in 5 yr does not match with orig?

fiveyr.found = ismember(fiveyr.name, orig.name);

% go through chunk by chunk. edit name typos etc. in data files to line it up
% fiveyr(100:end,{'record_id' 'name' 'found'}); 

% any five year Ps not found?
any(~fiveyr.found)

%% merge with original IDs, then de-identify
fiveyr_deid = join(fiveyr, orig, "Keys","name");
fiveyr_deid = removevars(fiveyr_deid, {'name', 'redcap_survey_identifier', 'first_name_fiveyr', 'first_name_orig', 'last_name_fiveyr', 'last_name_orig', 'email', 'found', 'compensation_timestamp'  'which_store_would_you_like' 'giftcard_by_email'  'compensation_complete' });

%% merge in group assignment

% -- WHEN YOU LOAD IN THE 12 MONTH DATA, SOME Ps ARE MISSING! NEED TO
% EXPLORE WHY, LATER ---

load('/Users/yoni/Repositories/OLP4CBP/data/final_acute_outcomes_wide.mat');
outcomes_wide = outcomes_wide(:, {'id', 'group'});
outcomes_wide.Properties.VariableNames{1} = 'record_id_orig';

fiveyr_deid = join(fiveyr_deid, outcomes_wide, "Keys","record_id_orig");

%% reorder, rename id, and save
fiveyr_deid = [fiveyr_deid(:,end-1:end) fiveyr_deid(:,1:end-2)];
fiveyr_deid.Properties.VariableNames{1} = 'id';

writetable(fiveyr_deid, fullfile(basedir, "5YearsFollowUp deidentified reIDed.csv"))
