% This script will merge the 5 yr data with the REDCap IDs from the
% original study. Once we have those IDs, we can merge with various other
% prepared datasets from that study

%% load both datasets

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu/data';

orig = readtable(fullfile(basedir,"OriginalProject Names RecordIds IDENTIFIED.csv"));
orig = orig(:, {'record_id', 'first_name', 'last_name'});
orig.name = lower(strcat(orig.first_name, '_', orig.last_name));


fiveyr = readtable(fullfile(basedir,"5YearsFollowUp IDENTIFIED.csv"));
fiveyr.name = lower(strcat(fiveyr.first_name, '_', fiveyr.last_name));

%% find duplicate entries and remove

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

% go through chunk by chunk and make edits to line it up
fiveyr(100:end,{'record_id' 'name' 'found'})


%%
join(fiveyr, orig, "Keys","name")

