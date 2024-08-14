%% set up
clear all

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu';
fiveyr = readtable(fullfile(basedir,'data',"5YearsFollowUp deidentified reIDed.csv"));

%% model outcomes -- need to think through!

model_olp4cbp_outcomes_at5yr('pain_avg', fiveyr)
