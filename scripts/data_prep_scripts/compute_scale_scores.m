%% set up
clear all

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu';
tabl = readtable(fullfile(basedir,'data',"5YearsFollowUp deidentified reIDed.csv"));

%% compute scale scores

% find the first field of each survey
timestamp_fields=regexp(tabl.Properties.VariableNames, '.*timestamp', 'match');
timestamp_fields = find(~cellfun(@isempty, timestamp_fields));

% compute scale scores
clc
for scale = 2:length(timestamp_fields)-1 % skip elig form at beginning and body map at end
    
    scalename = strsplit(tabl.Properties.VariableNames{timestamp_fields(scale)},'_timestamp');
    scalename = scalename{1}
    
    startind = timestamp_fields(scale)+1;
    endind = timestamp_fields(scale+1)-2;
    
    tabl.Properties.VariableNames(startind:endind)

    switch scalename
        
        case 'bpisf_lastweek_pain_ratings'
            tabl.bpi_intensity = mean(tabl{:, startind:startind+3},2);
            tabl.bpi_interference = mean(tabl{:, startind+4:endind},2);
        
        case 'olbpdq'
            tabl.odi = (sum(tabl{:, startind:endind},2) - 10) / 50 * 100; % #confirmed correct, Yoni 8.15.24
            
        case 'promis_4'
            tabl.promis_dep = sum(tabl{:, startind:startind+7},2);
            tabl.promis_sleep = (3*6 - sum(tabl{:, startind + [8 9 15]},2)) + sum(tabl{:, startind + [10:14]},2);
            tabl.promis_anger = sum(tabl{:, startind+16:startind+20},2);
            tabl.promis_anxiety = sum(tabl{:, startind+21:startind+28},2);      
            
        case 'panas_10'
            tabl.panas_pa = sum(tabl{:, startind+ [3 4 5 6 7]},2);
            tabl.panas_na = sum(tabl{:, startind+ [0 1 2 8 9]},2);
            
        case 'pcs'
            tabl.pcs = sum(tabl{:, startind:startind+12},2);
            
        case 'tsk11'
            tabl.tsk11 = sum(tabl{:, startind:startind+10},2);
            
        case 'chronic_pain_attribution_scale_v2'
            tabl.attribution_structural = tabl.chronic_pain_attribution_3a;
            tabl.attribution_mindbrain = tabl.chronic_pain_attribution_4;
            tabl.attribution_ratio = tabl.chronic_pain_attributio_5;
            
            warning('SKIPPING ATTRIBUTION FREE TEXT AND RANKING DATA FOR NOW')

        case 'sopacontrol_subscale'
            pos_sum = sum(tabl{:, startind:startind+5},2);
            rev_sum = (4*6) - sum(tabl{:, startind+6:endind},2); %scale ranges from 1 to 5. want 5s to become 1s. We have 4 items
            tabl.sopa_control = (pos_sum + rev_sum) / 10; % ten items

        case 'usual_care_measure_last_6_months'
             warning('SKIPPING USUAL CARE -- NO MEAINGFUL WAY TO SUM')
            
        otherwise
            warning('Skipping %s', scalename)
    end
end

%% save a full version with all columns, and a subset with scale scores


writetable(tabl, fullfile(basedir,'data',"5YearsFollowUp deidentified reIDed.csv"));

writetable(tabl(:,[1 2 10 214:end]), fullfile(basedir,'data',"5YearsFollowUp deidentified reIDed scale scores only.csv"));
