%% compute scale scores

orig_len = length(tabl.Properties.VariableNames);

% find the first field of each survey
timestamp_fields=regexp(tabl.Properties.VariableNames, '.*timestamp', 'match');
timestamp_fields = find(~cellfun(@isempty, timestamp_fields));

prev_injection_success_fields=regexp(tabl.Properties.VariableNames, 'treatmentsuccess.*', 'match');
prev_injection_success_fields = find(~cellfun(@isempty, prev_injection_success_fields));

%complete_fields=regexp(tabl.Properties.VariableNames, '.*complete', 'match');
%complete_fields = find(~cellfun(@isempty, complete_fields));
%complete_fields([2 3 end-1]) = []; % drop tlfb, admin, and false hit on "complete"

%tabl.Properties.VariableNames(timestamp_fields)'
%tabl.Properties.VariableNames(complete_fields)'

% compute scale scores
clc
for scale = 2:length(timestamp_fields)-2 % skip elig form, bpi daily, cliexa
    
    scalename = strsplit(tabl.Properties.VariableNames{timestamp_fields(scale)},'_timestamp');
    scalename = scalename{1}
    
    startind = timestamp_fields(scale)+1;
    endind = timestamp_fields(scale+1)-2;
    
    tabl.Properties.VariableNames(startind:endind)

    switch scalename
        
        case 'bpisf_last_week'
            tabl.bpi_intensity = mean(tabl{:, startind:startind+3},2);
            tabl.bpi_interference = mean(tabl{:, startind+4:endind},2);
        
        case 'olbpdq'
            tabl.odi = (sum(tabl{:, startind:endind},2) - 10) / 50 * 100;
            
        case 'promis_4_ceq'
            tabl.promis_dep = sum(tabl{:, startind:startind+7},2);
            tabl.promis_anger = sum(tabl{:, startind+16:startind+20},2);
            tabl.promis_anxiety = sum(tabl{:, startind+21:startind+28},2);
            tabl.promis_sleep = (3*6 - sum(tabl{:, startind + [8 9 15]},2)) + sum(tabl{:, startind + [10:14]},2);
            
            tabl.ceq = tabl{:,endind};
            
        case 'panas_10'
            tabl.panas_pa = sum(tabl{:, startind+ [3 4 5 6 7]},2);
            tabl.panas_na = sum(tabl{:, startind+ [0 1 2 8 9]},2);
            
        case 'pcslast_item_from_ipq'
            tabl.pcs = sum(tabl{:, startind:startind+12},2);
            
        case 'tsk11_sopa_emo_2_items'
            tabl.tsk11 = sum(tabl{:, startind:startind+10},2);
            tabl.sopa_emo = sum(tabl{:, startind+11:startind+12},2);
            
        case 'patient_global_impression_of_change_scale'
            tabl.pgic = tabl.pgics; % important so it will get found below
            
        case 'treatment_satisfaction_survey'
            tabl.tx_satisfaction = mean(tabl{:, startind:startind+1},2);
            
        case 'bmqspecific_modified'
            tabl.bmq = 4*6 - nanmean(tabl{:,startind:endind}, 2); % from 1 - 5 as in codebook. higher is stronger beliefs in potency (i am reverse scoring it).
            
        case 'back_pain_injection_treatment_history'
            tabl.prev_injections_number = tabl.injectionnumber;
            tabl.prev_injections_mean_success = nanmean(tabl{:,prev_injection_success_fields},2);
            tabl.prev_injections_max_success = nanmax(tabl{:,prev_injection_success_fields},[],2);
            
        case 'lotr'
            tabl.lotr = (7*6 - nansum(tabl{:, startind-1 + [1 2 4:6 8 10]},2)) + nansum(tabl{:, startind - 1 + [3 7 9]},2); % scoring is from 1 - 5.  items 3 7 9 reverse scored. higher is more optimism            
            
        case 'fear_of_pain_questionnaire'
            tabl.fear_of_pain = nanmean(tabl{:,startind:endind}, 2); % from 1 - 5 as in codebook and in original form. higher is more fear
            
        case 'generalized_selfefficacy_scale'
            tabl.ges = nanmean(tabl{:,startind:endind}, 2); % from 1 - 4 as in codebook and in original form. higher is more self-efficacy
            
        case 'mindful_attention_awareness_scale'
            tabl.maas = nanmean(tabl{:,startind:endind}, 2); % from 1 - 6 as in codebook and in original form. higher is more mindful
            
        case 'emotional_regulations_questionnaire'
            tabl.erq_cog_reapp = mean(tabl{:, startind-1+[1, 3, 5, 7, 8, 10]},2);
            tabl.erq_suppress = mean(tabl{:, startind-1+[2, 4, 6, 9]},2);            
        
        case 'aceclarke' % collected only at elig
            tabl.ace = nansum(tabl{:,startind:endind-1}, 2);
            tabl.ace_clarke = tabl{:,endind};
            
        otherwise
            warning('Skipping %s', scalename)
    end
end


