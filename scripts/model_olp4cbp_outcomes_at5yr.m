% Run a mixed-effect model for the OLP4CBP study with standard covariates
% and model set up.
%
% :Usage:
% ::
%     [lme_groupbytime, festats_groupbytime, festats_3groupsbytime, festats_prt_vs_combined, lme_p_vs_hc, festats_p_vs_hc, festats_simple_time1, festats_simple_time2, mm_groupbytime, mm_p_vs_hc] = model_olp4cbp_outcomes(outcome, covariates_table, varargin)
%
% ..
%
%     Copyright (C) 2020, Yoni Ashar, Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **outcome:**
%        A variable to model as the outcome of the mixed effects model, 1 x
%        n observations
%
%   **covariates_table**
%        A table including all the standard covariates. Can include other
%        columns as well. This MUST be lined up with outcome.
%
% :Optional Inputs:
%   **'pain', 'pain_avg'**
%        Will include pain (pain_avg) as a covariate, **CENTERED** within group
%        and within time and within healthies, for the group by time model and the patients vs.
%        healthies model. Pain_avg must be in the covariates_table
%
%   **'bladder_pain_mean'**
%        Will include average bladder pain rating as a covariate, **CENTERED** within group
%        and within time and within healthies, for the group by time model and the patients vs.
%        healthies model. 'bladder_pain_mean' must be in the covariates_table
%
%   **'sponpain_mean', 'sponpain_mean_rating'**
%        Will include average spon pain rating as a covariate, **CENTERED** within group
%        and within time and within healthies, for the group by time model and the patients vs.
%        healthies model. 'sponpain_mean_rating' must be in the covariates_table
%
%   **'backpaindur', 'backpain_length'**
%        Will include back pain duration (in years) as a covariate
%        Note this is strongly correlated with age
%
%   **'sponpain_motion'**
%        Will include head motion in sponpain run as a covariate (n_spike_regs_sponpain -- this variable must be covariates_table)
%
%   **'bladder_motion'**
%        Will include head motion in bladder run as a covariate (n_spike_regs -- this variable must be covariates_table)
%
%   **'noverbose'**
%        Do not show model output
%        
%   **'standardize'**
%        Standardize all continuous predictors and the outcome, providing
%        partial correlation coefficient estimates. For an important note
%        on interpretation, see http://www.daviddisabato.com/blog/2016/4/8/on-effect-sizes-in-multiple-regression
%
%   **'p_vs_hc_only'**
%        Do not run the group by time, or analyses of PRT vs. Control
%        diffs. This is useful for outcomes that are only of interest or
%        meaningful at baseline.
%
%   **'covs'**
%        Followed by a matrix of other covariates and a cell array of covariate names
%
%   **'doplots'**
%        Plot group by time interaction, and patients vs. healthies.  TO
%        IMPLEMENT
%
%
% :Outputs:
%
%        Results from the group by time model and results from the model comparing pre-tx
%        patients to healthies. For now, For now, both models drop placebo
%        patients. Also includes the marginal means tables for the group by
%        time and for the patients vs. healthy controls models
%
% :Examples:
% ::
%
%    % give examples of code here
%    model_olp4cbp_outcomes(hdi.betweenness, bs.subject_metadata)
%
%
% :See also:
%   - model_olp4cbp_wrapper.m, which calls this script for multiple
%   outcomes and collects the results in a usable format
%
% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..
function [lme_groupbytime, festats_groupbytime, festats_3groupsbytime, festats_prt_vs_combined, lme_p_vs_hc, festats_p_vs_hc, festats_simple_time1, festats_simple_time2, mm_groupbytime, mm_p_vs_hc, festats_OLP_vs_UC] = model_olp4cbp_outcomes_at5yr(outcome, covariates_table, varargin)

% parse input
% ----------------------------------------------------------------------
doplots = false;
verbose = true;
pain_covariate = false;
bladder_pain_covariate = false;
sponpain_covariate = false;
sponpain_motion = false;
bladder_motion = false;
backpaindur = false;
standardize = false;
p_vs_hc_only = false;

olp_vs_uc_only = false;

if any(strcmp(varargin, 'pain')) || any(strcmp(varargin, 'pain_avg')), pain_covariate = true; end
if any(strcmp(varargin, 'bladder_pain_mean')), bladder_pain_covariate = true; end
if any(strcmp(varargin, 'sponpain_mean')) || any(strcmp(varargin, 'sponpain_mean_rating')), sponpain_covariate = true; end
if any(strcmp(varargin, 'backpaindur')) || any(strcmp(varargin, 'backpain_length')), backpaindur = true; end
if any(strcmp(varargin, 'sponpain_motion')), sponpain_motion = true; end
if any(strcmp(varargin, 'bladder_motion')), bladder_motion = true; end
if any(strcmp(varargin, 'doplots')), doplots = true; end
if any(strcmp(varargin, 'noverbose')), verbose = false; end
if any(strcmp(varargin, 'standardize')), standardize = true; end
if any(strcmp(varargin, 'p_vs_hc_only')), p_vs_hc_only = true; end

if any(strcmp(varargin, 'olp_vs_uc_only')), olp_vs_uc_only = true; end

% initialize output
% ----------------------------------------------------------------------
[lme_groupbytime, festats_groupbytime, festats_3groupsbytime, festats_prt_vs_combined, lme_p_vs_hc, festats_p_vs_hc, festats_simple_time1, festats_simple_time2, mm_groupbytime, mm_p_vs_hc] = deal([]);


% covariates table
covs_matrix = [];
covs_names = {};
if any(strcmp(varargin, 'covs'))
    covs_matrix = varargin{ find(strcmp(varargin, 'covs')) + 1};
    covs_names = varargin{ find(strcmp(varargin, 'covs')) + 2};
    
    for i=1:length(covs_names)
        covariates_table.(covs_names{i}) = covs_matrix(:,i);
    end
end
    
if doplots
    warning('plotting is not yet implemented')
end

if numel(unique(covariates_table.time))==1
    warning('Cannot estimate Group by Time, data only has one timepoint')
end
    
if ~ismember('is_patient', covariates_table.Properties.VariableNames) || numel(unique(covariates_table.is_patient))==1
    warning('Cannot estimate Patient vs HC, data only has one group')
end

if isa(outcome, 'single'), outcome = double(outcome); end

X = covariates_table;
X.y = outcome; % add outcome to X

%% -- set up variables for modeling -- %%

% Centering matters when considering interactions (will influence main
% effects), so effects code. This way, main effects are uncorrelated with interaction term
% Although this has no effect on the interaction, it does have a substantial effect on the main effect estimates
% Time was [1 2] for Time 1 or 2, now -.5 (time 1) or .5 (time 2)
X.time = (X.time - 1); % assume time is passed in as 1 and 2

% Will not change P-values without any interactions with age specified.
X.age = X.age - nanmean(X.age);
X.gender = X.gender - 1.5;

% make group categorical; can do effects or reference coding below
X.grp = X.group; % save in integer format, for later
X.group = categorical(X.group, [1 2 3], {'PRT' 'PLA' 'TAU'}, 'Ordinal', false);

% back pain length: set to 0 for healthy controls; no back pain 
X.backpain_length(~X.is_patient) = 0; 

% mean-center pain within group, within time point. This way it is not
% correlated with group or time
if pain_covariate || bladder_pain_covariate
    for time=unique(X.time)'
        for grp = {'PRT' 'PLA' 'TAU'}
            wh = X.time==time & X.group == grp;
            X.pain_avg(wh) = X.pain_avg(wh) - nanmean(X.pain_avg(wh));
            X.bladder_pain_mean(wh) = X.bladder_pain_mean(wh) - nanmean(X.bladder_pain_mean(wh));
        end
    end

    % and mean center pain for the healthies
    X.pain_avg(~X.is_patient) = X.pain_avg(~X.is_patient) - nanmean(X.pain_avg(~X.is_patient));
    X.bladder_pain_mean(~X.is_patient) = X.bladder_pain_mean(~X.is_patient) - nanmean(X.bladder_pain_mean(~X.is_patient));
end


% standardize? if so, just standardize any var you might use, except for
% the categorical ones
if standardize
    X.age = scale(X.age);
    X.backpain_length = scale(X.backpain_length);
    X.pain_avg = scale(X.pain_avg);
    if isfield(X, 'n_spike_regs_sponpain'), X.n_spike_regs_sponpain = scale(X.n_spike_regs_sponpain); end
    if isfield(X, 'n_spike_regs'), X.n_spike_regs = scale(X.n_spike_regs); end
    X.y = scale(X.y);
end

%% create a dataset without placebo and without PRT

X_nopla = X;
X_nopla.group = categorical(X.grp, [1 3], {'PRT' 'TAU'}, 'Ordinal', false);

X_noprt = X;
%X_noprt = X(X.grp ~= 1,:);
%X_noprt.group = categorical(X_noprt.group, [2 3], {'PLA' 'TAU'}, 'Ordinal', false);
X_noprt.group = categorical(X.grp, [2 3], {'PLA' 'TAU'}, 'Ordinal', false);

%% -- run models of interest 

if olp_vs_uc_only

    % Group by time model: OLP vs TAU, no PRT
    % ----------------------------------------------
    % notes:    A*B         A + B + A:B, so time and group are
    % included
    % **********PLAY WITH THIS LINE OF CODE***************
    % analysis_str = 'y ~ 1 + age + gender + time*group + (1|id)';
    analysis_str = 'y ~ 1 + age + gender*time + gender*group + time*group + (1|id)';

    analysis_str = add_covariates(analysis_str);
    lme_groupbytime = fitlme(X_noprt, analysis_str, 'Exclude', X_noprt.is_patient == false, 'DummyVarCoding', 'reference', 'FitMethod', 'REML');
    [b, bnames, festats_groupbytime] = fixedEffects(lme_groupbytime, 'DFMethod', 'Satterthwaite');
    if verbose
        print_header('Group by time, OLP vs TAU, no PRT, gender moderation')
        festats_groupbytime
    end

%     analysis_str = 'y ~ 1 + gender + age*time*group + (1|id)';
% 
%     analysis_str = add_covariates(analysis_str);
%     lme_groupbytime = fitlme(X_noprt, analysis_str, 'Exclude', X_noprt.is_patient == false, 'DummyVarCoding', 'reference', 'FitMethod', 'REML');
%     [b, bnames, festats_groupbytime] = fixedEffects(lme_groupbytime, 'DFMethod', 'Satterthwaite');
%     if verbose
%         print_header('Group by time, OLP vs TAU, no PRT, age moderation')
%         festats_groupbytime
%     end
% 
%     analysis_str = 'y ~ 1 + gender + age*time + age*group + time*group + (1|id)';
% 
%     analysis_str = add_covariates(analysis_str);
%     lme_groupbytime = fitlme(X_noprt, analysis_str, 'Exclude', X_noprt.is_patient == false, 'DummyVarCoding', 'reference', 'FitMethod', 'REML');
%     [b, bnames, festats_groupbytime] = fixedEffects(lme_groupbytime, 'DFMethod', 'Satterthwaite');
%     if verbose
%         print_header('Group by time, OLP vs TAU, no PRT, age moderation decomposed')
%         festats_groupbytime
%     end


    % extract marginal means
    [margmeans, out_dat, levelnames] = table2marginalmeans(X_noprt, 'y', 'vars', {'group', 'time'}, 'exclude', X.is_patient == false | isundefined(X.group));
    mm_groupbytime.margmeans = margmeans;
    mm_groupbytime.out_dat = out_dat;
    mm_groupbytime.levelnames = levelnames;

end

% Group by time model: OLP vs TAU, no PRT
% ----------------------------------------------
% notes:    A*B         A + B + A:B, so time and group are included

% with gender interactions
analysis_str = 'y ~ 1 + age + time*group*gender + (1|id)';

% with gender as covariate
analysis_str = 'y ~ 1 + age + gender + time*group + (1|id)';

analysis_str = add_covariates(analysis_str);
lme_groupbytime_OLP_vs_UC = fitlme(X_noprt, analysis_str, 'Exclude', X_noprt.is_patient == false, 'DummyVarCoding', 'reference', 'FitMethod', 'REML');
[b, bnames, festats_OLP_vs_UC] = fixedEffects(lme_groupbytime_OLP_vs_UC, 'DFMethod', 'Satterthwaite');
if verbose
    print_header('Group by time, OLP vs TAU, no PRT')
    festats_OLP_vs_UC
end


if ~p_vs_hc_only & ~olp_vs_uc_only

    % ----------------------------------------------
    % Simple effectss Time 1
    % ----------------------------------------------
    analysis_str = 'y ~ 1 + age + gender + group';
    analysis_str = add_covariates(analysis_str);
    lme = fitlme(X, analysis_str, 'Exclude', X.time > 0, 'DummyVarCoding', 'reference', 'FitMethod', 'REML');
    % b, bnames, festats: Fixed effect stats using Satterthwaite df
    % confirmed: F in anova = t^2 for factors with 2 levels (ok), same P-values
    [b, bnames, festats_simple_time1] = fixedEffects(lme,'DFMethod', 'Satterthwaite');
    if verbose
        disp('Time 1 simple effects')
        disp(festats_simple_time1)
    end

    % Simple effects Time 2
    % ----------------------------------------------
    lme = fitlme(X, analysis_str, 'Exclude', X.time < 0, 'DummyVarCoding', 'effects', 'FitMethod', 'REML');
    [b, bnames, festats_simple_time2] = fixedEffects(lme,'DFMethod', 'Satterthwaite');
    if verbose
        disp('Time 2 simple effects')
        festats_simple_time2
    end

    % Group by time model: PRT vs TAU, no PLA
    % ----------------------------------------------
    % notes:    A*B         A + B + A:B, so time and group are included
    analysis_str = 'y ~ 1 + age + gender + time*group + (1|id)';
    analysis_str = add_covariates(analysis_str);
    lme_groupbytime = fitlme(X_nopla, analysis_str, 'Exclude', X_nopla.is_patient == false, 'DummyVarCoding', 'reference', 'FitMethod', 'REML');
    [b, bnames, festats_groupbytime] = fixedEffects(lme_groupbytime, 'DFMethod', 'Satterthwaite');
    if verbose
        print_header('Group by time, PRT vs TAU (no PLA)')
        festats_groupbytime
    end

    % Group by time model: OLP vs TAU, no PRT
    % ----------------------------------------------
    % notes:    A*B         A + B + A:B, so time and group are
    % included
    % **********PLAY WITH THIS LINE OF CODE***************
    % analysis_str = 'y ~ 1 + age + gender + time*group + (1|id)';
    analysis_str = 'y ~ 1 + age + gender*time + gender*group + time*group + (1|id)';

    analysis_str = add_covariates(analysis_str);
    lme_groupbytime = fitlme(X_noprt, analysis_str, 'Exclude', X_noprt.is_patient == false, 'DummyVarCoding', 'reference', 'FitMethod', 'REML');
    [b, bnames, festats_groupbytime] = fixedEffects(lme_groupbytime, 'DFMethod', 'Satterthwaite');
    if verbose
        print_header('Group by time, OLP vs TAU, no PRT')
        festats_groupbytime
    end


    % Group by time model: all 3 groups
    % ----------------------------------------------
    % notes:    A*B         A + B + A:B, so time and group are included
    analysis_str = 'y ~ 1 + age + gender + time*group + (1|id)';
    analysis_str = add_covariates(analysis_str);
    lme_groupbytime = fitlme(X, analysis_str, 'Exclude', X.is_patient == false, 'DummyVarCoding', 'reference', 'FitMethod', 'REML');
    [b, bnames, festats_3groupsbytime] = fixedEffects(lme_groupbytime, 'DFMethod', 'Satterthwaite');
    if verbose
        print_header('Group by time, all 3 groups')
        festats_3groupsbytime
    end


    % Group by time model: PRT vs combined controls
    % ----------------------------------------------
    % notes:    A*B         A + B + A:B, so time and group are included
    X.prt_vs_combined_controls = X.group=='PRT';
    analysis_str = 'y ~ 1 + age + gender + time*prt_vs_combined_controls + (1|id)';
    analysis_str = add_covariates(analysis_str);
    lme_groupbytime = fitlme(X, analysis_str, 'Exclude', X.is_patient == false, 'DummyVarCoding', 'reference', 'FitMethod', 'REML');
    [b, bnames, festats_prt_vs_combined] = fixedEffects(lme_groupbytime, 'DFMethod', 'Satterthwaite');
    if verbose
        print_header('Group by time, PRT vs combined controls')
        festats_prt_vs_combined
    end


    % extract marginal means
    [margmeans, out_dat, levelnames] = table2marginalmeans(X, 'y', 'vars', {'group', 'time'}, 'exclude', X.is_patient == false | isundefined(X.group));
    mm_groupbytime.margmeans = margmeans;
    mm_groupbytime.out_dat = out_dat;
    mm_groupbytime.levelnames = levelnames;

end % end ~p_vs_hc_only & ~olp_vs_uc_only

if numel(unique(X.is_patient))>1

    % Patient vs. healthy control, Time 1 only
    % ----------------------------------------------
    %X.is_patient = categorical(double(X.is_patient), 0:1, {'ZControl' 'Patient'},'Ordinal', false);
    X.is_healthy_cntrl = ~X.is_patient;
    analysis_str = 'y ~ 1 + age + gender + is_healthy_cntrl';
    analysis_str = add_covariates(analysis_str);
    % note: time has been recoded here, so time == 1 is time 2
    lme_p_vs_hc = fitlme(X, analysis_str, 'Exclude', X.time > 0, 'DummyVarCoding', 'reference', 'FitMethod', 'REML');
    [b, bnames, festats_p_vs_hc] = fixedEffects(lme_p_vs_hc, 'DFMethod', 'Satterthwaite');
    if verbose
        disp('Patients vs. Healthies')
        festats_p_vs_hc
    end

    % extract marginal means
    [margmeans, out_dat, levelnames] = table2marginalmeans(X, 'y', 'vars', {'is_patient'}, 'exclude', X.time > 0);
    mm_p_vs_hc.margmeans = margmeans;
    mm_p_vs_hc.out_dat = out_dat;
    mm_p_vs_hc.levelnames = levelnames;

end

% DESIGN MATRIX
% xnames = strrep(lme.CoefficientNames, '(Intercept)', 'Intercept');
% X_matrix = array2table(designMatrix(lme), 'VariableNames', strrep(xnames, ':', '_x_'));
% % display results
% X_matrix(1:5, :)



% if doplots
    % s1backT1 = get_wh_image(s1back_conn_masked, s1back_conn_masked.metadata_table.time==1);% & s1back_conn_masked.metadata_table.n_spike_regs_sponpain < 600);
    % figure; violinplot({s1backT1.metadata_table.pcl_somatosensory(s1backT1.metadata_table.is_patient), s1backT1.metadata_table.pcl_somatosensory(~s1backT1.metadata_table.is_patient)})
    % legend off
    % mes(s1backT1.metadata_table.pcl_somatosensory(s1backT1.metadata_table.is_patient), s1backT1.metadata_table.pcl_somatosensory(~s1backT1.metadata_table.is_patient), 'hedgesg')
% end



% *** return marginal means

    function analysis_str = add_covariates(analysis_str)
        
        % hack way to remove REs and add it back at end
        if analysis_str(end) == ')'
            [s1, s2] = split(analysis_str,'+ (');
            s2 = ['+ (' s1{2}];
            analysis_str = s1{1};
        else
            s2 = '';
        end
        
        if pain_covariate
            analysis_str = [analysis_str ' + pain_avg '];
        end
                
        if bladder_pain_covariate
            analysis_str = [analysis_str ' + bladder_pain_mean '];
        end
        
        if sponpain_covariate
            analysis_str = [analysis_str ' + spon_pain_rating_mean '];
        end
        
        if sponpain_motion
            analysis_str = [analysis_str ' + n_spike_regs_sponpain '];
        end
        
        if bladder_motion
            analysis_str = [analysis_str ' + n_spike_regs '];
        end
        
        if backpaindur
            analysis_str = [analysis_str ' + backpain_length '];
        end
        
        % add other covariates
        for c=1:length(covs_names)
            analysis_str = [analysis_str ' + ' covs_names{c} ' '];
        end
        
        analysis_str = [analysis_str s2];
    end
end


