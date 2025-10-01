%% set up
clear all

basedir =  '/Users/yoni/Repositories/Boulder_backpain_5_year_fu';
fiveyr = readtable(fullfile(basedir,'data',"5YearsFollowUp deidentified reIDed scale scores only.csv"));

wh1 = fiveyr.group==1;
wh2 = fiveyr.group==2;
wh3 = fiveyr.group==3;

%% follow up rates

obs = histc(fiveyr.group, 1:3)'

% Expected frequencies for an equal distribution
expected = height(fiveyr) / 3 * ones(1, 3);  % Since there are three categories

%% sig group diffs in follow up rates?
[h,p,stats] = chi2gof(fiveyr.group, 'expected', expected)

%% manual calculation of chi2

% Calculate the Chi-Squared statistic
chi2stat = sum((obs - expected).^2 ./ expected)

% Degrees of freedom (number of categories - 1)
df = height(fiveyr) - 1;

% Calculate p-value using Chi-Squared cumulative distribution function
p_value = 1 - chi2cdf(chi2stat, df)

%% age and sex of sample

fiveyr_long = readtable(fullfile(basedir,'data',"five_yr_long_final.csv"));
fiveyr_long_T1 = fiveyr_long(fiveyr_long.time==1,:);

mean(fiveyr_long_T1.age)
std(fiveyr_long_T1.age)

histc(fiveyr_long_T1.gender, 0:3)