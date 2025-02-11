set.seed(18913)
library(tidyverse)
## package for standard LMMs
library(lme4)
## package for fitting the GAULSS model
library(mgcv)

## path where data are stored
path <- file.path(".","data")
## path for saving results/figures
fig_out_path <- file.path(".","figs")
results_out_path <- file.path(".","figs")


## kead tge data
data <- read_csv(file.path(path, "five_yr_long_final.csv"))

outcome_vars <- c("pain_avg","promis_dep","promis_sleep","promis_anger","promis_anxiety",
                  "panas_pa","panas_na","tsk11","sopa_control",
                  "bpi_intensity","bpi_interference","odi","pcs")
## check distribution of variance across groups and time points for each outcome
data <- 
  data %>% 
  mutate(time_fac = factor(time, levels=c(-1,1), labels=c("baseline","5_years")),
         group_fac = factor(group, levels=c(1,2,3), labels=c("PRT","placebo","usual care")),
         gender_fac = factor(gender, levels=c(1,2), labels=c("male","female"))) %>%
  filter(!is.na(group)) ## remove participants without a group assignment

## create indicator for whether a participant has 5 year data
data <- 
  data %>% 
  group_by(id) %>% 
  mutate(has_5year = as.numeric(n() == 2))

## transform to long format for plotting, with 
## repeated outcomes per person (rowwise)
data_long <- 
  data %>% 
  pivot_longer(one_of(outcome_vars), names_to="outcome_variable", 
               values_to = "Y")
  
## plot distribution of sd of outcome by group, outcome, and time point 
## shows heterogeneous error
data_long %>% 
  group_by(time_fac, group_fac, outcome_variable) %>%
  summarize(mn=mean(Y),
            sd=sd(Y)) %>%
  ggplot() + 
  geom_point(aes(x=group_fac, y=sd, color=time_fac)) + 
  facet_wrap(~outcome_variable, ncol=4, scales="free_y")

  

## plot distribution of outcomes by group and time point
## can again see heterogeneous error
data_long %>% 
    ggplot() + geom_boxplot(aes(x=time_fac, y=Y,color=group_fac)) + 
    facet_wrap(~outcome_variable,ncol=4, scales="free_y")
  
## create factor variable for participant identifier (this is required by mgcv to fit the correct model we want)
data <- 
  data %>% 
  mutate(id_fac = factor(id))

## example of how to fit this model using one outcome
lmm_pain_avg <- gam(list(promis_dep ~ group_fac*time_fac + gender + age + s(id_fac,bs="re"),
                         ~ time_fac * group_fac), data=data, family=gaulss())

## specify the formula for the model
## we model outcomes (baseline and 5 year) as multivariate guassian
## with a random intercept per subject and gender, age, and group x time effects
right_fmla_mn <- 'group_fac*time_fac + gender_fac + age + s(id_fac,bs="re")'

## set up empty containers to store results
## coefs_ls: list of coefficients to store results from a standard LMM and GAULSS results
coefs_ls <- grp_mns_ls <- vector(mode="list",length=length(outcome_vars))
## empty array for storing bootstrap results
nboot <- 1000
cond_mn_arr <- 
        array(NA, 
              dim=c(length(outcome_vars), nboot, 3),
              dimnames=list("outcome"=outcome_vars,
                            "bootstrap"=1:nboot,
                            "group"=c("PRT","usual care","placebo")))
## empty array for storing Cohen's D estimates
d_mat <- data.frame("outcome"=outcome_vars, "PRT_vs_UC_cohens_d"=NA, "PRT_vs_placebo_cohens_d"=NA)

## set coefficient names (used to extract results from each model fit)
coef_names_output <- c("group_facplacebo:time_fac5_years","group_facusual care:time_fac5_years")

## set up data frame for obtaining model predictions 
## at the desired levels 
df_grp_mns <- expand.grid("group_fac" = levels(data$group_fac),
                          "time_fac" = levels(data$time_fac),
                          "gender_fac" = "female",
                          "age" = mean(data$age[data$time_fac == "baseline"]),
                          "id_fac" = data$id_fac[1])

## get a table for n's in each group at 5 years
## (used for variance estimates for conditional sample mean)
tab_n <- 
        data %>%
        filter(time_fac == "5_years") %>%
        group_by(group_fac) %>% 
        summarize(n_grp_time = n())

## set up a progress bar (pb) and index (inx) for keeping track of progreess
pb <- txtProgressBar(0, length(outcome_vars)*nboot, style=3)        
inx <- 1

## do the loop
for(p in seq_along(outcome_vars)){
  ## loop over outcomes
  outcome_p <- outcome_vars[[p]]
  ####### Model fitting
  ## set up the formula for the mean model for the current outcome
  fmla_p <- as.formula(paste0(outcome_p, "~", right_fmla_mn))
  ## fit the GAULSS model
  fit_gaulss_p <- gam(list(fmla_p,
                         ~ time_fac + group_fac), 
                      ## note this "b" constant needs to be modified from default to avoid convergence
                      ## issues
                    data=data, family=gaulss(b=0.00001))
  ## fit a standard LMM
  fit_lmm_p <- gam(fmla_p,
                   data=data)

  ####### Task 1: Group contrasts
  ## collect results for group contrasts, coefficients and p-values for each model
  df_p <- data.frame(contrast=rep(c("placebo_vs_PRT","UC_vs_PRT"),2),
                     model=rep(c("location-scale LMM","standard LMM"),each=2),
                     beta = c(fit_gaulss_p$coefficients[coef_names_output],
                              fit_lmm_p$coefficients[coef_names_output]),
                     p_val = c(summary(fit_gaulss_p)$p.table[coef_names_output,"Pr(>|z|)"],
                               summary(fit_lmm_p)$p.table[coef_names_output,"Pr(>|t|)"]
                                )
                     )
  ## save the results in the list container
  df_p$outcome <- outcome_p
  coefs_ls[[p]] <- df_p
  
  ####### Task 2: Get conditional mean estimates and corresponding uncertainty estimates
  #######         Estimand: E[Y_5 years| Y_baseline = x, other stuff, b_i]
  #######         Estimator: conditional mean with sample size of each group
  ## get group means/variances for calculating conditional mean and variance
  ## associated with distribution of 5 year outcome given baseline outcome
  preds_var <- 
          data.frame(df_grp_mns, 
                     "s2_e" = (1/predict(fit_gaulss_p, newdata=df_grp_mns, type='response')[,2]^2)
          ) %>% 
          pivot_wider(names_from="time_fac",values_from = "s2_e") %>% 
          mutate("s2_b"=(1/fit_gaulss_p$sp)^2) %>% 
          rename("s2_e_baseline" = "baseline",
                 "s2_e_5_years" = "5_years") %>% 
          mutate("rho_Y0_Y5" = s2_b/(sqrt(s2_e_baseline + s2_b)*sqrt(s2_e_5_years + s2_b)))
  preds_mns <- 
          data.frame(df_grp_mns, 
                     "mu_Y" = predict(fit_gaulss_p, newdata=df_grp_mns, type='response')[,1]) %>% 
          pivot_wider(names_from="time_fac",values_from = "mu_Y") %>% 
          rename("mu_Y_baseline" = "baseline",
                 "mu_Y_5_years" = "5_years")
  ## combine mean and variance results to get the ultimate quantity
  preds <- 
          inner_join(preds_var, preds_mns,by=c("group_fac","gender_fac","age","id_fac")) %>% 
          left_join(tab_n, by="group_fac") %>%
          mutate(outcome = outcome_p,
                 EY5_Y0 = mu_Y_5_years + sqrt(s2_e_5_years/s2_e_baseline)*rho_Y0_Y5*(mean(data[[outcome_p]][data$time_fac=="baseline"],na.rm=T) - mu_Y_baseline),
                 VY5_Y0 = (1-rho_Y0_Y5^2)*s2_e_5_years,
                 LB_Y5_Y0 = EY5_Y0 - qnorm(0.975)*sqrt(VY5_Y0/n_grp_time),
                 UB_Y5_Y0 = EY5_Y0 + qnorm(0.975)*sqrt(VY5_Y0/n_grp_time))
  grp_mns_ls[[p]] <- preds
  
  ## get cohen's d
  d_mat[p,"PRT_vs_UC_cohens_d"] <- (preds$mu_Y_5_years[preds$group_fac == "PRT"] - preds$mu_Y_5_years[preds$group_fac == "usual care"])/
          sqrt(mean(c(preds$s2_e_5_years[preds$group_fac == "PRT"] + preds$s2_b[preds$group_fac == "PRT"], preds$s2_e_5_years[preds$group_fac == "usual care"] + + preds$s2_b[preds$group_fac == "usual care"])))
  d_mat[p,"PRT_vs_placebo_cohens_d"] <- (preds$mu_Y_5_years[preds$group_fac == "PRT"] - preds$mu_Y_5_years[preds$group_fac == "placebo"])/
          sqrt(mean(c(preds$s2_e_5_years[preds$group_fac == "PRT"] + preds$s2_b[preds$group_fac == "PRT"], preds$s2_e_5_years[preds$group_fac == "placebo"]+ preds$s2_b[preds$group_fac == "placebo"])))

  ## get bootstrap estimates of uncertainty (this is a bit lazy coding -- copy and paste code above essentially, sorry!)
  uid <- unique(data$id)
  nid <- length(uid)
  for(b in 1:nboot){
          ## resample participants
          id_b <- sample(uid, size=nid, replace=T)
          ## create new data frame with new participant identifiers
          id_new_b <- data.frame("id_new"=1:nid, "id"=id_b)
          ## merge to get resampled data set
          data_b <-
                  data %>%
                  inner_join(id_new_b, by="id") %>%
                  mutate(id = id_new,
                         id_fac = factor(id))
          ## fit the model
          fit_gaulss_p_b <- try(gam(list(fmla_p,
                                   ~ time_fac + group_fac),
                              data=data_b, family=gaulss(b=0.00001)))
          if(inherits(fit_gaulss_p_b, "try-error")){
                  next
          }
          ## get conditional mean estimates for each group
          preds_var_b <-
                  data.frame(df_grp_mns,
                             "s2_e" = (1/predict(fit_gaulss_p_b, newdata=df_grp_mns, type='response')[,2]^2)
                  ) %>%
                  pivot_wider(names_from="time_fac",values_from = "s2_e") %>%
                  mutate("s2_b"=(1/fit_gaulss_p_b$sp)^2) %>%
                  rename("s2_e_baseline" = "baseline",
                         "s2_e_5_years" = "5_years") %>%
                  mutate("rho_Y0_Y5" = s2_b/(sqrt(s2_e_baseline + s2_b)*sqrt(s2_e_5_years + s2_b)))
          preds_mns_b <-
                  data.frame(df_grp_mns,
                             "mu_Y" = predict(fit_gaulss_p_b, newdata=df_grp_mns, type='response')[,1]) %>%
                  pivot_wider(names_from="time_fac",values_from = "mu_Y") %>%
                  rename("mu_Y_baseline" = "baseline",
                         "mu_Y_5_years" = "5_years")
          ## combine mean and variance results to get the ultimate quantity
          preds_b <-
                  inner_join(preds_var_b, preds_mns_b, by=c("group_fac","gender_fac","age","id_fac")) %>%
                  left_join(tab_n, by="group_fac") %>%
                  mutate(outcome = outcome_p,
                         EY5_Y0 = mu_Y_5_years + sqrt(s2_e_5_years/s2_e_baseline)*rho_Y0_Y5*(mean(data[[outcome_p]][data$time_fac=="baseline"],na.rm=T) - mu_Y_baseline)
                         )
          cond_mn_arr[p,b,"PRT"] <- preds_b$EY5_Y0[preds_b$group_fac=="PRT"]
          cond_mn_arr[p,b,"placebo"] <- preds_b$EY5_Y0[preds_b$group_fac=="placebo"]
          cond_mn_arr[p,b,"usual care"] <- preds_b$EY5_Y0[preds_b$group_fac=="usual care"]

          inx <- inx + 1
          setTxtProgressBar(pb, inx)
  }
  
}

## get coefficient estimates and conditional mean estimates across simulation scenarios and models
coefs_df <- bind_rows(coefs_ls)
grp_mns_df <- bind_rows(grp_mns_ls)
## combine bootstrap results across simulation scenarios
cond_mns_df <- as.data.frame.table(cond_mn_arr, responseName="cond_mn", stringsAsFactors = F)

## create matrix with Cohen's D results
d_mat <- 
        d_mat %>% 
        mutate_if(is.numeric, function(x) sprintf("%5.2f", round(abs(x), 3)))


## create plots of the results
plt_coef <- 
  coefs_df %>% 
  ggplot() + 
  geom_point(aes(x=contrast, y=beta, color=model),alpha=0.7) + 
  facet_wrap(~outcome, scales="free_y") + xlab("contrast") + ylab(expression(hat(beta))) + 
  theme_classic() + ggtitle("Coefficient Estimates")

plt_p <- 
  coefs_df %>% 
  ggplot() + 
  geom_point(aes(x=contrast, y=p_val, color=model),alpha=0.7) + 
  facet_wrap(~outcome, scales="free_y") + xlab("contrast") + ylab("p-value") + 
  theme_classic() + ggtitle("Statistical Significance") + 
  geom_hline(yintercept=0.05, col='grey',lty=2)

plt_var <- 
  data_long %>% 
  group_by(time_fac, group_fac, outcome_variable) %>%
  summarize(mn=mean(Y,na.rm=T),
            sd=sd(Y,na.rm=T)) %>%
  ggplot() + 
  geom_point(aes(x=group_fac, y=sd, color=time_fac)) + 
  facet_wrap(~outcome_variable, ncol=4, scales="free_y") + 
  theme_classic() + xlab("") + ylab("Standard Deviation of Outcome") + 
  ggtitle("Variability of Outcomes by Study Arm and Time")

## create table 1
table1 <- 
        grp_mns_df %>%
        mutate(mean_sd_conditional_mean = sprintf("%5.2f (%5.2f,%5.2f)", EY5_Y0, LB_Y5_Y0, UB_Y5_Y0)) %>% 
        dplyr::select(group_fac, outcome, mean_sd_conditional_mean) %>% 
        pivot_wider(names_from=group_fac, values_from =mean_sd_conditional_mean)
                

## write the results
write.csv(table1, file.path(results_out_path, "table1.csv"))
write.csv(d_mat, file.path(results_out_path, "cohens_d.csv"))

## save the figures
ggsave(file.path(fig_out_path, "coefficient_estimates_modelling_choice.jpeg"), plt_coef, height=5, width=11)
ggsave(file.path(fig_out_path, "p_values_modelling_choice.jpeg"), plt_p, height=5, width=11)
ggsave(file.path(fig_out_path, "plt_var.jpeg"), plt_var, height=5, width=11)


