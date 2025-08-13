rm(list=ls())
library(mgcv)
library(tidyverse)
library(lme4)
library(here)
library(ggplot2)
library(qgam)
## directory where data are saved locally (github will not have data uploaded)
path = here("data")
## read in the data
df <- read_csv(here(path,"behaviour_table_all.csv"))
## label study arm, assessment, and intensity
## create numeric variables for 
## Note: the code creating exposure sequence assumes that the first observations 
##       record: numeric counter for the exposure sequence 1,..., 20. Combines pre- and post-intervention in the counting
##       record_2: numeric counter for exposure sequence 1,...,10 for each participant-assessment
df <- 
    df %>% 
    mutate(group_fac = factor(Group, levels=c(1,2,3), labels=c("PRT","placebo","usual care")),
           time_fac = factor(Time, levels=c(0,1), labels=c("pre","post")),
           intensity_fac = factor(Intensity,levels=c(1,2), labels=c("low","high"))) %>% 
    group_by(Subject, time_fac) %>% 
    mutate(record = (1:n()) + Time*10,
           record_2 = 1:n()) %>% 
    ungroup()

## exploratory boxplots
df %>% 
    ggplot() + 
    geom_boxplot(aes(x=time_fac, y=Rating, fill=group_fac)) + facet_grid(~intensity_fac)


## fit a "naive" model with nested random intercepts/slopes for participants 
fit_naive <- lmer(Rating ~ group_fac*time_fac*intensity_fac + 
                      (1|Subject) + (record_2|Subject:time_fac), data=df)
summary(fit_naive)
## get residuals to visualized unaccoiunted for correlation present in the data
## given the random effects structure in our "naive" model
resids <- resid(fit_naive)
## add residuals back to the outcome data
df_resid <- 
    df %>% 
    mutate(resids = resids)
## calculate the MoM estimator for covariance/correlation of residuals
cov_mat <- 
    df_resid %>%
    dplyr::select(Subject, record, resids) %>%
    arrange(record) %>%
    pivot_wider(names_from=record, values_from = resids) %>%
    dplyr::select(-Subject) %>%
    as.matrix() %>%
    cov(use="pairwise.complete.obs")

## plot heatmats of residual correlation
par(mfrow=c(1,2))
fields::image.plot(1:20, 1:20, (cov_mat), xlab="Replicate", ylab="Replicate", main="Covariance Plot") 
contour(1:20,1:20, (cov_mat), levels=c(0), add=T)
fields::image.plot(1:20, 1:20, cov2cor(cov_mat), xlab="Replicate", ylab="Replicate", main="Correlation Plot") 
contour(1:20,1:20, cov2cor(cov_mat), levels=c(-0.1,0,0.1), add=T, method="edge")


## create some new factor variables needed to fit the model in mgcv
df <- 
    df %>% 
    mutate(Subject_fac = factor(Subject),
           Subject_time = (paste0(Subject,"_",time_fac)),
           Subject_time_fac = factor(paste0(Subject,"_",time_fac)))



## primary analysis: remove outliers
fit_gaulss <- gam(list(Rating ~ group_fac*time_fac + intensity_fac + record_2 + s(Subject_fac, bs="re") + 
                             s(Subject_time_fac, bs="re") + s(Subject_time_fac, by=record_2, bs="re"),
                         ~ time_fac + s(Subject_fac, bs="re")
                         ), data=df %>% filter(!Subject %in%  c(1231, 1029)), family=gaulss())
summary(fit_gaulss)


## sensitivity analysis: fit our proposed model to the whole dataset
fit_gaulss_all_data <- gam(list(Rating ~ group_fac*time_fac + intensity_fac + record_2 + s(Subject_fac, bs="re") + 
                             s(Subject_time_fac, bs="re") + s(Subject_time_fac, by=record_2, bs="re"),
                         ~ time_fac + s(Subject_fac, bs="re")
                         ), 
                    data=df, family=gaulss())
summary(fit_gaulss_all_data)


## sensitivity analysis: quantile additive mixed model 
fit_qgam <- qgam(Rating ~ group_fac:time_fac+intensity_fac + s(record_2) + 
                             s(Subject_fac, bs="re") + 
                             s(Subject_time_fac, bs="re") + 
                             s(Subject_fac, by=record_2, bs="re") + s(Subject_time_fac, by=record_2, bs="re"),
                             data=df,qu=0.5)
summary(fit_qgam)

