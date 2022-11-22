####################################################################################### 
# R Code for Heart Transplant Induction Project Using UNOS - Summer 2022
####################################################################################### 

## Load required R packages
library(haven)                  # Package to read in SAS files (sas7bdat)
library(Hmisc)                  # Package to describe data
library(gtsummary)              # Package to create Table 1
library(survival)               # Package for survival analyses
library(survminer)              # Package to create KM curves
library(rms)                    # Package to run regression analyses
library(arsenal)                # Package to create comparison tables
library(naniar)                 # Package to evaluate missing data patterns
library(mice)                   # Package to impute missing data
library(tidyverse)              # Package for data view & manipulation

## Load in SAS file from computer location and read into R
data_immuno_base <- read_sas('PATH')

## Look at number of patients & variables (n=7,845 & 612 variables)
dim(data_immuno_base)
colnames(data_immuno_base)

####################################################################################### 
# R Code for Baseline Characteristics
####################################################################################### 

## Select desired variables
vars_bl = c("AGE", "REC_MALE", "REC_RACE", "BMI_CALC", "REC_ABO", "HF_ETIOL", "REC_DM", 
            "REC_CIG", "CREAT_TRR", "REC_CLCR2", "REC_DIAL", "REC_ICD",
            "INOTROPES_TRR", "IABP_TRR", "ECMO_TRR", "REC_VAD", "VENTILATOR_TRR",
            "CPRA", "CPRA_PEAK", "REC_PRA", "MALIG", "INFECT_IV_DRUG_TRR",
            "ISCHTIME_MIN", "AGE_DON", "DON_ABO", "BMI_DON_CALC", "DON_HTN", "DON_DM", 
            "FU_time_yrs_TX", "GROUP", "TX_YEAR", "REC_CMV")

## Limit data to desired variables and look at number of patients & variables
data_immuno_bl = data_immuno_base[, vars_bl]

# Use the 'transform' function to change the type of each feature
data_immuno_bl <- transform(REC_MALE=as.factor(REC_MALE),
                         REC_RACE=as.factor(REC_RACE), HF_ETIOL=as.factor(HF_ETIOL),
                         REC_ABO=as.factor(REC_ABO), REC_DIAL=as.factor(REC_DIAL),
                         REC_DM=as.factor(REC_DM), REC_CIG=as.factor(REC_CIG),
                         REC_ICD=as.factor(REC_ICD),
                         ECMO_TRR=as.factor(ECMO_TRR), REC_VAD=as.factor(REC_VAD),
                         VENTILATOR_TRR=as.factor(VENTILATOR_TRR),
                         DON_ABO=as.factor(DON_ABO), DON_HTN=as.factor(DON_HTN),
                         DON_DM=as.factor(DON_DM), INOTROPES_TRR=as.factor(INOTROPES_TRR),
                         IABP_TRR=as.factor(IABP_TRR), 
                         GROUP=as.factor(GROUP), TX_YEAR=as.factor(TX_YEAR),
                         REC_CMV=as.factor(REC_CMV), REC_PRA=as.factor(REC_PRA),
                         MALIG=as.factor(MALIG),
                         INFECT_IV_DRUG_TRR=as.factor(INFECT_IV_DRUG_TRR),
                         data_immuno_bl)

## Look at selected data files
dim(data_immuno_bl)        # n = 7,845, 32 variables
colnames(data_immuno_bl)

describe(data_immuno_bl)

# Create baseline characteristics table using 'gtsummary' package
Table_Immuno <- tbl_summary(data_immuno_bl,
                          by = GROUP,
                          missing = "ifany"
) %>%
  add_p(simulate.p.value=TRUE)
Table_Immuno

####################################################################################### 
# R Code for Early Mortality After Successful Transplant
####################################################################################### 

## Load in SAS file from computer location and read into R
data_immuno_base <- read_sas('C:/Users/wbake/Desktop/UNOS/THORACIC/IMMUNO_DATA3.sas7bdat')

## Look at number of patients & variables (n=7,845 & 612 variables)
dim(data_immuno_base)

## Select desired variables
vars_immuno_surv = c("GROUP", "REC_ATGAM", "REC_CAMPATH", "REC_SIMULECT",
                "REC_THYMO", "REC_ZENAPAX", "TT_Out_Death_30d", "Out_Death_30d", "TT_Out_Death_1", 
                "Out_Death_1", "OUT_GRAFT_1", "TT_OUT_GRAFT_1", "OUT_TREAT_REJ_1",
                "OUT_HOSP_INF")

## Limit data to desired variables and look at number of patients & variables
data_immuno_surv = data_immuno_base[, vars_immuno_surv]

## Look at selected data files
dim(data_immuno_surv)        # n = 7,845, 14 variables

data_immuno_surv <- transform(GROUP=as.factor(GROUP), REC_ATGAM=as.factor(REC_ATGAM),
                         REC_CAMPATH=as.factor(REC_CAMPATH),
                         REC_SIMULECT=as.factor(REC_SIMULECT),
                         REC_THYMO=as.factor(REC_THYMO),
                         REC_ZENAPAX=as.factor(REC_ZENAPAX),
                         Out_Death_30d=as.numeric(Out_Death_30d),
                         Out_Death_1=as.numeric(Out_Death_1),
                         OUT_GRAFT_1=as.numeric(OUT_GRAFT_1),
                         OUT_HOSP_INF=as.factor(OUT_HOSP_INF),
                         data_immuno_surv)

describe(data_immuno_surv$OUT_HOSP_INF)

# Shows p-values as numbers vs. scientific
options(scipen = 999) 

tb_immuno <- tableby(GROUP ~ REC_ATGAM + REC_CAMPATH + 
                       REC_SIMULECT + REC_THYMO + REC_ZENAPAX, 
              data = data_immuno_surv)
summary(tb_immuno, text = TRUE)

tb_out <- tableby(GROUP ~ as.factor(Out_Death_30d) + as.factor(Out_Death_1) + 
                    as.factor(OUT_GRAFT_1) + as.factor(OUT_TREAT_REJ_1) +
                  as.factor(OUT_HOSP_INF),
                  data = data_immuno_surv)
summary(tb_out, text = TRUE)

################ 30-Day Survival Kaplan-Meier Curve ###############

Surv_Fit_30D <- survfit(Surv(TT_Out_Death_30d, Out_Death_30d) ~ GROUP, 
                        data=data_immuno_surv)
print(Surv_Fit_30D, print.rmean = TRUE)
survdiff(Surv(data_immuno_surv$TT_Out_Death_30d, data_immuno_surv$Out_Death_30d) ~ 
           data_immuno_surv$GROUP)

ggsurvplot(Surv_Fit_30D, data = data_immuno_surv,
           conf.int = TRUE,
           xlab = "Survival Time (days)",
           ylab = "Survival Probability",
           legend = c(0.2,0.2),
           ylim = c(0.9,1),
           xlim = c(0,30),
           break.time.by = 5,
           pval = TRUE,
           legend.labs = c("New System", "Old System"),
           risk.table = TRUE,
           risk.table.height = 0.25
)

dd <- datadist(data_immuno_surv)
options(datadist='dd')

cph(Surv(TT_Out_Death_30d, Out_Death_30d) ~ GROUP,
    x=TRUE, y=TRUE, surv=TRUE, data=data_immuno_surv)

################ 1-Year Survival Kaplan-Meier Curve ###############

Surv_Fit_1Y <- survfit(Surv(TT_Out_Death_1, Out_Death_1) ~ GROUP, 
                       data=data_immuno_surv)
print(Surv_Fit_1Y, print.rmean = TRUE)
survdiff(Surv(data_immuno_surv$TT_Out_Death_1, data_immuno_surv$Out_Death_1) ~ 
           data_immuno_surv$GROUP)

ggsurvplot(Surv_Fit_1Y, data = data_immuno_surv,
           conf.int = TRUE,
           xlab = "Survival Time (days)",
           ylab = "Survival Probability",
           legend = c(0.2,0.2),
           ylim = c(0.8,1),
           xlim = c(0,365),
           break.time.by = 90,
           pval = TRUE,
           legend.labs = c("New System", "Old System"),
           risk.table = TRUE,
           risk.table.height = 0.25,
)

dd <- datadist(data_immuno_surv)
options(datadist='dd')

cph_1Y <- cph(Surv(TT_Out_Death_1, Out_Death_1) ~ GROUP,
    x=TRUE, y=TRUE, surv=TRUE, data=data_immuno_surv)
cph_1Y

################ 1-Year Graft Failure Kaplan-Meier Curve ###############

Surv_Fit_1Y <- survfit(Surv(TT_OUT_GRAFT_1, OUT_GRAFT_1) ~ GROUP, 
                       data=data_immuno_surv)
print(Surv_Fit_1Y, print.rmean = TRUE)
survdiff(Surv(data_immuno_surv$TT_OUT_GRAFT_1, data_immuno_surv$OUT_GRAFT_1) ~ 
           data_immuno_surv$GROUP)

ggsurvplot(Surv_Fit_1Y, data = data_immuno_surv,
           conf.int = TRUE,
           xlab = "Survival Time (days)",
           ylab = "Survival Probability",
           legend = c(0.2,0.2),
           ylim = c(0.9,1),
           xlim = c(0,365),
           break.time.by = 90,
           pval = TRUE,
           legend.labs = c("New System", "Old System"),
           risk.table = TRUE,
           risk.table.height = 0.25,
)

dd <- datadist(data_immuno_surv)
options(datadist='dd')

cph(Surv(TT_OUT_GRAFT_1, OUT_GRAFT_1) ~ GROUP,
    x=TRUE, y=TRUE, surv=TRUE, data=data_immuno_surv)

####################################################################################### 
# R Code for Multiple Regression Post-Transplant Outcome Analysis
####################################################################################### 

## Load in SAS file from computer location and read into R
data_immuno_base <- read_sas('C:/Users/wbake/Desktop/UNOS/THORACIC/IMMUNO_DATA3.sas7bdat')

## Look at number of patients & variables (n=7,845 & 612 variables)
dim(data_immuno_base)

## Select desired variables
vars_immuno_model = c("GROUP", "TT_Out_Death_30d", "Out_Death_30d", "TT_Out_Death_1", 
                     "Out_Death_1", "OUT_GRAFT_1", "TT_OUT_GRAFT_1", "OUT_TREAT_REJ_1",
                     "HF_ETIOL", "SEX_MATCH", "DON_HTN", "DON_DM", "REC_DM", "ISCHTIME_MIN",
                     "AGE", "AGE_DON", "REC_RACE", "CPRA_PEAK", "REC_MALIG", "REC_CMV", "CREAT_TRR",
                     "OUT_HOSP_INF")

## Limit data to desired variables and look at number of patients & variables
data_immuno_model = data_immuno_base[, vars_immuno_model]

## Look at selected data files
dim(data_immuno_model)        # n = 7,845, 22 variables

data_immuno_model <- transform(GROUP=as.factor(GROUP), HF_ETIOL=as.factor(HF_ETIOL),
                              Out_Death_30d=as.numeric(Out_Death_30d),
                              Out_Death_1=as.numeric(Out_Death_1),
                              OUT_GRAFT_1=as.numeric(OUT_GRAFT_1),
                              SEX_MATCH=as.factor(SEX_MATCH), DON_HTN=as.factor(DON_HTN),
                              DON_DM=as.factor(DON_DM), REC_DM=as.factor(REC_DM),
                              REC_RACE=as.factor(REC_RACE), REC_MALIG=as.factor(REC_MALIG),
                              REC_CMV=as.factor(REC_CMV),
                              OUT_HOSP_INF=as.factor(OUT_HOSP_INF),
                              data_immuno_model)
colnames(data_immuno_model)

###### Impute missing data using 'MICE' package using predictive mean matching #####

## Get the number of missings per variable (n and %)
miss_var_summary(data_immuno_model)

immuno_model_mice <- mice(data_immuno_model, m = 20, seed = 0621, print = F, method = "pmm")

set.seed(0621)
dd <- datadist(data_immuno_model)
options(datadist="dd")

################ Cox Regression Using Missing Data ###############

# Convert death outcome variables to integer
data_immuno_model$Out_Death_30d = as.integer(data_immuno_model$Out_Death_30d)

# Fit Cox Regression for 30-Day Mortality Using Imputed Values
Surv_30d <- Surv(time = data_immuno_model$TT_Out_Death_30d, event = data_immuno_model$Out_Death_30d)
COX_30d <- fit.mult.impute(Surv_30d ~ GROUP + SEX_MATCH + HF_ETIOL + 
                             DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                             rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                             REC_CMV + rcs(CREAT_TRR,3),
                           fitter=cph, xtrans=immuno_model_mice, data=data_immuno_model, 
                           pr=FALSE, x = TRUE, y = TRUE)

print(COX_30d, coef=FALSE)
summary(COX_30d)
anova(COX_30d)

# Fit Cox Regression for 1-Year Mortality Using Imputed Values
data_immuno_model$Out_Death_1 = as.integer(data_immuno_model$Out_Death_1)

Surv_1Y <- Surv(time = data_immuno_model$TT_Out_Death_1, 
                event = data_immuno_model$Out_Death_1)
COX_1Y <- fit.mult.impute(Surv_1Y ~ GROUP + SEX_MATCH + HF_ETIOL + 
                            DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                            rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                            REC_CMV + rcs(CREAT_TRR,3),
                          fitter=cph, xtrans=immuno_model_mice, data=data_immuno_model, 
                          pr=FALSE, x = TRUE, y = TRUE)

print(COX_1Y, coef=FALSE)
summary(COX_1Y)
anova(COX_1Y)

# Fit Cox Regression for 1-Year Graft Failure Using Imputed Values
data_immuno_model$OUT_GRAFT_1 = as.integer(data_immuno_model$OUT_GRAFT_1)

Graft_1Y <- Surv(time = data_immuno_model$TT_OUT_GRAFT_1, 
                event = data_immuno_model$OUT_GRAFT_1)
GRAFT_1Y <- fit.mult.impute(Graft_1Y ~ GROUP + SEX_MATCH + HF_ETIOL + 
                            DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                            rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                            REC_CMV + rcs(CREAT_TRR,3),
                          fitter=cph, xtrans=immuno_model_mice, data=data_immuno_model, 
                          pr=FALSE, x = TRUE, y = TRUE)

print(GRAFT_1Y, coef=FALSE)
summary(GRAFT_1Y)
anova(GRAFT_1Y)

# Fit Cox Regression for 1-Year Drug-Treated Rejection Using Imputed Values
data_immuno_model$OUT_TREAT_REJ_1 = as.integer(data_immuno_model$OUT_TREAT_REJ_1)

INFECT_1Y <- fit.mult.impute(OUT_TREAT_REJ_1 ~ GROUP + SEX_MATCH + HF_ETIOL + 
                              DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                              rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                              REC_CMV + rcs(CREAT_TRR,3),
                            fitter=lrm, xtrans=immuno_model_mice, data=data_immuno_model,
                            pr=FALSE)

print(INFECT_1Y, coef=FALSE)
summary(INFECT_1Y)
anova(INFECT_1Y)

# Fit Logistic Regression for 1-Year Hospitalization for Infection Using Imputed Values
data_immuno_model$OUT_HOSP_INF = as.integer(data_immuno_model$OUT_HOSP_INF)

HOSP_1Y <- fit.mult.impute(OUT_HOSP_INF ~ GROUP + SEX_MATCH + HF_ETIOL + 
                               DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                               rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                               REC_CMV + rcs(CREAT_TRR,3),
                             fitter=lrm, xtrans=immuno_model_mice, data=data_immuno_model,
                             pr=FALSE)

print(HOSP_1Y, coef=FALSE)
summary(HOSP_1Y)
anova(HOSP_1Y)

####################################################################################### 
# R Code for Subgroup Analyses - Basiliximab Only
####################################################################################### 

# Limit dataset to only those who received basiliximab; n = 4,395
describe(as.factor(data_immuno_base$REC_SIMULECT))
data_immuno_basil <- subset(data_immuno_base, REC_SIMULECT == 1)
dim(data_immuno_basil)

## Select desired variables
vars_immuno_basil = c("GROUP", "TT_Out_Death_30d", "Out_Death_30d", "TT_Out_Death_1", 
                      "Out_Death_1", "OUT_GRAFT_1", "TT_OUT_GRAFT_1", "OUT_TREAT_REJ_1",
                      "HF_ETIOL", "SEX_MATCH", "DON_HTN", "DON_DM", "REC_DM", "ISCHTIME_MIN",
                      "AGE", "AGE_DON", "REC_RACE", "CPRA_PEAK", "REC_MALIG", "REC_CMV", "CREAT_TRR",
                      "OUT_HOSP_INF")

## Limit data to desired variables and look at number of patients & variables
data_immuno_basil = data_immuno_basil[, vars_immuno_basil]

## Look at selected data files
dim(data_immuno_basil)        # n = 4,395, 22 variables

data_immuno_basil <- transform(GROUP=as.factor(GROUP), HF_ETIOL=as.factor(HF_ETIOL),
                               Out_Death_30d=as.numeric(Out_Death_30d),
                               Out_Death_1=as.numeric(Out_Death_1),
                               OUT_GRAFT_1=as.numeric(OUT_GRAFT_1),
                               SEX_MATCH=as.factor(SEX_MATCH), DON_HTN=as.factor(DON_HTN),
                               DON_DM=as.factor(DON_DM), REC_DM=as.factor(REC_DM),
                               REC_RACE=as.factor(REC_RACE), REC_MALIG=as.factor(REC_MALIG),
                               REC_CMV=as.factor(REC_CMV), OUT_HOSP_INF=as.factor(OUT_HOSP_INF),
                               data_immuno_basil)

describe(data_immuno_basil$GROUP)

###### Impute missing data using 'MICE' package using predictive mean matching #####

## Get the number of missings per variable (n and %)
miss_var_summary(data_immuno_basil)

basil_model_mice <- mice(data_immuno_basil, m = 20, seed = 0621, print = F, method = "pmm")

set.seed(0621)
dd <- datadist(data_immuno_basil)
options(datadist="dd")

################ Cox Regression Using Missing Data ###############

# Convert death outcome variables to integer
data_immuno_basil$Out_Death_30d = as.integer(data_immuno_basil$Out_Death_30d)

# Fit Cox Regression for 30-Day Mortality Using Imputed Values
Surv_30d <- Surv(time = data_immuno_basil$TT_Out_Death_30d, event = data_immuno_basil$Out_Death_30d)
COX_30d <- fit.mult.impute(Surv_30d ~ GROUP + SEX_MATCH + HF_ETIOL + 
                             DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                             rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                             REC_CMV + rcs(CREAT_TRR,3),
                           fitter=cph, xtrans=basil_model_mice, data=data_immuno_basil, 
                           pr=FALSE, x = TRUE, y = TRUE)

print(COX_30d, coef=FALSE)
summary(COX_30d)
anova(COX_30d)

# Fit Cox Regression for 1-Year Mortality Using Imputed Values
data_immuno_basil$Out_Death_1 = as.integer(data_immuno_basil$Out_Death_1)

Surv_1Y <- Surv(time = data_immuno_basil$TT_Out_Death_1, 
                event = data_immuno_basil$Out_Death_1)
COX_1Y <- fit.mult.impute(Surv_1Y ~ GROUP + SEX_MATCH + HF_ETIOL + 
                            DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                            rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                            REC_CMV + rcs(CREAT_TRR,3),
                          fitter=cph, xtrans=basil_model_mice, data=data_immuno_basil, 
                          pr=FALSE, x = TRUE, y = TRUE)

print(COX_1Y, coef=FALSE)
summary(COX_1Y)
anova(COX_1Y)

# Fit Cox Regression for 1-Year Graft Failure Using Imputed Values
data_immuno_basil$OUT_GRAFT_1 = as.integer(data_immuno_basil$OUT_GRAFT_1)

Graft_1Y <- Surv(time = data_immuno_basil$TT_OUT_GRAFT_1, 
                 event = data_immuno_basil$OUT_GRAFT_1)
GRAFT_1Y <- fit.mult.impute(Graft_1Y ~ GROUP + SEX_MATCH + HF_ETIOL + 
                              DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                              rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                              REC_CMV + rcs(CREAT_TRR,3),
                            fitter=cph, xtrans=basil_model_mice, data=data_immuno_basil, 
                            pr=FALSE, x = TRUE, y = TRUE)

print(GRAFT_1Y, coef=FALSE)
summary(GRAFT_1Y)
anova(GRAFT_1Y)

# Fit Cox Regression for 1-Year Drug-Treated Rejection Using Imputed Values
data_immuno_basil$OUT_TREAT_REJ_1 = as.integer(data_immuno_basil$OUT_TREAT_REJ_1)

INFECT_1Y <- fit.mult.impute(OUT_TREAT_REJ_1 ~ GROUP + SEX_MATCH + HF_ETIOL + 
                               DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                               rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                               REC_CMV + rcs(CREAT_TRR,3),
                             fitter=lrm, xtrans=basil_model_mice, data=data_immuno_basil,
                             pr=FALSE)

print(INFECT_1Y, coef=FALSE)
summary(INFECT_1Y)
anova(INFECT_1Y)

# Fit Logistic Regression for 1-Year Hospitalization for Infection Using Imputed Values
data_immuno_basil$OUT_HOSP_INF = as.integer(data_immuno_basil$OUT_HOSP_INF)

HOSP_1Y <- fit.mult.impute(OUT_HOSP_INF ~ GROUP + SEX_MATCH + HF_ETIOL + 
                             DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                             rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                             REC_CMV + rcs(CREAT_TRR,3),
                           fitter=lrm, xtrans=basil_model_mice, data=data_immuno_basil,
                           pr=FALSE)

print(HOSP_1Y, coef=FALSE)
summary(HOSP_1Y)
anova(HOSP_1Y)

####################################################################################### 
# R Code for Subgroup Analyses - Thymo Only
####################################################################################### 

# Limit dataset to only those who received basiliximab; n = 3,083
describe(as.factor(data_immuno_base$REC_THYMO))
data_immuno_thymo <- subset(data_immuno_base, REC_THYMO == 1)
dim(data_immuno_thymo)

## Select desired variables
vars_immuno_thymo = c("GROUP", "TT_Out_Death_30d", "Out_Death_30d", "TT_Out_Death_1", 
                      "Out_Death_1", "OUT_GRAFT_1", "TT_OUT_GRAFT_1", "OUT_TREAT_REJ_1",
                      "HF_ETIOL", "SEX_MATCH", "DON_HTN", "DON_DM", "REC_DM", "ISCHTIME_MIN",
                      "AGE", "AGE_DON", "REC_RACE", "CPRA_PEAK", "REC_MALIG", "REC_CMV", "CREAT_TRR",
                      "OUT_HOSP_INF")

## Limit data to desired variables and look at number of patients & variables
data_immuno_thymo = data_immuno_thymo[, vars_immuno_thymo]

## Look at selected data files
dim(data_immuno_thymo)        # n = 3,083, 22 variables

data_immuno_thymo <- transform(GROUP=as.factor(GROUP), HF_ETIOL=as.factor(HF_ETIOL),
                               Out_Death_30d=as.numeric(Out_Death_30d),
                               Out_Death_1=as.numeric(Out_Death_1),
                               OUT_GRAFT_1=as.numeric(OUT_GRAFT_1),
                               SEX_MATCH=as.factor(SEX_MATCH), DON_HTN=as.factor(DON_HTN),
                               DON_DM=as.factor(DON_DM), REC_DM=as.factor(REC_DM),
                               REC_RACE=as.factor(REC_RACE), REC_MALIG=as.factor(REC_MALIG),
                               REC_CMV=as.factor(REC_CMV), OUT_HOSP_INF=as.factor(OUT_HOSP_INF),
                               data_immuno_thymo)

describe(data_immuno_thymo$GROUP)

###### Impute missing data using 'MICE' package using predictive mean matching #####

## Get the number of missings per variable (n and %)
miss_var_summary(data_immuno_thymo)

thymo_model_mice <- mice(data_immuno_thymo, m = 20, seed = 0621, print = F, method = "pmm")

set.seed(0621)
dd <- datadist(data_immuno_thymo)
options(datadist="dd")

################ Cox Regression Using Missing Data ###############

# Convert death outcome variables to integer
data_immuno_thymo$Out_Death_30d = as.integer(data_immuno_thymo$Out_Death_30d)

# Fit Cox Regression for 30-Day Mortality Using Imputed Values
Surv_30d <- Surv(time = data_immuno_thymo$TT_Out_Death_30d, event = data_immuno_thymo$Out_Death_30d)
COX_30d <- fit.mult.impute(Surv_30d ~ GROUP + SEX_MATCH + HF_ETIOL + 
                             DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                             rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                             REC_CMV + rcs(CREAT_TRR,3),
                           fitter=cph, xtrans=thymo_model_mice, data=data_immuno_thymo, 
                           pr=FALSE, x = TRUE, y = TRUE)

print(COX_30d, coef=FALSE)
summary(COX_30d)
anova(COX_30d)

# Fit Cox Regression for 1-Year Mortality Using Imputed Values
data_immuno_thymo$Out_Death_1 = as.integer(data_immuno_thymo$Out_Death_1)

Surv_1Y <- Surv(time = data_immuno_thymo$TT_Out_Death_1, 
                event = data_immuno_thymo$Out_Death_1)
COX_1Y <- fit.mult.impute(Surv_1Y ~ GROUP + SEX_MATCH + HF_ETIOL + 
                            DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                            rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                            REC_CMV + rcs(CREAT_TRR,3),
                          fitter=cph, xtrans=thymo_model_mice, data=data_immuno_thymo, 
                          pr=FALSE, x = TRUE, y = TRUE)

print(COX_1Y, coef=FALSE)
summary(COX_1Y)
anova(COX_1Y)

# Fit Cox Regression for 1-Year Graft Failure Using Imputed Values
data_immuno_thymo$OUT_GRAFT_1 = as.integer(data_immuno_thymo$OUT_GRAFT_1)

Graft_1Y <- Surv(time = data_immuno_thymo$TT_OUT_GRAFT_1, 
                 event = data_immuno_thymo$OUT_GRAFT_1)
GRAFT_1Y <- fit.mult.impute(Graft_1Y ~ GROUP + SEX_MATCH + HF_ETIOL + 
                              DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                              rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                              REC_CMV + rcs(CREAT_TRR,3),
                            fitter=cph, xtrans=thymo_model_mice, data=data_immuno_thymo, 
                            pr=FALSE, x = TRUE, y = TRUE)

print(GRAFT_1Y, coef=FALSE)
summary(GRAFT_1Y)
anova(GRAFT_1Y)

# Fit Cox Regression for 1-Year Drug-Treated Rejection Using Imputed Values
data_immuno_thymo$OUT_TREAT_REJ_1 = as.integer(data_immuno_thymo$OUT_TREAT_REJ_1)

INFECT_1Y <- fit.mult.impute(OUT_TREAT_REJ_1 ~ GROUP + SEX_MATCH + HF_ETIOL + 
                               DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                               rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                               REC_CMV + rcs(CREAT_TRR,3),
                             fitter=lrm, xtrans=thymo_model_mice, data=data_immuno_thymo,
                             pr=FALSE)

print(INFECT_1Y, coef=FALSE)
summary(INFECT_1Y)
anova(INFECT_1Y)

# Fit Logistic Regression for 1-Year Hospitalization for Infection Using Imputed Values
data_immuno_thymo$OUT_HOSP_INF = as.integer(data_immuno_thymo$OUT_HOSP_INF)

HOSP_1Y <- fit.mult.impute(OUT_HOSP_INF ~ GROUP + SEX_MATCH + HF_ETIOL + 
                             DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                             rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                             REC_CMV + rcs(CREAT_TRR,3),
                           fitter=lrm, xtrans=thymo_model_mice, data=data_immuno_thymo,
                           pr=FALSE)

print(HOSP_1Y, coef=FALSE)
summary(HOSP_1Y)
anova(HOSP_1Y)

####################################################################################### 
# R Code for Subgroup Analyses - Basiliximab Vs. Thymo in New and Old Separately
####################################################################################### 

# Create group variable for basil vs. thymo
data_immuno_base <- mutate(data_immuno_base, GROUP_BASIL = case_when(REC_SIMULECT==1 ~ "BASIL", REC_THYMO==1 ~ "THYMO"),
                   
)
describe(as.factor(data_immuno_base$GROUP_BASIL))



# Limit dataset to only those in the old era; n = 4,705
describe(as.factor(data_immuno_base$GROUP))
describe(as.factor(data_immuno_base$GROUP_BASIL))
data_immuno_old <- subset(data_immuno_base, GROUP == "OLD") # Limit to old system
data_immuno_old <- subset(data_immuno_old, GROUP_BASIL == "BASIL" | GROUP_BASIL == "THYMO") # Limit to basil or thymo
dim(data_immuno_old)

## Select desired variables
vars_immuno_old = c("GROUP_BASIL", "TT_Out_Death_30d", "Out_Death_30d", "TT_Out_Death_1", 
                      "Out_Death_1", "OUT_GRAFT_1", "TT_OUT_GRAFT_1", "OUT_TREAT_REJ_1",
                      "HF_ETIOL", "SEX_MATCH", "DON_HTN", "DON_DM", "REC_DM", "ISCHTIME_MIN",
                      "AGE", "AGE_DON", "REC_RACE", "CPRA_PEAK", "REC_MALIG", "REC_CMV", "CREAT_TRR",
                      "OUT_HOSP_INF")

## Limit data to desired variables and look at number of patients & variables
data_immuno_old = data_immuno_old[, vars_immuno_old]

## Look at selected data files
dim(data_immuno_old)        # n = 4,705, 22 variables

data_immuno_old <- transform(GROUP_BASIL=as.factor(GROUP_BASIL), HF_ETIOL=as.factor(HF_ETIOL),
                               Out_Death_30d=as.numeric(Out_Death_30d),
                               Out_Death_1=as.numeric(Out_Death_1),
                               OUT_GRAFT_1=as.numeric(OUT_GRAFT_1),
                               SEX_MATCH=as.factor(SEX_MATCH), DON_HTN=as.factor(DON_HTN),
                               DON_DM=as.factor(DON_DM), REC_DM=as.factor(REC_DM),
                               REC_RACE=as.factor(REC_RACE), REC_MALIG=as.factor(REC_MALIG),
                               REC_CMV=as.factor(REC_CMV), OUT_HOSP_INF=as.factor(OUT_HOSP_INF),
                             data_immuno_old)

describe(data_immuno_old$GROUP_BASIL)

###### Impute missing data using 'MICE' package using predictive mean matching #####

## Get the number of missings per variable (n and %)
miss_var_summary(data_immuno_old)

old_model_mice <- mice(data_immuno_old, m = 20, seed = 0621, print = F, method = "pmm")

set.seed(0621)
dd <- datadist(data_immuno_old)
options(datadist="dd")

################ Cox Regression Using Missing Data ###############

# Convert death outcome variables to integer
data_immuno_old$Out_Death_30d = as.integer(data_immuno_old$Out_Death_30d)

# Fit Cox Regression for 30-Day Mortality Using Imputed Values
Surv_30d <- Surv(time = data_immuno_old$TT_Out_Death_30d, event = data_immuno_old$Out_Death_30d)
COX_30d <- fit.mult.impute(Surv_30d ~ GROUP_BASIL + SEX_MATCH + HF_ETIOL + 
                             DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                             rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                             REC_CMV + rcs(CREAT_TRR,3),
                           fitter=cph, xtrans=old_model_mice, data=data_immuno_old, 
                           pr=FALSE, x = TRUE, y = TRUE)

print(COX_30d, coef=FALSE)
summary(COX_30d)
anova(COX_30d)

# Fit Cox Regression for 1-Year Mortality Using Imputed Values
data_immuno_old$Out_Death_1 = as.integer(data_immuno_old$Out_Death_1)

Surv_1Y <- Surv(time = data_immuno_old$TT_Out_Death_1, 
                event = data_immuno_old$Out_Death_1)
COX_1Y <- fit.mult.impute(Surv_1Y ~ GROUP_BASIL + SEX_MATCH + HF_ETIOL + 
                            DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                            rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                            REC_CMV + rcs(CREAT_TRR,3),
                          fitter=cph, xtrans=old_model_mice, data=data_immuno_old, 
                          pr=FALSE, x = TRUE, y = TRUE)

print(COX_1Y, coef=FALSE)
summary(COX_1Y)
anova(COX_1Y)

# Fit Cox Regression for 1-Year Graft Failure Using Imputed Values
data_immuno_old$OUT_GRAFT_1 = as.integer(data_immuno_old$OUT_GRAFT_1)

Graft_1Y <- Surv(time = data_immuno_old$TT_OUT_GRAFT_1, 
                 event = data_immuno_old$OUT_GRAFT_1)
GRAFT_1Y <- fit.mult.impute(Graft_1Y ~ GROUP_BASIL + SEX_MATCH + HF_ETIOL + 
                              DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                              rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                              REC_CMV + rcs(CREAT_TRR,3),
                            fitter=cph, xtrans=old_model_mice, data=data_immuno_old, 
                            pr=FALSE, x = TRUE, y = TRUE)

print(GRAFT_1Y, coef=FALSE)
summary(GRAFT_1Y)
anova(GRAFT_1Y)

# Fit Cox Regression for 1-Year Drug-Treated Rejection Using Imputed Values
data_immuno_old$OUT_TREAT_REJ_1 = as.integer(data_immuno_old$OUT_TREAT_REJ_1)

REJ_1Y <- fit.mult.impute(OUT_TREAT_REJ_1 ~ GROUP_BASIL + SEX_MATCH + HF_ETIOL + 
                               DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                               rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                               REC_CMV + rcs(CREAT_TRR,3),
                             fitter=lrm, xtrans=old_model_mice, data=data_immuno_old,
                             pr=FALSE)

print(REJ_1Y, coef=FALSE)
summary(REJ_1Y)
anova(REJ_1Y)

# Fit Logistic Regression for 1-Year Hospitalization for Infection Using Imputed Values
data_immuno_old$OUT_HOSP_INF = as.integer(data_immuno_old$OUT_HOSP_INF)

HOSP_1Y <- fit.mult.impute(OUT_HOSP_INF ~ GROUP_BASIL + SEX_MATCH + HF_ETIOL + 
                             DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                             rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                             REC_CMV + rcs(CREAT_TRR,3),
                           fitter=lrm, xtrans=old_model_mice, data=data_immuno_old,
                           pr=FALSE)

print(HOSP_1Y, coef=FALSE)
summary(HOSP_1Y)
anova(HOSP_1Y)

# Limit dataset to only those in the new era; n = 2,682
describe(as.factor(data_immuno_base$GROUP))
describe(as.factor(data_immuno_base$GROUP_BASIL))
data_immuno_new <- subset(data_immuno_base, GROUP == "NEW") # Limit to new system
data_immuno_new <- subset(data_immuno_new, GROUP_BASIL == "BASIL" | GROUP_BASIL == "THYMO") # Limit to basil or thymo
dim(data_immuno_new)

## Select desired variables
vars_immuno_new = c("GROUP_BASIL", "TT_Out_Death_30d", "Out_Death_30d", "TT_Out_Death_1", 
                    "Out_Death_1", "OUT_GRAFT_1", "TT_OUT_GRAFT_1", "OUT_TREAT_REJ_1",
                    "HF_ETIOL", "SEX_MATCH", "DON_HTN", "DON_DM", "REC_DM", "ISCHTIME_MIN",
                    "AGE", "AGE_DON", "REC_RACE", "CPRA_PEAK", "REC_MALIG", "REC_CMV", "CREAT_TRR",
                    "OUT_HOSP_INF")

## Limit data to desired variables and look at number of patients & variables
data_immuno_new = data_immuno_new[, vars_immuno_new]

## Look at selected data files
dim(data_immuno_new)        # n = 2,682, 22 variables

data_immuno_new <- transform(GROUP_BASIL=as.factor(GROUP_BASIL), HF_ETIOL=as.factor(HF_ETIOL),
                             Out_Death_30d=as.numeric(Out_Death_30d),
                             Out_Death_1=as.numeric(Out_Death_1),
                             OUT_GRAFT_1=as.numeric(OUT_GRAFT_1),
                             SEX_MATCH=as.factor(SEX_MATCH), DON_HTN=as.factor(DON_HTN),
                             DON_DM=as.factor(DON_DM), REC_DM=as.factor(REC_DM),
                             REC_RACE=as.factor(REC_RACE), REC_MALIG=as.factor(REC_MALIG),
                             REC_CMV=as.factor(REC_CMV), OUT_HOSP_INF=as.factor(OUT_HOSP_INF),
                             data_immuno_new)

describe(data_immuno_new$GROUP_BASIL)

###### Impute missing data using 'MICE' package using predictive mean matching #####

## Get the number of missings per variable (n and %)
miss_var_summary(data_immuno_new)

new_model_mice <- mice(data_immuno_new, m = 20, seed = 0621, print = F, method = "pmm")

set.seed(0621)
dd <- datadist(data_immuno_new)
options(datadist="dd")

################ Cox Regression Using Missing Data ###############

# Convert death outcome variables to integer
data_immuno_new$Out_Death_30d = as.integer(data_immuno_new$Out_Death_30d)

# Fit Cox Regression for 30-Day Mortality Using Imputed Values
Surv_30d <- Surv(time = data_immuno_new$TT_Out_Death_30d, event = data_immuno_new$Out_Death_30d)
COX_30d <- fit.mult.impute(Surv_30d ~ GROUP_BASIL + SEX_MATCH + HF_ETIOL + 
                             DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                             rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                             REC_CMV + rcs(CREAT_TRR,3),
                           fitter=cph, xtrans=new_model_mice, data=data_immuno_new, 
                           pr=FALSE, x = TRUE, y = TRUE)

print(COX_30d, coef=FALSE)
summary(COX_30d)
anova(COX_30d)

# Fit Cox Regression for 1-Year Mortality Using Imputed Values
data_immuno_new$Out_Death_1 = as.integer(data_immuno_new$Out_Death_1)

Surv_1Y <- Surv(time = data_immuno_new$TT_Out_Death_1, 
                event = data_immuno_new$Out_Death_1)
COX_1Y <- fit.mult.impute(Surv_1Y ~ GROUP_BASIL + SEX_MATCH + HF_ETIOL + 
                            DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                            rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                            REC_CMV + rcs(CREAT_TRR,3),
                          fitter=cph, xtrans=new_model_mice, data=data_immuno_new, 
                          pr=FALSE, x = TRUE, y = TRUE)

print(COX_1Y, coef=FALSE)
summary(COX_1Y)
anova(COX_1Y)

# Fit Cox Regression for 1-Year Graft Failure Using Imputed Values
data_immuno_new$OUT_GRAFT_1 = as.integer(data_immuno_new$OUT_GRAFT_1)

Graft_1Y <- Surv(time = data_immuno_new$TT_OUT_GRAFT_1, 
                 event = data_immuno_new$OUT_GRAFT_1)
GRAFT_1Y <- fit.mult.impute(Graft_1Y ~ GROUP_BASIL + SEX_MATCH + HF_ETIOL + 
                              DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                              rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                              REC_CMV + rcs(CREAT_TRR,3),
                            fitter=cph, xtrans=new_model_mice, data=data_immuno_new, 
                            pr=FALSE, x = TRUE, y = TRUE)

print(GRAFT_1Y, coef=FALSE)
summary(GRAFT_1Y)
anova(GRAFT_1Y)

# Fit Cox Regression for 1-Year Drug-Treated Rejection Using Imputed Values
data_immuno_new$OUT_TREAT_REJ_1 = as.integer(data_immuno_new$OUT_TREAT_REJ_1)

REJ_1Y <- fit.mult.impute(OUT_TREAT_REJ_1 ~ GROUP_BASIL + SEX_MATCH + HF_ETIOL + 
                            DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                            rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                            REC_CMV + rcs(CREAT_TRR,3),
                          fitter=lrm, xtrans=new_model_mice, data=data_immuno_new,
                          pr=FALSE)

print(REJ_1Y, coef=FALSE)
summary(REJ_1Y)
anova(REJ_1Y)

# Fit Logistic Regression for 1-Year Hospitalization for Infection Using Imputed Values
data_immuno_new$OUT_HOSP_INF = as.integer(data_immuno_new$OUT_HOSP_INF)

HOSP_1Y <- fit.mult.impute(OUT_HOSP_INF ~ GROUP_BASIL + SEX_MATCH + HF_ETIOL + 
                             DON_HTN + DON_DM + REC_DM + rcs(ISCHTIME_MIN,3) + rcs(AGE,3) + 
                             rcs(AGE_DON,3) + REC_RACE + rcs(CPRA_PEAK,3) + REC_MALIG + 
                             REC_CMV + rcs(CREAT_TRR,3),
                           fitter=lrm, xtrans=new_model_mice, data=data_immuno_new,
                           pr=FALSE)

print(HOSP_1Y, coef=FALSE)
summary(HOSP_1Y)
anova(HOSP_1Y)