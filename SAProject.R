#Load survival library
suppressMessages(library(survival))
#Load graphics library
suppressMessages(library(ggfortify))
#Load the dataset
#load("D:/DSTI/Survival Analysis/example_data/CRC_226_GSE14333.RData")
load("D:/moh/DSTI/Courses/Survival Analysis using R -S17/CRC_226_GSE14333.RData")
# Descriptive statistics
##Full statistics
summary(clinical_data)

# Analysis for time in months
## Kaplan-Meier non-parametric 
kmsurvival_month <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ 1)
summary(kmsurvival_month)
autoplot(kmsurvival_month, xlab="Time in Month", ylab="Survival Probability", surv.linetype = 'dashed', surv.colour = 'blue',
         conf.int.fill = 'dodgerblue3', conf.int.alpha = 0.5)
## Fleming-Harrington non-parametric analysis
fhsurvival_month <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ 1, type="fleming-harrington")
summary(fhsurvival_month)
autoplot(fhsurvival_month, xlab="Time in Month", ylab="Survival Probability", surv.linetype = 'dashed', surv.colour = 'blue',
         conf.int.fill = 'dodgerblue3', conf.int.alpha = 0.5)

# Analysis for time in years
## Kaplan-Meier non-parametric 
kmsurvival_year <- survfit(Surv(clinical_data$dfs_time/12, clinical_data$dfs_event) ~ 1)
summary(kmsurvival_year)
autoplot(kmsurvival_year, xlab="Time in Years", ylab="Survival Probability", surv.linetype = 'dashed', surv.colour = 'blue',
         conf.int.fill = 'dodgerblue3', conf.int.alpha = 0.5)
## Fleming-Harrington non-parametric analysis
fhsurvival_year <- survfit(Surv(clinical_data$dfs_time/12, clinical_data$dfs_event) ~ 1, type="fleming-harrington")
summary(fhsurvival_year)
autoplot(fhsurvival_year, xlab="Time in Years", ylab="Survival Probability", surv.linetype = 'dashed', surv.colour = 'blue',
         conf.int.fill = 'dodgerblue3', conf.int.alpha = 0.5)
# Analysis by gender
## Kaplan-Meier non-parametric
kmsurvival1_gender <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ clinical_data$gender)
summary(kmsurvival1_gender)
autoplot(kmsurvival1_gender,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time", ylab="Survival Probability")
## Fleming-Harrington non-parametric analysis
fhsurvival1_gender <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ clinical_data$gender, type="fleming-harrington")
summary(fhsurvival1_gender)
autoplot(fhsurvival1_gender,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time", ylab="Survival Probability")

# Analysis by tumor location
## Kaplan-Meier non-parametric
kmsurvival1_location <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ clinical_data$location)
summary(kmsurvival1_location)
autoplot(kmsurvival1_location,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time", ylab="Survival Probability")
## Fleming-Harrington non-parametric analysis
fhsurvival1_location <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ clinical_data$location, type="fleming-harrington")
summary(fhsurvival1_location)
autoplot(fhsurvival1_location,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time", ylab="Survival Probability")

# Analysis by cancer stage 
## Kaplan-Meier non-parametric
kmsurvival1_stage <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ clinical_data$dukes_stage)
summary(kmsurvival1_stage)
autoplot(kmsurvival1_stage, 
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time", ylab="Survival Probability")
## Fleming-Harrington non-parametric analysis
fhsurvival1_stage <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ clinical_data$dukes_stage, type="fleming-harrington")
summary(fhsurvival1_stage)
autoplot(fhsurvival1_stage,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time", ylab="Survival Probability")

# Analysis by adjuvant radiation therapy
## Kaplan-Meier non-parametric analysis
kmsurvival1_adjXRT <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ clinical_data$adjXRT)
summary(kmsurvival1_adjXRT)
autoplot(kmsurvival1_adjXRT,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time", ylab="Survival Probability")
## Fleming-Harrington non-parametric analysis
fhsurvival1_adjXRT <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ clinical_data$adjXRT, type="fleming-harrington")
summary(fhsurvival1_adjXRT)
autoplot(fhsurvival1_adjXRT, 
         censor.shape = '*', facets = TRUE, ncol = 2,xlab="Time", ylab="Survival Probability")

# Analysis by adjuvant chemotherapy
# Kaplan-Meier non-parametric analysis
kmsurvival1_adjCTX <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ clinical_data$adjXRT)
summary(kmsurvival1_adjCTX)
autoplot(kmsurvival1_adjCTX, 
         censor.shape = '*', facets = TRUE, ncol = 2,xlab="Time", ylab="Survival Probability")
## Fleming-Harrington non-parametric analysis
fhsurvival1_adjCTX <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ clinical_data$adjXRT, type="fleming-harrington")
summary(fhsurvival1_adjCTX)
autoplot(fhsurvival1_adjCTX,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time", ylab="Survival Probability")

# Analysis by age
## Crete age groups
clinical_data$age_groups <- ifelse(clinical_data$age_diag >= 20 & clinical_data$age_diag < 40, "20-40",
                                   ifelse(clinical_data$age_diag >= 40 & clinical_data$age_diag < 60, "40-60",
                                          ifelse(clinical_data$age_diag >= 60 & clinical_data$age_diag < 80, "60-80", "> 80")))
## Convert it to factor
clinical_data$age_groups <- factor(clinical_data$age_groups)
## Perform the estimations using Kaplan-Meier non-parametric 
kmsurvival1_age <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ clinical_data$age_groups)
summary(kmsurvival1_age)
autoplot(kmsurvival1_age, 
         censor.shape = '*', facets = TRUE, ncol = 2,xlab="Time", ylab="Survival Probability")
## Fleming-Harrington non-parametric analysis
fhsurvival1_age <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ clinical_data$age_groups, type="fleming-harrington")
summary(fhsurvival1_age)
autoplot(fhsurvival1_age,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time", ylab="Survival Probability")

# Cox proportional hazard model - coefficients and hazard rates
factors <- cbind(clinical_data$location, clinical_data$dukes_stage, clinical_data$age_diag,
                 clinical_data$gender, clinical_data$adjXRT, clinical_data$adjCTX)
coxph <- coxph(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ 
                 factors, method="breslow")
summary(coxph)


# Exponential, Weibull, and log-logistic parametric model coefficients
exponential <- survreg(Surv(clinical_data$dfs_time,clinical_data$dfs_event) ~ factors, dist="exponential")
summary(exponential)

weibull <- survreg(Surv(clinical_data$dfs_time,clinical_data$dfs_event) ~ factors, dist="weibull")
summary(weibull)

loglogistic <- survreg(Surv(clinical_data$dfs_time,clinical_data$dfs_event) ~ factors, dist="loglogistic")
summary(loglogistic)
## Plot exponential model prediction
survreg.curves(loglogistic, "red")
