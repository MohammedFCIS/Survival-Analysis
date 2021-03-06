---
title: "SA analysis for expression data from 290 primary colorectal cancers"
author: "Abdallah Bekhit, Mohammed Ali, Sami Emad"
date: "November 14, 2017"
output:
  pdf_document:
    toc: yes
    toc_depth: 2
  html_document:
    toc: yes
    toc_depth: 2
    number_sections: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
# Overview
This document explains a variety of survival analysis methods that performed on [Expression data from 290 primary colorectal cancers](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14333)

# Setup
This section aims to load needed libs to perform survival analysis functionalties.

### Load survival library
```{r}
suppressMessages(library(survival))
```

### Load graphics library
```{r}
suppressMessages(library(ggfortify))
```

### Load the dataset
```{r}
load("D:/moh/DSTI/Courses/Survival Analysis using R -S17/CRC_226_GSE14333.RData")
```

# Dataset Exploration
This section aims to explore dataset data and metadate before performing any kind of analysis.

### Metadata
```{r}
clinical_metadata
```

### Structure
```{r}
str(clinical_data)
```

### Sample Data
```{r}
head(clinical_data)
```

### Full statistics
```{r}
summary(clinical_data)
```

# Basic Non-Parametric
Here we are performing basic analysis using **_Kaplan-Meier_** and **_Fleming-Harrington_** methods using different time units in **Months** which is the default and in **Years**. 
The rest of analysis will use the default time unitin **Months**

##In months

### Kaplan-Meier non-parametric 
```{r}
kmsurvival_month <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ 1)
summary(kmsurvival_month)
autoplot(kmsurvival_month, xlab="Time in Month", ylab="Survival Probability",
         surv.linetype = 'dashed',
         surv.colour = 'blue',
         conf.int.fill = 'dodgerblue3', conf.int.alpha = 0.5)
```

### Fleming-Harrington non-parametric analysis
```{r}
fhsurvival_month <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ 1,
                            type="fleming-harrington")
summary(fhsurvival_month)
autoplot(fhsurvival_month, xlab="Time in Month", ylab="Survival Probability",
         surv.linetype = 'dashed',
         surv.colour = 'blue',
         conf.int.fill = 'dodgerblue3', conf.int.alpha = 0.5)
```

## In years

### Kaplan-Meier non-parametric 
```{r}
kmsurvival_year <- survfit(Surv(clinical_data$dfs_time/12, clinical_data$dfs_event) ~ 1)
summary(kmsurvival_year)
autoplot(kmsurvival_year, xlab="Time in Years", ylab="Survival Probability", 
         surv.linetype = 'dashed',
         surv.colour = 'blue',
         conf.int.fill = 'dodgerblue3', conf.int.alpha = 0.5)
```

### Fleming-Harrington non-parametric analysis
```{r}
fhsurvival_year <- survfit(Surv(clinical_data$dfs_time/12, clinical_data$dfs_event) ~ 1,
                           type="fleming-harrington")
summary(fhsurvival_year)
autoplot(fhsurvival_year, xlab="Time in Years", ylab="Survival Probability",
         surv.linetype = 'dashed',
         surv.colour = 'blue',
         conf.int.fill = 'dodgerblue3', conf.int.alpha = 0.5)
```

# Non-Parametric Groups Analysis
This section aims to perform different group analysis on the dataset.

## Gender Analysis
It seems from the analysis that males are living longer than females

### Kaplan-Meier non-parametric
```{r}
kmsurvival1_gender <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) 
                              ~ clinical_data$gender)
summary(kmsurvival1_gender)
autoplot(kmsurvival1_gender,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time",
         ylab="Survival Probability")
```

### Fleming-Harrington non-parametric analysis
```{r}
fhsurvival1_gender <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) 
                              ~ clinical_data$gender, type="fleming-harrington")
summary(fhsurvival1_gender)
autoplot(fhsurvival1_gender,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time",
         ylab="Survival Probability")
```

## Tumor Location
It seems from the analysis that Colon cancer is the most dengraous cancer location.

### Kaplan-Meier non-parametric
```{r}
kmsurvival1_location <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) 
                                ~ clinical_data$location)
summary(kmsurvival1_location)
autoplot(kmsurvival1_location,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time",
         ylab="Survival Probability")
```

### Fleming-Harrington non-parametric analysis
```{r}
fhsurvival1_location <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event)
                                ~ clinical_data$location, type="fleming-harrington")
summary(fhsurvival1_location)
autoplot(fhsurvival1_location,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time",
         ylab="Survival Probability")
```

## Cancer Stage
As expected, stage **A** is with the least death rate.

### Kaplan-Meier non-parametric
```{r}
kmsurvival1_stage <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) 
                             ~ clinical_data$dukes_stage)
summary(kmsurvival1_stage)
autoplot(kmsurvival1_stage, 
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time",
         ylab="Survival Probability")
```

### Fleming-Harrington non-parametric analysis
```{r}
fhsurvival1_stage <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event)
                             ~ clinical_data$dukes_stage, type="fleming-harrington")
summary(fhsurvival1_stage)
autoplot(fhsurvival1_stage,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time",
         ylab="Survival Probability")
```

## Adjuvant Radiation Therapy
It seems that who had the thetapy had a better chance to live.

### Kaplan-Meier non-parametric analysis
```{r}
kmsurvival1_adjXRT <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) 
                              ~ clinical_data$adjXRT)
summary(kmsurvival1_adjXRT)
autoplot(kmsurvival1_adjXRT,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time",
         ylab="Survival Probability")
```

### Fleming-Harrington non-parametric analysis
```{r}
fhsurvival1_adjXRT <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event)
                              ~ clinical_data$adjXRT, type="fleming-harrington")
summary(fhsurvival1_adjXRT)
autoplot(fhsurvival1_adjXRT, 
         censor.shape = '*', facets = TRUE, ncol = 2,xlab="Time",
         ylab="Survival Probability")
```

## Adjuvant Chemotherapy Analysis
Also the chemo therapy helps but not as much as radiation therapy.

### Kaplan-Meier non-parametric analysis
```{r}
kmsurvival1_adjCTX <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) 
                              ~ clinical_data$adjCTX)
summary(kmsurvival1_adjCTX)
autoplot(kmsurvival1_adjCTX, 
         censor.shape = '*', facets = TRUE, ncol = 2,xlab="Time", 
         ylab="Survival Probability")
```

### Fleming-Harrington non-parametric analysis
```{r}
fhsurvival1_adjCTX <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event)
                              ~ clinical_data$adjCTX, type="fleming-harrington")
summary(fhsurvival1_adjCTX)
autoplot(fhsurvival1_adjCTX,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time", 
         ylab="Survival Probability")
```

## Age
* People older than 80 are collapsing so fast.
* People between 20 and 40 death rate are normal butwith a huge amount of events at each time
* People between 60 and 80 are the strongest group.

### Create age groups
```{r}
clinical_data$age_groups <- ifelse(clinical_data$age_diag >= 20 &
                                     clinical_data$age_diag < 40, "20-40",
                                   ifelse(clinical_data$age_diag >= 40 &
                                            clinical_data$age_diag < 60, "40-60",
                                          ifelse(clinical_data$age_diag >= 60 & 
                                                   clinical_data$age_diag < 80, "60-80",
                                                 "> 80")))
## Convert it to factor
clinical_data$age_groups <- factor(clinical_data$age_groups)
## Perform the estimations using Kaplan-Meier non-parametric 
kmsurvival1_age <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ 
                             clinical_data$age_groups)
summary(kmsurvival1_age)
autoplot(kmsurvival1_age, 
         censor.shape = '*', facets = TRUE, ncol = 2,xlab="Time", 
         ylab="Survival Probability")
```

### Fleming-Harrington non-parametric analysis
```{r}
fhsurvival1_age <- survfit(Surv(clinical_data$dfs_time, clinical_data$dfs_event) 
                           ~ clinical_data$age_groups, type="fleming-harrington")
summary(fhsurvival1_age)
autoplot(fhsurvival1_age,
         censor.shape = '*', facets = TRUE, ncol = 2, xlab="Time", 
         ylab="Survival Probability")
```

# Semi-Parameter Analysis

### Cox proportional hazard model - coefficients and hazard rates
First we need to create a factor of the parameter to be used in **semi_paramteric** and **parametric** analysis
```{r}
factors <- cbind(clinical_data$location, clinical_data$dukes_stage,
                 clinical_data$age_diag,
                 clinical_data$gender, clinical_data$adjXRT, clinical_data$adjCTX)
coxph <- coxph(Surv(clinical_data$dfs_time, clinical_data$dfs_event) ~ 
                 factors, method="breslow")
summary(coxph)
```


# Paramteric

## Exponential
```{r}
exponential <- survreg(Surv(clinical_data$dfs_time,clinical_data$dfs_event) 
                       ~ factors, dist="exponential")
summary(exponential)
```

## Weibull
```{r}
weibull <- survreg(Surv(clinical_data$dfs_time,clinical_data$dfs_event) 
                   ~ factors, dist="weibull")
summary(weibull)
```

## log-logistic
```{r}
loglogistic <- survreg(Surv(clinical_data$dfs_time,clinical_data$dfs_event)
                       ~ factors, dist="loglogistic")
summary(loglogistic)
```

## Comparison between those methods
```{r}
plot(kmsurvival_month)
curve(plnorm(x, coef(weibull)[1], weibull$scale, lower.tail=FALSE), 
      from=0, to=140, col="blue", ylim=c(0,1), lwd=2, add=T, lty=5, xlab="", ylab="")
curve(plnorm(x, coef(loglogistic)[1], loglogistic$scale, lower.tail=FALSE),
      from=0, to=140, col="green", ylim=c(0,1), lwd=2, add=T, lty=5, xlab="", ylab="")
curve(plnorm(x, coef(exponential)[1], exponential$scale, lower.tail=FALSE), 
      from=0, to=140, col="red", ylim=c(0,1), lwd=2, add=T, lty=5, xlab="", ylab="")
legend('topright',c("None-Parametric", "weibull", "log", "exponential"), 
       col=c("black", "blue", "green", "red"), lty=c(1,1), lwd=c(5,5), cex=1.2,
       inset=c(.05,0), bty="n")
title(main="Parametric methods vs. Non Parametric",
        ylab = "Percent of Surviving Patients", xlab="Time")
```
