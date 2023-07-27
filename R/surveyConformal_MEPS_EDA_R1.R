## EDA code for "Design-based conformal prediction," Wieczorek (2023+)

## This file summarizes the initial EDA
## that led to the data filtering and modeling choices made in the paper.
## See also `surveyConformal_MEPS_Results_R1.R` for a streamlined analysis
## that simply reproduces results from the paper.




# We are using Medical Expenditure Panel Survey (MEPS) data,
# as in CQR paper by Y. Romano et al.:
# https://github.com/yromano/cqr/blob/master/get_meps_data/download_data.R

# MEPS requires the following Data Use Agreement when writing up results:
# https://meps.ahrq.gov/data_stats/data_use.jsp


library(foreign)
library(dplyr)
library(survey)
library(skimr)

# For illustrative purposes, we just chose one year of the MEPS data,
# starting with the 2015 full Year Consolidated Data File:
# This file contains MEPS survey data for calendar year 2015 obtained in
# rounds 3, 4, and 5 of Panel 19, and rounds 1, 2, and 3 of Panel 20.
df15 <- read.xport("h181.ssp")



#### What's in the data? And what were Romano et al. predicting? ####

# Romano et al. define a "utilization" variable as their response Y variable:
# https://github.com/yromano/cqr/blob/master/get_meps_data/meps_dataset_panel19_fy2015_reg.py

## def utilization(row):
##   return row['OBTOTV15'] + row['OPTOTV15'] + row['ERTOT15'] + row['IPNGTD15'] + row['HHTOTD15']


# Reading the MEPS documentation, from h181doc.pdf:

# OBTOTV15: total number of office-based visits reported for 2015
# OPTOTV15: total number of reported visits to hospital outpatient departments in 2015
#  ERTOT15: count of all emergency room visits reported for the survey year
# IPNGTD15: total number of nights associated with [hospital] discharges
# HHTOTD15: total number of days in 2015 where home health care was received ... from any type of paid or unpaid caregiver

## It appears that the first 3 are visits, but last 2 are days;
## but Romano et al. just sum them,
## so we'll do the same to be consistent with their examples.
## These all appear to exclude dental and vision visits,
##   as well as misc other services like chiro. or phys.therapy
##   (or are those included within the 1st?),
##   and don't deal with prescriptions at all.


# So we will need at least those variables.

# Also, what survey design features are available that we could use?
# There are approximate / anonymized PSUs and strata:
# VARPSU, VARSTR
# as well as person weights:
# PERWT15F

# Look at these variables:
df15 <- df15 |>
  mutate(utilization = OBTOTV15 + OPTOTV15 + ERTOT15 + IPNGTD15 + HHTOTD15)
df15 |>
  select(OBTOTV15, OPTOTV15, ERTOT15, IPNGTD15, HHTOTD15,
         utilization,
         VARPSU, VARSTR, PERWT15F) |>
  summary()

# HHTOTD15 == -9 is "not ascertained".
# Also, PERWT15F has some 0s. Why?

hist(df15$utilization)
length(unique(df15$utilization))





#### Subsetting to the SAQ respondents ####

## Are there any continuous predictor variables, so we can build a nontrivial model?
## https://github.com/yromano/cqr/blob/master/get_meps_data/meps_dataset_panel19_fy2015_reg.py
## Seems like only 'PCS42', 'MCS42', 'K6SUM42' are not simply categorical;
## and even those are indices, built to summarize answers to several yes/no Qs:
## h181doc.pdf explains these are from SAQ = Self Administered Questionnaire
## which was only fielded in SOME of the rounds of this panel survey,
## and only for people with certain criteria...

## "A special weight variable (SAQWT15F) has
## been designed to be used with the SAQ for persons who were age 18 and older at the interview
## date. This weight adjusts for SAQ non-response and weights to the U.S. civilian
## noninstitutionalized population..."

## So let's use THIS weight (and this subset of the data), instead of PERWT15F,
## so that we can focus on just the SAQ respondents,
## since they have at least SOME sort of "quant" predictor variable to be used.
summary(df15$SAQWT15F)
table(df15$SAQELIG)

df15 |>
  group_by(SAQELIG) |>
  skimr::skim(SAQWT15F)
## OK, so those SAQ weights really do drop the SAQ nonrespondents too,
## not only the ineligibles. Great!



# What's in this SAQ subset of the data?
df15saq <- df15 |>
  filter(SAQWT15F > 0)
summary(df15saq$SAQWT15F)
hist(df15saq$SAQWT15F)

df15saq |>
  skimr::skim(PCS42, MCS42, K6SUM42)
head(table(df15saq$PCS42))
head(table(df15saq$MCS42))
head(table(df15saq$K6SUM42))

## So: 'PCS42', 'MCS42', 'K6SUM42'
## The first two are based on Short-Form 12 Version 2
## which is 12 Qs that ask about healthy today, typical health, and past 4 weeks
## and become summarized as
## the Physical Component Summary (PCS) and the Mental Component Summary (MCS)
hist(df15saq$PCS42)
hist(df15saq$MCS42)
## and from the Kessler Index (K6) of non-specific psychological distress
## where 6 Qs ask about general mental health in past 30 days
## and you just sum these 6 answers, each on a 5-point scale from 0 (never) to 4 (always)
hist(df15saq$K6SUM42)
## BUT...
## many many people have values of -9 not ascertained or -1 inapplicable...
## Aha, in codebook for SAQELIG,
## around 10k were not eligible for these SAQ questions at all
## and another 3.5k just didn't answer them...
## So subsetting to only SAQWT15F>0 gets rid of all the inapplicables,
## but leaves a few not-ascertained (must be due to item NR rather than unit NR)
## ...
## Wait -- is high score or low score "better"?
# https://link.springer.com/referenceworkentry/10.1007/978-0-387-78665-0_6349
# For MCS and PCS, "Higher scores represent better physical health."
# For K6, based on the def'n it sounds like higher scores show "worse" mental health.




#### What to do with the remaining strata and PSUs? ####

table(df15saq$VARSTR, df15saq$VARPSU)
table(df15saq$VARPSU)

## Most strata have only exactly 2 PSUs (a few have 3)
## so we cannot really do a full example of k-fold CV appropriately.
## But for this paper, we just need to show a brief examples of how
##   there CAN be differences when doing conformal
##   when you do vs when you don't account for survey weights / svy design.
## So let's do a simple example...
## Let all the PSU==3 units be our test set.
## Show that if you use PSU==1 for training and PSU==2 for calibration,
## the MOE width for predictions on PSU==3 is different with vs without svy wts?
## In other words, if this WAS a one-off real data analysis,
##   we'd have ONE set of specific test cases,
##   and we'd do splitting ONCE
##   --and just within those constraints it matters whether or not
##     we weight to get conformal MOE.
## Finally, to look at conformal PI coverage, we can also compare:
##   "across many random shuffles of one-PSU-per-stratum each for test vs calib"
##   vs "across many shuffles of X% of the data each for test vs calib"

# For those "across many random shuffles" sims...
colSums(table(df15saq$VARSTR, df15saq$VARPSU) == 0)
# Very few strata have no PSU 1 or no PSU 2;
# so to keep the analysis simple, let's just drop those 3-4 rows.
# Then we can safely shuffle PSUs 1 vs 2 in each stratum separately
# to be in (proper) Training vs in Calibration sets
# as we repeat surveyCV-like approach for surveyConformal;
# and also we can repeatedly do SRS splits within PSUs 1+2 for iidConformal.

# Also,
colSums(table(df15saq$VARSTR, df15saq$VARPSU))
# Nearly all the ultimate units (indiv people) are in PSUs 1 and 2:
# around 10k each, while PSU 3 only have 1.7k total.
# Not ideal -- I prefer a larger test set -- but good enough for today's purposes.




#### Try an initial model on the SAQ respondents ####

# Do we use *all* of the SAQ data?

# We follow parts of this:
# https://github.com/yromano/cqr/blob/master/get_meps_data/meps_dataset_panel19_fy2015_reg.py
# They choose only panel 19 (data has both 19 and 20)
## filter(PANEL == 19)
# BUT the SAQ wts are for the whole 2015 year combined, NOT just for panel 19 alone
# so I will keep it as 2015 rather than split it by panel
# ...
# ALSO, when they chose the ___53 versions of many variables,
# that seems to be the version in round 5 of panel 19 or round 3 of panel 20
# (5th or 3rd time that person was in the survey)
# BUT since I am using SAQ, and that was ONLY used in round 4 of p19 or r2 of p20,
# we'll use the ___42 versions of those variables instead of ___53.
# ...
# ALSO, they used MANY categorical predictor variables.
# For today's paper, we will just pick a few of them for a simple predictive model.


# Redo data import, and subset to just the rows & variables we need:
df15raw = read.xport("h181.ssp")

df15 <- df15raw |>
  mutate(HHTOTD15 = na_if(HHTOTD15, -9),
         PCS42 = na_if(PCS42, -9),
         MCS42 = na_if(MCS42, -9),
         K6SUM42 = na_if(K6SUM42, -9),
         DIABDX = replace(DIABDX, DIABDX < 0, NA)) |>
  mutate(DIABDX01 = replace(DIABDX, DIABDX == 2, 0),   # 1 if has diabetes, 0 if not
         SEXM01 = replace(SEX, SEX == 2, 0),           # 1 if male, 0 if not
         INSCOV15PVT01 = replace(INSCOV15, INSCOV15 > 1, 0),          # 1 if private ins, 0 else
         INSCOV15PUB01 = replace(INSCOV15, INSCOV15 != 2, 1) - 1) |>  # 1 if public ins, 0 else
  mutate(utilization = OBTOTV15 + OPTOTV15 + ERTOT15 + IPNGTD15 + HHTOTD15)

# Subset just to respondents of the SAQ, self-administered questionnaire:
df15saq <- df15 |>
  filter(SAQWT15F > 0) |>
  select(OBTOTV15, OPTOTV15, ERTOT15, IPNGTD15, HHTOTD15,
         utilization,
         VARPSU, VARSTR, SAQWT15F,
         PCS42, MCS42, K6SUM42,
         AGE42X, DIABDX01, SEXM01,
         INSCOV15, INSCOV15PVT01, INSCOV15PUB01  ## INSCOV15 == 3 is uninsured
  )

df15saq |>
  summary()

svydes15saq <- svydesign(id = ~VARPSU, strata = ~VARSTR,
                         weights = ~SAQWT15F,
                         data = df15saq,
                         nest = TRUE)

# Drop the few strata with missing PSU 1 or 2.
strata_to_drop <- c(names(which(table(as.factor(df15saq$VARSTR)[df15saq$VARPSU == 1]) == 0)),
                    names(which(table(as.factor(df15saq$VARSTR)[df15saq$VARPSU == 2]) == 0)))
svydes15saq <- subset(svydes15saq, ! VARSTR %in% strata_to_drop)
df15saq <- df15saq |> filter(! VARSTR %in% strata_to_drop)




svyglm15saq <- svyglm(utilization ~ PCS42 + MCS42 + K6SUM42 +
                        AGE42X + DIABDX01 + SEXM01 + INSCOV15PVT01 + INSCOV15PUB01,
                      design = svydes15saq)
summary(svyglm15saq)

# At a first glance, higher utilization is positively assoc'ed with:
# lower PCS and MCS (worse health);
# lower K6 (better mental health?! --  but not stat signif despite massive n);
# older age, having diabetes, being female,
# and having private or public insurance.
#
# The only direction that surprises me was the negative sign for K6,
# but again, not signif;
# so it may just be a small & imprecise est (due to collinearity with MCS?)
plot(df15saq$MCS42, df15saq$K6SUM42)
cor(df15saq$MCS42, df15saq$K6SUM42, use = "complete.obs")




