## Data analysis code for Section 4.1 of
## "Design-based conformal prediction," Wieczorek (2023+)

## This file is streamlined to reproduce results from the paper.
## See also `surveyConformal_MEPS_EDA_R1.R` for some EDA
## and justification for some of the data filtering and modeling choices made below.




# We are using Medical Expenditure Panel Survey (MEPS) data,
# as in CQR paper by Y. Romano et al.:
# https://github.com/yromano/cqr/blob/master/get_meps_data/download_data.R

# MEPS requires the following Data Use Agreement when writing up results:
# https://meps.ahrq.gov/data_stats/data_use.jsp



#### Data prep ####


library(foreign)
library(dplyr)
library(survey)
library(matrixStats)

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


df15saqTest <- df15saq |>
  filter(VARPSU == 3)



#### Run survey-weighted simulations ####

## Replicate many times:
##   redefine df15saqTrain and df15saqCalib by *different* sets
##   of PSU 1 and PSU 2 in each;
##   train a model and calibrate the conformal quantiles of residuals;
##   find PI length and PI coverage on the test set.

set.seed(20230721)
B = 100

system.time({

  df15_sim_results_svyWtd <- replicate(B, {
    # Choose half the remaining strata.
    # Assign PSU 1 in those and PSU 2 in the others to Train,
    # and vice versa to Calib.
    strata_avail <- unique(df15saq$VARSTR[df15saq$VARPSU != 3])
    strata_Train_PSU_1 <- sample(strata_avail, floor(length(strata_avail)/2))
    strata_Train_PSU_2 <- setdiff(strata_avail, strata_Train_PSU_1)

    df15saqTrain <- df15saq |>
      filter((VARSTR %in% strata_Train_PSU_1 & VARPSU == 1) |
               (VARSTR %in% strata_Train_PSU_2 & VARPSU == 2))
    df15saqCalib <- df15saq |>
      filter((VARSTR %in% strata_Train_PSU_1 & VARPSU == 2) |
               (VARSTR %in% strata_Train_PSU_2 & VARPSU == 1))
    svydes15saqTrain <- subset(svydes15saq,
                               (VARSTR %in% strata_Train_PSU_1 & VARPSU == 1) |
                                 (VARSTR %in% strata_Train_PSU_2 & VARPSU == 2))
    svydes15saqCalib <- subset(svydes15saq,
                               (VARSTR %in% strata_Train_PSU_1 & VARPSU == 2) |
                                 (VARSTR %in% strata_Train_PSU_2 & VARPSU == 1))



    ## ?surveyoptions
    # When we use a dataset with one PSU per stratum,
    # there's no way to calculate variances.
    # That is OK for our purposes,
    # but {survey} will default to throwing errors unless we change options:
    # remove, adjust, or average all do diff things,
    # but all give the same point ests of betahats,
    # so any of them is fine for today's purposes.
    options(survey.lonely.psu = "remove")

    svyglm15saqTrain <- svyglm(utilization ~ PCS42 + MCS42 + K6SUM42 +
                                 AGE42X + DIABDX01 + SEXM01 + INSCOV15PVT01 + INSCOV15PUB01,
                               design = svydes15saqTrain, na.action = na.exclude)
    # Alas -- when we use predict.svyglm, it DOESN'T pad the result with NAs as needed!
    # When X's are NA and therefore Yhat is too (which is fine),
    # it just silently drops NAs (not fine!!!)
    # instead of returning a vector of same length as nrow(newdata).
    # BUT
    # if we FORCE it to use predict.glm() instead, we can fix that:
    df15saqCalib$yhat <- predict.glm(svyglm15saqTrain, newdata = df15saqCalib,
                                     na.action = na.pass)
    df15saqCalib$resid <- df15saqCalib$utilization - df15saqCalib$yhat

    df15saqTest$yhat <- predict.glm(svyglm15saqTrain, newdata = df15saqTest,
                                    na.action = na.pass)
    df15saqTest$resid <- df15saqTest$utilization - df15saqTest$yhat



    tmp <- sapply(1:nrow(df15saqTest), function(ii) {
      y <- abs(df15saqCalib$resid)
      w <- df15saqCalib$SAQWT15F
      isna_y <- is.na(y)
      y <- y[!isna_y]
      w <- w[!isna_y]
      length(y); length(w)
      # Reorder from smallest to largest y,
      # and keep weights in same order as their corresponding y's
      # (NOT from smallest to largest w)
      w <- w[order(y)]
      y <- y[order(y)]
      # Append new w_np1 at end, with "observed" y-value of Inf
      y <- c(y, Inf)
      w <- c(w, df15saqTest$SAQWT15F[ii])

      # Now use this padded version of w to find quantile
      q_n.80 <- y[which.max(cumsum(w)/sum(w) >= 0.8)]
      q_n.90 <- y[which.max(cumsum(w)/sum(w) >= 0.9)]
      q_n.95 <- y[which.max(cumsum(w)/sum(w) >= 0.95)]
      q_n.99 <- y[which.max(cumsum(w)/sum(w) >= 0.99)]

      out <- c(q_n.80, q_n.90, q_n.95, q_n.99,
               abs(df15saqTest$resid)[ii] <= c(q_n.80, q_n.90, q_n.95, q_n.99))
      return(out)
    })

    # rows 1:4 are the PI "MOEs", so their means will be avg PI half-width;
    # rows 5:8 are whether or not truth was in the MOE, so their means will be coverage.
    rownames(tmp) <- c("q_n.80", "q_n.90", "q_n.95", "q_n.99",
                       "In_q_n.80", "In_q_n.90", "In_q_n.95", "In_q_n.99")

    return(rowWeightedMeans(tmp, w = df15saqTest$SAQWT15F, na.rm = TRUE))
  })

}) # end system.time()


#  Means across all B reps
(svyWtd_means <- rowMeans(df15_sim_results_svyWtd))

# Approx 95% MOEs
(svyWtd_MOEs <- 2*rowSds(df15_sim_results_svyWtd, useNames = TRUE)/sqrt(B))






#### Run SRS simulations ####


# REPEAT the above analysis, except...
# * using iid splits instead of using PSUs within strata for the splits...
# * using SRS to fit the model, instead of the real survey design...
# * ignoring survey weights when finding the conformal quantiles...
# * and estimating test-set coverage two ways: with and without survey weights.

df15saq_NotPSU3 <- df15saq |> filter(VARPSU != 3)
svydes15saq_NotPSU3 <- subset(svydes15saq, VARPSU != 3)

system.time({

  df15_sim_results_srs <- replicate(B, {
    # Choose half the remaining data in a random partition.
    N <- nrow(df15saq_NotPSU3)
    i_train <- sample(N, floor(N/2))
    i_calib <- setdiff(1:N, i_train)

    df15saqTrain <- df15saq_NotPSU3[i_train, ]
    df15saqCalib <- df15saq_NotPSU3[i_calib, ]

    options(survey.lonely.psu = "remove")

    df15saqTrain_SRS <- svydesign(ids = ~1, weights = ~1, data = df15saqTrain)
    svyglm15saqTrain <- svyglm(utilization ~ PCS42 + MCS42 + K6SUM42 +
                                             AGE42X + DIABDX01 + SEXM01 +
                                             INSCOV15PVT01 + INSCOV15PUB01,
                               design = df15saqTrain_SRS, na.action = na.exclude)
    df15saqCalib$yhat <- predict.glm(svyglm15saqTrain, newdata = df15saqCalib,
                                     na.action = na.pass)
    df15saqCalib$resid <- df15saqCalib$utilization - df15saqCalib$yhat

    df15saqTest$yhat <- predict.glm(svyglm15saqTrain, newdata = df15saqTest,
                                    na.action = na.pass)
    df15saqTest$resid <- df15saqTest$utilization - df15saqTest$yhat



    tmp <- sapply(1:nrow(df15saqTest), function(ii) {
      y <- abs(df15saqCalib$resid)
      isna_y <- is.na(y)
      y <- y[!isna_y]
      # Reorder from smallest to largest y
      y <- y[order(y)]
      # Append at end a new "observed" y-value of Inf
      y <- c(y, Inf)

      # Now use this padded version to find conformal quantile, UNWEIGHTED
      Uw <- rep(1, length(y))
      Uq_n.80 <- y[which.max(cumsum(Uw)/sum(Uw) >= 0.8)]
      Uq_n.90 <- y[which.max(cumsum(Uw)/sum(Uw) >= 0.9)]
      Uq_n.95 <- y[which.max(cumsum(Uw)/sum(Uw) >= 0.95)]
      Uq_n.99 <- y[which.max(cumsum(Uw)/sum(Uw) >= 0.99)]

      out <- c(Uq_n.80, Uq_n.90, Uq_n.95, Uq_n.99,
               abs(df15saqTest$resid)[ii] <= c(Uq_n.80, Uq_n.90, Uq_n.95, Uq_n.99))
      return(out)
    })

    # rows 1:4 are the PI "MOEs", so their means will be avg PI half-width;
    # rows 5:8 are whether or not truth was in the MOE, so their means will be coverage.
    rownames(tmp) <- c("Uq_n.80", "Uq_n.90", "Uq_n.95", "Uq_n.99",
                       "In_Uq_n.80", "In_Uq_n.90", "In_Uq_n.95", "In_Uq_n.99")

    # Return a simple mean across the test set,
    # as well as a weighted mean (that should generalize to pop better)
    return(c(rowMeans(tmp, na.rm = TRUE),
             rowWeightedMeans(tmp, w = df15saqTest$SAQWT15F, na.rm = TRUE)))
  })

}) # end system.time()

df15_sim_results_srs_unwtdRowMeans <- df15_sim_results_srs[1:8, ]
df15_sim_results_srs_wtdRowMeans <- df15_sim_results_srs[9:16, ]


(svyWtd_means <- rowMeans(df15_sim_results_svyWtd))

# Approx 95% MOEs
(svyWtd_MOEs <- 2*rowSds(df15_sim_results_svyWtd, useNames = TRUE)/sqrt(B))


# Means across all B reps
(srs_unwtdTest_means <- rowMeans(df15_sim_results_srs_unwtdRowMeans))
(srs_wtdTest_means <- rowMeans(df15_sim_results_srs_wtdRowMeans))

# Approx 95% MOEs
(srs_unwtdTest_MOEs <- 2*rowSds(df15_sim_results_srs_unwtdRowMeans, useNames = TRUE)/sqrt(B))
(srs_wtdTest_MOEs <- 2*rowSds(df15_sim_results_srs_wtdRowMeans, useNames = TRUE)/sqrt(B))

# All the Uq_n.xx values SHOULD indeed be the same across wtd vs unwtd RowMeans,
# b/c they are based on unwtd conformal -- so every obs (in that rep) has same PI length,
# so it doesn't matter whether we take a wtd or unwtd mean of identical PI lengths.
# However, the In_Uq_n.xx values (coverage values) do change for wtd vs unwtd RowMeans.




#### Report results ####


# The coverages and PI lengths below are reported as approx 95% CI endpoints.



# (1) IF YOU'VE IGNORED WTS, YOU'LL THINK YOU'VE GOT THIS:
# NOT-correctly-gen'lizing test-set ests
# (i.e., what you'd THINK you're getting, based on unwtd means on test set)
# from when we use iid (unwtd) for all: lm(), splits, AND quantiles

# Coverage:
srs_unwtdTest_means[5:8] - srs_unwtdTest_MOEs[5:8]
srs_unwtdTest_means[5:8] + srs_unwtdTest_MOEs[5:8]

# PI lengths:
2*(srs_unwtdTest_means[1:4] - srs_unwtdTest_MOEs[1:4])
2*(srs_unwtdTest_means[1:4] + srs_unwtdTest_MOEs[1:4])


# (2) IF YOU'VE IGNORED WTS EVERYWHERE, INCL IN TEST SET,
#     YOU'LL ACTUALLY HAVE THIS W/O REALIZING IT:
# Correctly-gen'lizing test-set ests
# (i.e., what you're ACTUALLY likely getting, based on WTD means on test set)
# from when we use srs/iid (unwtd) for all training and calibration:
# lm(), splits, AND quantiles

# Coverage:
srs_wtdTest_means[5:8] - srs_wtdTest_MOEs[5:8]
srs_wtdTest_means[5:8] + srs_wtdTest_MOEs[5:8]

# PI lengths:
2*(srs_wtdTest_means[1:4] - srs_wtdTest_MOEs[1:4])
2*(srs_wtdTest_means[1:4] + srs_wtdTest_MOEs[1:4])

# Comparing (2) to (1), we notice that:
# (a) mean lengths won't change at all b/c they are constant across obs anyway,
# so wtd or unwtd mean of a constant is the same;
# (b) covg is actually a little higher than you thought in (1) (slight overcovg)
# although not much of a difference.


# (3) IF YOU DO IN FACT USE WTS, YOU WOULD ACTUALLY HAVE THIS:
# Correctly-gen'lizing test-set estimates of PI length and covg
# from when we use svy-wtd for all: lm(), splits, AND quantiles

# Coverage:
svyWtd_means[5:8] - svyWtd_MOEs[5:8]
svyWtd_means[5:8] + svyWtd_MOEs[5:8]

# PI lengths:
2*(svyWtd_means[1:4] - svyWtd_MOEs[1:4])
2*(svyWtd_means[1:4] + svyWtd_MOEs[1:4])

# Comparing (3) to (2), we notice that:
# (a) mean PI lengths are shorter now: by about 1-3 units for 80% to 95% PIs,
# and by a LOT for 99% PIs;
# (b) covg is very slightly lower again, typically closer to nominal than in (2).

