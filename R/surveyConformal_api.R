## Simulations, showing how conformal methods work
## in the case of repeated sampling from a finite population
## under several sampling designs:
## At least SRS and PPS; probably also stratified and cluster sampling


set.seed(20221119)

library(survey)
library(sampling)
library(matrixStats)
library(conformalInference)

# B_slow <- 1e2 # for sims where each rep is slow, use fewer reps when debugging
B_slow <- 1e3
B_fast <- 1e3 # for sims with fast reps, use more reps

## TODO: consistently use the same B=1e3 (nr of reps) everywhere,
## and rerun (though I know some of them are now 1e2 b/c they are so slow!)

## TODO: summarize coverage (means and CIs) visually,
## not just in tables and text alone

## TODO: add histograms of *achieved* probhat and length rep-by-rep,
## not just the average probhat and average length?
## At least find a FEW examples of where this is most informative:
## maybe pick one PI level, and compare hists for SRS vs PPS vs strat vs clus
## to show that clus has widest hists?

## TODO: better unequal-probs sampling WOR than just sample() ?
## I tried it two ways:
# y <- sample(apipop_clean$api00, 200, replace = FALSE,
#             prob = apipop_clean$probs)
# y <- apipop_clean$api00[sampling::UPtille(pik = 200 * apipop_clean$probs)]
## and the 2nd, Tille sampling, is what survey::election_pps apparently used,
## but it is VERY VERY SLOW. Way too slow for my purposes of 1000s of sim reps.
## So for now I will stick with sample() and hope it's OK.



#### Define survey-weighted quantile functions ####

## survey::svyquantile() seems to have a bug in it...
## https://stats.stackexchange.com/questions/576314/qrule-mathematic-interpolation-in-quantile-estimation-in-r-survey-package
## so let's define our own version,
## as well as a "conformal" counterpart
## which pads the n+1'th observation by its survey weight

## Both of these expect y and w to be vectors of same length,
## but p to be a scalar (and same for w_test).

## TODO: If I'm planning to use these for MANY test cases within each sim rep...
## or even just for two or more values of p...
## then I could probably write faster versions that sort (y,w) only once,
## instead of sapply()'ing same thing to many different p or w_test values.


q_wtd <- function(y, w, p) {
  w <- w[order(y)]
  y <- y[order(y)]
  q <- y[which.max(cumsum(w)/sum(w) >= p)]
}

q_wtd_conf <- function(y, w, p, w_test) {
  w <- w[order(y)]
  y <- y[order(y)]
  # Append new w_test at end, with "obs" value of Inf
  w <- c(w, w_test)
  y <- c(y, Inf)
  q <- y[which.max(cumsum(w)/sum(w) >= p)]
}




#### Clean up apipop ####


## Use the survey::apipop dataset as our finite pop,
## subset so that there are no missing values for the variables we'll use.
## For modeling, we'll use the svyglm() example:
##   svyglm(api00 ~ ell + meals + mobility)
## and for clusters and strata,
## we'll mimic apiclus1 and apistrat, so we need those IDs:
##   dnum, stype
## and for PPS, it makes sense to weight by enrollment:
##   enroll

data(api)

apipop_clean <- subset(apipop,
                       select = c(api00, ell, meals, mobility, dnum, stype, enroll))
summary(apipop_clean)

## Only `enroll` and `mobility` have any NAs, so drop those rows.
apipop_clean <- subset(apipop_clean,
                       !is.na(mobility) & !is.na(enroll))

## Add PPS weights, with prob proportional to enrollment
apipop_clean$probs <- apipop_clean$enroll / sum(apipop_clean$enroll)
apipop_clean$weights <- 1/apipop_clean$probs
## TODO: Should we rescale the weights so they all add to something specific?
sum(apipop_clean$probs)
sum(apipop_clean$weights)
sum(apipop_clean$enroll)
## Sum of probs is 1
## (these add up to 1, not add up to n,
##  so they are each prob of being included on *one* specific draw,
##  not overall prob of being included in *at least one* of n draws etc...)
## But sum of weights is neither N nor total enrollment.
## Think about whether we'd like to change that.

summary(apipop_clean)

## 6153 schools are left in the population.
(N <- nrow(apipop_clean))





#### The basic conformal quantile lemma works on SRS data ####

## Take a SRS of n=200.
## Use the standard Type 1 quantile.
## We get slight under-coverage.
B <- B_fast
system.time({
  probhats <- replicate(B, {
    y <- sample(apipop_clean$api00, 200)
    q_n.80 <- quantile(y, probs = 0.8, type = 1)
    q_n.95 <- quantile(y, probs = 0.95, type = 1)
    c(mean(apipop_clean$api00 <= q_n.80),
      mean(apipop_clean$api00 <= q_n.95))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## We are 95% confident that the true PI coverage is
## in the range (0.797, 0.798), less than our target of 0.8.
## 0.95: (0.946, 0.948)

SRS_sim_results <- data.frame(design = "SRS", padded = FALSE, weighted = FALSE, p = c(0.8, 0.95),
                                lo = rowMeans(probhats) - 2*SE,
                                hi = rowMeans(probhats) + 2*SE)

## Same thing, but now use the conformal-lemma corrected quantile.
## We get the desired coverage.
B <- B_fast
system.time({
  probhats <- replicate(B, {
    y <- sample(apipop_clean$api00, 200)
    q_np1.80 <- sort(y)[ceiling(0.8*201)]
    q_np1.95 <- sort(y)[ceiling(0.95*201)]
    c(mean(apipop_clean$api00 <= q_np1.80),
      mean(apipop_clean$api00 <= q_np1.95))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## Now we are 95% confident that the true PI coverage is
## in the range (0.801, 0.803), above our target of 0.8.

SRS_sim_results <- rbind(SRS_sim_results,
  data.frame(design = "SRS", padded = TRUE, weighted = FALSE, p = c(0.8,0.95),
             lo = rowMeans(probhats) - 2*SE,
             hi = rowMeans(probhats) + 2*SE))


## Summarize this section's results
SRS_sim_results$lo <- round(SRS_sim_results$lo, 3)
SRS_sim_results$hi <- round(SRS_sim_results$hi, 3)
SRS_sim_results




#### Basic conformal fails for PPS data, but weighted conformal works [api00] ####

## Take a PPS of n=200.
## Use the standard Type 1 quantile, IGNORING survey weights.
## We get slight under-coverage.
B <- B_fast
system.time({
  probhats <- replicate(B, {
    y <- sample(apipop_clean$api00, 200, replace = FALSE,
                prob = apipop_clean$probs)
    q_n.80 <- quantile(y, probs = 0.8, type = 1)
    q_n.95 <- quantile(y, probs = 0.95, type = 1)

    c(mean(apipop_clean$api00 <= q_n.80),
      mean(apipop_clean$api00 <= q_n.95))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## We are 95% confident that the true PI coverage is
## in the range (0.755, 0.759), less than our target of 0.8.
## For target 0.95, the range is (0.929, 0.931).

PPS_api00_sim_results <- data.frame(design = "PPS-api00",
                                     padded = FALSE, weighted = FALSE, p = c(0.8, 0.95),
                              lo = rowMeans(probhats) - 2*SE,
                              hi = rowMeans(probhats) + 2*SE)


## Same thing, but now use the conformal-lemma corrected quantile,
## but still IGNORING survey weights.
## We get under coverage again.
B <- B_fast
system.time({
  probhats <- replicate(B, {
    y <- sample(apipop_clean$api00, 200, replace = FALSE,
                prob = apipop_clean$probs)
    q_np1.80 <- sort(y)[ceiling(0.8*201)]
    q_np1.95 <- sort(y)[ceiling(0.95*201)]

    c(mean(apipop_clean$api00 <= q_np1.80),
      mean(apipop_clean$api00 <= q_np1.95))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## Now we are 95% confident that the true PI coverage is
## in the range (0.757, 0.761), still below our target of 0.8.
## For target 0.95, range is (0.935, 0.938).
## Changing from n to n+1 in the quantiles
## is NOT enough to fix coverage when we have unequal-probs sampling.

PPS_api00_sim_results <- rbind(PPS_api00_sim_results,
  data.frame(design = "PPS-api00", padded = TRUE, weighted = FALSE, p = c(0.8, 0.95),
             lo = rowMeans(probhats) - 2*SE,
             hi = rowMeans(probhats) + 2*SE))






## Take a PPS of n=200.
## Use the standard Type 1 quantile, but now SURVEY WEIGHTED.
## We get slight under-coverage, at least for 95% PIs.
B <- B_fast
system.time({
  probhats <- replicate(B, {
    i_samp <- sample(1:N, 200, replace = FALSE, prob = apipop_clean$probs)
    y <- apipop_clean$api00[i_samp]
    w <- apipop_clean$weights[i_samp]

    q_n.80 <- q_wtd(y = y, w = w, p = 0.8)
    q_n.95 <- q_wtd(y = y, w = w, p = 0.95)

    c(mean(apipop_clean$api00 <= q_n.80),
      mean(apipop_clean$api00 <= q_n.95))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## We are 95% confident that the true PI coverage is in the ranges
## (0.797, 0.801) for our target of 0.8, and
## (0.946, 0.949) for our target of 0.95.
## So we might be just barely on target for 0.8 but too low for 0.95.

PPS_api00_sim_results <- rbind(PPS_api00_sim_results,
   cbind(data.frame(design = "PPS-api00", padded = FALSE, weighted = TRUE, p = c(0.8, 0.95)),
         lo = rowMeans(probhats) - 2*SE,
         hi = rowMeans(probhats) + 2*SE))


## Now take a PPS of n=200.
## Use the standard Type 1 quantile, but SURVEY WEIGHTED, and CONFORMAL.
## We get slightly over the target coverage.
B <- B_slow  ## smaller B here b/c it takes longer, ~50 sec for B=1e2
system.time({
  probhats <- replicate(B, {
    i_samp <- sample(1:N, 200, replace = FALSE, prob = apipop_clean$probs)
    y <- apipop_clean$api00[i_samp]
    w <- apipop_clean$weights[i_samp]

    q_np1.80.vec <- sapply(apipop_clean$weights,
      function(x) q_wtd_conf(y = y, w = w, p = 0.80, w_test = x))
    q_np1.95.vec <- sapply(apipop_clean$weights,
      function(x) q_wtd_conf(y = y, w = w, p = 0.95, w_test = x))

    c(mean(apipop_clean$api00 <= q_np1.80.vec),
      mean(apipop_clean$api00 <= q_np1.95.vec))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## We are 95% confident that the true PI coverage is in the ranges
## (0.805, 0.817) for our target of 0.8, and
## (0.953, 0.959) for our target of 0.95.
## Now we have desired coverage for both targets.

PPS_api00_sim_results <- rbind(PPS_api00_sim_results,
   cbind(data.frame(design = "PPS-api00", padded = TRUE, weighted = TRUE, p = c(0.8, 0.95)),
         lo = rowMeans(probhats) - 2*SE,
         hi = rowMeans(probhats) + 2*SE))


## Summarize this section's results
PPS_api00_sim_results$lo <- round(PPS_api00_sim_results$lo, 3)
PPS_api00_sim_results$hi <- round(PPS_api00_sim_results$hi, 3)
PPS_api00_sim_results



#### Basic conformal fails for PPS data, but weighted conformal works [enroll] ####

## Same thing as above, but a more extreme example,
## where the samp probs are HIGHLY informative about the y-variable
## b/c the samp probs are directly proportional to `enroll`


## Take a PPS of n=200.
## Use the standard Type 1 quantile, IGNORING survey weights.
## We get major over-coverage.
B <- B_fast
system.time({
  probhats <- replicate(B, {
    y <- sample(apipop_clean$enroll, 200, replace = FALSE,
                prob = apipop_clean$probs)
    q_n.80 <- quantile(y, probs = 0.8, type = 1)
    q_n.95 <- quantile(y, probs = 0.95, type = 1)
    c(mean(apipop_clean$enroll <= q_n.80),
      mean(apipop_clean$enroll <= q_n.95))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## We are 95% confident that the true PI coverage is
## in the range (0.933, 0.934), WAY TOO HIGH for our target of 0.8.
## For 0.95, range is (0.987, 0.987), also too high.

PPS_enroll_sim_results <- data.frame(design = "PPS-enroll",
                                    padded = FALSE, weighted = FALSE, p = c(0.8, 0.95),
                                    lo = rowMeans(probhats) - 2*SE,
                                    hi = rowMeans(probhats) + 2*SE)



## Same thing, but now use the conformal-lemma corrected quantile,
## but still IGNORING survey weights.
## We get over-coverage again.
B <- B_fast
system.time({
  probhats <- replicate(B, {
    y <- sample(apipop_clean$enroll, 200, replace = FALSE,
                prob = apipop_clean$probs)
    q_np1.80 <- sort(y)[ceiling(0.8*201)]
    q_np1.95 <- sort(y)[ceiling(0.95*201)]
    c(mean(apipop_clean$enroll <= q_np1.80),
      mean(apipop_clean$enroll <= q_np1.95))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## Now we are 95% confident that the true PI coverage is
## in the range (0.934, 0.936), still WAY ABOVE our target of 0.8.
## For target 0.95, range is (0.988, 0.989).

PPS_enroll_sim_results <- rbind(PPS_enroll_sim_results,
   cbind(data.frame(design = "PPS-enroll", padded = TRUE, weighted = FALSE, p = c(0.8, 0.95)),
         lo = rowMeans(probhats) - 2*SE,
         hi = rowMeans(probhats) + 2*SE))



## Take a PPS of n=200.
## Use the standard Type 1 quantile, but now SURVEY WEIGHTED.
## We get slight under-coverage.
B <- B_fast
system.time({
  probhats <- replicate(B, {
    i_samp <- sample(1:N, 200, replace = FALSE, prob = apipop_clean$probs)
    y <- apipop_clean$enroll[i_samp]
    w <- apipop_clean$weights[i_samp]

    q_n.80 <- q_wtd(y = y, w = w, p = 0.8)
    q_n.95 <- q_wtd(y = y, w = w, p = 0.95)

    c(mean(apipop_clean$enroll <= q_n.80),
      mean(apipop_clean$enroll <= q_n.95))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## We are 95% confident that the true PI coverage is in the ranges
## (0.794, 0.797) for our target of 0.8, and
## (0.946, 0.948) for our target of 0.95.
## So we are just a bit too low for both targets.

PPS_enroll_sim_results <- rbind(PPS_enroll_sim_results,
  cbind(data.frame(design = "PPS-enroll", padded = FALSE, weighted = TRUE, p = c(0.8, 0.95)),
        lo = rowMeans(probhats) - 2*SE,
        hi = rowMeans(probhats) + 2*SE))


## Now take a PPS of n=200.
## Use the standard Type 1 quantile, but SURVEY WEIGHTED, and CONFORMAL.
## We get just-barely slight under-coverage.
B <- B_slow  ## smaller B here b/c it takes longer, ~50 sec for B=1e2
system.time({
  probhats <- replicate(B, {
    i_samp <- sample(1:N, 200, replace = FALSE, prob = apipop_clean$probs)
    y <- apipop_clean$enroll[i_samp]
    w <- apipop_clean$weights[i_samp]

    q_np1.80.vec <- sapply(apipop_clean$weights,
                           function(x) q_wtd_conf(y = y, w = w, p = 0.80, w_test = x))
    q_np1.95.vec <- sapply(apipop_clean$weights,
                           function(x) q_wtd_conf(y = y, w = w, p = 0.95, w_test = x))

    c(mean(apipop_clean$enroll <= q_np1.80.vec),
      mean(apipop_clean$enroll <= q_np1.95.vec))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## We are 95% confident that the true PI coverage is in the ranges
## (0.789, 0.799) for our target of 0.8, and
## (0.946, 0.950) for our target of 0.95.
## Now we have slight under-coverage for both targets,
## or possibly just-barely exact coverage for target 0.95.

PPS_enroll_sim_results <- rbind(PPS_enroll_sim_results,
  cbind(data.frame(design = "PPS-enroll", padded = TRUE, weighted = TRUE, p = c(0.8, 0.95)),
        lo = rowMeans(probhats) - 2*SE,
        hi = rowMeans(probhats) + 2*SE))



## Summarize this section's results
PPS_enroll_sim_results$lo <- round(PPS_enroll_sim_results$lo, 3)
PPS_enroll_sim_results$hi <- round(PPS_enroll_sim_results$hi, 3)
PPS_enroll_sim_results




#### Basic conformal fails for stratified data, but stratified conformal works ####

## Back to using api00
## But now instead of PPS, we stratify on stype
table(apistrat$stype)
## Use the same nr of obs per stratum: 100 E, 50 M, 50 H

## Note this is NOT proportional-allocation stratified sampling;
## if it were prop-alloc and we had 100 E schools,
## we'd only have 17 and 23 H and M schools
## instead of 50 of each
table(apipop$stype)/4421*100
## So this is a good example, b/c
## if we DID have prop alloc, we'd expect it to be very similar to SRS on avg,
## just slightly less variable from rep to rep.
## But here, we should get somewhat different results for SRS vs strat PIs.


## Using sampling::strata() which requires the data to be sorted by stratum...
apipop_clean_sort <- apipop_clean[order(apipop_clean$stype), ]

stratSizesSamp <- c(100, 50, 50)
stratSizesPop <- as.vector(table(apipop_clean_sort$stype))



## Take a stratified sample of n=200.
## Use the standard Type 1 quantile, IGNORING survey design.
## We get under-coverage.
B <- B_fast
system.time({
  probhats <- replicate(B, {
    y <- getdata(apipop_clean_sort,
                 sampling::strata(apipop_clean_sort, stratanames = "stype",
                        size = stratSizesSamp, method = "srswor"))$api00
    q_n.80 <- quantile(y, probs = 0.8, type = 1)
    q_n.95 <- quantile(y, probs = 0.95, type = 1)
    c(mean(apipop_clean$api00 <= q_n.80),
      mean(apipop_clean$api00 <= q_n.95))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## We are 95% confident that the true PI coverage is
## in the range (0.771, 0.774), too low for our target of 0.8.
## For 0.95, it's (0.935, 0.938).

strat_sim_results <- data.frame(design = "strat",
                                     padded = FALSE, weighted = FALSE, p = c(0.8, 0.95),
                                     lo = rowMeans(probhats) - 2*SE,
                                     hi = rowMeans(probhats) + 2*SE)



## Same thing, but now use the conformal-lemma corrected quantile,
## but still IGNORING survey design.
## We get under coverage again.
B <- B_fast
system.time({
  probhats <- replicate(B, {
    y <- getdata(apipop_clean_sort,
                 sampling::strata(apipop_clean_sort, stratanames = "stype",
                        size = stratSizesSamp, method = "srswor"))$api00
    q_np1.80 <- sort(y)[ceiling(0.8*201)]
    q_np1.95 <- sort(y)[ceiling(0.95*201)]
    c(mean(apipop_clean$api00 <= q_np1.80),
      mean(apipop_clean$api00 <= q_np1.95))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## Now we are 95% confident that the true PI coverage is
## in the range (0.757, 0.761), still below our target of 0.8.
## Changing from n to n+1 in the quantiles
## is NOT enough to fix coverage when we have unequal-probs sampling.

strat_sim_results <- rbind(strat_sim_results,
                                cbind(data.frame(design = "strat", padded = TRUE, weighted = FALSE, p = c(0.8, 0.95)),
                                      lo = rowMeans(probhats) - 2*SE,
                                      hi = rowMeans(probhats) + 2*SE))



## Take a stratified sample of n=200.
## Use the standard Type 1 quantile, but now STRATIFIED:
## calculate it separately in each stratum,
## and use the right quantile for the test case's stratum.
## We get slight under-coverage.
B <- B_slow
system.time({
  probhats <- replicate(B, {
    samp <- getdata(apipop_clean_sort,
                    sampling::strata(apipop_clean_sort, stratanames = "stype",
                        size = stratSizesSamp, method = "srswor"))
    y <- samp$api00
    st <- samp$stype
    q_n.80 <- sapply(unique(st), function(x) {
      tmp <- y[st == x]
      q <- sort(tmp)[ceiling(0.8*length(tmp))]
    })
    q_n.95 <- sapply(unique(st), function(x) {
      tmp <- y[st == x]
      q <- sort(tmp)[ceiling(0.95*length(tmp))]
    })

    c(mean(apipop_clean_sort$api00 <= rep(q_n.80, times = stratSizesPop)),
      mean(apipop_clean_sort$api00 <= rep(q_n.95, times = stratSizesPop)))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## We are 95% confident that the true PI coverage is in the ranges
## (0.781, 0.791) for our target of 0.8, and
## (0.935, 0.943) for our target of 0.95. So a bit too low overall.

strat_sim_results <- rbind(strat_sim_results,
                           cbind(data.frame(design = "strat", padded = FALSE, weighted = TRUE, p = c(0.8, 0.95)),
                                 lo = rowMeans(probhats) - 2*SE,
                                 hi = rowMeans(probhats) + 2*SE))



## Take a stratified sample of n=200.
## Use the standard Type 1 quantile, but now STRATIFIED and CONFORMAL.
## We get just over the target coverage!
B <- B_slow
system.time({
  probhats <- replicate(B, {
    samp <- getdata(apipop_clean_sort,
                    sampling::strata(apipop_clean_sort, stratanames = "stype",
                           size = stratSizesSamp, method = "srswor"))
    y <- samp$api00
    st <- samp$stype
    q_n.80 <- sapply(unique(st), function(x) {
      tmp <- y[st == x]
      q <- sort(tmp)[ceiling(0.8*(length(tmp)+1))]
    })
    q_n.95 <- sapply(unique(st), function(x) {
      tmp <- y[st == x]
      q <- sort(tmp)[ceiling(0.95*(length(tmp)+1))]
    })

    c(mean(apipop_clean_sort$api00 <= rep(q_n.80, times = stratSizesPop)),
      mean(apipop_clean_sort$api00 <= rep(q_n.95, times = stratSizesPop)))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## We are 95% confident that the true PI coverage is in the ranges
## (0.800, 0.811) for our target of 0.8, and
## (0.953, 0.960) for our target of 0.95. Great!

strat_sim_results <- rbind(strat_sim_results,
                           cbind(data.frame(design = "strat", padded = TRUE, weighted = TRUE, p = c(0.8, 0.95)),
                                 lo = rowMeans(probhats) - 2*SE,
                                 hi = rowMeans(probhats) + 2*SE))


## Summarize this section's results
strat_sim_results$lo <- round(strat_sim_results$lo, 3)
strat_sim_results$hi <- round(strat_sim_results$hi, 3)
strat_sim_results
## NOTE: "weighted" isn't quite right here --
## it's using strata but not using weights --
## but we'll leave it for now and possibly rename it later.




#### Basic conformal fails for cluster data, but cluster conformal works ####

## Now instead of PPS or strata, we cluster using dnum
table(apiclus1$dnum)
length(table(apiclus1$dnum))
## Use the same nr of clusters: 15
## (though it may correspond to VERY different nrs of schools per rep)

clusIDs <- unique(apipop_clean$dnum)

## Actually, it'd be great to take more than 15 clusters,
## b/c I'm doing "subsampling once"
## and 95% of 15 is 14.25 while 95% of 16 is 15.2,
## so it will ALWAYS round up to y=Inf...

## Can we simulate to figure out what nr of clusters
## tends to give n=200 total schools, on avg?
mean(replicate(1e4, sum(apipop_clean$dnum %in% sample(clusIDs, 15))))
mean(replicate(1e4, sum(apipop_clean$dnum %in% sample(clusIDs, 20))))
mean(replicate(1e4, sum(apipop_clean$dnum %in% sample(clusIDs, 24))))
## nclus=24 districts tends to give an average VERY close to n=200 schools,
## and .95*24 = 22.8 while .95*25 = 23.75
## so we'd end up using 23rd or 24th order statistic, NOT y=Inf.


## Take a SRSWOR cluster sample of nclus=24 school *districts*
##   (which may contain as few as 24 or as many as 100s of indiv schools total).
## Use the standard Type 1 quantile, IGNORING survey design.
## We get under-coverage.
B <- B_fast
system.time({
  probhats <- replicate(B, {
    i_samp <- apipop_clean$dnum %in% sample(clusIDs, 24)
    y <- apipop_clean[i_samp, ]$api00
    q_n.80 <- quantile(y, probs = 0.8, type = 1)
    q_n.95 <- quantile(y, probs = 0.95, type = 1)
    c(mean(apipop_clean$api00 <= q_n.80),
      mean(apipop_clean$api00 <= q_n.95))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## We are 95% confident that the true PI coverage is
## in the range (0.783, 0.796), too low for our target of 0.8.
## 0.95: (0.936, 0.941).

clus_sim_results <- data.frame(design = "clus",
                                padded = FALSE, weighted = FALSE, p = c(0.8, 0.95),
                                lo = rowMeans(probhats) - 2*SE,
                                hi = rowMeans(probhats) + 2*SE)




## Same thing, but now use the conformal-lemma corrected quantile,
## but still IGNORING survey design.
## We get under coverage again.
B <- B_fast
system.time({
  probhats <- replicate(B, {
    i_samp <- apipop_clean$dnum %in% sample(clusIDs, 24)
    y <- apipop_clean[i_samp, ]$api00
    q_np1.80 <- sort(y)[ceiling(0.8*(length(y)+1))]
    q_np1.95 <- sort(y)[ceiling(0.95*(length(y)+1))]
    c(mean(apipop_clean$api00 <= q_np1.80),
      mean(apipop_clean$api00 <= q_np1.95))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## Now we are 95% confident that the true PI coverage is
## in the range (0.787, 0.798),
## so it's still a bit too low for our target of 0.8.
## 0.95: (0.941, 0.946).

clus_sim_results <- rbind(clus_sim_results,
                           cbind(data.frame(design = "clus", padded = TRUE, weighted = FALSE, p = c(0.8, 0.95)),
                                 lo = rowMeans(probhats) - 2*SE,
                                 hi = rowMeans(probhats) + 2*SE))




## Take a cluster sample of nclus=24.
## Use the standard Type 1 quantile, but now using Dunn's "subsampling once":
## subsample to ONE school per cluster,
## treat that as an exchangeable sample...
## and use standard quantiles (w/o conformal correction!).
## We get slight under-coverage for 95%, but OK for 80%.
B <- B_fast
system.time({
  probhats <- replicate(B, {
    clusIDs_samp <- sample(clusIDs, 24)
    i_samp <- sapply(clusIDs_samp, function(x) {
      i_clus <- which(apipop_clean$dnum == x)
      if(length(i_clus) > 1) sample(i_clus, 1) else i_clus
    })
    y <- apipop_clean$api00[i_samp]
    q_n.80 <- sort(y)[ceiling(0.8*(length(y)))]
    q_n.95 <- sort(y)[ceiling(0.95*(length(y)))]

    c(mean(apipop_clean$api00 <= q_n.80),
      mean(apipop_clean$api00 <= q_n.95))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## We are 95% confident that the true PI coverage is in the ranges
## (0.804, 0.812) for our target of 0.8, and
## (0.920, 0.927) for our target of 0.95.
## So a bit too low for .95, but OK for .8.

clus_sim_results <- rbind(clus_sim_results,
                          cbind(data.frame(design = "clus", padded = FALSE, weighted = TRUE, p = c(0.8, 0.95)),
                                lo = rowMeans(probhats) - 2*SE,
                                hi = rowMeans(probhats) + 2*SE))



## Same as above, but WITH conformal correction:
## Take a cluster sample of nclus=24.
## Use the standard Type 1 quantile, but now using Dunn's "subsampling once":
## subsample to ONE school per cluster,
## treat that as an exchangeable sample...
## and use CONFORMAL quantiles.
## We get slight over-coverage.
B <- B_fast
system.time({
  probhats <- replicate(B, {
    clusIDs_samp <- sample(clusIDs, 24)
    i_samp <- sapply(clusIDs_samp, function(x) {
      i_clus <- which(apipop_clean$dnum == x)
      if(length(i_clus) > 1) sample(i_clus, 1) else i_clus
    })
    y <- apipop_clean$api00[i_samp]
    q_np1.80 <- sort(y)[ceiling(0.8*(length(y)+1))]
    q_np1.95 <- sort(y)[ceiling(0.95*(length(y)+1))]

    c(mean(apipop_clean$api00 <= q_np1.80),
      mean(apipop_clean$api00 <= q_np1.95))
  })
})
rowMeans(probhats)
(SE <- 2*rowSds(probhats)/sqrt(B))
t(sapply(1:2, function(i) rowMeans(probhats)[i] + c(-1, 1)*SE[i]))
## We are 95% confident that the true PI coverage is in the ranges
## (0.801, 0.810) for our target of 0.8, and
## (0.962, 0.966) for our target of 0.95.
## So OK but a bit too high for both.

clus_sim_results <- rbind(clus_sim_results,
                          cbind(data.frame(design = "clus", padded = TRUE, weighted = TRUE, p = c(0.8, 0.95)),
                                lo = rowMeans(probhats) - 2*SE,
                                hi = rowMeans(probhats) + 2*SE))


## Summarize this section's results
clus_sim_results$lo <- round(clus_sim_results$lo, 3)
clus_sim_results$hi <- round(clus_sim_results$hi, 3)
clus_sim_results
## NOTE: "weighted" isn't quite right here --
## it's using clusters but not using weights --
## but we'll leave it for now and possibly rename it later.








#### Basic conformal fails for MODELS with PPS data, but survey conformal works ####

## We COULD ask about effects of each of these 3 factors:
## (1) When splitting training data into "proper training" vs "calibration,"
##     should we use iid folds or surveyCV folds?
## (2) When fitting models on "proper training" set,
##     should we use weighted model-fitting or not?
## (3) When finding quantiles on calibration set and making PIs on test set,
##     should we use survey-weighted or iid conformal quantiles?

## To keep it simple, we'll only use PPS for now,
##   with wts prop. to enrollment at first (but later maybe to other things).
##   No strata or clusters here.
## That means (1) is off the table, since surveyCV for PPS data uses iid folds.

## Also, svy-weighted regression isn't built into Ryan's R package,
##   so let's wait with (2) for now as well.

## Just focus on (3): Does it matter whether or not we weight the conf quants?


## THIS TIME, ALSO COMPUTE PI LENGTHS.
## Clunky but simple way to get both the probhats
##   AND the PI lengths
## for every sim in one run of those sim settings:
## Each column is one rep of the sim, and
## the first N rows in each col are the 0/1s used for coverage,
## while the last N rows in each col are the indiv PI lengths.
ii_probs <- 1:N
ii_lens <- (N+1):(2*N)

## TODO: I guess that very rarely but not never...
## we'll get a set of weights
## for which conformal-padded 95th %ile is indeed Inf???
## and then we CAN'T calculate length directly;
## we'd have to talk about how often there's an Inf,
## and what's the avg length among the NON-Inf cases...???
## TODO: Has this been happening for the PI coverage above too,
## even when we were NOT looking at models?
## When we were JUST testing out the quantile lemma,
## we'd have ~200 obs, so it wouldn't be a problem to do 95%ile.
## ...
## AHA! Got it: We should do the same here.
## Use 400 obs for training overall,
## so that 200 each for "proper-training" and for calib.
## Then, calib set always has 200 and 95%ile won't be a problem.
## (Also make sure the resid-based weights aren't ridiculously extreme,
##  which could also cause a similar problem.)



## First,
## USE THE PPS PROBS that are proportional to `enroll` for now.


## Take a PPS of n=400.
## Use a 50/50 split of that into 200 each
## for "proper training" vs calibration sets.
## Use the ENTIRE finite pop as the test set.
## Now DO NOT use (survey)-weighted conformal quantiles for calib resids,
## just the vanilla unweighted quantiles.
## How often do the PIs capture true y-values on (full pop) test set?

## Once for 80% PIs...
B <- B_slow
runSim <- replicate(B, tryCatch({
  i_samp <- sample(1:N, 400, prob = apipop_clean$probs)
  apiPPS <- apipop_clean[i_samp, ]
  api_confpredsplit <- conformal.pred.split(
    x = as.matrix(apiPPS[c("ell", "meals", "mobility")]),
    y = apiPPS$api00,
    x0 = as.matrix(apipop_clean[c("ell", "meals", "mobility")]),
    train.fun = lm.funs()$train, predict.fun = lm.funs()$predict,
    alpha = 0.2,
    verb = FALSE)
  c(as.vector(api_confpredsplit$lo <= apipop_clean$api00 &
              apipop_clean$api00 <= api_confpredsplit$up),
    as.vector(api_confpredsplit$up - api_confpredsplit$lo))
}, error=function(err) NA))
hist(colMeans(runSim[ii_probs,], na.rm = TRUE))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_probs,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true PI coverage is in the range
## (0.802, 0.820) for our target of 0.8.
## Pretty good here despite lack of weighting.
model_covg_results <- data.frame(design = "PPS-enroll",
                                            weight_q = FALSE, weight_lm = FALSE,
                                            p = 0.8,
                                            lo = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) - 2*SE,
                                            hi = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + 2*SE)

hist(colMeans(runSim[ii_lens,], na.rm = TRUE))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_lens,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true mean PI length is in the range
## (193, 196) for our target of 0.8.
model_len_results <- data.frame(design = "PPS-enroll",
                                            weight_q = FALSE, weight_lm = FALSE,
                                            p = 0.8,
                                            lo = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) - 2*SE,
                                            hi = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + 2*SE)




## Now for 95% PIs...
B <- B_slow
runSim <- replicate(B, tryCatch({
  i_samp <- sample(1:N, 400, prob = apipop_clean$probs)
  apiPPS <- apipop_clean[i_samp, ]
  api_confpredsplit <- conformal.pred.split(
    x = as.matrix(apiPPS[c("ell", "meals", "mobility")]),
    y = apiPPS$api00,
    x0 = as.matrix(apipop_clean[c("ell", "meals", "mobility")]),
    train.fun = lm.funs()$train, predict.fun = lm.funs()$predict,
    alpha = 0.05,
    verb = FALSE)
  c(as.vector(api_confpredsplit$lo <= apipop_clean$api00 &
                apipop_clean$api00 <= api_confpredsplit$up),
    as.vector(api_confpredsplit$up - api_confpredsplit$lo))
}, error=function(err) NA))
hist(colMeans(runSim[ii_probs,], na.rm = TRUE))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_probs,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true PI coverage is in the range
## (0.954, 0.961) for our target of 0.95.
## Pretty good here despite lack of weighting.
model_covg_results <- rbind(model_covg_results,
                                       data.frame(design = "PPS-enroll",
                                            weight_q = FALSE, weight_lm = FALSE,
                                            p = 0.95,
                                            lo = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) - 2*SE,
                                            hi = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + 2*SE))
hist(colMeans(runSim[ii_lens,], na.rm = TRUE))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_lens,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true mean PI length is in the range
## (293, 301) for our target of 0.95.
model_len_results <- rbind(model_len_results,
                                       data.frame(design = "PPS-enroll",
                                                  weight_q = FALSE, weight_lm = FALSE,
                                                  p = 0.95,
                                                  lo = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) - 2*SE,
                                                  hi = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + 2*SE))




## Take a PPS of n=400.
## Use a 50/50 split of that into 200 each
## for "proper training" vs calibration sets.
## Use the ENTIRE finite pop as the test set.
## And DO use (survey)-WEIGHTED conformal quantiles for calib resids.
## How often do the PIs capture true y-values on (full pop) test set?

## Once for 80% PIs...
B <- B_slow
runSim <- replicate(B, tryCatch({
  i_samp <- sample(1:N, 400, prob = apipop_clean$probs)
  apiPPS <- apipop_clean[i_samp, ]
  api_confpredsplit <- conformal.pred.split(
    x = as.matrix(apiPPS[c("ell", "meals", "mobility")]),
    y = apiPPS$api00,
    x0 = as.matrix(apipop_clean[c("ell", "meals", "mobility")]),
    train.fun = lm.funs()$train, predict.fun = lm.funs()$predict,
    alpha = 0.2,
    verb = FALSE, w = apipop_clean$weights[c(i_samp, 1:N)])
  c(as.vector(api_confpredsplit$lo <= apipop_clean$api00 &
                apipop_clean$api00 <= api_confpredsplit$up),
    as.vector(api_confpredsplit$up - api_confpredsplit$lo))
}, error=function(err) NA))
hist(colMeans(runSim[ii_probs,], na.rm = TRUE))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_probs,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true PI coverage is in the range
## (0.792, 0.813) for our target of 0.8.
## Pretty good.
model_covg_results <- rbind(model_covg_results,
                                       data.frame(design = "PPS-enroll",
                                                  weight_q = TRUE, weight_lm = FALSE,
                                                  p = 0.8,
                                                  lo = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) - 2*SE,
                                                  hi = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + 2*SE))

hist(colMeans(runSim[ii_lens,], na.rm = TRUE))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_lens,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true mean PI length is in the range
## (187, 193) for our target of 0.8.
## Slightly NARROWER than for unweighted conformal quantiles!
model_len_results <- rbind(model_len_results,
                            data.frame(design = "PPS-enroll",
                                       weight_q = TRUE, weight_lm = FALSE,
                                       p = 0.8,
                                       lo = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) - 2*SE,
                                       hi = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + 2*SE))


## Now for 95% PIs...
B <- B_slow
runSim <- replicate(B, tryCatch({
  i_samp <- sample(1:N, 400, prob = apipop_clean$probs)
  apiPPS <- apipop_clean[i_samp, ]
  api_confpredsplit <- conformal.pred.split(
    x = as.matrix(apiPPS[c("ell", "meals", "mobility")]),
    y = apiPPS$api00,
    x0 = as.matrix(apipop_clean[c("ell", "meals", "mobility")]),
    train.fun = lm.funs()$train, predict.fun = lm.funs()$predict,
    alpha = 0.05,
    verb = FALSE, w = apipop_clean$weights[c(i_samp, 1:N)])
  c(as.vector(api_confpredsplit$lo <= apipop_clean$api00 &
                apipop_clean$api00 <= api_confpredsplit$up),
    as.vector(api_confpredsplit$up - api_confpredsplit$lo))
}, error=function(err) NA))
hist(colMeans(runSim[ii_probs,], na.rm = TRUE))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_probs,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true PI coverage is in the range
## (0.952, 0.962) for our target of 0.95.
## Pretty good!
model_covg_results <- rbind(model_covg_results,
                            data.frame(design = "PPS-enroll",
                                       weight_q = TRUE, weight_lm = FALSE,
                                       p = 0.95,
                                       lo = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) - 2*SE,
                                       hi = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + 2*SE))

hist(colMeans(runSim[ii_lens,], na.rm = TRUE))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_lens,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true mean PI length is in the range
## (289, 300) for our target of 0.95.
## Again, slightly narrower than for unweighted conformal quantiles!
model_len_results <- rbind(model_len_results,
                           data.frame(design = "PPS-enroll",
                                      weight_q = TRUE, weight_lm = FALSE,
                                      p = 0.95,
                                      lo = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) - 2*SE,
                                      hi = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + 2*SE))




## Try a DIFFERENT SET OF WEIGHTS.
## What if the samp probs were proportional to how hard it is to fit that point?
## That is, roughly proportional to residuals in the full-pop model?
## (actually use sqrt(abs(resids)) to fix huge resids having tiny weights,
##  and then add 1 to alleviate near-0 resids having massive weights)

apipop_clean_for_lm <- apipop_clean
apipop_clean_for_lm$resid <- resid(lm(api00 ~ ell + meals + mobility, data = apipop_clean_for_lm))
apipop_clean_for_lm$probs <- 1 + sqrt(abs(apipop_clean_for_lm$resid))
apipop_clean_for_lm$probs <- apipop_clean_for_lm$probs/sum(apipop_clean_for_lm$probs)
apipop_clean_for_lm$weights <- 1/apipop_clean_for_lm$probs



## Use vanilla UNWEIGHTED split conformal:

## Once for 80% PIs...
B <- B_slow
runSim <- replicate(B, tryCatch({
  i_samp <- sample(1:N, 400, prob = apipop_clean_for_lm$probs)
  apiPPS <- apipop_clean_for_lm[i_samp, ]
  api_confpredsplit <- conformal.pred.split(
    x = as.matrix(apiPPS[c("ell", "meals", "mobility")]),
    y = apiPPS$api00,
    x0 = as.matrix(apipop_clean_for_lm[c("ell", "meals", "mobility")]),
    train.fun = lm.funs()$train, predict.fun = lm.funs()$predict,
    alpha = 0.2,
    verb = FALSE)
  c(as.vector(api_confpredsplit$lo <= apipop_clean_for_lm$api00 &
                apipop_clean_for_lm$api00 <= api_confpredsplit$up),
    as.vector(api_confpredsplit$up - api_confpredsplit$lo))
}, error=function(err) NA))
hist(colMeans(runSim[ii_probs,], na.rm = TRUE))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_probs,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true PI coverage is in the range
## (0.872, 0.880) for our target of 0.8.
## WAY too high here.
model_covg_results <- rbind(model_covg_results,
                            data.frame(design = "PPS-resids",
                                       weight_q = FALSE, weight_lm = FALSE,
                                       p = 0.8,
                                       lo = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) - 2*SE,
                                       hi = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + 2*SE))
hist(colMeans(runSim[ii_lens,], na.rm = TRUE))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_lens,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true mean PI length is in the range
## (212, 217) for our target of 0.8.
## Wider than with probs based on enroll.
model_len_results <- rbind(model_len_results,
                           data.frame(design = "PPS-resids",
                                      weight_q = FALSE, weight_lm = FALSE,
                                      p = 0.8,
                                      lo = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) - 2*SE,
                                      hi = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + 2*SE))

## Now for 95% PIs...
B <- B_slow
runSim <- replicate(B, tryCatch({
  i_samp <- sample(1:N, 400, prob = apipop_clean_for_lm$probs)
  apiPPS <- apipop_clean_for_lm[i_samp, ]
  api_confpredsplit <- conformal.pred.split(
    x = as.matrix(apiPPS[c("ell", "meals", "mobility")]),
    y = apiPPS$api00,
    x0 = as.matrix(apipop_clean_for_lm[c("ell", "meals", "mobility")]),
    train.fun = lm.funs()$train, predict.fun = lm.funs()$predict,
    alpha = 0.05,
    verb = FALSE)
  c(as.vector(api_confpredsplit$lo <= apipop_clean_for_lm$api00 &
                apipop_clean_for_lm$api00 <= api_confpredsplit$up),
    as.vector(api_confpredsplit$up - api_confpredsplit$lo))
}, error=function(err) NA))
hist(colMeans(runSim[ii_probs,], na.rm = TRUE))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_probs,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true PI coverage is in the range
## (0.971, 0.975) for our target of 0.95.
## WAY too high here again.
model_covg_results <- rbind(model_covg_results,
                            data.frame(design = "PPS-resids",
                                       weight_q = FALSE, weight_lm = FALSE,
                                       p = 0.95,
                                       lo = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) - 2*SE,
                                       hi = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + 2*SE))
hist(colMeans(runSim[ii_lens,], na.rm = TRUE))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_lens,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true mean PI length is in the range
## (322, 333) for our target of 0.95.
## Again, MUCH wider than with probs based on enroll.
model_len_results <- rbind(model_len_results,
                           data.frame(design = "PPS-resids",
                                      weight_q = FALSE, weight_lm = FALSE,
                                      p = 0.95,
                                      lo = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) - 2*SE,
                                      hi = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + 2*SE))



## Now use WEIGHTED split conformal:

## Once for 80% PIs...
B <- B_slow
runSim <- replicate(B, tryCatch({
  i_samp <- sample(1:N, 400, prob = apipop_clean_for_lm$probs)
  apiPPS <- apipop_clean_for_lm[i_samp, ]
  api_confpredsplit <- conformal.pred.split(
    x = as.matrix(apiPPS[c("ell", "meals", "mobility")]),
    y = apiPPS$api00,
    x0 = as.matrix(apipop_clean_for_lm[c("ell", "meals", "mobility")]),
    train.fun = lm.funs()$train, predict.fun = lm.funs()$predict,
    alpha = 0.2,
    verb = FALSE, w = apipop_clean_for_lm$weights[c(i_samp, 1:N)])
  c(as.vector(api_confpredsplit$lo <= apipop_clean_for_lm$api00 &
                apipop_clean_for_lm$api00 <= api_confpredsplit$up),
    as.vector(api_confpredsplit$up - api_confpredsplit$lo))
}, error=function(err) NA))
hist(colMeans(runSim[ii_probs,], na.rm = TRUE))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_probs,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true PI coverage is in the range
## (0.794, 0.803) for our target of 0.8.
## Much better!
model_covg_results <- rbind(model_covg_results,
                            data.frame(design = "PPS-resids",
                                       weight_q = TRUE, weight_lm = FALSE,
                                       p = 0.8,
                                       lo = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) - 2*SE,
                                       hi = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + 2*SE))
hist(colMeans(runSim[ii_lens,], na.rm = TRUE))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_lens,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true mean PI length is in the range
## (178, 181) for our target of 0.8.
## Narrower than before!
model_len_results <- rbind(model_len_results,
                           data.frame(design = "PPS-resids",
                                      weight_q = TRUE, weight_lm = FALSE,
                                      p = 0.8,
                                      lo = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) - 2*SE,
                                      hi = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + 2*SE))


## Now for 95% PIs...
B <- B_slow
runSim <- replicate(B, tryCatch({
  i_samp <- sample(1:N, 400, prob = apipop_clean_for_lm$probs)
  apiPPS <- apipop_clean_for_lm[i_samp, ]
  api_confpredsplit <- conformal.pred.split(
    x = as.matrix(apiPPS[c("ell", "meals", "mobility")]),
    y = apiPPS$api00,
    x0 = as.matrix(apipop_clean_for_lm[c("ell", "meals", "mobility")]),
    train.fun = lm.funs()$train, predict.fun = lm.funs()$predict,
    alpha = 0.05,
    verb = FALSE, w = apipop_clean_for_lm$weights[c(i_samp, 1:N)])
  c(as.vector(api_confpredsplit$lo <= apipop_clean_for_lm$api00 &
                apipop_clean_for_lm$api00 <= api_confpredsplit$up),
    as.vector(api_confpredsplit$up - api_confpredsplit$lo))
}, error=function(err) NA))
hist(colMeans(runSim[ii_probs,], na.rm = TRUE))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_probs,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true PI coverage is in the range
## (0.948, 0.953) for our target of 0.95.
## Much better again!
model_covg_results <- rbind(model_covg_results,
                            data.frame(design = "PPS-resids",
                                       weight_q = TRUE, weight_lm = FALSE,
                                       p = 0.95,
                                       lo = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) - 2*SE,
                                       hi = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + 2*SE))
hist(colMeans(runSim[ii_lens,], na.rm = TRUE))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_lens,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true mean PI length is in the range
## (283, 289) for our target of 0.95.
## Again, narrower!
model_len_results <- rbind(model_len_results,
                           data.frame(design = "PPS-resids",
                                      weight_q = TRUE, weight_lm = FALSE,
                                      p = 0.95,
                                      lo = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) - 2*SE,
                                      hi = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + 2*SE))


## So... if we sampled with probs proportional to a variable
##   that's NOT super informative to the model,
##   then conformal works about as well with vs without weights.
## But if our probs WERE proportional to something INFORMATIVE,
##   and we can't use it as part of the model itself,
##   then weighted conformal DOES work substantially better.





## OK...
## So does it help if we DO use weighted lm() when the data are weighted?
## (and the weights are informative)

## Ran some checks writing this UNweighted lm function,
## to make sure we're getting the same results as earlier
train.fun.tmp <- function(x, y, out = NULL) {
  beta = lm.fit(x = cbind(1,x), y = y)$coefficients
  return(list(beta = beta, chol.R = NA))
}
## ...and indeed we were.

## So next I wrote this WEIGHTED lm function,
## so that we can fit weighted lm's
## & see how they relate to conformal PI covg and PI length.
train.fun.wtd.tmp <- function(x, y, out = NULL) {
  ## Hardcoded: first 3 cols are real X data, last col is the weights
  beta = lm.wfit(x = cbind(1,x[,1:3]), y = y, w = x[,4])$coefficients
  return(list(beta = c(beta,0), chol.R = NA))
}







## 80% PI, weighted lm, UNWEIGHTED conformal, props based on resids
B <- B_slow
runSim <- replicate(B, tryCatch({
  i_samp <- sample(1:N, 400, prob = apipop_clean_for_lm$probs)
  apiPPS <- apipop_clean_for_lm[i_samp, ]
  api_confpredsplit <- conformal.pred.split(
    x = as.matrix(apiPPS[c("ell", "meals", "mobility", "weights")]),
    y = apiPPS$api00,
    x0 = as.matrix(apipop_clean_for_lm[c("ell", "meals", "mobility", "weights")]),
    train.fun = train.fun.wtd.tmp,
    predict.fun = lm.funs()$predict,
    alpha = 0.2,
    verb = FALSE)
  c(as.vector(api_confpredsplit$lo <= apipop_clean_for_lm$api00 &
                apipop_clean_for_lm$api00 <= api_confpredsplit$up),
    as.vector(api_confpredsplit$up - api_confpredsplit$lo))
}, error=function(err) NA))
hist(colMeans(runSim[ii_probs,], na.rm = TRUE))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_probs,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true PI coverage is in the range
## (0.873, 0.880) for our target of 0.8.
## Way too high again!
model_covg_results <- rbind(model_covg_results,
                            data.frame(design = "PPS-resids",
                                       weight_q = FALSE, weight_lm = TRUE,
                                       p = 0.8,
                                       lo = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) - 2*SE,
                                       hi = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + 2*SE))
hist(colMeans(runSim[ii_lens,], na.rm = TRUE))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_lens,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true mean PI length is in the range
## (210, 214) for our target of 0.8.
## Wider than when weights were based on enroll;
## tiiiiny bit narrower than for unweighted lm.
model_len_results <- rbind(model_len_results,
                           data.frame(design = "PPS-resids",
                                      weight_q = FALSE, weight_lm = TRUE,
                                      p = 0.8,
                                      lo = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) - 2*SE,
                                      hi = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + 2*SE))

## 95% PI, weighted lm, UNWEIGHTED conformal, props based on resids
B <- B_slow
runSim <- replicate(B, tryCatch({
  i_samp <- sample(1:N, 400, prob = apipop_clean_for_lm$probs)
  apiPPS <- apipop_clean_for_lm[i_samp, ]
  api_confpredsplit <- conformal.pred.split(
    x = as.matrix(apiPPS[c("ell", "meals", "mobility", "weights")]),
    y = apiPPS$api00,
    x0 = as.matrix(apipop_clean_for_lm[c("ell", "meals", "mobility", "weights")]),
    train.fun = train.fun.wtd.tmp,
    predict.fun = lm.funs()$predict,
    alpha = 0.05,
    verb = FALSE)
  c(as.vector(api_confpredsplit$lo <= apipop_clean_for_lm$api00 &
                apipop_clean_for_lm$api00 <= api_confpredsplit$up),
    as.vector(api_confpredsplit$up - api_confpredsplit$lo))
}, error=function(err) NA))
hist(colMeans(runSim[ii_probs,], na.rm = TRUE))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_probs,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true PI coverage is in the range
## (0.972, 0.976) for our target of 0.95.
## Way too high again!
model_covg_results <- rbind(model_covg_results,
                            data.frame(design = "PPS-resids",
                                       weight_q = FALSE, weight_lm = TRUE,
                                       p = 0.95,
                                       lo = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) - 2*SE,
                                       hi = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + 2*SE))
hist(colMeans(runSim[ii_lens,], na.rm = TRUE))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_lens,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true mean PI length is in the range
## (328, 340) for our target of 0.95.
## Slightly wider than for unweighted conf with unweighted lm.
model_len_results <- rbind(model_len_results,
                           data.frame(design = "PPS-resids",
                                      weight_q = FALSE, weight_lm = TRUE,
                                      p = 0.95,
                                      lo = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) - 2*SE,
                                      hi = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + 2*SE))

## ...so using weighted lm() is NOT enough to fix issues.
## Even though it should give better-fitting models each time,
## it doesn't fix the distr. shift from calib-set to test-set.
## To fix THAT, we need to use weighted conformal quantiles.



## 80% PI, weighted lm, weighted conformal, props based on resids
B <- B_slow
runSim <- replicate(B, tryCatch({
  i_samp <- sample(1:N, 400, prob = apipop_clean_for_lm$probs)
  apiPPS <- apipop_clean_for_lm[i_samp, ]
  api_confpredsplit <- conformal.pred.split(
    x = as.matrix(apiPPS[c("ell", "meals", "mobility", "weights")]),
    y = apiPPS$api00,
    x0 = as.matrix(apipop_clean_for_lm[c("ell", "meals", "mobility", "weights")]),
    train.fun = train.fun.wtd.tmp,
    predict.fun = lm.funs()$predict,
    alpha = 0.2,
    verb = FALSE, w = apipop_clean_for_lm$weights[c(i_samp, 1:N)])
  c(as.vector(api_confpredsplit$lo <= apipop_clean_for_lm$api00 &
                apipop_clean_for_lm$api00 <= api_confpredsplit$up),
    as.vector(api_confpredsplit$up - api_confpredsplit$lo))
}, error=function(err) NA))
hist(colMeans(runSim[ii_probs,], na.rm = TRUE))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_probs,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true PI coverage is in the range
## (0.798, 0.808) for our target of 0.8.
## Much better again!
model_covg_results <- rbind(model_covg_results,
                            data.frame(design = "PPS-resids",
                                       weight_q = TRUE, weight_lm = TRUE,
                                       p = 0.8,
                                       lo = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) - 2*SE,
                                       hi = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + 2*SE))
hist(colMeans(runSim[ii_lens,], na.rm = TRUE))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_lens,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true mean PI length is in the range
## (175, 179) for our target of 0.8.
## A little narrower than when lm was unweighted!
model_len_results <- rbind(model_len_results,
                           data.frame(design = "PPS-resids",
                                      weight_q = TRUE, weight_lm = TRUE,
                                      p = 0.8,
                                      lo = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) - 2*SE,
                                      hi = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + 2*SE))


## 95% PI, weighted lm, weighted conformal, props based on resids
B <- B_slow
runSim <- replicate(B, tryCatch({
  i_samp <- sample(1:N, 400, prob = apipop_clean_for_lm$probs)
  apiPPS <- apipop_clean_for_lm[i_samp, ]
  api_confpredsplit <- conformal.pred.split(
    x = as.matrix(apiPPS[c("ell", "meals", "mobility", "weights")]),
    y = apiPPS$api00,
    x0 = as.matrix(apipop_clean_for_lm[c("ell", "meals", "mobility", "weights")]),
    train.fun = train.fun.wtd.tmp,
    predict.fun = lm.funs()$predict,
    alpha = 0.05,
    verb = FALSE, w = apipop_clean_for_lm$weights[c(i_samp, 1:N)])
  c(as.vector(api_confpredsplit$lo <= apipop_clean_for_lm$api00 &
                apipop_clean_for_lm$api00 <= api_confpredsplit$up),
    as.vector(api_confpredsplit$up - api_confpredsplit$lo))
}, error=function(err) NA))
hist(colMeans(runSim[ii_probs,], na.rm = TRUE))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_probs,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true PI coverage is in the range
## (0.949, 0.954) for our target of 0.95.
## Much better again!
model_covg_results <- rbind(model_covg_results,
                            data.frame(design = "PPS-resids",
                                       weight_q = TRUE, weight_lm = TRUE,
                                       p = 0.95,
                                       lo = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) - 2*SE,
                                       hi = mean(colMeans(runSim[ii_probs,], na.rm = TRUE)) + 2*SE))
hist(colMeans(runSim[ii_lens,], na.rm = TRUE))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE))
(SE = sd(colMeans(runSim[ii_lens,], na.rm = TRUE)) / sqrt(B))
mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + c(-1,1) * 2 * SE
## We are 95% confident that the true mean PI length is in the range
## (284, 292) for our target of 0.95.
## Similar (sliiightly higher) as to when lm was unweighted.
model_len_results <- rbind(model_len_results,
                           data.frame(design = "PPS-resids",
                                      weight_q = TRUE, weight_lm = TRUE,
                                      p = 0.95,
                                      lo = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) - 2*SE,
                                      hi = mean(colMeans(runSim[ii_lens,], na.rm = TRUE)) + 2*SE))


## Summarize this section's results
model_covg_results$lo <- round(model_covg_results$lo, 3)
model_covg_results$hi <- round(model_covg_results$hi, 3)
model_covg_results

model_len_results$lo <- round(model_len_results$lo, 1)
model_len_results$hi <- round(model_len_results$hi, 1)
model_len_results





#### Final summary of key results from all sections ####

# SRS:
# Without padding, we undercover,
# but with padding, we're on target.
SRS_sim_results

# PPS, predicting api00:
# With PPS, padding is not enough to fix undercoverage.
# Weighted quantiles help a lot,
# and padded wtd quantiles fix it completely.
PPS_api00_sim_results

# PPS, predicting enroll
#   (the same variable used to define PPS probs):
# Without weighting the quantiles,
#   we have extreme OVERcoverage,
#   and padding doesn't fix it at all.
# Weighted quantiles help a lot (turning into sliiight undercovg),
# and padded wtd quantiles fix it completely.
PPS_enroll_sim_results

# Strat:
# Ignoring strata undercovers, and padding isn't enough to fix it.
# Accounting for strata DOES fix it, esp with padding.
strat_sim_results

# Clus:
# Ignoring clusters undercovers (at least for 95% PIs),
# and padding isn't enough to fix it.
# Accounting for clusters DOES fix it,
# TODO: except for 95% PIs without padding???
#   why does that one's covg drop a little when we account for clusters,
#   instead of rising a little like all the others?
# I don't see any bugs in the code,
# and the effect persists even when I rerun with larger B.
# Must be something weird with clusters and "subsampling once"?
clus_sim_results

# Models and PIs based on regression residuals:
# Here we ALWAYS used padding; didn't bother trying without.
# For weights based on enroll (not one of the regr variables),
#   it makes very little difference whether or not we weight the quantiles,
#   whether to covg or to PI length.
# For more informative weights based on full-pop resids,
#   (a) there's big overcoverage with unwtd quantiles, but fixed by wtd quant,
#   (b) and PIs are much less wide with wtd quant than unwtd.
# Also in that setup,
#   weighting the LM model fit makes very little difference
#   whether to the covg or to the PI lengths.
#   So the benefit comes from wtd quantiles, NOT wtd models.
#   (I can imagine some scenarios where wtd models
#    could ALSO reduce the avg PI lengths, but this didn't seem to be one.)
model_covg_results
model_len_results



save(SRS_sim_results,
     PPS_api00_sim_results, PPS_enroll_sim_results,
     strat_sim_results, clus_sim_results,
     model_covg_results, model_len_results,
     file = "surveyConf_Api_SimResults_20221119.Rdata")
