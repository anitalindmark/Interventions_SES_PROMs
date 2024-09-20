################################################
# The function call used for the main analyses #
#  Exemplified for the outcome ADL-dependency  #
################################################

# Description of dataset used in the function call (note that the data are not publicly available):
# dataAnalysis = analysis dataset, with the following variables:
# ADL_3m = outcome, ADL-dependency at 3 month, with 1=dependent, 0=not dependent
# Mid = Indicator variable with 1=mid SES, 0=otherwise
# High = Indicator variable with 1=high SES, 0=otherwise
# smoker= mediator, 1=smoker, 0=non-smoker or unknown
# diabetes = mediator, 1=diabetes, 0=no diabetes
# antihypertensives = mediator, 1=prescribed antihypertensives, 0=no prescribed antihypertensives
# statins = mediator, 1=prescribed statins, 0=no prescribed statins
# atrial_fibrillation = mediator, 1=atrial fibrillation, 0=no atrial fibrillation
# hemorrhagic = mediator, 1=hemorrhagic stroke, 0=ischemic stroke
# lowered_consc = mediator, 1=lowered consciousness at hospital arrival, 0=not lowered consciousness
# sex = confounder, 1=male, 0=female
# age_centered = confounder, age in years centered around its mean
# age_centered_sq = age centered squared


# Load the required packages #
library(boot)
library(parallel)

# Since the estimation is time consuming parallelization was used #
# Here the clusters are prepared for parallelization #

cl <- makeCluster(7)
clusterExport(cl,"estimation_func")
clusterEvalQ(cl,library(boot))

# Effect estimation and bootstrap SEs #
# Run the estimation function through the boot function to obtain standard errors. The parallel option is used for parallelization. #

time <- proc.time() # Save start time

results <- boot(data = dataAnalysis, statistic = estimation_func, unexposed1="Mid", unexposed2="High", outcome="ADL_3m", 
                mediators=c("smoker","diabetes","antihypertensives","statins","atrial_fibrillation","hemorrhagic", "lowered_consc"),
                confounders=c("sex","age_centered","age_centered_sq"), 
                mcsim = 500, stype = "i", R = 1000, parallel = "snow", ncpus = 8,cl=cl) 

proc.time() - time # Time consumption

stopCluster(cl)

# Results #
Estimates <- results$t0
SE <- apply(results$t, 2, sd, na.rm=T)
# CI and p-values based on normal approximation. Remember to check the bootstrap distribution, if no approximately normal use e.g. percentile method.
CI_lower <- Estimates - qnorm(0.975)*SE
CI_upper <- Estimates + qnorm(0.975)*SE
p_value <- round(2*pnorm(abs(Estimates/SE), lower.tail = F),3)
prop_total_risk <- round(100*(c(Estimates[1:11]/Estimates[1],Estimates[12:22]/Estimates[12])), 1) 

cbind("Estimates"=round(Estimates*100,1),"CI lower"=round(100*CI_lower,1),"CI upper"=round(100*CI_upper,1),p_value, prop_total_risk)


