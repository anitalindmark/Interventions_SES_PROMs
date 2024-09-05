###############################################################
## Function for Monte Carlo simulation to estimate the       ##
##                interventional effects                     ##
##                                                           ##
## Note that the function is written specifically for the    ##
## application in the article and needs to be adapted if it  ##
## is to be used in other applications.                      ##
###############################################################

# Input parameters: 
# dat = data frame, 
# ind = indices for bootstrap, 
# unexposed1 = name of the indicator variable for the first unexposed contrast, 
# unexposed2 = name of the indicator variable for the second unexposed contrast, 
# outcome = name of the outcome variable,
# mediators = character vector with the names of the mediators,
# confounders = character vector with the names of the confounders,
# mcsim = number of Monte Carlo simulations,

estimation_func <- function(dat, ind = 1:nrow(dat), unexposed1, unexposed2, outcome, mediators, confounders, mcsim = 500)
{
  
  # Rename all variables and prepare dataset
  dat$Mid <- dat[, unexposed1]
  dat$High <- dat[, unexposed2]
  dat$Y <- dat[, outcome]
  dat$M1 <- dat[, mediators[1]]  # smoking
  dat$M21 <- dat[, mediators[2]] # diabetes
  dat$M22 <- dat[, mediators[3]] # antihypertensives
  dat$M23 <- dat[, mediators[4]] # statins
  dat$M3  <- dat[, mediators[5]] # atrial fibrillation
  dat$M41 <- dat[, mediators[6]] # hemorrhagic stroke
  dat$M42 <- dat[, mediators[7]] # lowered consciousness
  
  dat<-dat[, c("Mid","High", "M1", "M21", "M22", "M23","M3", "M41","M42", "Y", confounders)] 
  
  # Take bootstrap sample
  data <- dat[ind, ]
  
  # Data subsets, one for the low vs. mid SES contrast, one for the low vs. high contrast.
  dat2 <- data[which(data$High==0),] # low vs. mid
  n <- nrow(dat2)
  
  dat3 <- data[which(data$Mid==0),] # low vs. high
  n3 <- nrow(dat3)
  
  #Set flag to capure bootstrap samples to reject
  flag <- FALSE 
  
  # Prepare confounders for formulae
  confs <- paste(confounders, collapse = "+")
 
  # Pre-allocate vectors for risk contrasts 
  risk_joint_int <- risk_unexposed <- risk_exposed <- risk_M1_int <- risk_M2_int <- risk_M3_int <- risk_M4_int  <- rep(NA,mcsim)
  
  
  ### Step A: Model fitting ###
  
  ## Marginal models (not given other mediators (mediator groups)) ##
  
  # 1.a. M1 model
  fit_M1_marg <- glm(as.formula(paste("M1 ~ Mid + High + ", confs)), data = data, family = binomial,
                y = FALSE, model = FALSE) #
  if((!fit_M1_marg$converged)|any(is.na(fit_M1_marg$coefficients))) flag <- TRUE
  
  # 2.a. M2 models
  fit_M21_marg <- glm(as.formula(paste("M21 ~ Mid + High + ", confs)), data = data, family = binomial,
                      y = FALSE, model = FALSE) 
  if((!fit_M21_marg$converged)|any(is.na(fit_M21_marg$coefficients))) flag <- TRUE
  
  fit_M22_marg <- glm(as.formula(paste("M22 ~ (Mid + High + M21)^2 - Mid:High +", confs)), data = data, family = binomial,
                      y = FALSE, model = FALSE) #
  if((!fit_M22_marg$converged)|any(is.na(fit_M22_marg$coefficients))) flag <- TRUE
  
  fit_M23_marg <- glm(as.formula(paste("M23 ~ (Mid + High + M21 + M22)^2 - Mid:High +", confs)), data = data, family = binomial,
                      y = FALSE, model = FALSE) #
  if((!fit_M23_marg$converged)|any(is.na(fit_M23_marg$coefficients))) flag <- TRUE
  
  # 3.a. M3 model
  
  fit_M3_marg <- glm(as.formula(paste("M3 ~ Mid + High + ",confs)), data = data, family = binomial,
                     y = FALSE, model = FALSE)
  if((!fit_M3_marg$converged)|any(is.na(fit_M3_marg$coefficients))) flag <- TRUE
  
  # 4.a. M4 models
  
  fit_M41_marg <- glm(as.formula(paste("M41 ~ Mid + High + ",confs)), data = data, family = binomial,
                      y = FALSE, model = FALSE)
  if((!fit_M41_marg$converged)|any(is.na(fit_M41_marg$coefficients))) flag <- TRUE
  
  fit_M42_marg <- glm(as.formula(paste("M42 ~ (Mid + High + M41)^2 - Mid:High + ",confs)), data = data, family = binomial,
                     y = FALSE, model = FALSE)
  if((!fit_M42_marg$converged)|any(is.na(fit_M42_marg$coefficients))) flag <- TRUE
  
  
  ## Models conditional on one "preceding" mediator (mediator group) ##
  
  # 2.b. M2 models given M1
  fit_M21_cond <- glm(as.formula(paste("M21 ~ Mid + High + M1 + ", confs)), data = data, family = binomial,
                      y = FALSE, model = FALSE)
  if((!fit_M21_cond$converged)|any(is.na(fit_M21_cond$coefficients))) flag <- TRUE
  
  fit_M22_cond <- glm(as.formula(paste("M22 ~ (Mid + High + M1 + M21)^2 - Mid:High + ", confs)), data = data, family = binomial,
                      y = FALSE, model = FALSE)
  if((!fit_M22_cond$converged)|any(is.na(fit_M22_cond$coefficients))) flag <- TRUE
  
  fit_M23_cond <- glm(as.formula(paste("M23 ~ (Mid + High + M1 + M21 + M22)^2 - Mid:High + ", confs)), data = data, family = binomial,
                      y = FALSE, model = FALSE)
  if((!fit_M22_cond$converged)|any(is.na(fit_M22_cond$coefficients))) flag <- TRUE
  
  # 3.b. M3 model given M1
  fit_M3_condM1 <- glm(as.formula(paste("M3 ~ (Mid + High + M1)^2 - Mid:High  + ",confs)), data = data, family = binomial,
                        y = FALSE, model = FALSE)
  if((!fit_M3_condM1$converged)|any(is.na(fit_M3_condM1$coefficients))) flag <- TRUE
  
  # 3.c. M3 model given M2
  fit_M3_condM2 <- glm(as.formula(paste("M3 ~ (Mid + High + M21 + M22 + M23)^2 - Mid:High + ",confs)), data = data, family = binomial,
                      y = FALSE, model = FALSE)
  if((!fit_M3_condM2$converged)|any(is.na(fit_M3_condM2$coefficients))) flag <- TRUE
  
  
  ## Models conditional on two "preceding" mediators (mediator groups) ##
  
  # 3.d. M3 model given M1 and M2
  fit_M3_condM1M2 <- glm(as.formula(paste("M3 ~ (Mid + High + M1 + M21 + M22 + M23 )^2 - Mid:High +  ",confs)), 
                          data = data, family = binomial, y = FALSE, model = FALSE)
  if((!fit_M3_condM1M2$converged)|any(is.na(fit_M3_condM1M2$coefficients))) flag <- TRUE
  
  # 4.b. M4 models given M1 and M2
  fit_M41_condM1M2 <- glm(as.formula(paste("M41 ~ (Mid + High + M1 + M21 + M22 + M23 )^2 - Mid:High +  ",confs)), 
                         data = data, family = binomial, y = FALSE, model = FALSE)
  if((!fit_M41_condM1M2$converged)|any(is.na(fit_M41_condM1M2$coefficients))) flag <- TRUE
  
  fit_M42_condM1M2 <- glm(as.formula(paste("M42 ~ (Mid + High + M1 + M21 + M22 + M23 + M41)^2 - Mid:High + ",confs)), 
                         data = data, family = binomial, y = FALSE, model = FALSE)
  if((!fit_M42_condM1M2$converged)|any(is.na(fit_M42_condM1M2$coefficients))) flag <- TRUE
  
  # 4.c. M4 models given M1 and M3
  fit_M41_condM1M3 <- glm(as.formula(paste("M41 ~ (Mid + High + M1 + M3 )^2 - Mid:High +  ",confs)), 
                          data = data, family = binomial, y = FALSE, model = FALSE)
  if((!fit_M41_condM1M3$converged)|any(is.na(fit_M41_condM1M3$coefficients))) flag <- TRUE
  
  fit_M42_condM1M3 <- glm(as.formula(paste("M42 ~ (Mid + High + M1 + M3 + M41)^2 - Mid:High + ",confs)), 
                          data = data, family = binomial, y = FALSE, model = FALSE)
  if((!fit_M42_condM1M3$converged)|any(is.na(fit_M42_condM1M3$coefficients))) flag <- TRUE
  
  # 4.d. M4 models given M2 and M3
  fit_M41_condM2M3 <- glm(as.formula(paste("M41 ~ (Mid + High + M21 + M22 + M23 + M3 )^2 - Mid:High +  ",confs)), 
                          data = data, family = binomial, y = FALSE, model = FALSE)
  if((!fit_M41_condM2M3$converged)|any(is.na(fit_M41_condM2M3$coefficients))) flag <- TRUE
  
  fit_M42_condM2M3 <- glm(as.formula(paste("M42 ~ (Mid + High + M21 + M22 + M23 + M3 + M41)^2 - Mid:High + ",confs)), 
                          data = data, family = binomial, y = FALSE, model = FALSE)
  if((!fit_M42_condM2M3$converged)|any(is.na(fit_M42_condM2M3$coefficients))) flag <- TRUE
  
  
  ## Models conditional on three "preceding" mediators (mediator groups) ##
  
  # 4.e. M4 models given M1, M2 and M3
  fit_M41_condM1M2M3 <- glm(as.formula(paste("M41 ~ (Mid + High + M1 + M21 + M22 + M23 + M3 )^2 - Mid:High +  ",confs)), 
                          data = data, family = binomial, y = FALSE, model = FALSE)
  if((!fit_M41_condM1M2M3$converged)|any(is.na(fit_M41_condM1M2M3$coefficients))) flag <- TRUE
  
  fit_M42_condM1M2M3 <- glm(as.formula(paste("M42 ~ (Mid + High + M1 + M21 + M22 + M23 + M3 + M41)^2 - Mid:High + ",confs)), 
                          data = data, family = binomial, y = FALSE, model = FALSE)
  if((!fit_M42_condM1M2M3$converged)|any(is.na(fit_M42_condM1M2M3$coefficients))) flag <- TRUE
  
  
  ## Outcome model ##
  
  # 5.a. outcome model
  
  fit_Y <- glm(as.formula(paste("Y ~ (Mid + High + M1 + M21 + M22 + M23 + M3 + M41 + M42)^2 - Mid:High +",confs)), 
               data = data, family = binomial, y = FALSE, model = FALSE)
  if((!fit_Y$converged)|any(is.na(fit_Y$coefficients))) flag<-TRUE
   
 
  ### Step B: Randomly draw mediator values and predict risks of poor PROMs ###
  
  ## Predicted values that only depend on exposure and confounder values ##
  ##      and therefore do not need to be updated with new MC-draws      ##
  
  dat2$Mid <- 1 # Set Mid = 1, which corresponds to unexposed here
  M1_preds0_marg <- predict(fit_M1_marg, newdata = dat2, type = "response")
  M21_preds0_marg <- predict(fit_M21_marg, newdata = dat2, type = "response")
  M3_preds0_marg <- predict(fit_M3_marg, newdata = dat2, type = "response")
  M41_preds0_marg <- predict(fit_M41_marg, newdata = dat2, type = "response")
  
  dat2$Mid <-0 # Set Mid = 0, which corresponds to exposed (low ses) here
  M1_preds1_marg <- predict(fit_M1_marg, newdata = dat2, type = "response")
  M21_preds1_marg <- predict(fit_M21_marg, newdata = dat2, type = "response")
  M3_preds1_marg <- predict(fit_M3_marg, newdata = dat2, type = "response")
  M41_preds1_marg <- predict(fit_M41_marg, newdata = dat2, type = "response")
  
  dat3$High <- 1 # Set High = 1, which corresponds to unexposed here
  M1_preds0_marg3 <- predict(fit_M1_marg, newdata = dat3, type = "response")
  M21_preds0_marg3 <- predict(fit_M21_marg, newdata = dat3, type = "response")
  M3_preds0_marg3 <- predict(fit_M3_marg, newdata = dat3, type = "response")
  M41_preds0_marg3 <- predict(fit_M41_marg, newdata = dat3, type = "response")
  
  dat3$High <-0 # Set High = 0, which corresponds to exposed (low ses) here
  M1_preds1_marg3 <- predict(fit_M1_marg, newdata = dat3, type = "response")
  M21_preds1_marg3 <- predict(fit_M21_marg, newdata = dat3, type = "response")
  M3_preds1_marg3 <- predict(fit_M3_marg, newdata = dat3, type = "response")
  M41_preds1_marg3 <- predict(fit_M41_marg, newdata = dat3, type = "response")
  
  ## Loop for Monte Carlo simulations for the low vs. mid contrast ##
  
  for(i in 1:mcsim){ 
    print(i)
    
    
    ### Simulate from mediator distributions ###
    
    ## Marginal distributions ##
    
    # For the unexposed condition (mid SES) #
    
    dat2$Mid <- 1 # Set the unexposed condition
    M1_0_marg  <- rbinom(n, 1, M1_preds0_marg)
    
    M21_0_marg <- rbinom(n, 1, M21_preds0_marg)
    dat2$M21 <- M21_0_marg # Set the M21 variable to the simulated values since these are used to simulate M22.
    M22_0_marg  <- rbinom(n, 1, predict(fit_M22_marg, newdata = dat2, type = "response"))
    dat2$M22 <- M22_0_marg
    M23_0_marg  <- rbinom(n, 1, predict(fit_M23_marg, newdata = dat2, type = "response"))
    
    M3_0_marg  <- rbinom(n, 1, predict(fit_M3_marg, newdata = dat2, type = "response"))
    
    M41_0_marg <- rbinom(n, 1, M41_preds0_marg)
    dat2$M41 <- M41_0_marg
    M42_0_marg <- rbinom(n, 1, predict(fit_M42_marg, newdata = dat2, type = "response"))
    
    # For the exposed condition (low SES) #
    
    dat2$Mid<-0 # Set the exposed condition
    M1_1_marg  <- rbinom(n, 1, M1_preds1_marg)
    
    M21_1_marg <- rbinom(n, 1, M21_preds1_marg)
    dat2$M21 <- M21_1_marg
    M22_1_marg  <- rbinom(n, 1, predict(fit_M22_marg, newdata = dat2, type = "response"))
    dat2$M22 <- M22_1_marg
    M23_1_marg  <- rbinom(n, 1, predict(fit_M23_marg, newdata = dat2, type = "response"))
    
    M3_1_marg  <- rbinom(n, 1, predict(fit_M3_marg, newdata = dat2, type = "response"))
    
    M41_1_marg <- rbinom(n, 1, M41_preds1_marg)
    dat2$M41 <- M41_1_marg
    M42_1_marg <- rbinom(n, 1, predict(fit_M42_marg, newdata = dat2, type = "response"))
    
    
    ## Conditional distributions ##
    
    # M2 given M1 #
    
    # Unexposed
    dat2$Mid <- 1
    dat2$M1 <- M1_0_marg
    
    M21_0_cond <- rbinom(n, 1, predict(fit_M21_cond, newdata =dat2, type = "response"))
    dat2$M21 <- M21_0_cond
    M22_0_cond <- rbinom(n, 1, predict(fit_M22_cond, newdata =dat2, type = "response"))
    dat2$M22 <- M22_0_cond
    M23_0_cond <- rbinom(n, 1, predict(fit_M23_cond, newdata =dat2, type = "response"))
    
    # Exposed
    dat2$Mid <- 0
    dat2$M1 <- M1_1_marg
    
    M21_1_cond <- rbinom(n, 1, predict(fit_M21_cond, newdata =dat2, type = "response"))
    dat2$M21 <- M21_1_cond
    M22_1_cond <- rbinom(n, 1, predict(fit_M22_cond, newdata =dat2, type = "response"))
    dat2$M22 <- M22_1_cond
    M23_1_cond <- rbinom(n, 1, predict(fit_M23_cond, newdata =dat2, type = "response"))
    
    
    # M3 given M1 (only exposed)
    
    M3_1_condM1 <- rbinom(n, 1, predict(fit_M3_condM1, newdata = dat2, type = "response"))
    
    # M3 given M2 (only exposed)
    
    dat2$M21 <- M21_1_marg
    dat2$M22 <- M22_1_marg
    dat2$M23 <- M23_1_marg
    
    M3_1_condM2 <- rbinom(n, 1, predict(fit_M3_condM2, newdata = dat2, type = "response"))
    
    
    # M3 given M1 and M2
    
    # Unexposed
    dat2$Mid <- 1
    
    dat2$M1 <- M1_0_marg
    dat2$M21 <- M21_0_cond
    dat2$M22 <- M22_0_cond
    dat2$M23 <- M23_0_cond
    
    M3_0_condM1M2 <- rbinom(n, 1, predict(fit_M3_condM1M2, newdata = dat2, type = "response"))
    
    # Exposed
    dat2$Mid <- 0
    
    dat2$M1 <- M1_1_marg
    dat2$M21 <- M21_1_cond
    dat2$M22 <- M22_1_cond
    dat2$M23 <- M23_1_cond
    
    M3_1_condM1M2 <- rbinom(n, 1, predict(fit_M3_condM1M2, newdata = dat2, type = "response"))
    
    # Predicted risk of poor PROM under exposure, M4 under unexposed, M1,M2,M3 under exposure #
    dat2$M3 <- M3_1_condM1M2
    dat2$M41 <- M41_0_marg
    dat2$M42 <- M42_0_marg
    
    risk_M4_int[i] <- mean(predict(fit_Y, newdata = dat2, type = "response"))
    
    # M4 given M1 and M2 (only exposed)
    M41_1_condM1M2 <- rbinom(n, 1, predict(fit_M41_condM1M2, newdata = dat2, type = "response"))
    dat2$M41 <- M41_1_condM1M2
    M42_1_condM1M2 <- rbinom(n, 1, predict(fit_M42_condM1M2, newdata = dat2, type = "response"))
    
    # Predicted risk of poor PROM under exposure, M3 under unexposed, M1,M2,M4 under exposure #
    dat2$M3 <- M3_0_marg
    dat2$M42 <- M42_1_condM1M2
    
    risk_M3_int[i] <- mean(predict(fit_Y, newdata = dat2, type = "response"))
    
    # M4 given M1 and M3 (only exposed)
    dat2$M3 <- M3_1_condM1
    M41_1_condM1M3 <- rbinom(n, 1, predict(fit_M41_condM1M3, newdata = dat2, type = "response"))
    dat2$M41 <- M41_1_condM1M3
    M42_1_condM1M3 <- rbinom(n, 1, predict(fit_M42_condM1M3, newdata = dat2, type = "response"))
    
    # Predicted risk of poor PROM under exposure, M2 under unexposed, M1,M3,M4 under exposure #
    dat2$M21 <- M21_0_marg
    dat2$M22 <- M22_0_marg
    dat2$M23 <- M23_0_marg
    dat2$M42 <- M42_1_condM1M3
    
    risk_M2_int[i] <- mean(predict(fit_Y, newdata = dat2, type = "response"))
    
    # M4 given M2 and M3 (only exposed)
    dat2$M21 <- M21_1_marg
    dat2$M22 <- M22_1_marg
    dat2$M23 <- M23_1_marg
    
    dat2$M3 <- M3_1_condM2
    M41_1_condM2M3 <- rbinom(n, 1, predict(fit_M41_condM2M3, newdata = dat2, type = "response"))
    dat2$M41 <- M41_1_condM2M3
    M42_1_condM2M3 <- rbinom(n, 1, predict(fit_M42_condM2M3, newdata = dat2, type = "response"))
    
    # Predicted risk of poor PROM under exposure, M1 under unexposed, M2,M3,M4 under exposure #
    dat2$M1 <- M1_0_marg
    dat2$M42 <- M42_1_condM2M3
    
    risk_M1_int[i] <- mean(predict(fit_Y, newdata = dat2, type = "response"))
    
    # M4 given M1, M2 and M3 under exposure
    dat2$M1 <- M1_1_marg
    dat2$M21 <- M21_1_cond
    dat2$M22 <- M22_1_cond
    dat2$M23 <- M23_1_cond
    dat2$M3 <- M3_1_condM1M2
    
    M41_1_condM1M2M3 <- rbinom(n, 1, predict(fit_M41_condM1M2M3, newdata = dat2, type = "response"))
    dat2$M41 <- M41_1_condM1M2M3
    M42_1_condM1M2M3 <- rbinom(n, 1, predict(fit_M42_condM1M2M3, newdata = dat2, type = "response"))
    
    dat2$M42 <- M42_1_condM1M2M3
    
    # Predicted risk of poor PROM under exposure with joint distribution for M1, M2, M3 and M4 under exposure #
    risk_exposed[i] <- mean(predict(fit_Y, newdata = dat2, type = "response"))
    
    # M4 given M1, M2 and M3 under no exposure
    dat2$Mid <- 1
    dat2$M1 <- M1_0_marg
    dat2$M21 <- M21_0_cond
    dat2$M22 <- M22_0_cond
    dat2$M23 <- M23_0_cond
    dat2$M3 <- M3_0_condM1M2
    
    M41_0_condM1M2M3 <- rbinom(n, 1, predict(fit_M41_condM1M2M3, newdata = dat2, type = "response"))
    dat2$M41 <- M41_0_condM1M2M3
    M42_0_condM1M2M3 <- rbinom(n, 1, predict(fit_M42_condM1M2M3, newdata = dat2, type = "response"))
    
    
    # Predicted risk of poor PROM under exposure with joint distribution for M1, M2, M3 and M4 under no exposure #
    dat2$Mid <- 0
    dat2$M42 <- M42_0_condM1M2M3
    
    risk_joint_int[i] <- mean(predict(fit_Y, newdata = dat2, type = "response"))
    
    
    # Predicted risk of poor PROM under no exposure with joint distribution for M1, M2, M3 and M4 under no exposure #
    dat2$Mid <- 1
    
    risk_unexposed[i] <- mean(predict(fit_Y, newdata = dat2, type = "response"))
    
    
  }
  
  ### Step C: Estimate effects ###
  
  ## Total risk difference
  
  TRD <- mean(risk_exposed - risk_unexposed)
  
  
  ## Indirect effects
  
  reduction_M1 <- mean(risk_exposed - risk_M1_int)
  
  reduction_M2 <- mean(risk_exposed - risk_M2_int)
  
  reduction_M3 <- mean(risk_exposed - risk_M3_int)
  
  reduction_M4 <- mean(risk_exposed - risk_M4_int)
  
  reduction_joint <- mean(risk_exposed - risk_joint_int)
  
  ## Direct effects
  
  remaining_joint <- mean(risk_joint_int - risk_unexposed)
  
  remaining_M1 <- TRD - reduction_M1
  
  remaining_M2 <- TRD - reduction_M2
  
  remaining_M3 <- TRD - reduction_M3
  
  remaining_M4 <- TRD - reduction_M4
  
  ### Collect and return results
  
  resMid <- c("TRD_Mid"=TRD, "reduction_joint_Mid"=reduction_joint, "reduction_M1_Mid"=reduction_M1, "reduction_M2_Mid"=reduction_M2, "reduction_M3_Mid"=reduction_M3,
              "reduction_M4_Mid"=reduction_M4, "remaining_joint_Mid"=remaining_joint, "remaining_M1_Mid"=remaining_M1, "remaining_M2_Mid"=remaining_M2, "remaining_M3_Mid"=remaining_M3, "remaining_M4_Mid"=remaining_M4)
  
  # _______________________________________________________________________ #
  
  ### Step B: Randomly draw mediator values and predict risks of poor PROMs low vs. high contrast ###
  
  ## Loop for Monte Carlo simulations for the low vs. high contrast ##
  for(i in 1:mcsim){ 
    print(i)
    
    
    ### Simulate from mediator distributions ###
    
    ## Marginal distributions ##
    
    # For the unexposed condition (high SES) #
    
    dat3$High <- 1 # Set the unexposed condition
    M1_0_marg  <- rbinom(n3, 1, M1_preds0_marg3)
    
    M21_0_marg <- rbinom(n3, 1, M21_preds0_marg3)
    dat3$M21 <- M21_0_marg # Set the M21 variable to the simulated values since these are used to simulate M22.
    M22_0_marg  <- rbinom(n3, 1, predict(fit_M22_marg, newdata = dat3, type = "response"))
    dat3$M22 <- M22_0_marg
    M23_0_marg  <- rbinom(n3, 1, predict(fit_M23_marg, newdata = dat3, type = "response"))
    
    M3_0_marg  <- rbinom(n3, 1, predict(fit_M3_marg, newdata = dat3, type = "response"))
    
    M41_0_marg <- rbinom(n3, 1, M41_preds0_marg3)
    dat3$M41 <- M41_0_marg
    M42_0_marg <- rbinom(n3, 1, predict(fit_M42_marg, newdata = dat3, type = "response"))
    
    # For the exposed condition (low SES) #
    
    dat3$High<-0 # Set the exposed condition
    M1_1_marg  <- rbinom(n3, 1, M1_preds1_marg3)
    
    M21_1_marg <- rbinom(n3, 1, M21_preds1_marg3)
    dat3$M21 <- M21_1_marg
    M22_1_marg  <- rbinom(n3, 1, predict(fit_M22_marg, newdata = dat3, type = "response"))
    dat3$M22 <- M22_1_marg
    M23_1_marg  <- rbinom(n3, 1, predict(fit_M23_marg, newdata = dat3, type = "response"))
    
    M3_1_marg  <- rbinom(n3, 1, predict(fit_M3_marg, newdata = dat3, type = "response"))
    
    M41_1_marg <- rbinom(n3, 1, M41_preds1_marg3)
    dat3$M41 <- M41_1_marg
    M42_1_marg <- rbinom(n3, 1, predict(fit_M42_marg, newdata = dat3, type = "response"))
    
    
    ## Conditional distributions ##
    
    # M2 given M1 #
    
    # Unexposed
    dat3$High <- 1
    dat3$M1 <- M1_0_marg
    
    M21_0_cond <- rbinom(n3, 1, predict(fit_M21_cond, newdata =dat3, type = "response"))
    dat3$M21 <- M21_0_cond
    M22_0_cond <- rbinom(n3, 1, predict(fit_M22_cond, newdata =dat3, type = "response"))
    dat3$M22 <- M22_0_cond
    M23_0_cond <- rbinom(n3, 1, predict(fit_M23_cond, newdata =dat3, type = "response"))
    
    # Exposed
    dat3$High <- 0
    dat3$M1 <- M1_1_marg
    
    M21_1_cond <- rbinom(n3, 1, predict(fit_M21_cond, newdata =dat3, type = "response"))
    dat3$M21 <- M21_1_cond
    M22_1_cond <- rbinom(n3, 1, predict(fit_M22_cond, newdata =dat3, type = "response"))
    dat3$M22 <- M22_1_cond
    M23_1_cond <- rbinom(n3, 1, predict(fit_M23_cond, newdata =dat3, type = "response"))
    
    
    # M3 given M1 (only exposed)
    
    M3_1_condM1 <- rbinom(n3, 1, predict(fit_M3_condM1, newdata = dat3, type = "response"))
    
    # M3 given M2 (only exposed)
    
    dat3$M21 <- M21_1_marg
    dat3$M22 <- M22_1_marg
    dat3$M23 <- M23_1_marg
    
    M3_1_condM2 <- rbinom(n3, 1, predict(fit_M3_condM2, newdata = dat3, type = "response"))
    
    
    # M3 given M1 and M2
    
    # Unexposed
    dat3$High <- 1
    
    dat3$M1 <- M1_0_marg
    dat3$M21 <- M21_0_cond
    dat3$M22 <- M22_0_cond
    dat3$M23 <- M23_0_cond
    
    M3_0_condM1M2 <- rbinom(n3, 1, predict(fit_M3_condM1M2, newdata = dat3, type = "response"))
    
    # Exposed
    dat3$High <- 0
    
    dat3$M1 <- M1_1_marg
    dat3$M21 <- M21_1_cond
    dat3$M22 <- M22_1_cond
    dat3$M23 <- M23_1_cond
    
    M3_1_condM1M2 <- rbinom(n3, 1, predict(fit_M3_condM1M2, newdata = dat3, type = "response"))
    
    # Predicted risk of poor PROM under exposure, M4 under unexposed, M1,M2,M3 under exposure #
    dat3$M3 <- M3_1_condM1M2
    dat3$M41 <- M41_0_marg
    dat3$M42 <- M42_0_marg
    
    risk_M4_int[i] <- mean(predict(fit_Y, newdata = dat3, type = "response"))
    
    # M4 given M1 and M2 (only exposed)
    M41_1_condM1M2 <- rbinom(n3, 1, predict(fit_M41_condM1M2, newdata = dat3, type = "response"))
    dat3$M41 <- M41_1_condM1M2
    M42_1_condM1M2 <- rbinom(n3, 1, predict(fit_M42_condM1M2, newdata = dat3, type = "response"))
    
    # Predicted risk of poor PROM under exposure, M3 under unexposed, M1,M2,M4 under exposure #
    dat3$M3 <- M3_0_marg
    dat3$M42 <- M42_1_condM1M2
    
    risk_M3_int[i] <- mean(predict(fit_Y, newdata = dat3, type = "response"))
    
    # M4 given M1 and M3 (only exposed)
    dat3$M3 <- M3_1_condM1
    M41_1_condM1M3 <- rbinom(n3, 1, predict(fit_M41_condM1M3, newdata = dat3, type = "response"))
    dat3$M41 <- M41_1_condM1M3
    M42_1_condM1M3 <- rbinom(n3, 1, predict(fit_M42_condM1M3, newdata = dat3, type = "response"))
    
    # Predicted risk of poor PROM under exposure, M2 under unexposed, M1,M3,M4 under exposure #
    dat3$M21 <- M21_0_marg
    dat3$M22 <- M22_0_marg
    dat3$M23 <- M23_0_marg
    dat3$M42 <- M42_1_condM1M3
    
    risk_M2_int[i] <- mean(predict(fit_Y, newdata = dat3, type = "response"))
    
    # M4 given M2 and M3 (only exposed)
    dat3$M21 <- M21_1_marg
    dat3$M22 <- M22_1_marg
    dat3$M23 <- M23_1_marg
    
    dat3$M3 <- M3_1_condM2
    M41_1_condM2M3 <- rbinom(n3, 1, predict(fit_M41_condM2M3, newdata = dat3, type = "response"))
    dat3$M41 <- M41_1_condM2M3
    M42_1_condM2M3 <- rbinom(n3, 1, predict(fit_M42_condM2M3, newdata = dat3, type = "response"))
    
    # Predicted risk of poor PROM under exposure, M1 under unexposed, M2,M3,M4 under exposure #
    dat3$M1 <- M1_0_marg
    dat3$M42 <- M42_1_condM2M3
    
    risk_M1_int[i] <- mean(predict(fit_Y, newdata = dat3, type = "response"))
    
    # M4 given M1, M2 and M3 under exposure
    dat3$M1 <- M1_1_marg
    dat3$M21 <- M21_1_cond
    dat3$M22 <- M22_1_cond
    dat3$M23 <- M23_1_cond
    dat3$M3 <- M3_1_condM1M2
    
    M41_1_condM1M2M3 <- rbinom(n3, 1, predict(fit_M41_condM1M2M3, newdata = dat3, type = "response"))
    dat3$M41 <- M41_1_condM1M2M3
    M42_1_condM1M2M3 <- rbinom(n3, 1, predict(fit_M42_condM1M2M3, newdata = dat3, type = "response"))
    
    dat3$M42 <- M42_1_condM1M2M3
    
    # Predicted risk of poor PROM under exposure with joint distribution for M1, M2, M3 and M4 under exposure #
    risk_exposed[i] <- mean(predict(fit_Y, newdata = dat3, type = "response"))
    
    # M4 given M1, M2 and M3 under no exposure
    dat3$High <- 1
    dat3$M1 <- M1_0_marg
    dat3$M21 <- M21_0_cond
    dat3$M22 <- M22_0_cond
    dat3$M23 <- M23_0_cond
    dat3$M3 <- M3_0_condM1M2
    
    M41_0_condM1M2M3 <- rbinom(n3, 1, predict(fit_M41_condM1M2M3, newdata = dat3, type = "response"))
    dat3$M41 <- M41_0_condM1M2M3
    M42_0_condM1M2M3 <- rbinom(n3, 1, predict(fit_M42_condM1M2M3, newdata = dat3, type = "response"))
    
    
    # Predicted risk of poor PROM under exposure with joint distribution for M1, M2, M3 and M4 under no exposure #
    dat3$High <- 0
    dat3$M42 <- M42_0_condM1M2M3
    
    risk_joint_int[i] <- mean(predict(fit_Y, newdata = dat3, type = "response"))
    
    
    # Predicted risk of poor PROM under no exposure with joint distribution for M1, M2, M3 and M4 under no exposure #
    dat3$High <- 1
    
    risk_unexposed[i] <- mean(predict(fit_Y, newdata = dat3, type = "response"))
    
    
  }
  
  
  ### Step C: Estimate effects ###
  
  ## Total risk difference
  
  TRD <- mean(risk_exposed - risk_unexposed)
  
  
  ## Indirect effects
  
  reduction_M1 <- mean(risk_exposed - risk_M1_int)
  
  reduction_M2 <- mean(risk_exposed - risk_M2_int)
  
  reduction_M3 <- mean(risk_exposed - risk_M3_int)
  
  reduction_M4 <- mean(risk_exposed - risk_M4_int)
  
  reduction_joint <- mean(risk_exposed - risk_joint_int)
  
  ## Direct effects
  
  remaining_joint <- mean(risk_joint_int - risk_unexposed)
  
  remaining_M1 <- TRD - reduction_M1
  
  remaining_M2 <- TRD - reduction_M2
  
  remaining_M3 <- TRD - reduction_M3
  
  remaining_M4 <- TRD - reduction_M4
  
  ### Collect and return results
  
  resHigh <- c("TRD_High"=TRD, "reduction_joint_High"=reduction_joint, "reduction_M1_High"=reduction_M1, "reduction_M2_High"=reduction_M2, "reduction_M3_High"=reduction_M3,
              "reduction_M4_High"=reduction_M4, "remaining_joint_High"=remaining_joint, "remaining_M1_High"=remaining_M1, "remaining_M2_High"=remaining_M2, "remaining_M3_High"=remaining_M3, "remaining_M4_High"=remaining_M4)
  

  res <- c(resMid,resHigh)


  if(!flag)return(res) else
    return(rep(NA,length(res)))
  
}