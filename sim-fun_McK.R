##################
##function file cdsa
##################

gen.data <- function(numb.j = 100, 
                     model = 'beta', 
                     dist = "left"){
  if(model == 'beta'){
    if(dist == "left"){
      return(rbeta(numb.j, 2, 12))
    }
    if(dist == "unimodal"){
      return(rbeta(numb.j, 7, 7))
    }
    if(dist == "bimodal"){
      mode1	<- rbeta(0.5 * numb.j, 1, 6)
      mode2	<- rbeta(0.5 * numb.j, 6, 1)	
      return(c(mode1, mode2))
    }
  }
  if(model == 'probit'){
    if(dist == "left"){
      return(pnorm(rnorm(numb.j, -1.156, .415))) 
    }
    if(dist == 'unimodal'){
      return(pnorm(rnorm(numb.j, 0, 0.341)))
    }
    if(dist == 'bimodal'){
      mode1	<- pnorm(rnorm(0.5 * numb.j, -1.156, 1.35))
      mode2	<- pnorm(rnorm(0.5 * numb.j, 1.156, 1.35))
      return(c(mode1,mode2))
    }
  }
}


generate_drift <- function(strength = 'weak', true_prob = rbeta(500, 7,7), 
                           negative = FALSE, units = 500, 
                           lengthstream = 10000){
  modes <- seq(0+(1/(2*units)), 1 - ((1/(2*units))),1/units)
  upperlimit <- modes + 1/(2*units)
  lowerlimit <- modes - 1/(2*units)
  sort_prob <- sort(true_prob, decreasing = negative)
  abeta <- 1
  id <- c()
  if(strength == 'weak'){
    bbeta <- 1.27
  }
  if(strength == 'medium'){
    bbeta <- 1.85
  }
  if(strength == 'strong'){
    bbeta <- 2.8
  }
  dif <- bbeta - abeta 
  
  for(i in 1:lengthstream){
    concept_drift <- pbeta(upperlimit, abeta, bbeta) - 
      pbeta(lowerlimit, abeta, bbeta)
    if(sum(sort_prob < 0) > 0) {
      sort_prob[which(sort_prob < 0)] <- 0
    }
    id <- c(id,sample(1:units,1, replace = TRUE, prob = concept_drift ))
    abeta <- abeta + dif/(lengthstream)
    bbeta <- bbeta - dif/(lengthstream)
  }
  data <- data.frame(id = id, obs = rbinom(lengthstream, 1, sort_prob[id]))
  return(data)
}
update_average <- function(average = 0, obs = 1, numb.obs = 2){
  return(average + (obs - average) / numb.obs)
}

update_contribution <- function(statistic = 5, old_contr = 2, new_contr = 3){
  return( statistic - old_contr + new_contr)
}

transform_y <- function(n_bar = 2, 
                        average_j = 2, 
                        average = 3, 
                        n_j = 1, 
                        type ='morris'){
  if(type == 'morris'){
    return(sqrt(n_bar) * (asin(1 - 2 * average_j) - (asin(1 - 2 * average))) )
  }
  if(type == 'brown'){
    return(asin(sqrt( (average_j * n_j + .25) / (n_j + .5))))
  }
}

shrink_js <- function(units = 10, SumofSquares = 20){
  return(min((units - 2) / SumofSquares, 1))
}

#this one needs to be checked:
shrink_n <- function(numb.obs = 5, constant = 2){
  return(min(constant/(numb.obs + constant), 1))
}

#here there is a hardcoded .01 as learnrate, do we want to make this adaptable? 
# taumle <- function(tau = 1, 
#                    numb.obs = 1, 
#                    y_j = 1, 
#                    n_bar = 2, 
#                    n = 1, 
#                    average_y = 0){
# 	v_j	<- n_bar / numb.obs
# 	z_y_j	<- (y_j - average_y)**2 / v_j
# 	
# 	gradient <- (((z_y_j - 1) * v_j) - tau) / (2 * (v_j + tau)**2)  
# 	tau <- max(tau + .01 * gradient, 0)
# 	return(tau)
# }

shrink_ml <- function(varapp = 1, tauml = 1, n.id = 2){
  if(n.id < 1){
    B <- 1
  }
  else{
    B <- varapp / (varapp + tauml)
  }
  return(B)
}

shrink_mm <- function(sig_j = 1, tau_mm = 1, n.id = 1){
  if(n.id < 1){
    B <- 1
  }
  else{
    B <- sig_j / (sig_j + tau_mm)
  }
  return(B)
}

shrink_BB <- function(units = 10, 
                      numb.obs = 100, 
                      nj = 5, 
                      SumofSquares = 10, 
                      average = 5, 
                      inv_nj = 2){
  sigmasq	<- (units * SumofSquares) / ((units - 1) * numb.obs)
  M	<- max((average * (1 - average) - sigmasq) / 
             (sigmasq - ((average * (1 - average)) / units) * inv_nj), 0)
  b	<- M / (M + nj)
  if(M == 0 | sigmasq == 0) b <- 1
  return(b)
}

estimate_mu <- function(b = 1, 
                        average_j = 2, #transformed probability individual
                        average = 0,   #average probability
                        average_y = 0, #average transformed probabilities
                        average_n = 5, #average number of observations
                        n_j = 10,
                        transform = TRUE, 
                        type = 'morris'){
  if(length(b) == 1){
    if(is.na(b)) b <- 1
  }  
  if(transform == TRUE){
    y <- (1 - b) * average_j + b * average_y
    if(type == 'morris'){
      return(
        .5 * (1 - (sin((sqrt(average_n) * asin(1 - 2 * average) + y) / 
                         sqrt(average_n))))
      )
      # (sin((y / sqrt(average_n)) + asin(1 - 2 * average)) - 1) / (-2)    incorrect one       
    }
    if(type == 'brown'){
      return((sin(y)**2 * (n_j + .5) - .25) / n_j)
    }
  }
  if(transform == FALSE){
    return(b * average + (1 - b) * average_j)
  }
}

compute_error <- function(mu_j = c(1, 2, 3), 
                          true = c(3, 2, 1), 
                          units = 3){
  return(sum((mu_j - true)**2) / units )
}

offline_js <- function(average_n = 20, 
                       units = 5,
                       id_data = data.frame(id = 1:5, n_j = 18:22), 
                       y_j = c(1, 1.5, 2, 2.5, 3), 
                       average_y = 0, type = 'morris'){
  if(type == 'morris'){
    v_j	<- average_n / id_data$n_j
    z_y_j	<- (y_j - average_y)**2 / v_j
    SS	<- sum(z_y_j)
  }
  if(type == 'brown'){
    v_j <- 1 / (4 * id_data$n_j)
    SS <- sum((y_j - average_y)**2/v_j)
  }
  return(min((units - 2) / SS, 1))
}

offline_n <- function(numb.obs = c(1, 2), constant = 2){
  return(min(constant / (numb.obs + constant), 1))
}

#there is a hardcoded .01 in this code for learn rate should we change that?  
#for the formulation of tau I kept the same formulation as we did with the morris 
#tranformation but I do not know if this is actually correct 
offline_ml <- function(tau = 1, 
                       average_n = 10, 
                       numb.obs = c(1, 2, 3), 
                       y_j = c(1,2,3),
                       average_y = 0, 
                       type = 'morris'){
  if(type == 'morris'){
    v_j	<- average_n / numb.obs
    z_y_j	<- (y_j - average_y)**2
    for(s in 1:500){	
      gradient <- sum( (z_y_j - v_j - tau) / (2 * (v_j + tau)**2) )
      tau <- max(tau + 0.02 * gradient, 0)
      if(abs(gradient) < 0.0001){break}
    }
    return(v_j / (v_j + tau))  
  }
  if(type == 'brown'){
    v_j	<- 1 / (4 * numb.obs)			
    for(s in 1:1000){	
      mu <- sum( (y_j / (tau + v_j))) / sum(1 / (tau + v_j))
      z_y_j	<- (y_j - mu)**2    
      
      gradient <- sum( (z_y_j - v_j - tau) / (2 * (v_j + tau)**2) )
      tau <- max(tau + 0.000005 * gradient, 0)
      if(abs(gradient) < 0.0001){break}
    }
    return(c(v_j / (v_j + tau), mu))
  }#0.0001
}

offline_BB <- function(units = 10, 
                       id_data = id_df, 
                       average = 0.5, 
                       numb.obs = 2){
  sig	<- (units * sum(id_data$n_j * (id_data$p_j - average)**2)) / 
    ((units - 1) * sum(numb.obs))
  M	<- max((average * (1 - average) - sig) / 
             (sig - ((average * (1 - average)) / units) * sum(1 / id_data$n_j)),
           0)
  b	<- M / (M + id_data$n_j)
  if(M == 0){b <- 1}
  return(b)
}

compute_tauj <- function(tausq = 1, y = 1, bar_y = 0, sig_j = 1, type = 'morris'){
  if(type == 'morris'){
    return(0.02 * (-(tausq - (y-bar_y)**2 + sig_j))/(2 * (tausq + sig_j)**2))
  }
  if(type == 'brown'){
    return(0.00005 * (-(tausq - (y-bar_y)**2 + sig_j))/(2 * (tausq + sig_j)**2))
  }
}


classify	<- function(predict = .3, observed = 1){
  truepos	<- trueneg <- falseneg <- falsepos <- rep(0, 48)
  for(q in 1:48){
    testjs	<- predict[observed[, (q + 2)] != -1 & observed$stop > q, (q + 2)]
    testliss	<- observed[observed[, (q + 2)] != -1 & observed$stop > q, (q + 2)]
    for(i in 1:length(testliss)){
      if(testliss[i] == 1 & testjs[i] > .5) truepos[q] <- truepos[q] + 1		
      if(testliss[i] == 1 & testjs[i] < .5) falseneg[q] <- falseneg[q] + 1
      if(testliss[i] == 0 & testjs[i] > .5) falsepos[q] <- falsepos[q] + 1		
      if(testliss[i] == 0 & testjs[i] < .5) trueneg[q] <- trueneg[q] + 1		
    }
  }
  return(cbind(TP = truepos, TN = trueneg, FP = falsepos, FN = falseneg))
}

compute_mu_MM <- function(inverse_tot_variance = .2, 
                          sum_std_y = 4){
  return(sum_std_y/inverse_tot_variance)
}

compute_tau_MM <- function(units = 100, 
                           sumSquaredDif = 1, 
                           sumVariance_j = 4){
  tau <- max(sumSquaredDif / max(units - 1, 1) - sumVariance_j/units, 0)
  return(tau)
}

offline_MM <- function(vec_nj = c(5, 8), y_j, units = 2, sig_j = c(.05, .03)){
  mean_Sig_j <- sum(sig_j) / units
  mu_0 <- sum(y_j) / units 
  tau_0 <- 0
  for(s in 1:100){
    tau_1 <- (sum( (y_j - mu_0)**2 ) / (units - 1)) - mean_Sig_j
    mu_1 <- sum(y_j / (sig_j + tau_1))/sum(1 / (sig_j + tau_1))
    if(abs(tau_1-tau_0) < .0001 &
       abs(mu_1 - mu_0) < .0001){break}
    else{
      mu_0 <- mu_1
      tau_0 <- tau_1
    }
  }
  if(tau_0 < 0){tau_0 <- 0}
  return(c(sig_j / (sig_j + tau_0), mu_0))
}


compute_sig_j <- function(n_j){
  1/(4 * max(n_j, 1))
}

# McK: Check deze op basis van artikelen?
std_y_brown <- function(y_j = 0, sig_j = 0){
  std <- (y_j / sig_j)
  if(is.na(std)){
    std <- 0
  }
  return(std)
}

