#setwd("C:/Users/Lianne.Ippel/Dropbox/Lianne/Submissions/2.shrink_factors_Comput_stat/archive/Comput stat/R/liss")
#mydat<-read.table("liss_long.txt")
#library(lme4)
#setwd("I:/shrinkage factor paper")
N								<- nrow(mydat)

q <- 1
id_df <- data.frame(id = sort(unique(mydat$id)), stop=999,
                    p_j = 0, n_j = 0, v_j = 0, sig_j = 0, y_j_morris = 0, 
                    y_j_brown = 0, mu_j_js_morris = NA, mu_j_js_brown = NA,											  
                    mu_j_ml_morris = NA, mu_j_ml_brown = NA, 
                    mu_j_mm_morris = NA, mu_j_mm_brown = NA,	  
                    mu_j_nsf = NA, mu_j_BB = NA, std_y_brown = 0,
                    SSadd = 0, SSjsadd_morris = 0, SSjsadd_brown = 0,
                    bjs_morris = NA, bjs_brown = NA, bml_morris = NA,
                    bml_brown = NA, bmm_morris = NA, bmm_brown = NA, bn = NA, 
                    bBB = NA, tau_j_morris = 0, tau_j_brown = 0, 
                    dif_j_brown = 0, y_brown_div_tau_sig = 0, 
                    inv_tau_sig_j_brown_MM = 0,
                    inv_tau_sig_j_morris = 0, y_morris_div_tau_sig = 0,
                    dif_j_morris = 0, inv_tau_sig_j_brown = 0, 
                    y_brown_div_tau_sig_MM = 0)
#global statistics
pbar <- n <- j <- sbar <- 0
nbar <- 1

#james stein morris:
SSjs_morris <- 0
#james stein brown
SSjs_brown <- sum_std_y_brown <- sum_inverse_sig <- 0 
# method of moments morris 
sumdif_j_morris <- sum_y_morris_div_tau_sig <- sumsig_j_morris <- 
  sum_inverse_tot_var_morris <- 0
# method of moments brown 
sumdif_j_brown <- sum_y_brown_div_tau_sig_MM <- sumsig_j_brown <- 
  sum_inverse_tot_var_brown_MM <- 0 
# maximum likelihood brown
sum_inverse_tot_var_brown <- sum_y_brown_div_tau_sig <- 0
#beta binomial
SSbb <- inversnj <- 0

#variance starting values: 
tau_ml_hat_morris <- tau_mm_hat_morris <- tau_ml_hat_brown <- 
  tau_mm_hat_brown <- 0
#mu starting values:
mu_mm_brown <- ybar_morris <- mu_ml_brown <- mu_js_brown <- 0



dropouts							<- sort(unique(mydat$id[!is.na(mydat$stop)]))
for(t in 1:nrow(id_df)){
  if(sum(id_df$id[t] == dropouts) == 1){
    id_df$stop[t] <- mydat$stop[which(id_df$id[t] == mydat$id)][1]
  }
}
#m1								<- glmer(response~1+(1|id), data=mydat, family=binomial,nAGQ=20)
resJSonM_df <- resJSonB_df <- 
  resJSofM_df <- resJSofB_df <- 
  resMLonM_df <- resMLonB_df <- 
  resMLofM_df <- resMLofB_df <- 
  resMMonM_df <- resMMonB_df <- 
  resMMofM_df <- resMMofB_df <- 
  resBBon_df <- resBBof_df <-
  resHNon_df	<- resHNof_df	<- resglmof_df <- data.frame(id=sort(unique(mydat$id)))

resJSonM_df[,2:50] <- resJSonB_df[,2:50] <- 
  resJSofM_df[,2:50] <- resJSofB_df[,2:50] <- 
  resMLonM_df[,2:50] <- resMLonB_df[,2:50] <- 
  resMLofM_df[,2:50] <- resMLofB_df[,2:50] <- 
  resMMonM_df[,2:50] <- resMMonB_df[,2:50] <- 
  resMMofM_df[,2:50] <- resMMofB_df[,2:50] <- 
  resBBon_df[,2:50] <- resBBof_df[,2:50] <-
  resHNon_df[,2:50]	<- resHNof_df[,2:50]	<- NA


for(i in 1:N){
  id				<- which(mydat$id[i]==id_df$id)
  id_row			<- id_df[id,] # select working row
  
  if(id_row$n_j == 0 | (j < 3 | n == j)){
    id_row$mu_j_js_morris = pbar 
    id_row$mu_j_js_brown = pbar											  
    id_row$mu_j_ml_morris = pbar 
    id_row$mu_j_ml_brown = pbar 
    id_row$mu_j_mm_morris = pbar 
    id_row$mu_j_mm_brown = pbar	  
    id_row$mu_j_nsf = pbar 
    id_row$mu_j_BB = pbar
  }
  
  
  else{
    #update mean and variance transformed probabilities: morris 
    id_row$y_j_morris <- transform_y(n_bar = nbar, average_j = id_row$p_j,
                                     average = pbar, type = 'morris')
    id_row$v_j <- nbar / id_row$n_j	
    
    #update mean and variance transformed probabilities: brown 
    id_row$y_j_brown <- transform_y(average_j = id_row$p_j,
                                    n_j = id_row$n_j, type = 'brown')
    
    #update shrinkage factors
    #heuristic shrinkage factor 
    id_row$bn <- shrink_n(constant = 2, numb.obs = id_row$n_j)
    #Beta Binomial 
    id_row$bBB	<- shrink_BB(units = j, 
                            numb.obs = n, 
                            nj = id_row$n_j,
                            SumofSquares = SSbb, 
                            average = pbar, 
                            inv_nj = inversnj)	 
    #####################################################
    #Maximum likelihood Morris transformation
    id_row$bml_morris <- shrink_ml(varapp = id_row$v_j,
                                   tauml = tau_ml_hat_morris,
                                   n.id = id_row$n_j)
    
    #James Stein morris transformation
    id_row$bjs_morris	<- shrink_js(units = j, SumofSquares = SSjs_morris)
    
    #Metod of moments Morris transformation: this one should be checked
    id_row$bmm_morris <- shrink_mm(sig_j = id_row$v_j, 
                                   tau_mm = tau_mm_hat_morris)
    
    ######################################################
    #Maximum likelihood brown transformation: this one should be checked
    id_row$bml_brown <- shrink_ml(varapp = id_row$sig_j,
                                  tauml = tau_ml_hat_brown,
                                  n.id = id_row$n_j)
    #James Stein brown transformation
    id_row$bjs_brown	<- shrink_js(units = j, SumofSquares = SSjs_brown)
    
    #Metod of moments brown transformation
    
    id_row$bmm_brown <- shrink_mm(sig_j = id_row$sig_j, 
                                  tau_mm = tau_mm_hat_brown)
    ################################################################
    #predictions
    #heuristic shrinkage factor 
    id_row$mu_j_nsf <- estimate_mu(b = id_row$bn, average_j = id_row$p_j, 
                                   average = pbar, average_n = nbar, 
                                   transform = FALSE)
    #beta binomial
    id_row$mu_j_BB <- estimate_mu(b = id_row$bBB, average_j = id_row$p_j, 
                                  average = pbar, average_n = nbar, 
                                  transform = FALSE)		
    #maximum likelihood - morris
    id_row$mu_j_ml_morris	<- estimate_mu(b = id_row$bml_morris, 
                                         average_j = id_row$y_j_morris, 
                                         average_y = ybar_morris, 
                                         average = pbar, average_n = nbar, 
                                         transform = TRUE, 
                                         type = 'morris')		
    #james stein - morris
    id_row$mu_j_js_morris	<- estimate_mu(b = id_row$bjs_morris, 
                                         average_j = id_row$y_j_morris,
                                         average_y = ybar_morris, 
                                         average = pbar, average_n = nbar, 
                                         transform = TRUE, 
                                         type = 'morris')	
    #method of moments - morris: this needs to be checked
    id_row$mu_j_mm_morris	<- estimate_mu(b = id_row$bmm_morris, 
                                         average_j = id_row$y_j_morris,
                                         average_y = ybar_morris, 
                                         average = pbar, average_n = nbar, 
                                         transform = TRUE, 
                                         type = 'morris')	
    ###############################################################
    #maximum likelihood - brown : this needs to be checked
    id_row$mu_j_ml_brown	<- estimate_mu(b = id_row$bml_brown, 
                                        average_j = id_row$y_j_brown, 
                                        average_y = mu_ml_brown, 
                                        average = pbar, average_n = nbar, 
                                        transform = TRUE, 
                                        type = 'brown')		
    #james stein - brown
    id_row$mu_j_js_brown	<- estimate_mu(b = id_row$bjs_brown, 
                                        average_j = id_row$y_j_brown,
                                        average_y = mu_js_brown, 
                                        average = pbar, average_n = nbar, 
                                        transform = TRUE, 
                                        type = 'brown')
    #method of moments - brown
    id_row$mu_j_mm_brown    <- estimate_mu(b = id_row$bmm_brown, 
                                           average_j = id_row$y_j_brown,
                                           average_y = mu_mm_brown,
                                           average = pbar, average_n = nbar, 
                                           transform = TRUE, 
                                           type = 'brown')	
    
  }
  
  obs <- mydat[i, 2]
  n <- n + 1
  pbar <- update_average(average = pbar, obs = obs, numb.obs = n)
  
  if(id_row$n_j == 0){
    j <- j + 1
    id_row$v_j <- nbar
    inversnj <- update_contribution(statistic = inversnj, old_contr = 0,
                                    new_contr = 1)
    ## this is for james stein brown transformation 
    sum_inverse_sig <- update_contribution(statistic = sum_inverse_sig, 
                                           old_contr = 0, 
                                           new_contr = 4)#1/(1/(1*4))
    
  }
  id_row$n_j <- id_row$n_j + 1
  id_row$p_j <- update_average(average = id_row$p_j, obs = obs, 
                               numb.obs = id_row$n_j)
  id_row$sig_j <- 1 / (4 * id_row$n_j)
  
  
  nbar <- n / j	
  
  ###################################################################
  #update transformed mean morris 
  sbar <- update_contribution(statistic = sbar, 
                              old_contr = id_df[id, ]$y_j_morris,
                              new_contr = id_row$y_j_morris)
  ybar_morris <- sbar / j
  
  # update Beta binomial parameters
  id_row$SSadd <- id_row$n_j * (id_row$p_j - pbar)**2
  if(id_df[id, ]$n_j > 0){
    ##beta binomial 
    inversnj <- update_contribution(statistic = inversnj, 
                                    old_contr = 1 / id_df[id, ]$n_j,
                                    new_contr = 1 / id_row$n_j)
    ## this is for james stein brown transformation 
    sum_inverse_sig <- update_contribution(statistic = sum_inverse_sig, 
                                           old_contr = 1 / id_df[id, ]$sig_j, 
                                           new_contr = 1 / id_row$sig_j)
    
  }
  SSbb <- update_contribution(statistic = SSbb, 
                              old_contr = id_df[id, ]$SSadd,
                              new_contr = id_row$SSadd)
  # update Maximum likelihood parameters - Morris
  id_row$tau_j_morris <- compute_tauj(tausq = tau_ml_hat_morris,
                                      bar_y = ybar_morris,
                                      y = id_row$y_j_morris, 
                                      sig_j = id_row$v_j)
  tau_ml_hat_morris <- max(update_contribution(statistic = tau_ml_hat_morris, 
                                               old_contr = id_df[id, ]$tau_j_morris,
                                               new_contr = id_row$tau_j_morris), 0)
  # update james stein shrinkage factor - Morris
  id_row$SSjsadd_morris <- (id_row$y_j_morris - ybar_morris)**2 / id_row$v_j
  SSjs_morris <- update_contribution(statistic = SSjs_morris, 
                                     old_contr = id_df[id, ]$SSjsadd_morris, 
                                     new_contr = id_row$SSjsadd_morris)
  # update method moments parameters - Morris: this should be checked
  # update contribution numerator mu:
  
  id_row$y_morris_div_tau_sig <- id_row$y_j_morris / 
    (tau_mm_hat_morris + id_row$v_j)
  # update contribution denominator mu: 
  id_row$inv_tau_sig_j_morris <- 1/(tau_mm_hat_morris + id_row$v_j)
  
  # compute numerator mu 
  sum_y_morris_div_tau_sig <- update_contribution(
    statistic = sum_y_morris_div_tau_sig, 
    old_contr = id_df[id, ]$y_morris_div_tau_sig, 
    new_contr = id_row$y_morris_div_tau_sig) 
  # compute denominator mu
  sum_inverse_tot_var_morris <- update_contribution(
    statistic = sum_inverse_tot_var_morris, 
    old_contr = id_df[id, ]$inv_tau_sig_j_morris, 
    new_contr = id_row$inv_tau_sig_j_morris)
  # compute mu 
  mu_mm_morris <- compute_mu_MM(
    sum_std_y = sum_y_morris_div_tau_sig, 
    inverse_tot_variance = sum_inverse_tot_var_morris)
  
  # update contribution numerator first part tau
  id_row$dif_j_morris <- (id_row$y_j_morris - mu_mm_morris)**2
  # compute numerator first part tau
  sumdif_j_morris <- update_contribution(statistic = sumdif_j_morris, 
                                         old_contr = id_df[id, ]$dif_j_morris, 
                                         new_contr = id_row$dif_j_morris)
  # numerator second part tau
  sumsig_j_morris <- update_contribution(statistic = sumsig_j_morris, 
                                         old_contr = id_df[id, ]$v_j, 
                                         new_contr = id_row$v_j)
  
  tau_mm_hat_morris <- compute_tau_MM(units = j, 
                                      sumSquaredDif = sumdif_j_morris,
                                      sumVariance_j = sumsig_j_morris )
  ##########################################################################
  # update maximum likelihood parameters brown transformation: 
  #this should be checked: i used the formulae of the morris transformation 
  # to estimate tau and i dont know whether that is correct
  id_row$tau_j_brown <- compute_tauj(tausq = tau_ml_hat_brown, 
                                     bar_y = mu_ml_brown, 
                                     y = id_row$y_j_brown, 
                                     sig_j = id_row$sig_j, 
                                     type = 'brown')
  tau_ml_hat_brown <- max(update_contribution(statistic = tau_ml_hat_brown, 
                                              old_contr = id_df[id, ]$tau_j_brown,
                                              new_contr = id_row$tau_j_brown), 0)
  
  id_row$inv_tau_sig_j_brown <- 1 / (tau_ml_hat_brown + id_row$sig_j)
  
  sum_inverse_tot_var_brown <- update_contribution(
    statistic = sum_inverse_tot_var_brown, 
    old_contr = id_df[id, ]$inv_tau_sig_j_brown,
    new_contr = id_row$inv_tau_sig_j_brown)
  
  id_row$y_brown_div_tau_sig <- id_row$y_j_brown / 
    (tau_ml_hat_brown + id_row$sig_j)
  
  sum_y_brown_div_tau_sig <- update_contribution(
    statistic = sum_y_brown_div_tau_sig, 
    old_contr = id_df[id, ]$y_brown_div_tau_sig, 
    new_contr = id_row$y_brown_div_tau_sig) 
  
  mu_ml_brown <- compute_mu_MM(
    sum_std_y = sum_y_brown_div_tau_sig, 
    inverse_tot_variance = sum_inverse_tot_var_brown)
  
  # update james stein brown tranformation
  id_row$std_y_brown <- std_y_brown(y_j = id_row$y_j_brown, 
                                    sig_j = id_row$sig_j)
  
  sum_std_y_brown <- update_contribution(statistic = sum_std_y_brown, 
                                         old_contr = id_df[id, ]$std_y_brown, 
                                         new_contr = id_row$std_y_brown)
  
  
  mu_js_brown <- sum_std_y_brown / sum_inverse_sig
  id_row$SSjsadd_brown <- (id_row$y_j_brown - mu_js_brown)**2 / id_row$sig_j
  SSjs_brown <- update_contribution(statistic = SSjs_brown, 
                                    old_contr = id_df[id, ]$SSjsadd_brown,
                                    new_contr = id_row$SSjsadd_brown)
  
  
  #update method of moments - brown 
  # update contribution denominator mu: 
  id_row$inv_tau_sig_j_brown_MM <- 1 / (tau_mm_hat_brown + id_row$sig_j)
  id_row$y_brown_div_tau_sig_MM <- id_row$y_j_brown/(tau_mm_hat_brown + id_row$sig_j)
  # compute numerator mu 
  sum_y_brown_div_tau_sig_MM <- update_contribution(
    statistic = sum_y_brown_div_tau_sig_MM, 
    old_contr = id_df[id, ]$y_brown_div_tau_sig_MM, 
    new_contr = id_row$y_brown_div_tau_sig_MM) 
  # compute denominator mu
  sum_inverse_tot_var_brown_MM <- update_contribution(
    statistic = sum_inverse_tot_var_brown_MM, 
    old_contr = id_df[id, ]$inv_tau_sig_j_brown_MM, 
    new_contr = id_row$inv_tau_sig_j_brown_MM)
  # compute mu 
  mu_mm_brown <- compute_mu_MM(
    sum_std_y = sum_y_brown_div_tau_sig_MM, 
    inverse_tot_variance = sum_inverse_tot_var_brown_MM)
  
  # update contribution numerator first part tau
  id_row$dif_j_brown <- (id_row$y_j_brown - mu_mm_brown)**2
  # compute numerator first part tau
  sumdif_j_brown <- update_contribution(statistic = sumdif_j_brown, 
                                        old_contr = id_df[id, ]$dif_j_brown, 
                                        new_contr = id_row$dif_j_brown)
  # numerator second part tau
  sumsig_j_brown <- update_contribution(statistic = sumsig_j_brown, 
                                        old_contr = id_df[id, ]$sig_j, 
                                        new_contr = id_row$sig_j)
  
  tau_mm_hat_brown <- compute_tau_MM(units = j, 
                                     sumSquaredDif = sumdif_j_brown,
                                     sumVariance_j = sumsig_j_brown )
  
  #update storage dataframe 
  id_df[id, ] <- id_row
  
  if(mydat$vragenlijst_nummer[i]!=mydat$vragenlijst_nummer[i+1]| i==nrow(mydat))
  {
    
    temp_df		<- id_df[id_df$n_j > 0 & (id_df$stop>=q | is.na(id_df$stop)),]
    offline_y_j_morris	<-transform_y(n_bar = nbar, 
                                      average_j = temp_df$p_j, 
                                      average = pbar, 
                                      type = 'morris')	
    offline_y_j_brown	<-transform_y(average_j = temp_df$p_j, 
                                    n_j = temp_df$n_j,
                                    type = 'brown')	
    
    
    ##############################################################33
    ##shrinkage factors
    #james stein shrinkage factor - morris
    off_js_morris	<- offline_js(average_n = nbar, 
                                units = j,
                                id_data = temp_df, 
                                y_j = offline_y_j_morris,
                                average_y = mean(offline_y_j_morris), 
                                type = 'morris')
    #james stein shrinkage factor - brown: there seems to be a mistake here, 
    ## shrinkage factor doesnt become smaller than 1 (it returns 1 of the estimated
    ## shrinkage factor is larger than 1) the morris transformed james stein sf 
    ## does work and the online morris and online brown transformed also seem to work 
    # check equation 4.17 page 20 of brown(2008)
    mu_js_brown_offline <- sum(offline_y_j_brown / (1 / (4 * temp_df$n_j)))/
      sum(1 / (1 / (4 * temp_df$n_j)))
    
    off_js_brown	<- offline_js(id_data = temp_df,
                               units = j,
                               y_j = offline_y_j_brown,
                               average_y = mu_js_brown_offline, 
                               type = 'brown')
    #maximum likelihood shrinkage factor - morris 
    off_ml_morris	<- offline_ml(tau = tau_ml_hat_morris, 
                                average_n = nbar, 
                                numb.obs = temp_df$n_j,
                                average_y = mean(offline_y_j_morris), 
                                y_j = offline_y_j_morris, 
                                type = 'morris')
    #maximum likelihood shrinkage factor - brown, i dont know whether this is correct 
    off_ml_brown	<- offline_ml(tau = tau_ml_hat_brown,
                               numb.obs = temp_df$n_j,
                               y_j = offline_y_j_brown, 
                               type = 'brown')
    #method of moments - morris
    off_mm_morris <- offline_MM(vec_nj = temp_df$n_j, 
                                y_j = offline_y_j_morris, 
                                units = j, sig_j = nbar/temp_df$n_j)
    
    #method of moments - brown
    off_mm_brown <- offline_MM(vec_nj = temp_df$n_j, 
                               y_j = offline_y_j_brown, 
                               units = j, sig_j = 1/(4*temp_df$n_j))
    
    #beta binomial, shrinkage factor
    off_BB <- offline_BB(units = j, id_data = temp_df, average = pbar, 
                         numb.obs = temp_df$n_j)		
    #heurist shrinkage factor
    off_n		<- offline_n(constant = 2, numb.obs = temp_df$n_j)
    
    test_dat		<- mydat[1:i,]
    test_dat		<- test_dat[test_dat$stop>=q | is.na(test_dat$stop),]
    
    m1			<- glmer(response~1+(1|id), data=test_dat[1:i,], family=binomial,nAGQ=20)
    
    #################################################
    #predictions 
    ##james stein morris
    temp_df$offline_js_mu_morris <- estimate_mu(b = off_js_morris, 
                                        average_j = offline_y_j_morris, 
                                        average_y = mean(offline_y_j_morris),
                                        average = pbar, average_n = nbar, 
                                        transform = TRUE, 
                                        type = 'morris')
    ##james stein brown
    temp_df$offline_js_mu_brown <- estimate_mu(b = off_js_brown, 
                                       average_j = offline_y_j_brown, 
                                       average_y = mu_js_brown_offline, 
                                       transform = TRUE, 
                                       type = 'brown')
    ##maximum likelihood morris
    temp_df$offline_ml_mu_morris <- estimate_mu(b = off_ml_morris, 
                                        average_j = offline_y_j_morris, 
                                        average_y = mean(offline_y_j_morris),
                                        average = pbar, average_n = nbar, 
                                        transform = TRUE, 
                                        type = 'morris')
    ##maximum likelihood brown
    temp_df$offline_ml_mu_brown <- estimate_mu(b = off_ml_brown[1:j], 
                                       average_j = offline_y_j_brown, 
                                       average_y = off_ml_brown[j+1],
                                       average = pbar, average_n = nbar, 
                                       transform = TRUE, 
                                       type = 'brown')
    ##method of moments morris
    temp_df$offline_mm_mu_morris <- estimate_mu(b = off_mm_morris[1:j], 
                                        average_j = offline_y_j_morris,
                                        average_y = off_mm_morris[j+1], 
                                        average = pbar, average_n = nbar, 
                                        transform = TRUE, 
                                        type = 'morris')
    ##method of moments brown
    temp_df$offline_mm_mu_brown <- estimate_mu(b = off_mm_brown[1:j], 
                                       average_j = offline_y_j_brown,
                                       average_y = off_mm_brown[j+1], 
                                       average = pbar, average_n = nbar, 
                                       transform = TRUE, 
                                       type = 'brown')
    ##beta binomial
    temp_df$offline_BB_mu <- estimate_mu(b = off_BB, average_j = temp_df$p_j, 
                                 average = pbar, average_n = nbar, 
                                 transform = FALSE)
    ## heuristic
    temp_df$offline_n_mu  <- estimate_mu(b = off_n, average_j = temp_df$p_j, 
                                 average = pbar, average_n = nbar, 
                                 transform = FALSE)
    ##glm 
    temp_df$offline_glm_mu <- exp(as.matrix(ranef(m1)$id) + fixef(m1)) / 
      (1 + exp(fixef(m1) + as.matrix(ranef(m1)$id)))
    
    
    identify_full	<- merge(x=id_df, y=temp_df, all.x=T)
      resJSofM_df[,q+1] <- identify_full$offline_js_mu_morris 
      resJSofB_df[,q+1] <- identify_full$offline_js_mu_brown 
      resMLofM_df[,q+1] <- identify_full$offline_ml_mu_morris
      resMLofB_df[,q+1] <- identify_full$offline_ml_mu_brown
      resMMofM_df[,q+1] <- identify_full$offline_mm_mu_morris
      resMMofB_df[,q+1] <- identify_full$offline_mm_mu_brown
      resBBof_df[,q+1] <- identify_full$offline_BB_mu
      resHNof_df[,q+1] <- identify_full$offline_n_mu
      resglmof_df[,q+1] <- identify_full$offline_glm_mu
    
      resJSonM_df[,q+1] <- id_df$mu_j_js_morris
      resJSonB_df[,q+1] <- id_df$mu_j_js_brown 
      resMLonM_df[,q+1] <- id_df$mu_j_ml_morris
      resMLonB_df[,q+1] <- id_df$mu_j_ml_brown
      resMMonM_df[,q+1] <- id_df$mu_j_mm_morris
      resMMonB_df[,q+1] <- id_df$mu_j_mm_brown
      resBBon_df[,q+1] <- id_df$mu_j_BB
      resHNon_df[,q+1] <- id_df$mu_j_nsf
      

    temp_data		<- id_df[id_df$stop==(q),]
    temp_data$sumj	<- temp_data$p_j*temp_data$n_j		
    j			<- j-nrow(temp_data)
    responsesum		<- (pbar*n)- sum(temp_data$sumj)
    n			<- n-sum(temp_data$n_j)
    pbar			<- responsesum/n
    sbar			<- sbar - sum(temp_data$y_j)
    nbar			<- n/j	
    
    id_df$mu_j_js_morris[id_df$stop<=(q)] <- NA 
    id_df$mu_j_js_brown[id_df$stop<=(q)] <- NA											  
    id_df$mu_j_ml_morris[id_df$stop<=(q)] <- NA 
    id_df$mu_j_ml_brown[id_df$stop<=(q)] <- NA 
    id_df$mu_j_mm_morris[id_df$stop<=(q)] <- NA 
    id_df$mu_j_mm_brown[id_df$stop<=(q)] <- NA	  
    id_df$mu_j_nsf[id_df$stop<=(q)] <- NA 
    id_df$mu_j_BB[id_df$stop<=(q)] <- NA
    q			<- q+1
    print(c(i,q))
    
  }
}

write.table(resJSonM_df, "resJSonM.txt", sep=" ")
write.table(resJSonB_df, "resJSonB.txt", sep=" ")

write.table(resJSofM_df, "resJSofM.txt", sep=" ")
write.table(resJSofB_df, "resJSofB.txt", sep=" ")

write.table(resMLonM_df, "resMLonM.txt", sep=" ")
write.table(resMLonB_df, "resMLonB.txt", sep=" ")

write.table(resMLofM_df, "resMLofM.txt", sep=" ")
write.table(resMLofB_df, "resMLofB.txt", sep=" ")

write.table(resMMonM_df, "resMMonM.txt", sep=" ")
write.table(resMMonB_df, "resMMonB.txt", sep=" ")

write.table(resMMofM_df, "resMMofM.txt", sep=" ")
write.table(resMMofB_df, "resMMofB.txt", sep=" ")

write.table(resBBon_df ,"resBBon.txt", sep=" ")
write.table(resBBof_df ,"resBBof.txt", sep=" ")
write.table(resHNon_df ,"resHNon.txt", sep=" ")
write.table(resHNof_df ,"resHNof.txt", sep=" ")
write.table(resglmof_df, "resglmof.txt", sep=" ")
