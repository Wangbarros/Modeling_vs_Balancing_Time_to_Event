rm(list=ls())
source("functions.R")
set.seed(5324)

n_rep=1000
n = 1500
m_model = 1

# This is just to choose which algorithms will be run. CBPS is by far the slowest and sometimes you may want faster results
do_ps = TRUE
do_raw = TRUE
do_glm = TRUE
do_calibration = TRUE
do_entropy = TRUE # entropy and calibration give the same results
do_sbw = TRUE
do_CBPS = TRUE # This should be false when doing overspecification

#This is just related to the name of the file that will be generated.
state = '_misspec_'

type_censoring = 'uniform'
overspecification = FALSE

k_list = c(1)
targets = c(0.8)
sizes = c(1500)
censoring_ps = c(0.4)
scenario_list = c('B')
#scenario_list = c('G','F','E','D','C','B','A')


for (k_model in k_list) {
for (scenario in scenario_list){
   for (censoring_p in censoring_ps){
   for (n in sizes){
   for (target_beta in targets) {
  
  print(n)
  
  results <- matrix(NA, nrow = n_rep, ncol = 38)
  colnames(results) = c('beta_PS', 'exp_beta_PS', 'ciL_PS','ciU_PS','pval_PS',
                        'beta_GLM', 'exp_beta_GLM', 'ciL_GLM','ciU_GLM','pval_GLM',
                        'beta_raw', 'exp_beta_raw', 'ciL_raw','ciU_raw','pval_raw',
                        'beta_ebal', 'exp_beta_ebal', 'ciL_ebal','ciU_ebal','pval_ebal',
                        'beta_cal', 'exp_beta_cal', 'ciL_cal','ciU_cal','pval_cal',
                        'beta_CBPS', 'exp_beta_CBPS', 'ciL_CBPS','ciU_CBPS','pval_CBPS',
                        'beta_SBW', 'exp_beta_SBW', 'ciL_SBW','ciU_SBW','pval_SBW','sbw_tol',
                        'beta_cond','target'
  )
  
  n_errors = 0
  
  i = 1
  while(i <= n_rep){
    
        if (scenario == 'A'){
          non_linear = 0
          non_additive = 0
        }
        if (scenario == 'B'){
          non_linear = 1
          non_additive = 0
        }
        if (scenario == 'C'){
          non_linear = 3
          non_additive = 0
        }
        if (scenario == 'D'){
          non_linear = 0
          non_additive = 4
        }
        if (scenario == 'E'){
          non_linear = 1
          non_additive = 4
        }
        if (scenario == 'F'){
          non_linear = 0
          non_additive = 10
        }
        if (scenario == 'G'){
          non_linear = 3
          non_additive = 10
        }
      print(paste("i =", i))
      
      sim_res = misspec_surv_data(n = n,
                                    scenario = scenario, m = 1, k = k_model, target_beta = target_beta,
                                    censoring_p = censoring_p, type_censoring = type_censoring)
      
      data = sim_res$data
      sum(data$event)
      sum(data$Treatment)
      
      # This is done for overspecification 
      # The idea is to generate all possible square and interactions
      if (overspecification){
        new_matrix = model.matrix( ~(.-1-event-Y-Treatment)^2, data=data)
        new_matrix = cbind(data$Y, data$event, data$Treatment, new_matrix)
        new_matrix = as.data.frame(new_matrix)
        new_matrix$'X2:X2' = new_matrix$X2*new_matrix$X2
        new_matrix$'X4:X4' = new_matrix$X4*new_matrix$X4
        new_matrix$'X7:X7' = new_matrix$X7*new_matrix$X7
        new_matrix$'X10:X10' = new_matrix$X10*new_matrix$X10
        colnames(new_matrix)[1] = 'Y'
        colnames(new_matrix)[2] = 'event'
        colnames(new_matrix)[3] = 'Treatment'
        data = new_matrix
      }
      
      #tapply(data$Y, data$event, summary)
      
      betaT = sim_res$betaT
      results[i, 37] <- betaT
      results[i, 38] <- target_beta
      t_id <- which(data$Treatment == 1)
      c_id <- which(data$Treatment == 0)
      data <- data[c(t_id, c_id), ]
      
      # For all algorithms the steps are the same:
      # 1) Generate weights
      # 2) Use the weights to try to estimate the MHR with a Cox regression
      
      #####################
      #IPTW-LASSO
      #####################
      if(do_ps){
        glmnet_ps <- ps_est_func(data)
        weights <- ifelse(data$Treatment == 1, 1, glmnet_ps/(1-glmnet_ps))
        ps_weighting_data <- cbind(data, weights)
        result_line = weights_cox_reg(weighting_data = ps_weighting_data)
        result_line = unlist(result_line)
        results[i, 1:5] <- c(result_line[1:5])
      }
      #####################
      #IPTW - GLM
      #####################
      if(do_glm){
        glm_res = glm(Treatment ~. -Y -event, data = data, family = 'binomial')
        weights = ifelse(data$Treatment == 1, 1, glm_res$fitted.values/(1-glm_res$fitted.values))
        glm_weighting_data <- cbind(data, weights)
        result_line = weights_cox_reg(glm_weighting_data)
        result_line = unlist(result_line)
        results[i, 6:10] <- c(result_line[1:5])
      }
      #####################
      #RAW
      #####################
      if(do_raw){
      data_raw = data
      data_raw$weights = rep(1,nrow(data_raw)) # all weights = 1, same as no weights
      result_line = weights_cox_reg(weighting_data = data_raw)
      result_line = unlist(result_line)
      results[i, 11:15] <- c(result_line[1:5])
      }
      #####################
      #ENTROPY
      #####################
      if(do_entropy){
      col_notkeep = which(colnames(data)  %in% c('Y', 'event', 'Treatment'))
      
      ebal_res = ebalance(Treatment = data$Treatment, X = data[, -col_notkeep], )
      ebal_weights = c(rep(1, sum(data$Treatment)), ebal_res$w)
      ebal_data = data
      ebal_data$weights = ebal_weights
      result_line = weights_cox_reg(weighting_data = ebal_data)
      result_line = unlist(result_line)
      results[i, 16:20] <- c(result_line[1:5])
      }
      #####################
      #CALIBRATION
      #####################      
      if(do_calibration){
      # Calibration uses only the covariates, so we remove the rest
      col_notkeep = which(colnames(data)  %in% c('Y', 'event', 'Treatment'))
      
      res_cal = ATE::ATE(Y = data$Y, Ti = data$Treatment, X = data[, -col_notkeep], ATT = TRUE)
      cal_weights = ifelse(data$Treatment == 1, 1/length(t_id), res_cal$weights.q)
      cal_data = data
      cal_data$weights = cal_weights
      result_line = weights_cox_reg(weighting_data = cal_data)
      result_line = unlist(result_line)
      results[i, 21:25] <- c(result_line[1:5])
      }
      #####################
      #CBPS - non parametric
      #####################  
      if(do_CBPS){
      # CBPS doesn't want to receive time to event and event markers, so we remove
      col_notkeep = which(colnames(data)  %in% c('Y', 'event'))
      CBPS_data = data[,-col_notkeep]
      
      CBPS_data$Treatment = as.factor(CBPS_data$Treatment)
      CBPS_res = npCBPS(Treatment ~ ., data = CBPS_data) #note the use of npCBPS instead of just CBPS
      CBPS_data_weighted = data
      CBPS_data_weighted$weights = CBPS_res$weights
      
      result_line = weights_cox_reg(weighting_data = CBPS_data_weighted)
      result_line = unlist(result_line)
      results[i, 26:30] <- c(result_line[1:5])
      }
      #####################
      #SBW
      #####################
      if(do_sbw){
        ###### WITH ZUBIZARRETA ALG
        
        # Just getting the covariates names and getting a dataset with only Treatment and covariates to throw into the function
        col_remove <- which(colnames(data)  %in% c('Y', 'event', 'Treatment'))
        mom_covs <- data[ , -c(col_remove)]
        col_remove <- which(colnames(data)  %in% c('Y', 'event'))
        
        # This is a bunch of parameters, but it is all in fact the "standard" parameters here
        sbw_weighting_object <- sbw::sbw(dat = data[ , -c(col_remove)],
                                         ind = 'Treatment', out = NULL,
                                         bal = list(bal_cov = names(mom_covs), bal_alg = TRUE,
                                                    bal_std = 'group',
                                                    bal_sam = 1000), wei = list(wei_sum = TRUE, wei_pos =TRUE),
                                         sol = list(sol_nam = "quadprog", sol_dis = FALSE),
                                         par = list(par_est = "att", par_tar = NULL))
        
        col_keep = which(colnames(data)  %in% c('Y', 'event'))
        sbw_weighting_data <- cbind(data[, col_keep], sbw_weighting_object$dat_weights)
        names(sbw_weighting_data)[which(names(sbw_weighting_data) == "sbw_weights")] = 'weights'
        
        result_line = weights_cox_reg(sbw_weighting_data)
        result_line = unlist(result_line)
        results[i, 31:35] <- c(result_line[1:5])
        results[i, 36] = sbw_weighting_object$balance_parameters$bal_tol_ori
      }
    print(paste0(state,'sbw_results_n', n, '_nrep', n_rep,
                 '_MHR', target_beta, 
                 '_nadditive', non_additive, '_nlinear', non_linear,
                 '_censoring', censoring_p,'_overspecification', overspecification,
                 '_m', m_model, '_k',k_model))
    
    i = i+1
  }
  
  results = as.data.frame(results)

  results$relative_ps = ((results$exp_beta_PS - results$target)/results$target)
  results$relative_sbw = ((results$exp_beta_SBW - results$target)/results$target)
  results$relative_glm = ((results$exp_beta_GLM - results$target)/results$target)
  results$relative_raw = ((results$exp_beta_raw - results$target)/results$target)
  results$relative_ebal = ((results$exp_beta_ebal - results$target)/results$target)
  results$relative_cal = ((results$exp_beta_cal - results$target)/results$target)
  results$relative_cbps = ((results$exp_beta_CBPS - results$target)/results$target)
  
  results$se_ps = (log(results$ciU_PS) - results$beta_PS)/1.96
  results$se_sbw = (log(results$ciU_SBW) - results$beta_SBW)/1.96
  results$se_glm = (log(results$ciU_GLM) - results$beta_GLM)/1.96
  results$se_raw = (log(results$ciU_raw) - results$beta_raw)/1.96
  results$se_ebal = (log(results$ciU_ebal) - results$beta_ebal)/1.96
  results$se_cal = (log(results$ciU_cal) - results$beta_cal)/1.96
  results$se_cbps = (log(results$ciU_CBPS) - results$beta_CBPS)/1.96
  
  #results_tols$rob_se = (log(results_tols$ciU_SBW) - results_tols$beta_SBW)/1.96
  save(results, file = paste0('Results/',state,'sbw_results_n', n, '_nrep', n_rep,
                              '_MHR', target_beta, 
                              '_nadditive', non_additive, '_nlinear', non_linear,
                              '_censoring', censoring_p, '_overspecification', overspecification,
                              '_m', m_model, '_k',k_model,
                              '_.RData'))
  
 
  
   }
   }
   }
   }
}
