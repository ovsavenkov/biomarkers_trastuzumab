###################################################################
###################################################################
###                                                             ###
###                Main function                                ###
###                                                             ###
###################################################################

### combined data set contains all information for each patient
### we split it randomly into five non-overlaped subsets
### fit model on 3/4 of the data and then use it on other 1/4 of the data

### randomly permute row numbers and then split


set.seed(321)
N = nrow(combined)
id_row = 1:N
row_perm = sample(id_row, N, replace = F)
row_split = split(row_perm, row_perm%%5)

### create five test data frames
test_df = list()
for(i in 1:5){
  test_df[[i]] = combined[row_split[[i]], ]
}

### create five train data frames
train_df = list()
for(i in 1:5){
  train_df[[i]] = combined[-row_split[[i]], ]
}

## we run a super learner five times on training data
## then output rusult of the model on the independent data 1/5 of
## the whole data set

## create function that takes training data set and testing data set
## and returns predicted values for training set

library(SuperLearner)
source("/Users/savenkov/Dropbox/projects/precision/code/long_df.R")
 ## function long_data returns long data frame for RMST estimation
 ## the input data frame should contain following columns:
 ## - id
 ## - time
 ## - censored
 ## - arm
 source("/Users/savenkov/Dropbox/projects/precision/code/SL_function.R")
 create.SL.xgboost = function(ntrees = 5:20*10) {
   for(mm in seq(length(ntrees))){
     eval(parse(text = paste('SL.xgboost.', ntrees[mm], '<- function(..., ntrees = ', ntrees[mm], ') SL.xgboost(..., ntrees = ntrees)', sep = '')), envir = .GlobalEnv)
   }
   invisible(TRUE)
   }

create.SL.xgboost()

predicted_sl = function(training_df, testing_df){

  train_df =  training_df
  test_df =   testing_df

  signif_genes = data.frame(gene_name = genes_endo, magnitude = NA,
                            p_value = NA, stringsAsFactors = FALSE)

  ### create list of formulas
  formulas = sapply(genes_endo,
                    function(x) as.formula(paste('Surv(rfs_days, RFSc) ~ nodalstatus +
                                                 age + tumorgrade + tumorsizecat +
                                                 hrihcstatus + arm + ',x,' + arm*',x)))

  ### create list of models based on formulas defined above
  cp_models = lapply(formulas, function(x){
                     cp_fit = coxph(x, data = train_df)
                     coef(summary(cp_fit))%>%as.matrix()
  })

  for(i in 1:dim(signif_genes)[1]){
    signif_genes[i, c(2,3)] = cp_models[[i]][9, c(2, 5)]
  }

  ### extract genes with significant interaction term
  signif_genes = signif_genes[signif_genes$p_value < 0.05, ]

  ## create data short df with baseline variables
  ## for training data
  base_vars = c("subjid", "rfs_months", "RFSc", "arm")
  base_df  = train_df[base_vars]
  colnames(base_df)  = c("id",  "time", "censored", "arm")
  row.names(base_df) = NULL

  ## create long data frame
  base_long = long_data(base_df)

  ## add important genes and baseline variables
  cov_vars  = c("subjid", "age", "arm"  ,"nodalstatus", "tumorgrade", "tumorsizecat", "hrihcstatus", signif_genes$gene_name)
  cov_df   = train_df[cov_vars]
  names(cov_df)[names(cov_df) == 'subjid'] = 'id'

  ## create full data set in long format for estimation
  full_long  = merge(base_long, cov_df)

  atrisk_df = full_long[full_long$Im == 1, ]
  atrisk_df = subset(atrisk_df, select = -c( Im, Rm, Jm)) ## id?
  ## create data frames for two arms (trt & control)
  atrisk0 = atrisk_df[atrisk_df$arm == 0, ]
  atrisk1 = atrisk_df[atrisk_df$arm == 1, ]

  ### cross-validation for Super Learner

  y_train0 = atrisk0$Lm  # response
  x_train0 = subset(atrisk0, select = -c(id, Lm, arm, tumorgrade)) # covariates

  x_train1 = subset(atrisk1, select = -c(id, Lm, arm, tumorgrade))
  y_train1 = atrisk1$Lm
  SL.library = c("SL.xgboost.50", "SL.xgboost.60",
                 "SL.xgboost.70", "SL.xgboost.80", "SL.xgboost.90",
                 "SL.xgboost.100", "SL.xgboost.110", "SL.xgboost.120",
                 "SL.xgboost.130", "SL.xgboost.140", "SL.xgboost.150",
                 "SL.xgboost.160", "SL.xgboost.170", "SL.xgboost.180",
                 "SL.xgboost.190", "SL.xgboost.200", #"SL.glmnet", "SL.glm")

   method = "method.CC_nloglik"
   family = "binomial"
   id0 = atrisk0$id
   id1 = atrisk1$id
   set.seed(1001)

   fit_SL0 = SuperLearner(Y = y_train0, X = x_train0, SL.library = SL.library,
                                                         cvControl = list(V = 5),
                                                         family = family, method = method, id = id0, verbose = TRUE)

   fit_SL1 = SuperLearner(Y = y_train1, X = x_train1, SL.library = SL.library,
                                                         cvControl = list(V = 5),
                                                         family = family, method = method, id = id1, verbose = TRUE)




   ## predict based on the test data
   base_test  = test_df[base_vars]
   colnames(base_test)  = c("id",  "time", "censored", "arm")
   row.names(base_test) = NULL

   long_test = long_data(base_test)
   cov_vars  = c("subjid", "age", "arm"  ,"nodalstatus", "tumorgrade", "tumorsizecat", "hrihcstatus", signif_genes$gene_name)
   cov_test   = test_df[cov_vars]
   names(cov_test)[names(cov_test) == 'subjid'] = 'id'

   ## create full data set in long format for estimation
   full_test  = merge(long_test, cov_test)

   f_test = full_test[, -c(1, 3, 4, 5, 6, 8, 10)]

   pred_SL0_test = predict(fit_SL0, data.matrix(f_test))

   pred0_test = pred_SL0_test$pred
   pred_SL1_test = predict(fit_SL1, data.matrix(f_test))
   pred1_test = pred_SL1_test$pred
   haz_test = data.frame(id = full_test$id, time = full_test$m,
                          h0 = pred0_test, h1 = pred1_test)
   haz_test$h0_neg = 1 - haz_test$h0
   haz_test$h1_neg = 1 - haz_test$h1
   haz_test$S0  = ave(haz_test$h0_neg, by = haz_test$id, FUN = cumprod)
   haz_test$S1 = ave(haz_test$h1_neg, by = haz_test$id, FUN = cumprod)

   haz_test$S_diff = haz_test$S1 - haz_test$S0
   return(haz_test)
  }

responders_pred = list()
for (i in 1:5){
    responders_pred[[i]] = predicted_sl(train_df[[i]], test_df[[i]])
}

resps = do.call("rbind", responders_pred)

# find responders and nonresponders in the test data
resp_nonresp = aggregate(S_diff ~ id, resps, sum)
resp_nonresp$responder = ifelse(resp_nonresp$S_diff > 0, 1, 0)

resp_nonresp$S_diff = NULL

resp_id = resp_nonresp[resp_nonresp$responder == 1, ]$id
nonresp_id = resp_nonresp[resp_nonresp$responder == 0, ]$id

resp_test = resps[resps$id%in% resp_id, ]
nonresp_test = resps[resps$id%in% nonresp_id, ]

resp_avg  =  aggregate(cbind(S0, S1) ~ time, resp_test, mean)
nonresp_avg = aggregate(cbind(S0, S1) ~ time, nonresp_test, mean)

## create long data frames for ggplot2

resp_long = melt(resp_avg, id.vars=c("time"), variable.name = "trt",
                 value.name = "surv_est")
nonresp_long = melt(nonresp_avg, id.vars=c("time"), variable.name = "trt",
                    value.name = "surv_est")


## Create Figure 2 in the paper

resp_p  = ggplot(resp_long, aes(time, surv_est, color = trt)) + geom_line() +
          ylim(c(0,1)) +
          coord_cartesian(xlim = c(0, 150)) +
          ggtitle("Responders") +
          ylab("Estimated Survival Probability") +
          xlab("Time in months") +
          scale_color_discrete(name="Treatment\nGroup",
                               breaks=c("S0", "S1"),
                               labels=c("Control", "Trastuzumab")) +
          theme_bw() +
          theme(legend.position = c(0.25, 0.2))



nonresp_p = ggplot(nonresp_long, aes(time, surv_est, color = trt)) + geom_line() +
            ylim(c(0,1)) +
            coord_cartesian(xlim = c(0, 150)) +
            ggtitle("Non-responders") +
            ylab("Estimated Survival Probability") +
            xlab("Time in months") +
            scale_color_discrete(name="Treatment\nGroup",
                       breaks=c("S0", "S1"),
                       labels=c("Control", "Trastuzumab")) +
            theme_bw() +
            theme(legend.position = c(0.25, 0.2))

gene_plot = grid.arrange(nonresp_p, resp_p, ncol = 2)





