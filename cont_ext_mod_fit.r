### Load libraries
source("load_packages.r")

### Load simulation data
if(!exists("sim_res")) sim_res = readRDS("sim_res.RDS")

### Summary statistics for replicates
ext.s7.df = data.frame(
  sel = sapply(sim_res, with, selstrength), 
  ext = sapply(sim_res, with, extinct), 
  sdopt =  sapply(sim_res, with, sdopt), 
  meanopt = sapply(sim_res, with, meanopt),
  mean.pheno = sapply(sim_res, with, mean(phenomean[(250-10):250])),
  ext_year = sapply(sim_res, with, yearextinct),
  age.sex.mat = sapply(sim_res, with, age.sex.mat),
  avg.surv = sapply(sim_res, with,avg.surv),
  sons.mean = sapply(sim_res, with,sons.mean),
  cor.all.est = sapply(sim_res, with, cor.all.est),
  cor.all.p = sapply(sim_res, with, cor.all.p), 
  cat.freq = sapply(sim_res, with, cat.freq),
  catastr.mort = sapply(sim_res, with,catastr.mort)
)
ext.s7.df$sim_rep = 1:nrow(ext.s7.df)


ext.s7.df$year_bc = NA

for (i in 1:nrow(ext.s7.df)) {
  
  if (ext.s7.df$ext[i] == 1) {
    
    cat_happ = which(sim_res[[i]]$catastr.vett == 1) 
    cat_happ_bef = cat_happ[which(cat_happ < ext.s7.df$ext_year[i])]
    if (length(cat_happ) > 0) {
      cat_close = min(ext.s7.df$ext_year[i] - cat_happ_bef)}
    if (length(cat_happ) == 0) {
      cat_close = NA}
    
    ext.s7.df$year_bc[i] = cat_close 
    
  }
  
  
}


### I select the extinct replicates for the test data set

seed_num = 50


set.seed(seed_num)

ext_df = filter(ext.s7.df, ext == 1)

s_w = 10 # length of the sample window (between 5 and 20)

ext_df = filter(ext_df, ext_year > (s_w + 2))

num_ext_out = 2300


set.seed(seed_num)

oft = sample(1:nrow(ext_df), size = num_ext_out, replace = F)
out_for_test = ext_df[oft,]
ext_df = ext_df[-oft,]


ext_df_3 = rbind(ext_df,ext_df,ext_df)



ext_df_3$sampl_bef_ext = NA

set.seed(seed_num)

ext_df_3$sampl_bef_ext = sample(x = 1:s_w ,size = nrow(ext_df_3), replace = T)

ext_df_3$sampl_end = ext_df_3$ext_year - ext_df_3$sampl_bef_ext
ext_df_3$sampl_end = ifelse(ext_df_3$sampl_end>=s_w ,ext_df_3$sampl_end,s_w)
ext_df_3$sampl_beg = ext_df_3$sampl_end - (s_w -1)

ext_df_3$min_pop_win = NA
ext_df_3$max_pop_win = NA
ext_df_3$mean_pop_win = NA
ext_df_3$max_min_pop_win = NA
ext_df_3$cat_win = NA
ext_df_3$max_opt_win = NA
ext_df_3$min_opt_win = NA
ext_df_3$max_min_opt_win = NA
ext_df_3$max_pheno_dist_win = NA
ext_df_3$mean_pheno_dist_win = NA


for (i in 1:nrow(ext_df_3)) {
  if (i%%100 == 0){ print(i)}
  targ = ext_df_3$sim_rep[i]
  if(ext_df_3$ext_year[i] > s_w) {
    ext_df_3$min_pop_win[i] = min(sim_res[[targ]]$popsize.post[ext_df_3$sampl_beg[i]:ext_df_3$sampl_end[i]])
    ext_df_3$max_pop_win[i] = max(sim_res[[targ]]$popsize.post[ext_df_3$sampl_beg[i]:ext_df_3$sampl_end[i]])
    ext_df_3$max_min_pop_win[i] = ext_df_3$max_pop_win[i] - ext_df_3$min_pop_win[i]
    ext_df_3$mean_pop_win[i] = mean(sim_res[[targ]]$popsize.post[ext_df_3$sampl_beg[i]:ext_df_3$sampl_end[i]])
    ext_df_3$cat_win[i] = ifelse(sum(which(sim_res[[targ]]$catastr.vett == 1) %in% (ext_df_3$sampl_beg[i]:ext_df_3$sampl_end[i]) > 0),"y","n")
    ext_df_3$max_opt_win[i] = max(sim_res[[targ]]$optimum[ext_df_3$sampl_beg[i]:ext_df_3$sampl_end[i]])
    ext_df_3$min_opt_win[i] = min(sim_res[[targ]]$optimum[ext_df_3$sampl_beg[i]:ext_df_3$sampl_end[i]])
    ext_df_3$max_min_opt_win[i] = ext_df_3$max_opt_win[i] - ext_df_3$min_opt_win[i] 
    ext_df_3$max_pheno_dist_win[i] = max(abs(sim_res[[targ]]$phenomean[ext_df_3$sampl_beg[i]:ext_df_3$sampl_end[i]] - sim_res[[targ]]$optimum[ext_df_3$sampl_beg[i]:ext_df_3$sampl_end[i]]))
    ext_df_3$mean_pheno_dist_win[i] = mean(abs(sim_res[[targ]]$phenomean[ext_df_3$sampl_beg[i]:ext_df_3$sampl_end[i]] - sim_res[[targ]]$optimum[ext_df_3$sampl_beg[i]:ext_df_3$sampl_end[i]]))
  }
}


### For out for test ####

out_for_test$sampl_bef_ext = NA

set.seed(seed_num)

out_for_test$sampl_bef_ext = sample(x = 1:s_w ,size = nrow(out_for_test), replace = T)

out_for_test$sampl_end = out_for_test$ext_year - out_for_test$sampl_bef_ext
out_for_test$sampl_end = ifelse(out_for_test$sampl_end>=s_w ,out_for_test$sampl_end,s_w)
out_for_test$sampl_beg = out_for_test$sampl_end - (s_w -1)

out_for_test$min_pop_win = NA
out_for_test$max_pop_win = NA
out_for_test$mean_pop_win = NA
out_for_test$max_min_pop_win = NA
out_for_test$cat_win = NA
out_for_test$max_opt_win = NA
out_for_test$min_opt_win = NA
out_for_test$max_min_opt_win = NA
out_for_test$max_pheno_dist_win = NA
out_for_test$mean_pheno_dist_win = NA


for (i in 1:nrow(out_for_test)) {
  if (i%%100 == 0){ print(i)}
  targ = out_for_test$sim_rep[i]
  if(out_for_test$ext_year[i] > s_w) {
    out_for_test$min_pop_win[i] = min(sim_res[[targ]]$popsize.post[out_for_test$sampl_beg[i]:out_for_test$sampl_end[i]])
    out_for_test$max_pop_win[i] = max(sim_res[[targ]]$popsize.post[out_for_test$sampl_beg[i]:out_for_test$sampl_end[i]])
    out_for_test$max_min_pop_win[i] = out_for_test$max_pop_win[i] - out_for_test$min_pop_win[i]
    out_for_test$mean_pop_win[i] = mean(sim_res[[targ]]$popsize.post[out_for_test$sampl_beg[i]:out_for_test$sampl_end[i]])
    out_for_test$cat_win[i] = ifelse(sum(which(sim_res[[targ]]$catastr.vett == 1) %in% (out_for_test$sampl_beg[i]:out_for_test$sampl_end[i]) > 0),"y","n")
    out_for_test$max_opt_win[i] = max(sim_res[[targ]]$optimum[out_for_test$sampl_beg[i]:out_for_test$sampl_end[i]])
    out_for_test$min_opt_win[i] = min(sim_res[[targ]]$optimum[out_for_test$sampl_beg[i]:out_for_test$sampl_end[i]])
    out_for_test$max_min_opt_win[i] = out_for_test$max_opt_win[i] - out_for_test$min_opt_win[i] 
    out_for_test$max_pheno_dist_win[i] = max(abs(sim_res[[targ]]$phenomean[out_for_test$sampl_beg[i]:out_for_test$sampl_end[i]] - sim_res[[targ]]$optimum[out_for_test$sampl_beg[i]:out_for_test$sampl_end[i]]))
    out_for_test$mean_pheno_dist_win[i] = mean(abs(sim_res[[targ]]$phenomean[out_for_test$sampl_beg[i]:out_for_test$sampl_end[i]] - sim_res[[targ]]$optimum[out_for_test$sampl_beg[i]:out_for_test$sampl_end[i]]))
  }
}




surv_df = filter(ext.s7.df, ext == 0)

surv_df$sampl_bef_ext = NA
surv_df$sampl_bef_ext = NA

set.seed(seed_num)

surv_df$sampl_end = sample(x = 20:240,size = nrow(surv_df), replace = T)
surv_df$sampl_beg = surv_df$sampl_end - (s_w - 1)

surv_df$min_pop_win = NA
surv_df$max_pop_win = NA
surv_df$mean_pop_win = NA
surv_df$max_min_pop_win = NA
surv_df$cat_win = NA
surv_df$max_opt_win = NA
surv_df$min_opt_win = NA
surv_df$max_min_opt_win = NA


for (i in 1:nrow(surv_df)) {
  if (i%%100 == 0){ print(i)}
  targ = surv_df$sim_rep[i]
  surv_df$min_pop_win[i] = min(sim_res[[targ]]$popsize.post[surv_df$sampl_beg[i]:surv_df$sampl_end[i]])
  surv_df$max_pop_win[i] = max(sim_res[[targ]]$popsize.post[surv_df$sampl_beg[i]:surv_df$sampl_end[i]])
  surv_df$max_min_pop_win[i] = surv_df$max_pop_win[i] - surv_df$min_pop_win[i]
  surv_df$mean_pop_win[i] = mean(sim_res[[targ]]$popsize.post[surv_df$sampl_beg[i]:surv_df$sampl_end[i]])
  surv_df$cat_win[i] = ifelse(sum(which(sim_res[[targ]]$catastr.vett == 1) %in% (surv_df$sampl_beg[i]:surv_df$sampl_end[i]) > 0),"y","n")
  surv_df$max_opt_win[i] = max(sim_res[[targ]]$optimum[surv_df$sampl_beg[i]:surv_df$sampl_end[i]])
  surv_df$min_opt_win[i] = min(sim_res[[targ]]$optimum[surv_df$sampl_beg[i]:surv_df$sampl_end[i]])
  surv_df$max_min_opt_win[i] = surv_df$max_opt_win[i] - surv_df$min_opt_win[i] 
  surv_df$max_pheno_dist_win[i] = max(abs(sim_res[[targ]]$phenomean[surv_df$sampl_beg[i]:surv_df$sampl_end[i]] - sim_res[[targ]]$optimum[surv_df$sampl_beg[i]:surv_df$sampl_end[i]]))
  surv_df$mean_pheno_dist_win[i] = mean(abs(sim_res[[targ]]$phenomean[surv_df$sampl_beg[i]:surv_df$sampl_end[i]] - sim_res[[targ]]$optimum[surv_df$sampl_beg[i]:surv_df$sampl_end[i]]))
}


pred_df = rbind(ext_df_3,surv_df)


######## start preparation of dataset for model training and testing 
######## 
require(rms)
require(sampling)

pred_df$test = 0

pos_0 = which(pred_df$ext == 0)

set.seed(seed_num)  #make the analysis repeatable

pos_oft = sample(pos_0, size = num_ext_out, replace = F)
pred_df$test[pos_oft] = 1

out_for_test$test = 1

test_pred_df = rbind(pred_df, out_for_test)
test_pred_df$cont_cons = 1:nrow(test_pred_df)

saveRDS(test_pred_df,"test_pred_df.RDS")


### standardize the continuous variables

ext.stand.df = scale(test_pred_df[,c("min_pop_win","max_pop_win", "max_min_pop_win", "mean_pop_win", "max_opt_win", "max_min_opt_win","sons.mean","age.sex.mat", "max_pheno_dist_win","mean_pheno_dist_win")], center = TRUE, scale = TRUE) #standardize predictors

ext.stand.df = as.data.frame(ext.stand.df) #transform in data frame

ext.stand.df = cbind(ext.stand.df,test_pred_df$ext, test_pred_df$cat_win,test_pred_df$sel)  #add a column with extinction (1/0) coming from test_pred_df

colnames(ext.stand.df)[ncol(ext.stand.df)-2] = "ext"
colnames(ext.stand.df)[ncol(ext.stand.df)-1] = "cat_win"#rename the column
colnames(ext.stand.df)[ncol(ext.stand.df)] = "sel"
ext.stand.df$sel = as.factor(ext.stand.df$sel)

# ext.stand.df %>% left_join(., select(test_pred_df, cont_cons, test,sim_rep,year_bc, e))


ext.stand.df$test = test_pred_df$test

ext.stand.df$sim_rep = test_pred_df$sim_rep
ext.stand.df$year_bc = test_pred_df$year_bc
ext.stand.df$sampl_bef_ext = test_pred_df$sampl_bef_ext
ext.stand.df$sampl_end = test_pred_df$sampl_end
ext.stand.df$sampl_beg = test_pred_df$sampl_beg 
ext.stand.df$cont_cons = test_pred_df$cont_cons
ext.stand.df$ext_year = test_pred_df$ext_year


######## dataset without standardization #####

# ext.nostand.df = test_pred_df[,c("sel","cat.freq","pop.ext10.rand","addvar.mean.10.rand")]
# 
# ext.nostand.df = as.data.frame(ext.nostand.df) #transform in data frame
# 
# ext.nostand.df = cbind(ext.nostand.df,test_pred_df$ext)  #add a column with extinction (1/0) coming from test_pred_df
# 
# colnames(ext.nostand.df)[ncol(ext.nostand.df)] = "ext"  #rename the column

########

set.seed(seed_num)  #make the analysis repeatable



ext.stand.df.val = filter(ext.stand.df,test == 1)  #validation dataset
ext.stand.df.cal = filter(ext.stand.df,test == 0) #calibration dataset

saveRDS(ext.stand.df.val,"ext.stand.df.val.RDS")
saveRDS(ext.stand.df.cal,"ext.stand.df.cal.RDS")


############## end of preparation of dataset ###############

#### Full GLM
ext.lrm.full =  with(ext.stand.df.cal, lrm(ext ~  sel + min_pop_win +  mean_pop_win + max_opt_win + sons.mean + age.sex.mat + cat_win + max_pheno_dist_win + mean_pheno_dist_win))  #full glm (using Harrell's rms package) model
saveRDS(ext.lrm.full,"ext.lrm.full.RDS")

dd <- datadist(ext.stand.df.cal)
options(datadist="dd")

### Reduced GLM
### 
ext.lrm.base =  with(ext.stand.df.cal, lrm(ext ~   sel + max_opt_win + sons.mean + age.sex.mat + cat_win + max_pheno_dist_win + mean_pheno_dist_win, x = TRUE, y = TRUE))  #full glm (using Harrell's rms package) model
# 
saveRDS(ext.lrm.base,"ext.lrm.base.RDS")

### Full GAM
ext.gam.full = gam(as.factor(ext) ~ as.factor(sel)  + s(min_pop_win) + s(mean_pop_win) + s(max_opt_win)  + sons.mean + age.sex.mat + as.factor(cat_win) + 
                     s(max_pheno_dist_win) + s(mean_pheno_dist_win), family = "binomial",data = ext.stand.df.cal,na.action = "na.fail", select = TRUE)

saveRDS(ext.gam.full,"ext.gam.full.RDS")

### Reduce GAM
ext.gam.base = gam(as.factor(ext) ~ as.factor(sel) + s(max_opt_win)  + sons.mean + age.sex.mat + as.factor(cat_win) + 
                     s(max_pheno_dist_win) + s(mean_pheno_dist_win), family = "binomial",data = ext.stand.df.cal,na.action = "na.fail", select = TRUE)

saveRDS(ext.gam.base,"ext.gam.base.RDS")


# Full RF in parallel (better to run in the shell)
set.seed(seed_num)
library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1

fitControl <- trainControl(
  ## Repeated 5–fold CV
  method = "repeatedcv",
  number = 2,
  ## repeated 10 times
  repeats = 2,
  allowParallel = TRUE,
  #verboseIter = TRUE,
  returnResamp = "all",search = "random")

library(doParallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)

ext.rf.full <- train(as.factor(ext) ~  sel + min_pop_win +  mean_pop_win + max_opt_win + sons.mean + age.sex.mat + cat_win + max_pheno_dist_win + mean_pheno_dist_win,
                     data = ext.stand.df.cal,
                     method = 'ranger',
                     
                     # should be set high at least p/3
                     tuneLength = 10,
                     trControl = fitControl,
                     
                     ## parameters passed onto the ranger function
                     # the bigger the better.
                     num.trees = 1000,
                     importance = "permutation")


saveRDS(ext.rf.full,"ext.rf.full.RDS")
stopCluster(cl)
registerDoSEQ()

# Reduced RF in parallel (better to run in the shell)

set.seed(seed_num)
library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1

library(doParallel)
# create the cluster for caret to use
cl <- makePSOCKcluster(no_cores)
registerDoParallel(cl)


ext.rf.base <- train(as.factor(ext) ~  sel + max_opt_win + sons.mean + age.sex.mat + cat_win + max_pheno_dist_win + mean_pheno_dist_win,
                    data = ext.stand.df.cal,
                    method = 'ranger',
                    
                    # should be set high at least p/3
                    tuneLength = 10,
                    trControl = fitControl,
                    
                    ## parameters passed onto the ranger function
                    # the bigger the better.
                    num.trees = 1000,
                    importance = "permutation")

saveRDS(ext.rf.base,"ext.rf.base.RDS")
stopCluster(cl)
registerDoSEQ()

