### Load libraries
source("load_packages.r")


#Load data (or source("all_count.r"), but it takes some time to run the routine)
pers_df = readRDS("pers_df.RDS")
list_tb = readRDS("list_tb.RDS")

seed_num = 110

set.seed(seed_num)

### change to factors some of the predictors
pers_df$sdopt = as.factor(pers_df$sdopt)
pers_df$cat.freq = as.factor(pers_df$cat.freq)
pers_df$catastr.mort = as.factor(pers_df$catastr.mort)
pers_df$age.sex.mat = as.factor(pers_df$age.sex.mat)
pers_df$sons.mean = as.factor(pers_df$sons.mean)
pers_df$meanopt = as.factor(pers_df$meanopt)
pers_df$sel = as.factor(pers_df$sel)

### 20% randomly chosen for test data set
val_row = sample(1:nrow(pers_df), size = round(0.2 * nrow(pers_df)), replace = F)

all_cal_df = pers_df[-val_row,]
all_val_df = pers_df[val_row,]

all_cal_df = filter(all_cal_df, !is.na(tot_all)) # training data
all_val_df = filter(all_val_df, !is.na(tot_all)) # test data

saveRDS(all_cal_df,"all_cal_df.RDS")
saveRDS(all_val_df,"all_val_df.RDS")

### Full GAM model

all.gam.full = mgcv::gam(tot_all ~ meanopt + sdopt + sel +  s(tot_pop_un_2, k = 4) + cat.freq + catastr.mort + age.sex.mat + sons.mean, data = all_cal_df)

saveRDS(all.gam.full,"all.gam.full.RDS")

### Full GLM model 
all.lm.full = lm(tot_all ~ meanopt + sdopt + sel + tot_pop_un_1 + cat.freq + catastr.mort + age.sex.mat + sons.mean, data = all_cal_df)

saveRDS(all.lm.full, "all.lm.full.RDS")


### Reduced GLM model
all.lm.red = lm(tot_all ~ meanopt + sdopt + sel + cat.freq + catastr.mort + age.sex.mat + sons.mean, data = all_cal_df)

saveRDS(all.lm.red, "all.lm.red.RDS")

### Full RF model
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


data_y = all_cal_df[,"tot_all"]
data_x = all_cal_df[, c("meanopt","sdopt", "sel","cat.freq","catastr.mort","age.sex.mat","sons.mean","tot_pop_un_1")]

all.rf.full <- train(x = data_x, y = data_y,
                     method = 'ranger',
                     # should be set high at least p/3
                     tuneLength = 10,
                     trControl = fitControl,
                     
                     ## parameters passed onto the ranger function
                     # the bigger the better.
                     num.trees = 1000,
                     importance = "permutation")


saveRDS(all.rf.full,"all.rf.full.RDS")
stopCluster(cl)
registerDoSEQ()


### Reduced RF model
set.seed(seed_num)
library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1

library(doParallel)
# create the cluster for caret to use
cl <- makePSOCKcluster(no_cores)
registerDoParallel(cl)


data_y = all_cal_df[,"tot_all"]
data_x = all_cal_df[, c("meanopt","sdopt","sel", "cat.freq","catastr.mort","age.sex.mat","sons.mean")]
all.rf.red <- train(x = data_x, y = data_y,method = 'ranger',
                    
                    # should be set high at least p/3
                    tuneLength = 10,
                    trControl = fitControl,
                    
                    ## parameters passed onto the ranger function
                    # the bigger the better.
                    num.trees = 1000,
                    importance = "permutation")


saveRDS(all.rf.red,"all.rf.red.RDS")
stopCluster(cl)
registerDoSEQ()


all.rf.full = readRDS("all.rf.full.RDS")
all.rf.red = readRDS("all.rf.red.RDS")
all.gam.full = readRDS("all.gam.full.RDS")
all.lm.full = readRDS("all.lm.full.RDS")
all.lm.red = readRDS("all.lm.red.RDS")

### Predictions for each model for the test data set

all_val_df$pred_rf_full = predict(all.rf.full, newdata = all_val_df)
all_val_df$pred_rf_red = predict(all.rf.red, newdata = all_val_df)
all_val_df$pred_gam_full = predict(all.gam.full, newdata = all_val_df)
all_val_df$pred_lm_full = predict(all.lm.full, newdata = all_val_df)
all_val_df$pred_lm_red = predict(all.lm.red, newdata = all_val_df)



all_mod_rsq_df = data.frame(mod = rep(NA,5), r_sq = NA, mean_err = NA) ## data frame with r_sq wrt 1:1 line and mean absolute error over the test data set

SSE.rf_full = sum((all_val_df$tot_all-all_val_df$pred_rf_full)^2, na.rm = T)

SST.base = sum((all_val_df$tot_all-mean(all_val_df$tot_all))^2, na.rm = T)

all_mod_rsq_df$mod[1] = "rf_full"
all_mod_rsq_df$r_sq[1] = 1-(SSE.rf_full/SST.base)
all_mod_rsq_df$mean_err[1] = mean(abs(all_val_df$tot_all-all_val_df$pred_rf_full))


SSE.rf_red = sum((all_val_df$tot_all-all_val_df$pred_rf_red)^2, na.rm = T)

SST.base = sum((all_val_df$tot_all-mean(all_val_df$tot_all))^2, na.rm = T)

all_mod_rsq_df$mod[2] = "rf_red"
all_mod_rsq_df$r_sq[2] = 1-(SSE.rf_red/SST.base)
all_mod_rsq_df$mean_err[2] = mean(abs(all_val_df$tot_all-all_val_df$pred_rf_red))


SSE.gam_full = sum((all_val_df$tot_all-all_val_df$pred_gam_full)^2, na.rm = T)

SST.base = sum((all_val_df$tot_all-mean(all_val_df$tot_all))^2, na.rm = T)

all_mod_rsq_df$mod[3] = "gam_full"
all_mod_rsq_df$r_sq[3] = 1-(SSE.gam_full/SST.base)
all_mod_rsq_df$mean_err[3] = mean(abs(all_val_df$tot_all-all_val_df$pred_gam_full))


SSE.lm_full = sum((all_val_df$tot_all-all_val_df$pred_lm_full)^2, na.rm = T)

SST.base = sum((all_val_df$tot_all-mean(all_val_df$tot_all))^2, na.rm = T)

all_mod_rsq_df$mod[4] = "lm_full"
all_mod_rsq_df$r_sq[4] = 1-(SSE.lm_full/SST.base)
all_mod_rsq_df$mean_err[4] = mean(abs(all_val_df$tot_all-all_val_df$pred_lm_full))



SSE.lm_red = sum((all_val_df$tot_all-all_val_df$pred_lm_red)^2, na.rm = T)

SST.base = sum((all_val_df$tot_all-mean(all_val_df$tot_all))^2, na.rm = T)

all_mod_rsq_df$mod[5] = "lm_red"
all_mod_rsq_df$r_sq[5] = 1-(SSE.lm_red/SST.base)
all_mod_rsq_df$mean_err[5] = mean(abs(all_val_df$tot_all-all_val_df$pred_lm_red))


save

### Plot the GAM smooth for n_low

plt_gam_list = plot(all.gam.full)[[1]]

gam_all_df = data.frame(x = plt_gam_list$x, y = plt_gam_list$fit,
                        se_up = plt_gam_list$fit + 2 * plt_gam_list$se, se_down = plt_gam_list$fit - 2 * plt_gam_list$se)




size.title = 15
line.lwd = 1
size.label.x = 18
size.text.x = 16
size.point = 2
size.label.y = 18
size.text.y = 16
size.legend.text = 15
size.legend.title = 20
unit.legend.h = 1.8
unit.legend.w = 1.8
size.ann = 10
colour.axis = "gray20"
colour.theme = "black"
colour.axis.line = "gray20"
colour.line = "gray50"
label.T = "Heterozygosity"
max_size_dot = 8
leg.x = 0.15
leg.y = 0.15

## Theme to be used for all plots

theme.pop =  theme(plot.title = element_text(lineheight=.8, face="bold", size = size.title,hjust = 0.5), 
                   plot.background = element_blank()
                   ,panel.grid.major = element_blank()
                   ,panel.grid.minor = element_blank()
                   ,panel.border = element_blank()
                   ,panel.background = element_blank(),
                   axis.line = element_line(color = 'black'),
                   plot.margin = unit(c(t = 1,r = 1,b = 2, l = 1), "cm"),
                   axis.title.x = element_text(size=size.label.x,vjust=-2),
                   axis.text.x  = element_text(size=size.text.x, vjust = 0.5),
                   axis.title.y = element_text(size=size.label.x, vjust = 2),
                   axis.text.y  = element_text(size=size.text.x),
                   legend.title = element_blank(),
                   legend.text = element_text(size = size.legend.text),
                   legend.spacing.y = unit(5,"cm"),
                   legend.position = c(leg.x, leg.y),
                   legend.key = element_rect(fill = "white", size = 5),
                   legend.key.size = unit(1.5,"lines")) 




plot_gam_all = ggplot(data = gam_all_df, aes(x = x, y = y)) +
  geom_line() +
  geom_line(data = gam_all_df, aes(x = x, y = se_up), lty = 2) +
  geom_line(data = gam_all_df, aes(x = x, y = se_down), lty = 2) + 
  theme.pop + 
  labs(x = "n_low (year)", y = "s(n_low, k = 3)") +
  scale_y_continuous(limits = c(-140,20), breaks =seq(-140, 20, 20))

plot_gam_all

# uncomment below to save the plot

# ggsave(plot_gam_all, device = "pdf", filename = "Plot_gam_all.pdf",width = 10, height = 10)  
  