all_val_df = readRDS("all_val_df.RDS")
all_cal_df = readRDS("all_cal_df.RDS")

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


saveRDS(all_mod_rsq_df,"all_mod_rsq_df.RDS")

### Plot the GAM smooth for n_low


plot.data <- {
  dev.new()
  plt_gam_list = plot(all.gam.full)[[1]]
  dev.off()
  plt_gam_list
}


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




plot_nonlinear_n_low = ggplot(data = gam_all_df, aes(x = x, y = y)) +
  geom_line() +
  geom_line(data = gam_all_df, aes(x = x, y = se_up), lty = 2) +
  geom_line(data = gam_all_df, aes(x = x, y = se_down), lty = 2) + 
  theme.pop + 
  labs(x = "n_low (year)", y = "s(n_low, k = 3)") +
  scale_y_continuous(limits = c(-140,20), breaks =seq(-140, 20, 20))

plot_nonlinear_n_low

# uncomment below to save the plot

# ggsave(plot_gam_all, device = "pdf", filename = "Plot_gam_all.pdf",width = 10, height = 10)  
