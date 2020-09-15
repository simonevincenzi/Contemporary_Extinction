source("load_packages.r")


plot_sim.f = function (targ = 10, res_list = sim_res) {
  
  size.title = 15
  line.lwd = 1
  size.label.x = 18
  size.text.x = 14
  size.point = 2
  size.label.y = 18
  size.text.y = 14
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
                     plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
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
  
  
  
  ext_time = max(which(sim_res[[targ]]$popsize.post>0))
  mult_cat = -5
  pheno_ts = sim_res[[targ]]$phenomean *20
  pheno_ts[which(pheno_ts == 0)] = NA
  
  data_ts = data.frame(year = 1:250, time_s = sim_res[[targ]]$popsize.post, cat = sim_res[[targ]]$catastr.vett, opt = sim_res[[targ]]$optimum *50 + 200, pheno = pheno_ts + 50)
  
  data_ts[(ext_time+1) : nrow(data_ts),c("time_s","cat","opt")] = NA
  
  data_ts$cat = ifelse(data_ts$cat == 1,data_ts$cat *  mult_cat,NA)
  
  pp = ggplot(data = data_ts, aes(x = year, y = time_s)) +
    geom_point() +
    geom_line() +
    theme.pop +
    geom_point(data = data_ts, aes(x = year, y = cat), col = "black", size = 2.5, shape = 1) +
    scale_x_continuous(limits = c(0,250)) +
    scale_y_continuous(limits = c(-15,550), breaks = seq(0,500,100)) +
    labs(x = "Year", y = "Population size") + 
    #geom_point(data = data_ts, aes(x = year, y = opt), size = 0.5) +
    geom_line(data = data_ts, aes(x = year, y = opt), lwd = 0.65, lty = 2, col = "gray") +
    geom_line(data = data_ts, aes(x = year, y = pheno), lwd = 0.8, lty = 1, col = "gray")
    
  
  return(pp)
  

}



plot_sampl_pred.f = function (cont = 10, val_df = ext.stand.df.val, res_list = sim_res) {
  
  #  cont is a consecutive number of the replicate from ext.stand.df.val
  #  val_df = ext.stand.df.val is the data frame used for validation
  #  res_list = sim_res it is the list with all simulation results
  
  
  size.title = 15
  line.lwd = 1
  size.label.x = 18
  size.text.x = 14
  size.point = 2
  size.label.y = 18
  size.text.y = 14
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
                     plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
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
  
  
  
  targ = val_df$sim_rep[val_df$cont_cons == cont]
  pheno_ts = sim_res[[targ]]$phenomean *20
  pheno_ts[which(pheno_ts == 0)] = NA
  
  ext_time = max(which(sim_res[[targ]]$popsize.post>0))
  mult_cat = -5
  
  data_ts = data.frame(year = 1:250, time_s = sim_res[[targ]]$popsize.post, cat = sim_res[[targ]]$catastr.vett, opt = sim_res[[targ]]$optimum *50 + 200, pheno = pheno_ts + 50)
  
  data_ts[(ext_time+1) : nrow(data_ts),c("time_s","cat","opt")] = NA
  
  data_ts$cat = ifelse(data_ts$cat == 1,data_ts$cat *  mult_cat,NA)
  
  data_ts$sampl_in = ifelse(data_ts$year %in% seq(from = val_df$sampl_beg[val_df$cont_cons == cont],to = val_df$sampl_end[val_df$cont_cons == cont], by = 1), "y","n")
  
  pp = ggplot(data = data_ts, aes(x = year, y = time_s)) +
    geom_point() +
    geom_line() +
    theme.pop +
    geom_point(data = data_ts, aes(x = year, y = cat), col = "black", size = 2.5, shape = 1) +
    scale_x_continuous(limits = c(0,250)) +
    scale_y_continuous(limits = c(-15,600), breaks = seq(0,500,100)) +
    labs(x = "Year", y = "Population size") + 
    geom_line(data = data_ts, aes(x = year, y = opt), lwd = 0.65, lty = 2, col = "gray") + 
    geom_line(data = data_ts, aes(x = year, y = pheno), lwd = 0.8, lty = 1, col = "gray") +
    #geom_point(data = data_ts, aes(x = year, y = opt), size = 0.5) +
    # geom_line(data = filter(data_ts, sampl_in == "y"), aes(x = year, y = time_s), lwd = 2) +
    # geom_line(data = filter(data_ts, sampl_in == "y"), aes(x = year, y = opt), lwd = 1.5, lty = 2, col = "gray20") +
    # geom_point(data = filter(data_ts, sampl_in == "y"), aes(x = year, y = cat), col = "black", size = 2.5, shape = 21, stroke = 2) +
    # geom_rect(aes(xmin = min(filter(data_ts, sampl_in == "y")$year), xmax = max(filter(data_ts, sampl_in == "y")$year), ymin = -Inf, ymax = Inf), fill = "lightgray", alpha = 0.01) + 
    geom_vline(xintercept = (max(filter(data_ts, sampl_in == "y")$year) + 10), linetype = 1,  lwd = 0.5, col = "gray50") +
    geom_vline(xintercept = (min(filter(data_ts, sampl_in == "y")$year)), linetype = 1,  lwd = 0.5, col = "gray50") +
    geom_vline(xintercept = (max(filter(data_ts, sampl_in == "y")$year)), linetype = 1,  lwd = 0.5, col = "gray50") +
    annotate("text", x = (min(filter(data_ts, sampl_in == "y")$year) + max(filter(data_ts, sampl_in == "y")$year))/2, y = 540, 
             label = "O",
             color = "black", fontface = 2, size = 4) +
    annotate("text", x = (max(filter(data_ts, sampl_in == "y")$year) + (max(filter(data_ts, sampl_in == "y")$year)+10))/2, y = 540, 
             label = "P",
             color = "black", fontface = 2, size =4)
  
  return(pp)
  
  
}



plot_all.f = function (targ = 11, list_top_bt = list_tb, surv_sim_df = pers_df, res_list = sim_res) {
  
  # ggplot(data = filter(ext.s7.df, ext_year<249 & ext == 1), aes(x = ext_year)) + geom_histogram(fill = "white", col = "red") + theme.pop
  
  library("tidyverse")
  library("Metrics")
  library("caret")
  library("MuMIn")
  library("lubridate")
  library("mgcv")
  library("parallel")
  library("brms")
  library("rlist")
  library("e1071")
  library("ranger")
  library("pdp")
  library("tidymv")
  library("cowplot")
  library("rlist")
  
  
  
  size.title = 15
  line.lwd = 0.5
  size.label.x = 18
  size.text.x = 14
  size.point = 2
  size.label.y = 18
  size.text.y = 14
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
                     plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
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
  
  
  
  
  
  # cont = surv_sim_df$cont_cons[surv_sim_df$sim_rep == targ]
  
  cont = list.which(list_top_bt, sim_rep == targ)
  top_bt_df = list_top_bt[[cont]]$t_b_df
  agg_top_bt = list_top_bt[[cont]]$agg_t_b
  
  col_top = "gray10"
  col_bot = "gray70"
  cols <- c("t" = col_top, "b" = col_bot)
  lwds <- c("t" = 1, "b" = 2)
  
  pp = ggplot(data = top_bt_df, aes(x = Time, y = Freq, group = as.factor(Allelic_value), col = as.factor(t_b), lty = as.factor(t_b))) +
    geom_point(lwd = line.lwd, show.legend = FALSE) +
    geom_line(lwd = line.lwd,show.legend = FALSE) +
    scale_colour_manual(values=cols) +
    scale_linetype_manual(values=lwds) +
    theme.pop +
    geom_point(data = filter(agg_top_bt, t_b == "t"), aes(x = Time, y = freq_mean, group = t_b), col = col_top, size = 4) +
    geom_line(data = filter(agg_top_bt, t_b == "t"), aes(x = Time, y = freq_mean, group = t_b), col = col_top, lwd = 3) +
    geom_point(data = filter(agg_top_bt, t_b == "b"), aes(x = Time, y = freq_mean, group = t_b), col = col_bot, size = 4) +
    geom_line(data = filter(agg_top_bt, t_b == "b"), aes(x = Time, y = freq_mean, group = t_b), col = col_bot, lwd = 3) +
    scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1)) +
    theme(legend.position = "none") +
    labs(x = "Year", y = "Allelic frequency")
  
  
  
  return(pp)
  
  
}
