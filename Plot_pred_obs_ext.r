### Load simulation data
if(!exists("sim_res")) sim_res = readRDS("sim_res.RDS")
### Load test data (needed to see which models has predicted what)
if(!exists("ext.stand.df.val")) ext.stand.df.val = readRDS("ext.stand.df.val.RDS")

### Load plotting functions
source("Plot_functions_for_simulations.r")



y1 = 570
y2 = 540
y3 = 510

pred_1 = plot_sampl_pred.f(cont = 42319, val_df = ext.stand.df.val, res_list = sim_res)
pred_1 = pred_1 + annotate("text", x = 230, y = y1, 
                           label = "GLM = 0",
                           color = "black", fontface = 2, size =4,hjust = 0) +
  annotate("text", x = 230, y = y2, 
           label = "GAM = 0",
           color = "black", fontface = 2, size =4,hjust = 0) +
  annotate("text", x = 230, y = y3, 
           label = "RF = 1",
           color = "black", fontface = 2, size =4,hjust = 0)


pred_2 = plot_sampl_pred.f(cont = 41881, val_df = ext.stand.df.val, res_list = sim_res)
pred_2 = pred_2 + annotate("text", x = 230, y = y1, 
                           label = "GLM = 1",
                           color = "black", fontface = 2, size =4,hjust = 0) +
  annotate("text", x = 230, y = y2, 
           label = "GAM = 0",
           color = "black", fontface = 2, size =4,hjust = 0) +
  annotate("text", x = 230, y = y3, 
           label = "RF = 0",
           color = "black", fontface = 2, size =4,hjust = 0)

pred_3 = plot_sampl_pred.f(cont = 19922, val_df = ext.stand.df.val, res_list = sim_res)
pred_3 = pred_3 + annotate("text", x = 230, y = y1, 
                           label = "GLM = 0",
                           color = "black", fontface = 2, size =4,hjust = 0) +
  annotate("text", x = 230, y = y2, 
           label = "GAM = 1",
           color = "black", fontface = 2, size =4,hjust = 0) +
  annotate("text", x = 230, y = y3, 
           label = "RF = 0",
           color = "black", fontface = 2, size =4,hjust = 0)

pred_4 = plot_sampl_pred.f(cont = 41868, val_df = ext.stand.df.val, res_list = sim_res)
pred_4 = pred_4 + annotate("text", x = 230, y = y1, 
                           label = "GLM = 0",
                           color = "black", fontface = 2, size =4,hjust = 0) +
  annotate("text", x = 230, y = y2, 
           label = "GAM = 0",
           color = "black", fontface = 2, size =4,hjust = 0) +
  annotate("text", x = 230, y = y3, 
           label = "RF = 0",
           color = "black", fontface = 2, size =4,hjust = 0)

plot_pred_obs_ext = plot_grid(pred_1,
                          pred_2,
                          pred_3,
                          pred_4,
                          labels = c("(a)", "(b)", "(c)", "(d)"),
                          nrow = 2, align = "v",hjust = -4)

plot_pred_obs_ext

# Plot in pdf, width = 16, height = 10