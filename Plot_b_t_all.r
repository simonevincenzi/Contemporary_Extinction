source("Plot_functions_for_simulations.r")

if (!exists("list_bt")) list_bt = readRDS("list_tb.RDS")
if (!exists("pers_df")) list_bt = readRDS("pers_df.RDS")
if (!exists("sim_res")) list_bt = readRDS("sim_res.RDS")

all_p1 = plot_all.f(targ = 5466)
sim_p1 = plot_sim.f(targ = 5466)

all_p2 = plot_all.f(targ = 7952)
sim_p2 = plot_sim.f(targ = 7952)

plot_b_t_all = plot_grid(
                         sim_p1,
                         sim_p2,
                         all_p1,
                         all_p2,
                         labels = c("(a)", "(b)", "(c)", "(d)"),
                         ncol = 2, align = "v",hjust = -4)

plot_b_t_all

# 16 9