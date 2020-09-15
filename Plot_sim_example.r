### Plot simulation replicates. The replicates used here are the same plotted in the paper. Targ is the number of the simulation in sim_res.RDS


source("Plot_functions_for_simulations.r")

### Load simulation data
if(!exists("sim_res")) sim_res = readRDS("sim_res.RDS")


sim_3 = plot_sim.f(targ = 18425) + labs(y="")
sim_2 = plot_sim.f(targ = 1979) + labs(y="")
sim_1 = plot_sim.f(targ = 16400)


# 20124 another interesting replicate

Plot_sim_example = plot_grid(sim_1,
                         sim_2,
                         sim_3,
                         labels = c("(a)", "(b)", "(c)"),
                         nrow = 1, align = "v",hjust = -4)

Plot_sim_example

# 20 8