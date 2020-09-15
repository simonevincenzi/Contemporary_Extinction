### Load libraries

source('load_packages.r')


### Source main script and data

source("phenoplast.r")
source("data_for_par.r")



### Run simulations

logFile = "log_file.txt"
cat("This is a log file for simulation models", file=logFile, append=FALSE, sep = "\n")

sim_row = 1:200

dataforpar = dataforpar[sim_row,]  # this is as an example, dataforpar has 34560 sets of parameters, if you want to run simulations for all sets of parameters (34560), comment this out

max <- 50 ## number of simulations before re-building the cluster. Simulations are joined in a unique list 
x <- 1:nrow(dataforpar)
d1 <- split(1:nrow(dataforpar), ceiling(x/max))

### the routine below runs the simulations in parallel (they still take hours for thousands of them! The simulations in the paper have been saved as sim_res.RDS)

sim_res_list = list()
for (i in 1:length(d1)) {
  sim_res_list_chunk = mclapply(d1[[i]],function (x) do.call(Pheno.Plast,as.list(dataforpar[x,])),
                                mc.cores = detectCores(), mc.preschedule = F)
  sim_res_list = c(sim_res_list,sim_res_list_chunk)
  rm(sim_res_list_chunk)
  # saveRDS(ll_list_temp,"sim_res.RDS") # uncomment if you want to over-write the 34560 simulation replicates already saved 
  cat(cat(as.character(Sys.time()), file=logFile, append=TRUE, sep = "\n"), file=logFile, append=TRUE, sep = "\n")
}


# sim_res_list is the list with all the results (each simulation can be accessed via sim_res_list[[number]]). Using the `rlist` package, filtering can be easily applied to lists. For example, if you want to find which simulations (number as above) have gone extinct -- sim_ext = list.which(sim_res_list, extinct == 1). 
# data.frame("items" = names(sim_res_list[[1]])) to get the names of items that are saved for each simulations. For example, to access the frequency of point extremes for simulation number 1, sim_res_list[[1]]$cat.freq
# 
# 
# to run a single simulation, you can use sim_list_single = do.call(Pheno.Plast,as.list(dataforpar[x,])), where x is the number of the simulation you want to run