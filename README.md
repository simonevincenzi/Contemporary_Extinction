<strong>September 2020</strong>

# Data and code for the manuscript "Contemporary risk of extinction in an extreme environment"


<strong>Here is the abstract of the paper, which gives context to the modeling done.</strong>

> The increased frequency and intensity of extreme events are recognized among the most worrisome aspects of climate change. However, despite increased attention from scientists and conservationists, developing and testing general theories and hypotheses on the effects of extreme events on natural populations remains intrinsically challenging.
Using numerical simulations, I tested some of the hypotheses on risk of extinction and population and genetic dynamics in an environment in which both climate (e.g., temperature, rainfall) and point (e.g., fires, floods) extremes occur. A quantitative trait is selected for by a climate variable, but point extremes cause trait-independent massive mortalities.
I found additive effects between age at first reproduction and fecundity on risk of extinction for the range of values I simulated. The extent of population bottlenecks (operationally, the number of years in which a population was at low numbers) was a good predictor of allelic richness for the quantitative trait selected for by the climate. A simple model including basic demographic and vital rates information, along with climate/environmental measures, provided excellent predictions of the contemporary risk of population extinction. Mean and minimum population size measured in a 10-year “observation window” were largely the most important predictors of risk of population extinction in the following 10-year “extinction window”.



## 1. Simulations

I ran the scripts and models with R version 4.0.0 (2020-04-24), Platform: x86_64-apple-darwin17.0 (64-bit), Running under: macOS Catalina 10.15.6. The parallel routine uses mclapply, which works only on Mac and Linux machines. I recommend to use the terminal for running the model fitting routines in parallel (RStudio is more likely to crash, according to my experience).

The script `load_packages.r` installs and loads the libraries needed. The script `phenoplast.r` (quite well commented) contains the function `Pheno.Plast()`, which is the main function for running the simulations. The script `data_for_par.r` creates the data frame with the input for the simulations as run in the manuscript. The command `source("run_model.r")` runs the first 200 simulations, which are saved in the list `sim_res_list`. To run all the simulations (>30,000), follow the instructions in the script (i.e., uncomment a line). All simulation results are saved in the RDS object (a list) `sim_res.RDS`.


## 2. Analysis of simulation results    

First, download some pre-saved data and model files from [https://figshare.com/articles/dataset/Contemporary_Extinction/12959948](https://figshare.com/articles/dataset/Contemporary_Extinction/12959948) and move them in the same folder as the other scripts.

Then, the command `source("Plot_sim_example.r")` produces a panel of three simulation replicates (details are in the code). The plot object within R is saved as `plot_sim_example`.

![Plot_sim_example](https://github.com/simonevincenzi/Contemporary_Extinction/blob/master/Plots/Plot_sim_example.png)

The command `source("Plot_cat_ext.r")` (1) creates the data.frame (`ext_10_5_df`) and (2) plots simulation results on the number or point extreme events 5 or 10 years before extinction and on the proportion of replicates that went extinct (`plot_cat_ext`).

![Plot_cat_ext](https://github.com/simonevincenzi/Contemporary_Extinction/blob/master/Plots/Plot_repro_cat.png)

The script `Plot_repro_cat.r` creates the plot (saved within R as `plot_repro_cat`) and data frames (`s_mat_all_df` and `cat_fm_df`) with the summary information on how the frequency of extinction in simulations that are run with the same parameter values varies with age at first reproduction, expected number of offspring, and intensity and severity of points extremes.  

![Plot_repro_cat](https://github.com/simonevincenzi/Contemporary_Extinction/blob/master/Plots/Plot_repro_cat.png)

The command `source("all_count.r")` creates the data.frame (`pers_df` within R) and list (`list_tb`) with data and statistics for populaiton that persisted, including frequencies of alleles in the top and bottom 10% of allelic values (the more detailed description of the objects created by this script are in the file). The correlation object `all_cor_test` is also created. The plot shown here below, which can be obtained with the command `source("Plot_b_t_all.r")` is an example of alleles frequencies changing over time.  

![Plot_b_t_all](https://github.com/simonevincenzi/Contemporary_Extinction/blob/master/Plots/Plot_b_t_all.png)


The command `source("all_mod_fit.r")` fits the models (GAMs, GLMs, and RFs) for predicting the number of alleles at the end of simulation time (models are already saved— if you run the script, they will be overwritten, with names `all.gam.full`, `all.rf.full` etc.). The command `source("all.pred.r")` makes predictions for each model and saves the summary results (r squared, MAE) in the object `all_mod_rsq_df` (also saved as `all_mod_rsq_df.RDS`). The plot for the smooth function for n_low is also created (`plot_nonlinear_n_low`).  

![Plot_nonlinear_n_low](https://github.com/simonevincenzi/Contemporary_Extinction/blob/master/Plots/Plot_nonlinear_n_low.png)

The command `source("cont_ext_mod_fit.r")` fits the logistic regressions (GLMs and GAMs) and classification Random Forest for prediction of extinction (or persistence) in the 10-year "extinction window". It takes some time to fit the models — I already saved the data sets for training (`ext.stand.df.cal`) and testing (`ext.stand.df.val`), along with the fitted models (`ext.lrm.full.RDS`, `ext.lrm.base.RDS`, `ext.gam.full.RDS`, `ext.gam.base.RDS`, `ext.rf.full.RDS`, `ext.rf.base.RDS`).  

The command `source("test_pred_full.r")` makes predictions for the full models; summary results are saved in `acc_table_full_df` (predictions are also added to the data frame `ext.stand.df.val`). The command `source("cont_ext_mod_fit.r")` makes predictions for the base models, summary results are saved in `acc_table_base_df` (predictions are also added to the data frame `ext.stand.df.val`).  

The command `source("Plot_pred_obs_ext.r")` plots 4 examples of matching/mis-matching predictions of the models (saved in object `plot_pred_obs_ext`).



![Plot_pred_obs_ext](https://github.com/simonevincenzi/Contemporary_Extinction/blob/master/Plots/Plot_pred_obs_ext.png)


The command `source("Plot_imp")` produces the importance plot for the RF (full and base) models for prediction of contemporary extinction.


![Plot_imp](https://github.com/simonevincenzi/Contemporary_Extinction/blob/master/Plots/Plot_imp.png)


## 4. Manuscript

A pre-print can be found at https://www.biorxiv.org/content/10.1101/2020.09.15.298919v1
