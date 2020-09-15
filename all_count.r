### Load simulation data
if (!exists("sim_res")) sim_res = readRDS("sim_res.RDS")

### Summary statistics in a data frame for all populations
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


### Compute additional summary statistics only for surviving replicates and assigns them to replicates in a data frame (names should be self-explanatory) -- it takes some time to run the routine below, I saved the objects in pers_df.RDS and list_tb.RDS


min_pop_1 = 150
min_pop_2 = 100

set.seed(50)
list_tb = list() # list with the frequencies of alleles in top and bottom 10%
pers_df = filter(ext.s7.df, ext == 0)   # subset of replicates that survived

pers_df$pheno_mean = NA
pers_df$pheno_sd = NA
pers_df$hetero_mean = NA
pers_df$hetero_sd = NA
pers_df$pop_mean = NA
pers_df$pop_sd = NA
pers_df$pop_cv = NA
pers_df$test_p = NA
pers_df$diff_all = NA
pers_df$l_seq_low_1 = 0 # maximum continuous number of years under min_pop_1
pers_df$l_seq_low_2 = 0 # maximum continuous number of years under min_pop_2
pers_df$tot_pop_un_1 = 0 # total number of years under min_pop_1
pers_df$tot_pop_un_2 = 0 # total number of years under min_pop_2
pers_df$tot_all = NA
pers_df$cont_cons = 1:nrow(pers_df)

for (i in 1:nrow(pers_df)) {

  if (i%%100 == 0) {print(i)}
  
  r_sim = pers_df$sim_rep[i]
  
  if(last(sim_res[[r_sim]]$popsize.post)>10) { 
  
    pers_df$pheno_mean[i] = mean(sim_res[[r_sim]]$phenomean[240:250])
    pers_df$pheno_sd[i] = sd(sim_res[[r_sim]]$phenomean[240:250])
    pers_df$hetero_mean[i] = mean(sim_res[[r_sim]]$heter.mean.year[240:250])
    pers_df$hetero_sd[i] = sd(sim_res[[r_sim]]$heter.mean.year[240:250])
    
    
  bottom_10_df = filter(sim_res[[r_sim]]$all.freq, sim_res[[r_sim]]$all.freq$Allelic_value < as.numeric(quantile(sim_res[[r_sim]]$all.freq$Allelic_value,probs = c(0.1,0.9))[1]))
  
  top_10_df = filter(sim_res[[r_sim]]$all.freq, sim_res[[r_sim]]$all.freq$Allelic_value > as.numeric(quantile(sim_res[[r_sim]]$all.freq$Allelic_value,probs = c(0.1,0.9))[2]))

  
  if (sum(c(bottom_10_df$Freq_time_250,top_10_df$Freq_time_250)) > 0) { ###  check if at year 250 there are remaining alleles in the top and bottom 10%
  
  
  if(min(sim_res[[r_sim]]$popsize.post) < min_pop_1) {
  
    tt = ifelse(as.numeric(sim_res[[r_sim]]$popsize.post) < min_pop_1, 1, 0)
  
  pers_df$tot_pop_un_1[i] = sum(tt)
    
  pers_df$l_seq_low_1[i] = data.frame(l_seq = rle(tt)$length, val = rle(tt)$values) %>% filter(., val == 1) %>% arrange(., desc(l_seq)) %>% top_n(., wt = l_seq, n = 1) %>%  head(.,1) %>% .$l_seq   
  
  }
    
    if(min(sim_res[[r_sim]]$popsize.post) < min_pop_2) {
      
      tt = ifelse(as.numeric(sim_res[[r_sim]]$popsize.post) < min_pop_2, 1, 0)
      
      pers_df$tot_pop_un_2[i] = sum(tt)
      
      pers_df$l_seq_low_2[i] = data.frame(l_seq = rle(tt)$length, val = rle(tt)$values) %>% filter(., val == 1) %>% arrange(., desc(l_seq)) %>% top_n(., wt = l_seq, n = 1) %>%  head(.,1) %>% .$l_seq   
      
    }
    
    
  ### mean, sd, and cv of population size over simulation time
    
  pers_df$pop_mean[i] = mean(sim_res[[pers_df$sim_rep[i]]]$popsize.post)   
  pers_df$pop_sd[i] = sd(sim_res[[pers_df$sim_rep[i]]]$popsize.post) 
  pers_df$pop_cv[i] = pers_df$pop_sd[i]/pers_df$pop_mean[i]
    
  pers_df$tot_all[i] = sum(sim_res[[r_sim]]$all.freq$Freq_time_250 > 0)
    
diff_test = t.test(bottom_10_df$Freq_time_250,top_10_df$Freq_time_250)
  
pers_df$test_p[i] = ifelse(diff_test$p.value > 0.05, 0, 1)
  

list_prov = list()
list_prov[[1]] = as.list(data.frame(sim_rep = NA, t_b_df = NA, agg_t_b = NA))
list_prov[[1]]$sim_rep = pers_df$sim_rep[i]
bt_10 = rbind(top_10_df,bottom_10_df)
bt_10$t_b = c(rep("t", nrow(top_10_df)),rep("b",nrow(bottom_10_df)))

bt_10 = bt_10 %>% select(., -Gene, -t_b) %>% 
  pivot_longer(!Allelic_value, names_to = "Time", values_to = "Freq") %>% 
  left_join(., select(bt_10, Gene, t_b, Allelic_value)) %>% filter(., Time!="Freq_diff") %>%
  mutate(., sim_rep = pers_df$sim_rep[i], cont_cons = pers_df$cont_cons[i])

bt_10$Time = rep(c(1,50,100,150,200,250),(nrow(bt_10)/6))

agg_t_b = bt_10 %>% group_by(t_b, Time) %>% summarise(n = n(),
                                            freq_mean = mean(Freq),
                                            freq_sd = sd(Freq),
                                            all_mean = mean(Allelic_value)) %>%
  mutate(., sim_rep = pers_df$sim_rep[i], cont_cons = pers_df$cont_cons[i])
                                            

list_prov[[1]]$t_b_df = bt_10
list_prov[[1]]$agg_t_b = agg_t_b

pers_df$diff_all[i] = filter(list_prov[[1]]$agg_t_b, Time == 250, t_b == "t")$freq_mean - filter(list_prov[[1]]$agg_t_b, Time == 250,t_b == "b")$freq_mean 

list_tb = c(list_tb,list_prov)
  
  }  
  
  
  }
  
  
}

# saveRDS(pers_df, "pers_df.RDS")
# saveRDS(list_tb, "list_tb.RDS")

pers_df = readRDS("pers_df.RDS") ## data.frame with the main statistics for the surviving populations (diff_all is the difference in frequency between the top and bottom 10% of alleles, tot_pop_un_1 is the total number of years below 150 individuals, tot_pop_un_2 under 100 individuals. l_seq_low_1 is the maximum continuous number of years under 150 ind and l_seq_low_2 is the maximum continuous number of years under 100 ind)
list_tb = readRDS("list_tb.RDS") # list with number of replicate (correspondence in sim_rep), t_b_df has the allelic value of gene in the top (t_b = t) or bottom (t_b = b) 10% of alleles, and Freq the allelic frequencies of single alleles every 25 years during simulation. agg_t_b has the aggregate values (mean and sd of frequency) for the top and bottom 10% of alleles.

### Correlation between number of alleles at the end of simulation time and mean phenotype in the population

all_cor_test = cor.test(pers_df$diff_all, pers_df$pheno_mean)

#plot(pers_df$tot_all, pers_df$pheno_mean)




