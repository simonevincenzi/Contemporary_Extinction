library(pbapply)
library(plyr)
library(dplyr)
library(rlist)

end_ = 1:6
#paste("phenoBIG",i,".rds",sep="")
for(i in end_){
  if(i == 1) {
    ris.tot = readRDS(paste("phenoBIG",i,".rds",sep=""))}
  else{
    nex = readRDS(paste("phenoBIG",i,".rds",sep=""))
    ris.tot = c(ris.tot,nex)
  }

  }
rm(nex)
#nex1 = readRDS("phenoBIG_bis1.rds")
#nex2 = readRDS("phenoBIG_bis2.rds")
#ris.tot = c(ris.tot,nex1,nex2)

#ris.tot = prev

##### This is for surv = 0.7

end_ = 1:6
#paste("phenoBIG",i,".rds",sep="")
for(i in end_){
  if(i == 1) {
    ris.tot.s7 = readRDS(paste("phenoBIG_S07",i,".rds",sep=""))}
  else{
    nex = readRDS(paste("phenoBIG_S07",i,".rds",sep=""))
    ris.tot.s7 = c(ris.tot.s7,nex)
  }
}
rm(nex)
#ris.tot.s7 = prev

ris.tot.s7.filter = list.filter(ris.tot.s7)

ris.tot.s7.filter = sim_res

ext.s7.df = data.frame(
  sel = sapply(ris.tot.s7.filter, with, selstrength), 
  ext = sapply(ris.tot.s7.filter, with, extinct), 
  sdopt =  sapply(ris.tot.s7.filter, with, sdopt), 
  meanopt = sapply(ris.tot.s7.filter, with, meanopt),
  mean.pheno = sapply(ris.tot.s7.filter, with, mean(phenomean[(250-10):250])),
  ext_year = sapply(ris.tot.s7.filter, with, yearextinct),
  age.sex.mat = sapply(ris.tot.s7.filter, with, age.sex.mat),
  avg.surv = sapply(ris.tot.s7.filter, with,avg.surv),
  sons.mean = sapply(ris.tot.s7.filter, with,sons.mean),
  cor.all.est = sapply(ris.tot.s7.filter, with, cor.all.est),
  cor.all.p = sapply(ris.tot.s7.filter, with, cor.all.p), 
  cat.freq = sapply(ris.tot.s7.filter, with, cat.freq),
  catastr.mort = sapply(ris.tot.s7.filter, with,catastr.mort)
)
ext.s7.df$sim_rep = 1:nrow(ext.s7.df)


ext.s7.df$year_bc = NA

for (i in 1:nrow(ext.s7.df)) {
  
  if (ext.s7.df$ext[i] == 1) {
    
    cat_happ = which(ris.tot.s7.filter[[i]]$catastr.vett == 1) 
    cat_happ_bef = cat_happ[which(cat_happ < ext.s7.df$ext_year[i])]
    if (length(cat_happ) > 0) {
    cat_close = min(ext.s7.df$ext_year[i] - cat_happ_bef)}
    if (length(cat_happ) == 0) {
    cat_close = NA}
    
    ext.s7.df$year_bc[i] = cat_close 
    
    }
  
  
}


#####


ris.tot.filter = list.filter(ris.tot, sons.mean == 2)

ext.s1.df = data.frame(
  sel = sapply(ris.tot.filter, with, selstrength), 
  ext = sapply(ris.tot.filter, with, extinct), 
  sdopt =  sapply(ris.tot.filter, with, sdopt), 
  meanopt = sapply(ris.tot.filter, with, meanopt),
  mean.pheno = sapply(ris.tot.filter, with, mean(phenomean[(250-10):250])),
  ext_year = sapply(ris.tot.filter, with, yearextinct),
  age.sex.mat = sapply(ris.tot.filter, with, age.sex.mat),
  avg.surv = sapply(ris.tot.filter, with,avg.surv),
  sons.mean = sapply(ris.tot.filter, with,sons.mean),
  cor.all.est = sapply(ris.tot.filter, with, cor.all.est),
  cor.all.p = sapply(ris.tot.filter, with, cor.all.p), 
  cat.freq = sapply(ris.tot.filter, with, cat.freq),
  catastr.mort = sapply(ris.tot.filter, with,catastr.mort)
)
ext.s1.df$sim_rep = 1:nrow(ext.s1.df)

years_wind = sample(125:229,length(ris.tot),replace = T)

sim_pos = which(ext.s1.df$ext_year>125)


sim_pos = which(ext.s7.df$ext_year>125)

sim_rep = sample(sim_pos, length(ris.tot), replace = T)


sim_years = as.data.frame(cbind(sim_rep, years_wind))
sim_years$ext_year = ext.s1.df$ext_year[sim_years$sim_rep]
sim_years$ext = ext.s1.df$ext[sim_years$sim_rep]
sim_years$with_ = ifelse((sim_years$years_wind+10-sim_years$ext_year)>0,"yes","no" )
sim_years$foll_ext = ifelse((sim_years$years_wind+20-sim_years$ext_year)>0,"yes","no" )
sim_years = sim_years[which(sim_years$with_ == "no"),]
set.seed(10)
sim_years_red = sim_years[sample(1:nrow(sim_years),nrow(sim_years),replace = F),]

### extract information from list of results
### 10 rows per each replicate

ris_df = data.frame(sim_rep = rep(sim_years_red$sim_rep, each = 10), years_init_wind = rep(sim_years_red$years_wind, each = 10), 
                    years_wind= 0, ext_year = rep(sim_years_red$ext_year,each = 10),
                    ext = rep(sim_years_red$ext, each = 10), heter_year = 0, pop_size = 0,pop_size_min = 0, pop_size_wind_end = 0,
                    cat_freq = 0, cat_mort = 0, catastr_vett = 0, 
                    age_sex = 0, sons_mean = 0, sd_opt = 0, sel = 0, meanopt = 0,foll_ext = 0) 

#for (i in 1:1000){
for (i in 1:(nrow(sim_years_red))){
  if(i == 1) {
    range_df = 1:10} else {range_df = ((i-1)*10+1):((i-1)*10+10)}
  print(i)
  range_years = ris_df$years_init_wind[min(range_df)]:(ris_df$years_init_wind[min(range_df)]+9)
  
  ris_df$heter_year[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$heter.mean.year[range_years]
  ris_df$pop_size[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$popsize.post[range_years]
  ris_df$pop_size_min[range_df] = min(ris.tot[[sim_years_red$sim_rep[i]]]$popsize.post[range_years])
  ris_df$pop_size_wind_end[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$popsize.post[range_years[length(range_years)]]
  ris_df$catastr_vett[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$catastr.vett[range_years]
  ris_df$years_wind[range_df] = ris_df$years_init_wind[min(range_df)]:(ris_df$years_init_wind[min(range_df)]+9)
  ris_df$cat_freq[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$cat.freq
  ris_df$age_sex[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$age.sex.mat
  ris_df$sons_mean[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$sons.mean
  ris_df$sd_opt[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$sdopt
  ris_df$sel[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$selstrength
  ris_df$meanopt[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$meanopt
  ris_df$cat_mort[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$catastr.mort
  ris_df$foll_ext[range_df] = sim_years_red$foll_ext[i]
}
saveRDS(ris_df,"ris_df.RDS")








cont_l = 1:nrow(sim_years_red)
#cont_l = 1:10000
tt = pblapply(cont_l, function (x) {
 
  ris_df = data.frame(sim_rep = rep(sim_years_red$sim_rep[x], 10), 
                      years_init_wind = rep(sim_years_red$years_wind[x], 10), 
                      years_wind= 0, ext_year = rep(sim_years_red$ext_year[x], 10),
                      ext = rep(sim_years_red$ext[x],10), heter_year = 0, pop_size = 0,pop_size_min = 0,pop_size_wind_end = 0,
                      cat_freq = 0, cat_mort = 0, catastr_vett = 0, 
                      age_sex = 0, sons_mean = 0, sd_opt = 0, sel = 0, meanopt = 0,foll_ext = 0,sum_cat_bef = 0, sum_cat_wind = 0, sum_cat_ext = 0) 
  

  range_years = ris_df$years_init_wind[1]:(ris_df$years_init_wind[1]+9)

  
  ris_df$heter_year = ris.tot[[sim_years_red$sim_rep[x]]]$heter.mean.year[range_years]
  ris_df$pop_size = ris.tot[[sim_years_red$sim_rep[x]]]$popsize.post[range_years]
  ris_df$pop_size_min = min(ris.tot[[sim_years_red$sim_rep[x]]]$popsize.post[range_years])
  ris_df$pop_size_wind_end = ris.tot[[sim_years_red$sim_rep[x]]]$popsize.post[range_years[length(range_years)]]
  ris_df$catastr_vett = ris.tot[[sim_years_red$sim_rep[x]]]$catastr.vett[range_years]
  ris_df$years_wind = range_years
  ris_df$cat_freq = ris.tot[[sim_years_red$sim_rep[x]]]$cat.freq
  ris_df$age_sex = ris.tot[[sim_years_red$sim_rep[x]]]$age.sex.mat
  ris_df$sons_mean = ris.tot[[sim_years_red$sim_rep[x]]]$sons.mean
  ris_df$sd_opt = ris.tot[[sim_years_red$sim_rep[x]]]$sdopt
  ris_df$sel = ris.tot[[sim_years_red$sim_rep[x]]]$selstrength
  ris_df$meanopt = ris.tot[[sim_years_red$sim_rep[x]]]$meanopt
  ris_df$cat_mort = ris.tot[[sim_years_red$sim_rep[x]]]$catastr.mort
  ris_df$foll_ext = sim_years_red$foll_ext[x] 
  ris_df$sum_cat_bef = sum(ris.tot[[sim_years_red$sim_rep[x]]]$catastr.vett[100:ris_df$years_init_wind[1]])
  ris_df$sum_cat_wind = sum(ris.tot[[sim_years_red$sim_rep[x]]]$catastr.vett[range_years])
  ris_df$sum_cat_ext = sum(ris.tot[[sim_years_red$sim_rep[x]]]$catastr.vett[(range_years[1]+1):(range_years[1]+10)])
  return(ris_df)
})
tt = rbind.fill(tt)













#mc.cores=18











for (i in 1:(nrow(sim_years_red))){
  if(i == 1) {
    range_df = 1:10} else {range_df = ((i-1)*10+1):((i-1)*10+10)}
  print(i)
  range_years = ris_df$years_init_wind[min(range_df)]:(ris_df$years_init_wind[min(range_df)]+9)
  
  ris_df$heter_year[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$heter.mean.year[range_years]
  ris_df$pop_size[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$popsize.post[range_years]
  ris_df$catastr_vett[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$catastr.vett[range_years]
  ris_df$years_wind[range_df] = ris_df$years_init_wind[min(range_df)]:(ris_df$years_init_wind[min(range_df)]+9)
  ris_df$cat_freq[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$cat.freq
  ris_df$age_sex[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$age.sex.mat
  ris_df$sons_mean[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$sons.mean
  ris_df$sd_opt[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$sdopt
  ris_df$sel[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$selstrength
  ris_df$meanopt[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$meanopt
  ris_df$cat_mort[range_df] = ris.tot[[sim_years_red$sim_rep[i]]]$catastr.mort
  ris_df$foll_ext[range_df] = sim_years_red$foll_ext[i]
}








tt$repl = rep(seq(1:(nrow(tt)/10)), each = 10)
tt$foll_ext = ifelse(tt$foll_ext == "yes",1,0)

pred_df = tt %>% 
  group_by(repl) %>% 
  summarise_each(funs(mean))

test = with(prova, glm(foll_ext ~ heter_year + pop_size + cat_freq +  age_sex * sons_mean + sd_opt + meanopt + cat_mort, family = 'binomial'))

with(prova, cor.test(heter_year, pop_size))
with(prova, cor.test(heter_year, sum_cat_bef))


####


test = ext.s1.df %>%
  group_by(age.sex.mat,sdopt,cat.freq) %>%
  summarize(mean_ext = mean(ext), se_ext = sqrt((mean(ext)*(1-mean(ext)))/n())) %>%
  ggplot(aes(x = age.sex.mat, y = mean_ext, group = cat.freq)) +
  geom_line()

with(ext.s1.df, table(ext,age.sex.mat)/colSums(table(ext,age.sex.mat)))


test = filter(ext.s1.df,ext == 0) %>%
  group_by(sdopt,cat.freq) %>%
  summarize(mean_cor = mean(cor.all.est), sd_cor = (sd(cor.all.est))) %>%
  ggplot(aes(x = age.sex.mat, y = mean_ext, group = cat.freq)) +
  geom_line()




test_lm = with(filter(test,meanopt!=0), lm(mean_ext ~ age.sex.mat*cat.freq + age.sex.mat*sdopt))