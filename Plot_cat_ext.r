### Load libraries

source("load_packages.r")

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


set.seed(50)

ext_df = filter(ext.s7.df, ext == 1) # only extinct populations


ext_df$cat_5 = NA   # yes/no point extreme in the last 5 years before extinction
ext_df$cat_10 = NA # yes/no point extreme in the last 5 years before extinction
ext_df$cat_5_tot = NA  # total number of point extremes in the last 5 years before extinction
ext_df$cat_10_tot = NA # total number of point extremes in the last 10 years before extinction
ext_df$max_pheno_dist_win = NA # maximum distance mean phenotype - optimal phenotype in the 5 years before extinction
ext_df$mean_pheno_dist_win = NA # mean distance mean phenotype - optimal phenotype in the 5 years before extinction

### Assign to simulations

for (i in 1:nrow(ext_df)) {
  
  if (i%%100 == 0){ print(i)}
  targ = ext_df$sim_rep[i]
  ext_df$cat_5[i] = ifelse(sum(which(sim_res[[targ]]$catastr.vett == 1) %in% ((ext_df$ext_year[i]-5):(ext_df$ext_year[i]-1)) > 0),"y","n")
  ext_df$cat_10[i] = ifelse(sum(which(sim_res[[targ]]$catastr.vett == 1) %in% ((ext_df$ext_year[i]-10):(ext_df$ext_year[i]-1)) > 0),"y","n")
  
  ext_df$cat_5_tot[i] = sum(which(sim_res[[targ]]$catastr.vett == 1) %in% ((ext_df$ext_year[i]-5):(ext_df$ext_year[i]-1)))
  ext_df$cat_10_tot[i] = sum(which(sim_res[[targ]]$catastr.vett == 1) %in% ((ext_df$ext_year[i]-10):(ext_df$ext_year[i]-1)))
  
  ext_df$max_pheno_dist_win[i] = max(abs(sim_res[[targ]]$phenomean[((ext_df$ext_year[i]-5):(ext_df$ext_year[i]-1))] - sim_res[[targ]]$optimum[((ext_df$ext_year[i]-5):(ext_df$ext_year[i]-1))]))
  ext_df$mean_pheno_dist_win[i] = mean(abs(sim_res[[targ]]$phenomean[((ext_df$ext_year[i]-5):(ext_df$ext_year[i]-1))] - sim_res[[targ]]$optimum[((ext_df$ext_year[i]-5):(ext_df$ext_year[i]-1))]))
}


### Ratios (proportion of replicates going extinct with ane extreme in the 5 or 10 years before extinction)
### 
ext_10 = as.data.frame(table(ext_df$cat_10_tot)/sum(table(ext_df$cat_10_tot))) %>% rename(., n_cat = Var1)
ext_5 = as.data.frame(table(ext_df$cat_5_tot)/sum(table(ext_df$cat_5_tot))) %>% rename(., n_cat = Var1)

### Join the 5 and 10 years data set together
### 
ext_10_5_df = rbind(ext_10, ext_5) %>% mutate(., y_bef = c(rep(10,nrow(ext_10)),rep(5,nrow(ext_5))) )


size.title = 15
line.lwd = 1
size.label.x = 18
size.text.x = 14
size.point = 4
size.label.y = 18
size.text.y = 14
size.legend.text = 15
size.legend.title = 15
unit.legend.h = 1.8
unit.legend.w = 1.8
size.ann = 10
colour.axis = "gray20"
colour.theme = "black"
colour.axis.line = "gray20"
colour.line = "gray50"
label.T = "Heterozygosity"
max_size_dot = 8
leg.x = 0.85
leg.y = 0.85

## Theme to be used for all plots

theme.pop =  theme(plot.title = element_text(lineheight=.8, face="bold", size = size.title,hjust = 0.5), 
                   plot.background = element_blank()
                   ,panel.grid.major = element_blank()
                   ,panel.grid.minor = element_blank()
                   ,panel.border = element_blank()
                   ,panel.background = element_blank(),
                   axis.line = element_line(color = 'black'),
                   plot.margin = unit(c(1.2,1.2,1.2,1.2), "cm"),
                   axis.title.x = element_text(size=size.label.x,vjust=-3),
                   axis.text.x  = element_text(size=size.text.x, vjust = 0.5),
                   axis.title.y = element_text(size=size.label.x, vjust = 3),
                   axis.text.y  = element_text(size=size.text.x),
                   #legend.title = element_blank(),
                   legend.text = element_text(size = size.legend.text),
                   legend.spacing.y = unit(0.25,"cm"),
                   legend.position = c(leg.x, leg.y),
                   legend.key = element_rect(fill = "white", size = 5),
                   legend.key.size = unit(1.5,"lines")
) 



### Plot with proportion of replicates going extinct given number of extreme point event occurring 5 or 10 years before the extinction window

plot_cat_ext = ggplot(data = ext_10_5_df, aes(x = n_cat, y = Freq, shape = as_factor(y_bef))) +
  geom_point(size = size.point) +
  theme.pop +
  labs(x = "Number of extreme events", y = "Proportion of replicates") +
  theme(legend.title = element_text(size = size.legend.title)) + 
  scale_shape_manual(name = "Years before extinction", labels = c("5", "10"),values=c(1,8)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0,0.6,0.1))

plot_cat_ext

# to save as pdf, width = 12, height = 8