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



size.title = 15
line.lwd = 1
size.label.x = 18
size.text.x = 14
size.point = 5
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
leg.x = 0.15
leg.y = 0.9

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

### data frame with proportion of replicates going extinct grouped by fecundity, age at first reproduction, frequency and intensity of extreme events

s_mat_df = ext.s7.df %>%  group_by(sons.mean, age.sex.mat) %>%  
  summarise(n = n(),
            prop.ext = sum(ext == 1)/n) %>% arrange(., prop.ext)

s_mat_cat_df = ext.s7.df %>% filter(., cat.freq == max(cat.freq), catastr.mort == max(catastr.mort),
                                                       meanopt == max(meanopt), sdopt == max(sdopt), sel == max(sel)) %>% group_by(sons.mean, age.sex.mat) %>%  
  summarise(n = n(),
            prop.ext = sum(ext == 1)/n) %>% arrange(., prop.ext)


s_mat_all_df = rbind(s_mat_df, s_mat_cat_df)
s_mat_all_df$repl = c(rep("all", nrow(s_mat_df)), rep("extr", nrow(s_mat_cat_df)))

### same but group only by intensity and severity of extreme events


cat_fm_df = ext.s7.df %>%  group_by(cat.freq, catastr.mort) %>%  
  summarise(n = n(),
            prop.ext = sum(ext == 1)/n) 


### Plots

plot_repro = ggplot(data = filter(s_mat_all_df, repl == "all"), aes(x = age.sex.mat, y = prop.ext, group = sons.mean,  shape = as.factor(sons.mean)), linetype = 1) + 
  geom_point(size = size.point) +
  geom_line(lwd = line.lwd) +
  theme.pop +
  labs(x = "Age at reproduction (year)", y = "Frequency of extinction") +
  theme(legend.title = element_text(size = size.legend.title)) + 
  scale_shape_manual(name = "# offspring", labels = c("1.0", "1.5", "2.0", "2.5"),values=c(0,1,2,3))  +
  scale_linetype_manual(name = "replicates", labels = c("all", "most extreme"),values=c(1,2)) +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(from = 0, to = 1, by = 0.2)) +
geom_point(data = filter(s_mat_all_df, repl == "extr"), aes(x = age.sex.mat, y = prop.ext, group = sons.mean,  shape = as.factor(sons.mean)), size = size.point) +
  geom_line(data = filter(s_mat_all_df, repl == "extr"), aes(x = age.sex.mat, y = prop.ext, group = sons.mean),lwd = line.lwd, linetype = 2) +
  theme(legend.position = c(0.15, 0.91))
  

plot_cat = ggplot(data = cat_fm_df, aes(x = cat.freq, y = prop.ext, group = catastr.mort, shape = as.factor(catastr.mort))) + 
  geom_point(size = size.point) +
  geom_line(lwd = line.lwd) +
  theme.pop +
  labs(x = "Annual probability of extreme events", y = "Frequency of extinction") +
  theme(legend.title = element_text(size = size.legend.title)) + 
  scale_shape_manual(name = "Mortality", labels = c("0.3", "0.5", "0.7"),values=c(4,5,6)) +
  scale_x_continuous(breaks = c(0.05, 0.10, 0.15)) +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(from = 0, to = 1, by = 0.2))  +
  theme(legend.position = c(0.15, 0.94))


plot_repro_cat = plot_grid(plot_repro,
                         plot_cat,
                          labels = c("(a)", "(b)"),
                          nrow = 1, align = "v",hjust = -5)

plot_repro_cat




# Device = 15,8

