
### Importance plot for Random Forest

size.title = 15
line.lwd = 1
size.label.x = 16
size.text.x = 16
size.point = 4
size.label.y = 16
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


ext.rf.full = readRDS("ext.rf.full.RDS")
imp = (varImp(ext.rf.full,)$importance)

imp = data.frame(var = rownames(imp), importance = imp$Overall)


imp$var = as.factor(imp$var)

imp = imp %>% mutate(var = recode(var, "min_pop_win" = "min(pop)",
                            "mean_pop_win" = "mean(pop)",
                            "max_pheno_dist_win" = "max(pheno_d)",
                            "mean_pheno_dist_win" = "mean(pheno_d)",
                            "max_opt_win"  = "max(opt)",
                            "cat_winy" = "extreme",
                            "age.sex.mat" = "repro",
                            "sons.mean" = "fec",
                            "sel" = "sel"))


min_n_lab = bquote(paste("min(", italic('N'),")"))
mean_n_lab = bquote(bar(italic('N')))
max_pheno_d_lab= bquote(paste("max(", italic(Theta), " - ", italic('z'),")")) 
mean_pheno_d_lab = bquote(paste("mean(", italic(Theta), " - ", italic('z'),")")) 
max_opt_lab = bquote(paste("max(", italic(Theta), ")"))
extreme_lab = bquote(italic('E'))
repro_lab = bquote(italic(a)[f])
fec_lab = bquote(italic(lambda)[0])
sel_lab = bquote(italic('s'))
                    
tt = data.frame(x = as.factor(1:2), y = 1:10)                    
ggplot(data = tt, aes(x= x, y = y)) +
  scale_x_discrete(labels = c(min_n_lab,fec_lab))



plot_imp_full = ggplot(data = imp, aes(x=reorder(var,importance), y=importance)) +
geom_point(size = size.point) +
  geom_segment(aes(x=var,xend=var,y=0,yend=importance)) +
  ylab("Importance") +
  xlab("") +
  scale_x_discrete(labels = c(sel_lab,fec_lab,repro_lab,extreme_lab,max_opt_lab,mean_pheno_d_lab,max_pheno_d_lab,mean_n_lab,min_n_lab)) +                  
    theme.pop + 
  coord_flip()


ext.rf.full = readRDS("ext.rf.base.RDS")
imp = (varImp(ext.rf.full,)$importance)
imp = data.frame(var = rownames(imp), importance = imp$Overall)
imp$var = as.factor(imp$var)

imp = imp %>% mutate(var = recode(var, "min_pop_win" = "min(pop)",
                            "mean_pop_win" = "mean(pop)",
                            "max_pheno_dist_win" = "max(pheno_d)",
                            "mean_pheno_dist_win" = "mean(pheno_d)",
                            "max_opt_win"  = "max(opt)",
                            "cat_winy" = "extreme",
                            "age.sex.mat" = "repro",
                            "sons.mean" = "fec",
                            "sel" = "sel"))

plot_imp_red = ggplot(data = imp, aes(x=reorder(var,importance), y=importance)) +
  geom_point(size = size.point) +
  geom_segment(aes(x=var,xend=var,y=0,yend=importance)) +
  ylab("Importance") +
  xlab("") +
  scale_x_discrete(labels = c(sel_lab,max_opt_lab,extreme_lab,max_pheno_d_lab, mean_pheno_d_lab,repro_lab,fec_lab)) +                  
  theme.pop + 
  coord_flip()




# ext.gam.full = readRDS("ext.gam.full.RDS")
# 
# plot(ext.gam.full)


Plot_imp = plot_grid(plot_imp_full,
                         plot_imp_red,
                         labels = c("(a)", "(b)"),
                         nrow = 1, align = "v",hjust = -5)

Plot_imp

# 15 7
