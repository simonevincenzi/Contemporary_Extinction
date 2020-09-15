source("load_packages.r")

ext.stand.df.val = readRDS("ext.stand.df.val.RDS")
ext.stand.df.cal = readRDS("ext.stand.df.cal.RDS")

acc_table_base_df = data.frame(model = c("GLM", "GAM", "RF"), accuracy = NA, false_pos = NA, false_neg = NA) # table with accuracy measures

ext.lrm.full = readRDS("ext.lrm.base.RDS")

y = ext.stand.df.cal$ext   ###vector of extinct/persisted

### this function calculate specific and sensitivity for each cut-off point

perf = function(cut, mod, y)
{
  
  if("lrm" %in% class(mod)) {
    fitted = predict(mod,type="fitted")}
  if("gam" %in% class(mod)) {
    fitted = predict(mod,type="response")}
  
  yhat = (fitted>cut)
  w = which(y==1)
  sensitivity = mean( yhat[w] == 1 ) 
  specificity = mean( yhat[-w] == 0 ) 
  c.rate = mean( y==yhat ) 
  d = cbind(sensitivity,specificity)-c(1,1)
  d = sqrt( d[1]^2 + d[2]^2 ) 
  out = t(as.matrix(c(sensitivity, specificity, c.rate,d)))
  colnames(out) = c("sensitivity", "specificity", "c.rate", "distance")
  return(out)
}

####

s = seq(.01,.99,length=1000)
OUT = matrix(0,1000,4)

par(oma =c(1,1,1,1))
for(i in 1:1000) OUT[i,]=perf(s[i],ext.lrm.full,ext.stand.df.cal$ext)
plot(s,OUT[,1],xlab="Cutoff value",ylab="Response",cex.lab=1.5,cex.axis=1.5,ylim=c(0,1),type="l",lwd=3,axes=FALSE,col="gray30") #plotting sensitivity, specificity, classification rate and distance
axis(1,seq(0,1,length=5),seq(0,1,length=5),cex.lab=1.5)
axis(2,seq(0,1,length=5),seq(0,1,length=5),cex.lab=1.5)
lines(s,OUT[,2],col="gray50",lwd=3)
lines(s,OUT[,3],col="gray70",lwd=3)
lines(s,OUT[,4],col="gray90",lwd=3)
box()

cutoff = s[which(OUT[,4]==min(OUT[,4]))]   #### cutoff 
points(cutoff,OUT[which(OUT[,4]==min(OUT[,4])),4],pch=8,cex=1.5) ###plot start to identify the selected cutoff
text(0.38,0.62,bquote(paste("Optimal cutoff = ", .(cutoff[1]))),cex=1.05)
legend(0.2,.5,col=c("gray30","gray50","gray70","gray90"),lwd=rep(3,4),c("Sensitivity","Specificity","Classification Rate","Distance"),box.lwd = 0,box.col = "white",bg = "white")


### predictions with the cutoff with equal sensitivity and specificity (pred_lrm)
### 
ext.stand.df.val$pred_lrm = ifelse(predict(ext.lrm.full,type="fitted", newdata=ext.stand.df.val)>cutoff[1],1,0)

val.table.lrm = with(ext.stand.df.val,table(pred_lrm,ext)) # 2X2 contingency table

accuracy.lrm = sum(val.table.lrm[1,1],val.table.lrm[2,2])/sum(val.table.lrm)

acc_table_base_df$accuracy[which(acc_table_base_df$model == "GLM")] = accuracy.lrm

false.pos.val.lrm = val.table.lrm[2,1]/(val.table.lrm[2,1]+val.table.lrm[1,1]) #false positives for the test dataset

acc_table_base_df$false_pos[which(acc_table_base_df$model == "GLM")] = false.pos.val.lrm

false.neg.val.lrm = val.table.lrm[1,2]/(val.table.lrm[1,2]+val.table.lrm[2,2])

acc_table_base_df$false_neg[which(acc_table_base_df$model == "GLM")] = false.neg.val.lrm


### load gam model 
ext.gam.full = readRDS("ext.gam.base.RDS")

### it takes quite a bit of time to do the cutoff routine for GAM models, so I already saved the results in cutoff. Otherwise, uncomment the section to calculate and save the cutoffs

# s = seq(.01,.99,length=100)
# OUT = matrix(0,length(s),4)
# 
# par(oma =c(1,1,1,1))
# for(i in 1:length(s)) {
#   print(i)
#   OUT[i,]=perf(s[i],ext.gam.full,ext.stand.df.cal$ext)}
# 
# plot(s,OUT[,1],xlab="Cutoff value",ylab="Response",cex.lab=1.5,cex.axis=1.5,ylim=c(0,1),type="l",lwd=3,axes=FALSE,col="gray30") #plotting sensitivity, specificity, classification rate and distance
# axis(1,seq(0,1,length=5),seq(0,1,length=5),cex.lab=1.5)
# axis(2,seq(0,1,length=5),seq(0,1,length=5),cex.lab=1.5)
# lines(s,OUT[,2],col="gray50",lwd=3)
# lines(s,OUT[,3],col="gray70",lwd=3)
# lines(s,OUT[,4],col="gray90",lwd=3)
# box()
# 
# cutoff = s[which(OUT[,4]==min(OUT[,4]))]   #### cutoff identified as minimum of the distance as explained in Supplementary Information #####
# points(cutoff,OUT[which(OUT[,4]==min(OUT[,4])),4],pch=8,cex=1.5) ###plot start to identify the selected cutoff
# text(0.38,0.62,bquote(paste("Optimal cutoff = ", .(cutoff[1]))),cex=1.05)
# legend(0.2,.5,col=c("gray30","gray50","gray70","gray90"),lwd=rep(3,4),c("Sensitivity","Specificity","Classification Rate","Distance"),box.lwd = 0,box.col = "white",bg = "white")
# 
# saveRDS(cutoff, "cutoff_gam_base.RDS")

cutoff = readRDS("cutoff_gam_base.RDS")

### predictions for the GAM model
### 
ext.stand.df.val$pred_gam = ifelse(predict(ext.gam.full,type="response", newdata=ext.stand.df.val)>cutoff[1],1,0)

val.table.gam = with(ext.stand.df.val,table(pred_gam,ext)) # 2X2 contingency table

accuracy.gam = sum(val.table.gam[1,1],val.table.gam[2,2])/sum(val.table.gam)

acc_table_base_df$accuracy[which(acc_table_base_df$model == "GAM")] = accuracy.gam

false.pos.val.gam = val.table.gam[2,1]/(val.table.gam[2,1]+val.table.gam[1,1]) #false positives for the test dataset

acc_table_base_df$false_pos[which(acc_table_base_df$model == "GAM")] = false.pos.val.gam

false.neg.val.gam = val.table.gam[1,2]/(val.table.gam[1,2]+val.table.gam[2,2])

acc_table_base_df$false_neg[which(acc_table_base_df$model == "GAM")] = false.neg.val.gam


### Load the RF model, we don't need cutoff because we use the binary prediction directly
ext.rf.full = readRDS("ext.rf.base.RDS")

### predictions for the RF model
ext.stand.df.val$pred_rf = as.numeric(as.character(predict(ext.rf.full,newdata=ext.stand.df.val)))

val.table.rf = with(ext.stand.df.val,table(pred_rf,ext)) # 2X2 contingency table

accuracy.rf = sum(val.table.rf[1,1],val.table.rf[2,2])/sum(val.table.rf)

acc_table_base_df$accuracy[which(acc_table_base_df$model == "RF")] = accuracy.rf

false.pos.val.rf = val.table.rf[2,1]/(val.table.rf[2,1]+val.table.rf[1,1]) #false positives for the validation dataset

acc_table_base_df$false_pos[which(acc_table_base_df$model == "RF")] = false.pos.val.rf

false.neg.val.rf = val.table.rf[1,2]/(val.table.rf[1,2]+val.table.rf[2,2])

acc_table_base_df$false_neg[which(acc_table_base_df$model == "RF")] = false.neg.val.rf


saveRDS(acc_table_base_df,"acc_table_base_df.RDS")


### Number of predictions that are the same among models
cont = 0

for (i in 1:nrow(ext.stand.df.val)){
  
  if(sum(ext.stand.df.val$pred_gam[i],
         ext.stand.df.val$pred_lrm[i],
         ext.stand.df.val$pred_rf[i]) %in% c(0,3)) {
    
    cont = cont+1
  }
  
}

ratio_const_pred = cont/nrow(ext.stand.df.val)

### some data frames saved with inconsistent or consistent predictions
### 

# test_f1 = filter(ext.stand.df.val, pred_rf == 0, pred_gam == 1, pred_lrm == 1, ext == 1)
# 
# test_f2 = filter(ext.stand.df.val, pred_rf == 0, pred_gam == 1, pred_lrm == 0, ext == 0)
# 
# test_f3 = filter(ext.stand.df.val, pred_rf == 1, pred_gam == 1, pred_lrm == 1, ext == 1)
# 
# test_f4 = filter(ext.stand.df.val, pred_rf == 0, pred_gam == 0, pred_lrm == 0, ext == 1)
# 
# test_f5 = filter(ext.stand.df.val, pred_rf == 1, pred_gam == 1, pred_lrm == 1, ext == 0)