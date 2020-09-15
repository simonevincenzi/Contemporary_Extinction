S.vett = 500
N.vett = 500
iter.vett = 250
selection.vett = c(0.08,0.11)
sdopt.vett = c(0,0.01,0.015)
num.loci.vett = 10
plotorno.vett = 0
meanopt.vett = c(0,0.015)
major.vett = 0
num.alleles.vett = 10
sd.alleles.vett = 0.05
recomb.vett = 0
mutation.vett = 0
mut.alfa.vett = 0.3
catastr.vett = 1
p.pre.vett = 0.05
p.post.vett = c(0.05,0.1,0.15)
catastr.mort.vett = c(0.3,0.5,0.7)
age.sex.mat.vett = c(1,2,3,4)
avg.surv.vett = c(1,1.5) #added 1.5
sons.mean.vett = c(1,1.5,2,2.5)    ##### with 1 every population goes extinct
first.phase.opt.vett = 100
var.amb.init.vett = 1


###### remember that the columns must be named!!!


dataforpar = expand.grid(S = S.vett, N = N.vett ,iter = iter.vett , selection = selection.vett ,sdopt = sdopt.vett ,num.loci = num.loci.vett,plotorno = plotorno.vett,meanopt = meanopt.vett,major = major.vett,num.alleles = num.alleles.vett,sd.alleles = sd.alleles.vett ,recomb = recomb.vett ,mutation = mutation.vett,mut.alfa = mut.alfa.vett,catastr = catastr.vett,p.pre = p.pre.vett,p.post = p.post.vett,catastr.mort = catastr.mort.vett,
                         age.sex.mat = age.sex.mat.vett,avg.surv = avg.surv.vett,sons.mean = sons.mean.vett,first.phase = first.phase.opt.vett,var.amb = var.amb.init.vett) 


replicate.sim = 10

dataforpar = do.call(rbind, replicate(replicate.sim, dataforpar, simplify=FALSE))