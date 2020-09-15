#setwd("/Users/mrag/Dropbox/Articoli/Gene/PhenoPlast")

#library("Rlab")
#library("MASS")
#library('hierfstat')

library(Rlab)
library(MASS)
library(parallel)
library(dplyr)

################################

## CYCLE FOR A SINGLE POPULATION  

Pheno.Plast <- function (S = 500,N = 500,iter = 300,selection = 0.05,sdopt = 0.005,num.loci =10,
                         plotorno = 1,meanopt = 0.01,major=0,num.alleles = 5,sd.alleles = 0.1,recomb=1,mutation=1,
                         mu.mut=0.0002,mut.alfa=0.05,
                         catastr = 0, p.pre = 0.05, p.post = 0.07, catastr.mort = 0.3, age.sex.mat = 2,
                         avg.surv = 0.9,sons.mean = 2,
                         first.phase.opt = 150, var.amb.init = 0)
  
  # S = maximum (and initial number of individuals in the population)
  # N = not used
  # iter = total number of simulation steps/year
  # selection = selection strength
  # sdopt = yearly increase of the variance of the climate variable
  # num.loci = number of genes
  # plotorno = 1 if plot of phenotype/optimum/population
  # meanopt = yearly increase of the mean of the distribution of the climate variable
  # major = 1 if there is a major gene coding for the quantitative trait, 0 otherwise
  # num.alleles = num alleles for each locus (randomly drawn at the start of the simulation)
# sd.alleles = standard deviation of the allelic effects (not used, as it is defined in the function)
# recomb = 1 for full recombination, 0 for no recombination
# mutation = 1 if mutations occur, 0 otherwise
# mu.mut = mutation rate per locus
# mut.alfa = mutational variance (variance of the distribution of mutational effects)
# catastr = 1 if point extremes occur, 0 otherwise
# p.pre = annual probability of point extreme before climate change
# p.post = annual probability of point extreme after climate change
# catastr.mort = probability of dying due to the impact of a point extreme
# age.sex.mat = age at sexual maturity (fixed for every individual)
# avg.surv = average maximum probability of annual survival (average because maximum survival each year
# is drawn from a uniform distribution with min and max = avg.surv - or + 1
# sons.mean = mean (and variance) of the Poisson distribution of number of offspring per pair
# first.phase.opt = number of years with no climate change
# var.amb.init = variance of the distribution of the climate variable

#### example of run
# source("phenoplast.r")

# ris.pp <- Pheno.Plast(S = 500,N = 500,iter = 250,selection = 0.09,sdopt = 0.01,
#                       num.loci = 10,plotorno = 1,
#                       meanopt = 0.01,major=0,num.alleles=5,
#                       sd.alleles=0.05,recomb=0,mutation=1, mut.alfa = 0.3,catastr = 1, p.pre = 0.0, 
#                       p.post = 0.05, 
#                       catastr.mort = 0,age.sex.mat = 1,avg.surv = 0.7,sons.mean = 4,first.phase.opt = 100,
#                       var.amb.init = 1)




{
  
  extinct = 0
  options(warn = 0) ###stop when there is warning when warn = 2
  
  ### I define the standard deviation of distribution of single alleles, which depends on heritability at time 1
  ### and environmental variance
  
  ### this can be adapted to other traits that have means != 0
  # the correct way is to proceed from a phenotypic variance and heritability and then
  # obtain the starting additive genetic variance
  var.env = 1
  heritability = 0.3
  starting_genetic_variance = (heritability * var.env)/(1-heritability)
  
  #starting_genetic_variance = 0.2
  
  
  
  #sd.alleles = (heritability * var.env)/(2*num.loci)
  # this is the standard deviation of the normal distribution of allelic effects 
  sd.alleles = sqrt(starting_genetic_variance/(2*num.loci))
  
  ############ Define values of some model parameters not included in the function
  
  first.phase.opt = first.phase.opt ## starting years with no climate change
  
  second.phase.opt = iter-first.phase.opt ## second part with climate change (iter is the total number
  # of simulation years)
  
  
  ## create a vector of random annual survival probabiity - avg.surv is given as input
  avg.surv.vett = runif(iter,min = avg.surv - 0.1,max=avg.surv + 0.1)
  
  catastr.vett = c(rbern(first.phase.opt,p.pre),rbern(second.phase.opt,p.post))
  
  end_increase_variance = 25 #####years of increase variance of the climate variable after climate change
  # , then it remains constant up to the end of the simulation
  
  var.amb.init = var.amb.init ####### should be 1variance of the climate variable up to climate change
  ##### define properties of individuals (position on tehe individual vector) - sector, age, sex ######
  
  ## the population stays on a matrix (area.pop), each column is an individual. Here below I define
  # the position of entities in the column
  
  pheno.pos = 1  ## this means that the first row is the phenotype of the individuals
  age = 2  ## second row is the age of the individual
  #sex = 3
  #sector = 4
  #plastic = 3
  
  length_prop = 2 #  length of entities of individual vector excluding loci 
  # (phenotipic value and age, first 2 top rows)
  
  
  heteroz.mat=matrix(0,num.loci*num.alleles,(iter%/%50))  ##### define the matrix for heterozigosyis, same as all.freq100 but with a col less
  
  area.pop<-matrix(0,length_prop+(num.loci*2),S) ##  matrix space area.pop of S elements, 
  # for each col the first row is the pheno value, second is age,  
  #the other num.loci*2 rows are the values for each gene
  #the first num.loci (2-num.loci+1)for the first chromosome, the second num.loci (num.loci+2-(2*num.loci)+2) for the second chromosome.
  
  fchrom = seq((length_prop+1),num.loci+length_prop,1) ##identifier of rows for the 1st chromosome
  schrom = seq(num.loci+(length_prop+1),num.loci*2+(length_prop),1) ##identifier of rows for the 2nd chromosome
  
  ### discard
  #######################assign age,sex,sector##########
  #sectors = c(1,2,3,4,5,6,7)  ##sectors in the stream
  #firstsect.size = round(rep(S/(length(sectors)),(length(sectors)-1)))
  
  #maxsize.sect = c(firstsect.size,(S-(sum(firstsect.size))))
  
  #initial_ages = rep(1,S) ### start with all fish of age 1
  ## end discard 
  
  ## define initial ages for the S individuals in the population at time = 1
  # I defined some criteria
  initial_ages = round(runif(S,min = 1, max = max(age.sex.mat,1)*1.3))
  
  # initial_plastic = rnorm(S,0,1)
  
  
  #initial_sex = sample(c(0,1),S,replace = T) ###0 for females, 1 for females
  #initial_sector = sample(sectors,S,replace=T) #### randomly draw sectors
  #initial_sector = rep(6,S) #### randomly draw sectors
  
  area.pop[age,] = initial_ages # assign ages to individuals
  #area.pop[plastic,] = initial_plastic
  #area.pop[sex,] = initial_sex
  #area.pop[sector,] = initial_sector
  
  
  
  
  #####CREATE THE MATRIX OF ALLELES IN THE POPULATION##########
  alleles.mat = matrix(rnorm(num.loci*num.alleles,0,sd.alleles), 
                       num.loci,num.alleles)  
  
  #in matrix alleles.mat I extract randomly, for each of the num.loci, num.alleles number of alleles 
  
  alleles.mat.keep = alleles.mat  
  
  
  alleles.mat.gen = matrix(0, num.loci,num.alleles)
  
  for (nl in 1:num.loci) {
    
    alleles.mat.gen[nl,] = seq(1:num.alleles)
    
  } ###I assign a number from 1 to num.alleles to each allele (for Fst since I need the genotype)
  
  #####MAJOR GENE####
  
  if (major==1) {
    alleles.mat[1,] = rnorm(num.alleles,1.5,0.2) ##major genes is in the first locus and it has 1.5 as a mean value
  }
  
  #####CLOSE MAJOR GENE####
  #double the matrix alleles.mat for ease of extraction of alleles for the diploid organism 
  doublealleles.mat = rbind(alleles.mat,alleles.mat)
  
  ##########assign alleles to individuals at the start of simulation
  for (a in 1:(num.loci*2)) { #assign alleles to loci in the population matrix area.pop
    
    area.pop[length_prop+a,] = sample(doublealleles.mat[a,],dim(area.pop)[2],replace=T) } #end assignment 
  
  area.pop[1,] = colSums(area.pop[(length_prop+1):dim(area.pop)[1],]) #the first row of the population matrix is the sum of the genetic values of each alleles 
  #of an individual
  ############################################ 
  pheno = area.pop[1,] + rnorm(dim(area.pop)[2],0,sqrt(var.env)) 
  # phenotypic value = genotypic value + environmental deviate from N(0,sqrt(var.env)) (all stored in pheno)
  
  # print(area.pop)
  
  f = rep(0,S) ## vector to be filled with  individuals fitness (each year is recycled)
  p = rep(0,S) ## vector to be filled with survival probability (each year is recycled)
  
  
  #####create vectors and matrix for later use
  vettopt = rep(0,iter)  # vector with realized values of the climate variable
  mediaphen = rep(0,iter) # vector with mean of the phenotype (one value each year)
  sdphen = rep(0,iter)  # vector with standard deviation of the phenotype (one value each year)
  optimum = rep(0,iter)  # vector with realized values of the climate variable (same as vettopt)
  varphen = rep(0,iter) # vector with variance of the phenotype
  mean.age = rep(0,iter) # vector with mean age at each time step
  
  mutant_surviving = 100  ####initialize to 100 the number of mutant alleles surviving, so I can differentiate between a 
  #value of 0 when there are no mutants surviving and when there was a break due to extinction
  
  vettad = rep(0,iter)    #######population size post_selection at each time step
  vettad_pre = rep(0,iter)  ###### population size pre_selection at each time step
  addgenvar = rep(0,iter)   ####additive genetic variance of the whole population at each time step
  addgenmean = rep(0,iter)  ### mean genetic value in the whole population at each time step
  phenomean = rep(0,iter)  #### mean of the phenotype pre_selection, including the newborns of the year i-1
  mean.var.age = matrix(0,iter,2)  ## matrix with iter row and 2 cols (1 = mean age, 2 = sd of age every year)
  
  
  w_gen_mean = rep(0,iter)  ###### mean of the phenotype in a new generation at each time step
  w_gen_var = rep(0,iter)		###### sampling variance of the phenotype in a new generation at each time step
  
  heritability_paroff = rep(100,iter) #####vector to record heritabilities as parent offspring slope of regression (cov/var)
  num.sons.check = matrix(0,2,iter) #check the number of offspring recruited in the population
  herit.coef = rep(0,iter) ###vector initialized to 0 with slope of the linear model 
  
  fitness.mean = rep(0,iter) # vector with mean of fitness at each time step
  fitness.variance = rep(0,iter) # vector with variance of fitness at each time step
  #Fst.sect = rep(100,100)
  #Fst.time = list()
  
  ###########
  
  
  all.freq100 = matrix(0,num.loci*num.alleles,(iter%/%50)+3) ## initialize matrix with 
  #frequencies of alleles during simulation, every 50 years I record the value
  dim(alleles.mat) = c(num.loci*num.alleles,1) ##unravel matrix with alleles (alleles.mat) in order to have a vector
  
  heter_mat_year = matrix(0,num.loci*num.alleles,iter) # saves allelic frequency each year
  heter_vect_mean_year = rep(0,iter) #saves mean heterozygosity each year
  heter_vect_sd_year = rep(0,iter) #saves sd of heterozygosity each year
  
  
  col.freq = 2 # counter for column number when checking for frequencies of alleles at different time steps  
  
  mutants.mat = matrix(0,num.loci*num.alleles*(iter/5),(iter%/%50)+3)  # initialize matrix to be filled with mutant alleles, 
  # year 1 + every 50 years + difference in frequency at the end
  
  mutants.check = matrix(0,num.loci*num.alleles*(iter/5),1) #initialize vector to be filled with alleles that mutate
  row.mutants=0 #counter for rows when filling the matrix of mutants
  
  
  a = 0
  #t10=1
  
  ##################### START FOR TIME LOOP ########################
  
  
  for(i in 1:iter)  ##time-steps or cohorts or year
    
  {
    print(i)  # prints time-steps
    
    ###### I compute mean and variance of the phenotype for the first generation
    if (i==1) {
      w_gen_mean[i] = mean(pheno)
      w_gen_var[i]  = var(pheno)
    }
    #####################################
    
    #########check if the population went extinct####
    
    #print("a")
    
    if(length(area.pop[1,area.pop[1,]!=0])<2){  ##is population extinct?  #open for to check for extinction
      print(paste("Extinct at year", i)) #prints the year of extinction
      extinct = 1 # when an individual is dead the first column is 0
      
      
      if(iter>=50) {
        
        
        
        all.freq100[,dim(all.freq100)[2]] = all.freq100[,col.freq] - all.freq100[,2]  #in case of population going extinct, it fills the 
        #last col of all.freq100 with the difference between 
        heteroz = 0		
      }#frequency at time 50 and frequency at the last year to be divided
      #by 50 before extinction
      break}   #close for to check for extinction
    
    ######################################  
    
    ### here below I create the matrix of allele frequencies every 50 years starting from year 1
    
    if (i==1) { #if first year of the simulation
      
      
      all.freq100[,1] = alleles.mat  #first column is the vector of alleles
      
      for (zz in (1:(num.loci*num.alleles))) {  #for each alleles in the population (present at time 1)
        
        
        all.freq100[zz,col.freq] = length(area.pop[area.pop==alleles.mat[zz]])/
          (2*length(area.pop[1,area.pop[1,]!=0])) #divide the number of alleles in the 
        #population by 2 times the number of individuals alive
        # (diploid organism)
      }
      
    }
    
    
    #print(1)
    if (i%%50==0) {  #if year can be divided by 50
      
      col.freq=col.freq+1  #update counter of rows
      #print(alleles.mat)
      for (zz in (1:(num.loci*num.alleles))) {  #for each alleles in the population (present at time 1)
        
        
        all.freq100[zz,col.freq] = length(area.pop[area.pop==alleles.mat[zz]])/(2*length(area.pop[1,area.pop[1,]!=0])) #divide the number of alleles in the 
        #population by 2 times the number of individuals alive
        #print(i)																							# (diploid organism)
      }
      #print(all.freq100)  
      
    }
    
    
    
    
    if (i==iter & iter>=50){ #if last year of the simulation
      
      prova = area.pop
      all.freq100[,dim(all.freq100)[2]] = all.freq100[,(dim(all.freq100)[2]-1)] - all.freq100[,2] #it fills the 
      #last col of all.freq100 with the difference between 
      #allelic frequency at time 50 and frequency at the last year to be divided
      #by 50 before extinction
      
      heteroz.mat = all.freq100[,(dim(all.freq100)[2]-1)]   #matrix for heterozigosity (no mutants included)
      
      if (length(mutants.mat[mutants.mat[,1]!=0,1])>0) {  #if there are mutants (otherwise mutants.mat is all zero)
        
        mutants.mat = mutants.mat[mutants.mat[,1]!=0,]  #strip away all zeros from the matrix of mutants (in the first column)
        
        for (zz in (1:dim(mutants.mat)[1])) { # cycle for all mutants  
          mutants.mat[zz,(dim(all.freq100)[2]-1)] = length(area.pop[area.pop==mutants.mat[zz,1]])/(2*length(area.pop[1,area.pop[1,]!=0]))
        }
        
        
        all.freq100 = rbind(all.freq100,mutants.mat[mutants.mat[,(dim(all.freq100)[2]-1)]!=0,])  ### check this
      }
      
      
    }
    
    a = a+1   ## counter for the creation of matrix of allelic values
    
    
    #### allelic frequency each year
    
    for (zz in (1:(num.loci*num.alleles))) {  #for each alleles in the population (present at time 1)
      
      
      heter_mat_year[zz,i] = length(area.pop[area.pop==alleles.mat[zz]])/(2*length(area.pop[1,area.pop[1,]!=0])) #divide the number of alleles in the 
      #population by 2 times the number of individuals alive
      #print(i)																							# (diploid organism)
    }
    
    
    
    
    ########### optimum before and after climate change #######
    # remember that optimum = value of the climate variable
    if  (i<first.phase.opt)
    {optimum[i] = rnorm(1,0,sqrt(var.amb.init))}  # optimum of the environment up to year first.phase.opt
    ## the mean is 0, but it can be changed if we want to simulate a real climate variable
    
    else
    {
      #optimum[i] <- rnorm(1,0+(i-(iter/2))*meanopt, 0.5+(i-(iter/2))*sdopt) # original is 0.01, try with double increase in direction
      
      meanopt.increase = 0 + (i-first.phase.opt)*meanopt # mean of the distribution of the climate variable
      # after climate change
      sdopt.increase = sqrt(var.amb.init) + min((i-(first.phase.opt)), end_increase_variance)*sdopt
      # standard deviation of the mean of the distribution of the climate variable after climate changed
      
      optimum[i] = rnorm(1,meanopt.increase, sdopt.increase)
      
      #varianceopt <- 0.5+(i-(iter/2))*sdopt
      
    }  #open and close else
    
    
    ##### close before and after climate change #######
    
    
    ##### compute mean age excluding newborns ######
    
    mean.age.who = which(area.pop[2,]>=0 & area.pop[1,]!=0)
    
    mean.age[i] = mean(area.pop[2,mean.age.who],na.rm = T)
    
    
    
    
    ############increase age by 1 if t>1 ##########
    
    if (i>1) {
      area.pop[age,] = area.pop[age,] + 1 #### since eggs are -1, juveniles are 0 and so on
    }
    ###########
    
    # print("b")
    
    
    ####vector of optimum and mean and sd of phenotype in the population ######
    
    vettopt[i] = optimum[i]
    
    mediaphen[i] = mean(area.pop[1,area.pop[1,]!=0])
    
    sdphen[i] = sd(area.pop[1,area.pop[1,]!=0])
    
    #### close vector of optimum and mean and sd of phenotype in the population ######
    
    #print(2)
    
    #######PLOT GRAPHS######################
    
    ## plots are every 50 years starting from year 150
    
    if (i%%50 == 0 & plotorno==1 & i>=150) {   # open if for plot
      
      if (i==150) { # open if to check for first plot
        
        #quartz()
        #png("Sim.png",700,350)
        
        ## 6 panels, this can be modified according to goals  
        
        par(mfrow=c(2,3),mar=c(5,5,3,2)) 
        title <- bquote(paste("Year ", .(a))) 
        
        p1 <- hist(pheno[area.pop[1,]!=0], plot=F)
        
        #p2 <- hist(area.pop[1,area.pop[1,]!=0], plot=F, breaks= seq(-3,3,0.4))
        
        
        plot(p1, col="white", xlim=c(-6,6), xlab="z", main = title, ylab="Frequency", cex.main=1.5,cex.lab=1.5,cex.axis=1.6,freq=F,ylim=c(0,1))  # first histogram
        #plot(p2, col=rgb(1,0,0,1/4), add=T, xlim=c(-6,6), cex.main=2,cex.lab=1.6,cex.axis=1.6,freq=F)
        
        abline(v = mean(pheno[area.pop[1,]!=0]),lty=2,lwd=2, col="black")
        
        #abline(v = mean(area.pop[area.pop!=0]),lty=2,lwd=3, col="black")
        
        distp = -3.4 - 0.4*(p.post>0.09)
        distb = -3.9 + 0.2*((sdopt==0.005 | sdopt==0.015))
        dists = -4.4 - 0.1*(selection>0.08)
        
        text(distp,1,bquote(paste("p(",E[a], ") = ",.(p.post))),cex=1.5)
        text(dists,0.88,bquote(paste("s = ",.(selection))),cex = 1.5)
        text(distb,0.76,bquote(paste(beta[sigma],"," [Theta], " = ",.(sdopt))),cex = 1.5)
        text(-3.9,0.64,bquote(paste(beta[mu],"," [Theta], " = ",.(meanopt))),cex = 1.5)
        
        
        
        
        
        
        
        mutsel.distr <- fitdistr(pheno[area.pop[1,]!=0],"normal")
        
        
        
      }  # close if to check for first plot
      
      else {title <- bquote(paste(.(a))) 
      
      
      p1 <- hist(pheno[area.pop[1,]!=0], plot=F)
      
      #p2 <- hist(area.pop[area.pop[1,]!=0],plot=F, breaks= seq(-3,3,0.4))
      
      
      plot(p1, col="white", xlim=c(-6,6), xlab="z", main = title, ylab="", cex.main=1.5,cex.lab=1.4,cex.axis=1.6,freq=F,ylim=c(0,1))  # first histogram
      #plot( p2, col=rgb(1,0,0,1/4),add=T, xlim=c(-6,6) , cex.main=2,cex.lab=1.6,cex.axis=1.6,freq=F)
      
      abline(v = mean(pheno[area.pop[1,]!=0]),lty=2,lwd=2, col="black")
      
      #abline(v = mean(area.pop[area.pop!=0]),lty=2,lwd=3, col="black")
      
      }
      
    }  # close if for plot
    
    
    ############ CLOSE plots ############
    
    
    
    # print("c")
    
    addgenvar[i] <- var(area.pop[1,area.pop[1,]!=0])
    addgenmean[i] <- mean(area.pop[1,area.pop[1,]!=0])
    phenomean[i] <- mean(pheno[area.pop[1,]!=0])
    varphen[i] = var(pheno[pheno!=0])
    
    
    vettad_pre[i] = length(area.pop[1,area.pop[1,]!=0]) ########population size pre_selection
    
    ############### MEAN AND SD OF AGE #####
    
    who = which(area.pop[1,]!=0) 
    # where the individuals are (because empty columns have zero in first position)
    mean.var.age[i,1] = mean(area.pop[age,who]) #mean age
    mean.var.age[i,2] = sd(area.pop[age,who]) #sd of age
    
    
    
    ############mortality##############
    
    f = rep(0,S) 
    p = rep(0,S)
    
    who = which(area.pop[1,]!=0)
    # where the individuals are (because empty columns have zero in first position) 
    
    f[who] = optimum[i]-pheno[who] # distance between optimum and individual phenotype
    
    
    
    p[who] = (exp(-(selection*(f[who])^2))) * avg.surv.vett[i]  ########important!! 
    #multiply x average (max) survival to obtain the individual survival
    
    fitness.mean[i] = mean(exp(-(selection*(f[who])^2)))
    
    fitness.variance[i] = var(exp(-(selection*(f[who])^2)))
    
    r = runif(S,0,1)
    
    vettad_pre[i] = length(area.pop[1,area.pop[pheno.pos,]!=0]) ########population size pre_selection including juveniles
    
    dead <- which(r>p)   #### since the p of absent columns are also empty, they are always "dead"
    area.pop[,dead] <- 0 ##remove who dies, ie all columns gets 0
    pheno[dead] <- 0  #phenotypes of the dead (thus empty cols) are 0
    vettad[i] = length(area.pop[1,area.pop[pheno.pos,]!=0])  ########population size post_selection includig juveniles  
    
    
    
    
    
    #### in case of scenario with mass mortality and mass mortality happening
    #### I draw individuals alive that will die, fraction given by catastr.mort
    if (catastr != 0 & catastr.vett[i]>0) {  #scenario with catastrophes and catastrophes happening
      catastr.dead = sample(which(area.pop[1,]!=0),size = catastr.mort * length(which(area.pop[1,]!=0)),replace=F) # who of the guys alive die
      area.pop[,catastr.dead]  = 0 ##remove who dies (all columns of the dead go to 0)
      pheno[catastr.dead] = 0 #### 
      
    }
    
    ########## close mortality ##############
    
    
    ####################   REPRODUCTION ##########################
    
    
    # print("d")
    nsonscont = 0
    
    cosons = 0
    
    vettsons = rep(100,S*2)  # define with 100s and then delete those 100s when manipulating the vector
    matsons = rbind(rep(100,S*2),matrix(0,(length_prop+(num.loci*2)-1),S*2))# define with 100s the first and then delete the remaining 100s 
    #when manipulating the vector
    pheno_parents = rep(100,S*2)  ##empty vector that will contain the mean phenotype of parents
    #sex.from.parents = sample(c(0,1),S*2,replace=T)
    #matsons[(sex),] = sex.from.parents 
    matsons[(age),] = -1  # age of newborns is -1, so next year they will be 0 (it is fish related, -1 is basically the egg,
    # but it can be changed and the age of newborns can be 0)
    
    reproductor = which(area.pop[1,]!=0 & area.pop[age,]>=age.sex.mat)       ##### vector of reproductor
    
    
    ####this implies that every mature individual reproduces, it can be modified by
    ## randomly select a proportion of the potential reproductors actually reproducing
    ## think that pheno value has also an effect on reproduction
    ## think about another trait related to either fecundity or prob of mating   
    
    gene.from.parents = rbind(rep(100,S*2),matrix(0,(num.loci*2*2)-1,S*2))
    number.of.sons = rep(100,S*2)
    
    if (length(reproductor) >1)   { #open if for reproduction
      
      
      #########open control for odd number of reproductors
      
      if (length(reproductor)%%2 != 0) {		#open and close control for odd number of reproductors, since I assume that a female mates with only one male
        
        reproductor <- reproductor[-sample(1:length(reproductor),1)]
      }  # close for odd number of reproductors 
      
      
      ############## close control for odd number of reproductors #########
      
      repro.matrix <- matrix(sample(reproductor,length(reproductor)),2,length(reproductor)/2) 
      ## produce the matrix of reproductors, each column is a couple 
      
      #print(4)
      
      
      ######### begin loop for reproductors ###########
      
      
      for(rr in 1:length(repro.matrix[1,])) {   #begin loop for of reproduction, go on for all the couples
        
        nsonscouple <- rpois(1,sons.mean)   ### this is the number of newborns for a couple, 
        #it is a deviate from a poisson distribution with sons.mean mean and variance
        
        
        
        if(nsonscouple>0) {  # open if for couples with at least one offspring
          
          
          
          #matsons[sector,(nsonscont+1):(nsonscont + nsonscouple)] = sector.from.parents #assign the sector to offspring
          
          
          #nsonscont = nsonscont + nsonscouple
          
          gene.from.father = matrix(0,num.loci*2,nsonscouple)  # alleles coming from father (empty) used in case of recombination
          gene.from.mother = matrix(0,num.loci*2,nsonscouple) # alleles coming from mother (empty) used in case of recombination
          
          gene.sons = matrix(0,num.loci*2,nsonscouple)  #create a matrix with cols = number of offspring of the couple and rows = num.loci * 2
          
          ####### open recombination####
          
          if (recomb==0) {  #whether to recombine or not
            for (aa in 1:nsonscouple) {   #cicle for number of offspring
              
              sampleseq = sample(c(1,2),2,replace=T)    #draw the chromosome to be passed to offspring (First for father, second for mother)
              
              if (sampleseq[1]==1) {gene.sons[1:num.loci,aa] = area.pop[fchrom,repro.matrix[1,rr]]}  ## if 1 for father I pick the first chromosome (fchrom = 2:num.loci+1, schrom = num.loci + 2:max) 
              else if (sampleseq[1]==2) {gene.sons[1:num.loci,aa] = area.pop[schrom,repro.matrix[1,rr]]} # if 2 for father I pick the second chromosome
              
              if (sampleseq[2]==1) {gene.sons[(num.loci+1):(num.loci*2),aa] = area.pop[fchrom,repro.matrix[2,rr]]} ## if 1 for mother I pick the first chromosome (fchrom = 2:num.loci+1, schrom = num.loci + 2:max) 
              else if (sampleseq[2]==2) {gene.sons[(num.loci+1):(num.loci*2),aa] = area.pop[schrom,repro.matrix[2,rr]]} # if 2 for mother I pick the second chromosome
              
            }
          }  #close for recombination = 0
          
          #### here starts for recombination (recombination ==1)
          
          else if (recomb==1) { #### recombination ==1
            
            
            for (aa in 1:nsonscouple) {   #cicle for number of offspring
              
              rchrom.f = 0
              rchrom.m = 0
              
              
              sampleseq = sample(c(1,2),num.loci,replace=T)    #draw the chromosome to be passed to offspring (First for father, second for mother)
              rchrom.f = diag(cbind(fchrom,schrom)[,sampleseq])
              gene.sons[1:num.loci,aa] = area.pop[rchrom.f,repro.matrix[1,rr]]  ## if 1 for father I pick the first chromosome (fchrom = 2:num.loci+1, schrom = num.loci + 2:max) 
              
              
              gene.from.father[,aa] = area.pop[rchrom.f,repro.matrix[1,rr]] 
              
              
              sampleseq = sample(c(1,2),num.loci,replace=T)    #draw the chromosome to be passed to offspring (First for father, second for mother)
              rchrom.m = diag(cbind(fchrom,schrom)[,sampleseq])
              gene.sons[(num.loci+1):(num.loci*2),aa] = area.pop[rchrom.m,repro.matrix[2,rr]] ## if 1 for mother I pick the first chromosome (fchrom = 2:num.loci+1, schrom = num.loci + 2:max) 
              # if 2 for mother I pick the second chromosome
              
              gene.from.mother[,aa] = area.pop[rchrom.m,repro.matrix[2,rr]]
              
            }
          } ##close for recombination == 1
          
          ##### close recombination options #########
          
          
          matsons[((length_prop+1):nrow(matsons)),(cosons+1):(cosons + nsonscouple)] <- gene.sons
          
          matsons[1,(cosons+1):(cosons + nsonscouple)] <- 0
          pheno_parents[(cosons+1):(cosons + nsonscouple)] = mean(pheno[repro.matrix[1,rr]],pheno[repro.matrix[1,rr]])
          
          gene.from.parents[,(cosons+1):(cosons + nsonscouple)] = rbind(gene.from.father,gene.from.mother)
          
          number.of.sons[(cosons+1):(cosons + nsonscouple)] = nsonscouple
          
          cosons <- length(matsons[1,matsons[1,]!=100])      
          
        }  # close if  for couples with at least one offspring
      }   #close loop for reproduction in one sector
      
    } ### close for if (length(reproductor) >1)
    
    empty<-which(area.pop[1,]==0)  ### identify the available spots for newborns
    
    #print(i)
    
    ##### open if there are offspring produced ########
    
    
    if (length(matsons[1,matsons[1,]!=100]) >1) ##if there are offspring produced
    { 
      
      #print(1.5)
      matsons = matsons[,matsons[1,]!=100]  ## delete the empty places in matsons
      pheno_parents = pheno_parents[pheno_parents!=100]  ##delete the empty places in pheno_parents
      number.of.sons = number.of.sons[number.of.sons!=100]  ##delete the empty places in number.of.sons
      gene.from.parents = gene.from.parents[,gene.from.parents[1,]!=100] ##delete the empty places in gene.from.parents
      
      #########MUTATION######################
      
      #print("e")
      
      if (mutation == 1 && i<iter) {     ### if the options is for mutation and if simulation time is not the last year
        
        mutants = 0
        sample.locus = 0
        mut.off = 0
        num.mutants = 0
        confr.mut = 0
        mutants.prov = 0
        
        confr.mut = runif(dim(matsons)[2])
        if (length(confr.mut[confr.mut<(mu.mut*num.loci)])>0) {
          
          num.mutants = length(confr.mut[confr.mut<(mu.mut*num.loci)])
          
          mut.off = which(confr.mut<(mu.mut*num.loci))
          
          sample.locus = sample(((length_prop+1):nrow(matsons)),num.mutants,replace=F)
          if (num.mutants == 1) {
            
            
            mutants.prov = matsons[sample.locus, mut.off]
            matsons[sample.locus, mut.off] = matsons[sample.locus, mut.off] + rnorm(1,0,sqrt(mut.alfa))
            
            mutants = matsons[sample.locus, mut.off]
            
          }
          
          
          else {
            
            for (yy in 1:num.mutants) {
              
              mutants.prov[yy] = matsons[sample.locus[yy],mut.off[yy]]
              matsons[sample.locus[yy],mut.off[yy]] = matsons[sample.locus[yy],mut.off[yy]] + rnorm(1,0,mut.alfa)
              
              mutants[yy] = matsons[sample.locus[yy],mut.off[yy]]
              
            }
          }
          
          
          #print("f")
          
          mutants.mat[(row.mutants+1):(row.mutants + num.mutants),1] = mutants
          
          mutants.mat[(row.mutants+1):(row.mutants + num.mutants),2] = i
          
          
          mutants.mat[(row.mutants+1):(row.mutants + num.mutants ),3] =num.mutants 
          
          
          mutants.check [(row.mutants+1):(row.mutants + num.mutants )]=mutants.prov
          
          
          row.mutants <- length(mutants.mat[mutants.mat[,1]!=0,1])
          
          
        }
        
        
        
      }  ### this closes the mutation
      
      ###########END MUTATION###########
      
      sons = 1  ## this is related to if (length(matsons[1,matsons[1,]!=100]) >1)
    }
    else # if there are no newborns produced
    { matsons[,] == 0  
      sons=0
    }
    
    ####### end of if (length(matsons[1,matsons[1,]!=100]) >1)  (including the else)
    
    
    
    
    vettpheno = rep(50,S*2) #### I prepare the empty vector of newborn phenotypes (50 means there is nothing)
    
    
    #### vector of offspring is empty if sons == 0 or no space is available #####
    
    if(sons==0 | length(empty)==0){
      vettpheno=0
      
      #print(9)
    }  #open and close  ############sistemare
    
    ####### close of offspring is empty if sons == 0 or no space is available #####
    
    
    
    ######## open if there are sons ######
    
    else if (sons==1){
      
      if(length(empty)>dim(matsons)[2]){
        area.pop[1:dim(area.pop)[1],empty[1:dim(matsons)[2]]] = matsons
        
        area.pop[1,empty[1:dim(matsons)[2]]] = colSums(area.pop[(length_prop+1):dim(area.pop)[1],
                                                                empty[1:dim(matsons)[2]]])
        vettpheno[1:dim(matsons)[2]] <- area.pop[1,empty[1:dim(matsons)[2]]] + 
          rnorm(dim(matsons)[2],0,sqrt(var.env)) 
        
        vettpheno = vettpheno[vettpheno!=50]
        
        pheno[empty[1:dim(matsons)[2]]] <- vettpheno
        
        pheno_parents_unique = unique(pheno_parents)
        
        off.mean = rep(0,length(pheno_parents_unique))
        
        sons.for.couple = rep(0,length(pheno_parents_unique))
        
        
        
        for (hh in 1:length(pheno_parents_unique))
          
        {off.mean[hh] = mean(vettpheno[pheno_parents == pheno_parents_unique[hh]]) 
        
        sons.for.couple[hh] = min(number.of.sons[pheno_parents == pheno_parents_unique[hh]])
        
        
        if (i==50) {
          
          check.parents  = rbind(gene.from.parents,matsons)
          
          
        }
        
        
        }
        
        
        
        if (length(pheno_parents)>1 & length(vettpheno)>1){
          #heritability_paroff[i] = cov(vettpheno,pheno_parents)/var(pheno_parents)
          heritability_paroff[i] = cov(off.mean[sons.for.couple>1],pheno_parents_unique[sons.for.couple>1])/var(pheno_parents_unique[sons.for.couple>1])
          
        }
        
        num.sons.check[1,i] = 1   ##this tells me that there were more empty spots than offspring (when 0 is more offspring than spots)
        num.sons.check[2,i] = length(vettpheno) ##this tells me how many offspring were introduced
        #herit.coef[i] = summary(lm(vettpheno  ~ pheno_parents))$coef[2]
        #print(10)
        
        
        
        ####plot of parent-offspring in case heritability is < 0 
        #if (heritability_paroff[i]<0) {
        #dev.new()
        #plot(off.mean[sons.for.couple>1]  ~ pheno_parents_unique[sons.for.couple>1],pch=16,main = bquote(.(i)))
        
        #}
        ######### 
        
        
      } ## space is more than the number of newborns, puts all sons     #open and close
      
      
      
      else if(length(empty)==1){   #this is special case because I cannot use colSums with only one column
        
        
        
        area.pop[1:dim(area.pop)[1],empty]<-matsons[,1:length(empty)]
        area.pop[1,empty] = sum(area.pop[(length_prop+1):dim(area.pop)[1],empty])
        
        vettpheno[1:length(empty)] = area.pop[1,empty] + rnorm(length(empty),0,sqrt(var.env))
        
        vettpheno = vettpheno[vettpheno!=50]
        
        pheno[empty]<-vettpheno
        
      }
      
      
      
      
      else if(length(empty)<= dim(matsons)[2]){
        
        
        
        area.pop[1:dim(area.pop)[1],empty]<-matsons[,1:length(empty)]
        area.pop[1,empty] = colSums(area.pop[(length_prop+1):dim(area.pop)[1],empty])
        
        vettpheno[1:length(empty)] = area.pop[1,empty] + rnorm(length(empty),0,sqrt(var.env))
        
        vettpheno = vettpheno[vettpheno!=50]
        
        pheno[empty]<-vettpheno
        
        pheno_parents_unique = unique(pheno_parents[1:length(empty)])
        
        off.mean = rep(0,length(pheno_parents_unique))
        sons.for.couple = rep(0,length(pheno_parents_unique))
        
        
        
        for (hh in 1:length(pheno_parents_unique))     ### compute the mean phenotype of newborns from a couple
          
        {off.mean[hh] = mean(vettpheno[pheno_parents[1:length(empty)] == pheno_parents_unique[hh]]) 
        
        sons.for.couple[hh] = min(number.of.sons[pheno_parents[1:length(empty)] == pheno_parents_unique[hh]])
        
        }
        
        
        if (length(pheno_parents)>1 & length(vettpheno)>1){
          #heritability_paroff[i] = cov(vettpheno,pheno_parents[1:length(empty)])/var(pheno_parents[1:length(empty)])
          heritability_paroff[i] = cov(off.mean[sons.for.couple>1],pheno_parents_unique[sons.for.couple>1])/var(pheno_parents_unique[sons.for.couple>1])
        }
        
        
        ####plot of parent-offspring in case heritability is < 0
        #if (heritability_paroff[i]<0) {
        #dev.new()
        #plot(off.mean[sons.for.couple>1]  ~ pheno_parents_unique[sons.for.couple>1],pch=16,main = bquote(.(i)))
        
        #}
        ##############
        
        num.sons.check[2,i] = length(vettpheno) ##this tells me how many offspring were introduced
        #herit.coef[i] = summary(lm(vettpheno  ~ pheno_parents[1:length(empty)]))$coef[2]
        #print(13)
        
        
        if (i==50) {
          
          check.parents  = rbind(gene.from.parents[,1:length(empty)],matsons[,1:length(empty)])
          
          
        }
        
        
      }
      
    }
    
    
    ####here I compute the mean and the variance of the new generation  phenotypes####
    if (i>1 & sons==1 & length(empty)>0) {
      w_gen_mean[i] = mean(vettpheno)
      w_gen_var[i]  = var(vettpheno)
    }
    
    
    
    
    ####################   CLOSE REPRODUCTION ##########################
    
  }  ##close time
  
  ### heteroz  vector for each locus at the end of simulation time
  if (i == iter & extinct ==0) {
    
    dim(heteroz.mat) =  c(num.loci,num.alleles)
    heteroz = matrix(0,num.loci,1)
    heteroz = 1- rowSums(heteroz.mat^2)
  } 
  
  ### mean and sd of heteroz vector for each year of simulation time
  
  for (hh in 1:i) {
    
    all_freq_year = heter_mat_year[,hh]
    dim(all_freq_year) =  c(num.loci,num.alleles)
    heter_vect_mean_year[hh] = mean((1- rowSums(all_freq_year^2)), na.rm =T)
    heter_vect_sd_year[hh] = sd((1- rowSums(all_freq_year^2)), na.rm =T)
  } 
  
  ##################### CLOSE FOR TIME LOOP ########################
  
  if (plotorno >0) {
    
    plot(vettad,type="b",pch=16,xlab="Year",ylab="Population size",box=NULL,
         ylim=c(0,S+50),cex.lab=1.5,bty="n",cex.axis=1.6)
    
    if (extinct==1) {
      text(250,450,"extinct",cex = 1.5)
    }
    
    points((fitness.mean * 100),type="l",pch=1, col="blue",lwd=1.3,cex=0.7)
    
    points(catastr.vett[which(catastr.vett>0)]*500 ~ which(catastr.vett>0),pch="|",cex=1.5)
    
    plot(vettopt,type="b",pch=16,xlab="Year",ylab="Optimum",box=NULL,
         ylim=c(-4,8),cex.lab=1.5,bty="n",cex.axis=1.6,col="gray45", cex=0.3)
    
    meanott = c(rep(0,first.phase.opt),meanopt*(1:second.phase.opt))
    
    
    points(phenomean,type="b",cex=0.5)
    
    points(meanott,type="l",lty=2,lwd=3,col="red")
    
    
    
    #plot(mediaphen,type="b",pch=16,xlab="Year",ylab="meanPhen",box=NUL#L,ylim=c(-0.5,0.5))
    
    #plot(sdphen,type="b",pch=16,xlab="Year",ylab="sdPhen",box=NULL,yli#m=c(-0.5,0.5))
  }
  
  varianceopt = 0
  meanpopend <- mean(addgenmean[(iter-5):iter])
  
  meanphenoend <- mean(phenomean[(iter-5):iter])
  
  
  #cat("iteration = ", 50, "\n")
  
  ### Pearson's correlation
  cor.all.est = as.numeric(cor.test(all.freq100[1:(num.loci*num.alleles),1],
                                    all.freq100[1:(num.loci*num.alleles),ncol(all.freq100)])$estimate)
  
  ### p-value of Pearson's correlation
  cor.all.p = as.numeric(cor.test(all.freq100[1:(num.loci*num.alleles),1],
                                  all.freq100[1:(num.loci*num.alleles),ncol(all.freq100)])$p.value)
  
  n.mut = nrow(all.freq100) - (num.loci*num.alleles)
  
  
  ### change all.freq100 to data.frame #####
  
  all.freq100 = as.data.frame(all.freq100)
  steps = c(1,seq(0,iter,50)[-1])
  #steps[length(steps)] = min(i,steps[length(steps)])  ### if the population goes extinct
  colnames(all.freq100)[1] = "Allelic_value"
  colnames(all.freq100)[2:(length(steps)+1)] = paste("Freq_time",steps,sep ="_")
  colnames(all.freq100)[ncol(all.freq100)] = "Freq_diff"
  all.freq100$Gene = rep(0,nrow(all.freq100))
  
  for (j in 1:num.loci) {
    all.freq100$Gene[seq(j,(num.loci*num.alleles),num.loci)] = j
    
    
  }
  ### Gene 0 is the mutant alleles
  
  all.freq100 = arrange(all.freq100, Gene)
  
  ris.list = list("extinct"=extinct, 
                  #"addvar" = addgenvar, 
                  #"phenvar" = varphen,
                  "area.pop" = area.pop,
                  "yearextinct" = i-1,
                  "phenomean" = phenomean,
                  "heteroz" = heteroz,
                  "sdopt" = sdopt, 
                  "meanopt"= meanopt,
                  "selstrength"=selection,
                  "numloci"=num.loci,
                  "num.alleles"=num.alleles,
                  "sd.alleles" = sd.alleles, 
                  "recomb"= recomb, 
                  "major" = major, 
                  "mumut" = mu.mut,
                  "mutalfa" = mut.alfa, 
                  "fitness.mean" = fitness.mean, 
                  #"fitness.variance" = fitness.variance,
                  "optimum" = optimum, 
                  "catastr.vett" = catastr.vett,
                  "popsize.post"=vettad,
                  "cat.freq" = p.post, 
                  "age.sex.mat" = age.sex.mat,
                  "avg.surv" = avg.surv, 
                  "sons.mean" = sons.mean,
                  "cor.all.est" = cor.all.est, 
                  "cor.all.p" = cor.all.p,
                  "catastr.mort" = catastr.mort, 
                  #"n.mut" = n.mut,
                  "all.freq" = all.freq100,
                  #"heteroz.mat" = heteroz.mat,
                  #"heter.mat" = heter_mat_year,
                  "heter.mean.year" = heter_vect_mean_year,
                  "heter.sd.year" = heter_vect_sd_year)
  
  return(ris.list)
  
  dev.off()
  
  #########other potential items to be returned######
  # "check.parents" = check.parents  		at year 50 matrix with alleles for parents from row 1 to row num.loci*2 and remaining rows for alleles of offspring
  
  # "herit" = heritability_paroff  heritability at each year with cov(p,o)/var(o)
  # "num.sons.check" = num.sons.check number of offspring recruiting in the population each year
  #"herit.coef" = herit.coef  heritability estimated with the linear model (regression parent-offspring, same as "herit")
  
  
}# close all

#### The function returns a list of results that can be accessed with $
# extinct = 1 if a population wen extinct, 0 otherwise
# addvar = vector of additive genetic variance, one value for each year
# phenvar = vector of phenotypic variance, one value for each year 
# yearextinct = year of exinction if the population went extinct, iter-1 otherwise
# phenomean = vector of mean phenotype, one value for each year
# heteroz = vector of mean heterozigosity, one value for each year
# sdopt = input sdopt
# meanopt = input meanopt
# selstrength = input selection,
# numloci = input num.loci
# num.alleles = input num.alleles,
# sd.alleles = input sd.alleles
# recomb = input recomb 
# major = input major
# mumut = input mu.mut
# mutalfa = input mut.alfa, 
# fitness.mean = vector of mean fitness, one value for each year 
# fitness.variance = vector of fitness variance, one value for each year 
# optimum = vector of realized climate variable values, one value for each year
# catastr.vett = vector of point extremes, 1 if the point extreme occurred, 0 otherwise
# popsize.post = population size (after mortality, both normal and due to point extreme), one value for eacy year
# cat.freq = input p.post 
# age.sex.mat = input age.sex.mat
# avg.surv = input avg.surv
# sons.mean = input sons.mean
# cor.all.est = correlation (Pearson's r) between allele frequency of alleles at time = 1 and allele frequency 
##    at the end of simulation time
# cor.all.p" = p-value of the correlation
# n.mut = number of mutant alleles in the population at the end of simulation time


