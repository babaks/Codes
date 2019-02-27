###############################################################################################
#  This R code is to calculate the posterior probability of gene being differentially
#  expressed based on the Bayesian mixture model. The Bayesian ODP is also computed.
#
#  Please acknowledge your use of BODP in publications by referencing:
#  Cao J, Xie X, Zhang S, Whitehurst A, and White MA. Bayesian optimal discovery procedure    
#  for simultaneous significance testing. BMC Bioinformatics, 2009. 
#
#  This code is free for academic use only. It is distributed in the hope that it will 
#  be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#   
#  Copyright (C) 2008 Jing Cao (jcao@smu.edu)
#  Department of Statistical Science, Southern Methodist University, Dallas, TX
###############################################################################################

 
###############################################################################################
#
# Notations:
# BODP: function which calculates the posterior probability and the Bayesian ODP
# x: expression measurements under control
# y: expression measurements under treatment
# n0: number of burn-in iterations in MCMC
# ntotal: number of total iterations in MCMC
#
###############################################################################################



 BODP = function(x,y, n0, ntotal)
 { 
  ngene = nrow(x)
  nx = ncol(x)
  ny = ncol(y)
  nxy = nx + ny 

##  centering the gene expressions
  for (i in 1:ngene) 
  { 
    temp= mean(x[i,])
    x[i,] = x[i,]-temp
    y[i,] = y[i,]-temp
   }

## set the starting values

  ux = uy = var = rep(NA,ngene)

  for (i in 1:ngene) 
  {
    ux[i] = mean(x[i,])
    uy[i] = mean(y[i,])
    var[i] = (var(x[i,])*(nx-1)+var(y[i,])*(ny-1))/(nxy-2)
    if(var[i] < 0.0000001) var[i]=0.0000001
   }

   var0 = mean(var)

   d = rd = rv = BODP = rep(0,ngene) 

## set the hyperparameters of the priors

   pd = pv = 1

   ag = 2.0; bg = var0; ab = 1.0; bb = 1.0; 
 
   su=10; sdiff=10       ## may need adjustment based on the expression level of genes

   pi2 = sqrt(2*pi); invsu2=1/su^2; invsd2=1/sdiff^2

   sumux = sumvar = sumd = sumrd = rep(0,ngene)

## start MCMC

   for (j in 1:ntotal)
   { print(j)

##  update mean expression level under control     
    for (i in 1:ngene)
    {
      temp1 = (sum(x[i,]) + sum(y[i,]-d[i]))/var[i] 
      temp2 = 1/((nxy)/var[i]+invsu2)
      ux[i] = rnorm(1)*sqrt(temp2)+temp1*temp2
     }

##  update indicator of whether gene is differentially expressed
    rd = rep(0,ngene)
    temp1 = 1 - pd
    for (i in 1:ngene)
    {
      temp =  ny/var[i]+invsd2
      temp2 = pd/sdiff/sqrt(temp)
      temp2 = temp2*exp(0.5*(sum(y[i,]-ux[i])/var[i])^2/temp)
      prob = temp2/(temp1+temp2)
      if (temp2 > 100000000*temp1)  prob = 1            # to avoid overflow in calculation
      if(runif(1) < prob)   rd[i]=1
    }

##  update the difference in the expression under two conditions    
    d = rep(0,ngene)
    for (i in 1:ngene)
    { 
     if (rd[i] > 0)
     {
       temp1 = sum(y[i,]-ux[i])/var[i]  
       temp2 = 1/(ny/var[i]+invsd2)
       d[i] = rnorm(1)*sqrt(temp2)+temp1*temp2
      }
     }


##  update indicator of whether gene shares a common variance
    rv = rep(0,ngene)
    temp3 = pv/(pi2)^nxy*bg^ag*gamma(nxy/2+ag)/gamma(ag)
    for (i in 1:ngene)
    {
      temp4 = (sum((x[i,]-ux[i])^2)+sum((y[i,]-ux[i]-d[i])^2))/2

      temp1 = (1-pv)/(pi2*sqrt(var0))^nxy*exp(-temp4/var0)
      temp2 = temp3/(temp4+bg)^(nxy/2+ag)
      prob = temp2/(temp1+temp2)
      if (temp2 > 100000000*temp1)  prob = 1
      if(runif(1) < prob)   rv[i]=1
    }


##  update gene-specific variance
    var = rep(var0,ngene)
    for (i in 1:ngene)
    { 
     if (rv[i] > 0)
     {
       shape = ag + nxy/2
       rate = bg + (sum((x[i,]-ux[i])^2)+sum((y[i,]-ux[i]-d[i])^2))/2
       var[i] = rate/rgamma(1,shape=shape)
      }
     }

##  update the common variance
     shape = ag + nxy/2*sum(1-rv)
     rate = 0
     for (i in 1:ngene) if(rv[i]==0) rate = rate + sum((x[i,]-ux[i])^2)+sum((y[i,]-ux[i]-d[i])^2)
     rate = bg + rate/2
     var0 = rate/rgamma(1,shape=shape)

##  update prior parameters
     temp = sum(rd)
     a = ab + temp
     b = bb + ngene - temp
     pd = rbeta(1,a,b)

     temp = sum(rv)
     a = ab + temp
     b = bb + ngene - temp
     pv = rbeta(1,a,b)

     if(j > n0)
     {
       sumux = sumux + ux
       sumd = sumd + d
       sumvar = sumvar + var
       sumrd = sumrd + rd
      }

   }

## end of MCMC
   

## calculate posterior means
   nuse = ntotal-n0
   ux = sumux/nuse
   d = sumd/nuse
   var = sumvar/nuse
   rd = sumrd/nuse
  

## calculate Bayesian ODP   
   for(j in 1:ngene)
   {print(j)
    temp0 = temp1 = 0
    for (i in 1:ngene)
    {
    temp0 = temp0 + (1-rd[i])*prod(dnorm(y[j,],ux[i],sqrt(var[i])))
    temp1 = temp1 + prod(dnorm(y[j,],(ux[i]+d[i]),sqrt(var[i])))
    }
    BODP[j] = temp1/temp0
   }

## return with output of posterior probability and Bayesian ODP
   out = list(rd = rd, BODP = BODP)
   return(out)

  }

   
    

 








