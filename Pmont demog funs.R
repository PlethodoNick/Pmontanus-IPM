
library(doBy)
library(tidyverse)


### Define the demographic functions ###############

#####  Growth PDF, VB growth, gaussian errors
# num = replicate size (first argument of apply)
# y = size
# y1 = size in year t+1
# growth.par = data.frame of growth parameter posterior distributions (9,000)
# enviro.par = list of environmental variables
# atemp.num = active temp number (Low, Mid, High - site-specific quantiles)
# aprecip.num = active precip number (Low, Mid, High - site-specific quantiles)
# elev.num = elevation number (SPG, IMG, HG, BBT, CG)
# gp.num = growth parameter number - stored outside function, random sample of pd

growth.fnc <- function(y, y1, num, growth.par, enviro.par, 
                       atemp.num, aprecip.num, 
                       gp.num=gp.num){
  
  logk.active <- growth.par[['mean.k.active']][gp.num[num]] + 
    growth.par[['beta.k.active.temp']][gp.num[num]]*
      enviro.par[['active.temp.scaled']][atemp.num] +
    growth.par[['beta.k.active.precip']][gp.num[num]]*
      enviro.par[['active.precip.scaled']][aprecip.num] + 
    growth.par[['beta.k.active.elev']][gp.num[num]]*enviro.par[['elev.scaled']] +
    growth.par[['beta.k.active.elev.temp']][gp.num[num]]*enviro.par[['elev.scaled']]*
      enviro.par[['active.temp.scaled']][atemp.num]
  logk.inactive <- growth.par[['mean.k.inactive']][gp.num[num]]
  
  k.active <- exp(logk.active)
  k.inactive <- exp(logk.inactive)
  
  k <- k.active*(138/365) - k.inactive*(227/365)
  
  mu.size.t1 <- y + (growth.par[['a']][gp.num[num]]-y)*(1-exp(-k))
  
  p.den.grow <- dnorm(y1, mean=mu.size.t1, sd=growth.par[['growth.sd']][gp.num[num]])
  return(p.den.grow)
}


##### Survival function, logit link 
# num = replicate size (first argument of apply)
# y = size (scaled inside function)
# surv.par = data.frame of survival parameter posterior distributions (18,000)
# enviro.par = list of  environmental variables 
# sv.num = survival parameter number - stored outside function, random sample of pd
# atemp.num = active temp number (Low, Mid, High - site-specific quantiles)
# aprecip.num = active precip number (Low, Mid, High - site-specific quantiles)
# itemp.num = inactive temp number (Low, Mid, High - site-specific quantiles)
# iswe.num = inactive precip number (Low, Mid, High - site-specific quantiles)
# elev.num = elevation number (SPG, IMG, HG, BBT, CG)


surv.fnc <- function(y, num, mean.size=mean.size, sd.size=sd.size, 
                     surv.par, enviro.par, sv.num=sv.num,
                     atemp.num, aprecip.num,
                     itemp.num, iswe.num,
                     small.phi, min.phi){
  
  y.scale <- (y-mean.size)/sd.size
  
  mu.phi.a <- surv.par[['mean.phi']][sv.num[num]] +
    surv.par[['beta.size.phi']][sv.num[num]]*y.scale +
    surv.par[['beta.elev.phi']][sv.num[num]]*enviro.par[['elev.scaled']] +
    surv.par[['beta.temp.active.phi']][sv.num[num]]*
    enviro.par[['active.temp.scaled']][atemp.num] +
    surv.par[['beta.precip.active.phi']][sv.num[num]]*
    enviro.par[['active.precip.scaled']][aprecip.num]


  mu.phi.i <- surv.par[['mean.phi']][sv.num[num]] +
    surv.par[['beta.size.phi']][sv.num[num]]*y.scale +
    surv.par[['beta.elev.phi']][sv.num[num]]*enviro.par[['elev.scaled']] +
    surv.par[['beta.season.phi']][sv.num[num]] +
    surv.par[['beta.elev.season.phi']][sv.num[num]]*enviro.par[['elev.scaled']] +
    surv.par[['beta.temp.inactive.phi']][sv.num[num]]*
      enviro.par[['inactive.temp.scaled']][itemp.num] +
  surv.par[['beta.elev.temp.inactive.phi']][sv.num[num]]*
      enviro.par[['inactive.temp.scaled']][itemp.num]*enviro.par[['elev.scaled']] +
    surv.par[['beta.swe.inactive.phi']][sv.num[num]]*
      enviro.par[['inactive.swe.scaled']][iswe.num] +
    surv.par[['beta.elev.swe.inactive.phi']][sv.num[num]]*
      enviro.par[['inactive.swe.scaled']][iswe.num]*enviro.par[['elev.scaled']] +
    surv.par[['beta.temp.swe.inactive.phi']][sv.num[num]]*
      enviro.par[['inactive.temp.scaled']][itemp.num]*
      enviro.par[['inactive.swe.scaled']][iswe.num] +
    surv.par[['beta.elev.temp.swe.inactive.phi']][sv.num[num]]*
      enviro.par[['inactive.temp.scaled']][itemp.num]*
      enviro.par[['inactive.swe.scaled']][iswe.num]*enviro.par[['elev.scaled']]

  phi.a <- 1/(1+exp(-mu.phi.a))
  phi.i <- 1/(1+exp(-mu.phi.i))
  phi <- (phi.a*138 + phi.i*227)/365
  phi2 <- ifelse(phi<small.phi, min.phi, phi)
  return(phi2)
}


## Reproduction prob, logit link
# num = replicate size (first argument of apply)
# y = size
# repro.prob.par = data frame of reproduction probability parameters pd (10000)
# rp.num = reproduction probability number - stored outside function, random sample of pd

repro.prob.fnc <- function(y, num, repro.prob.par, rp.num, 
                           mean.size.repro, sd.size.repro, bio4.scale){
  y.scale <- (y-mean.size.repro)/sd.size.repro
  mu.repro.prob <- repro.prob.par[['int']][rp.num[num]] +
    repro.prob.par[['beta.size']][rp.num[num]]*y.scale +
    repro.prob.par[['beta.bio4']][rp.num[num]]*bio4.scale
  
  repro <- 1/(1+exp(-mu.repro.prob))
  return(repro)
}



## Fecundity, log link
# mean.eggs = log link of mean eggs

fecund.fnc <- function(num, mean.eggs, y){
  N <- exp(mean.eggs)    # Mean number of eggs per female
  return(rep(N, times=length(y)))
}



## Recruit size PDF
# y1 = size in year t+1
# recruit.size.par = list of elevation specific means  (SPG, IMG, HG70, BBT70, CG)
# SD is 1 

recruit.size.fnc <- function(num, y1, recruit.size.par)
{
  p.deRecr <- dnorm(y1, mean=recruit.size.par[['mean.svl']], 
                    sd=recruit.size.par[['sd.svl']]) # pdf of a size z1 recruit
  return(p.deRecr)
}


######### FUNCTIONS TO BUILD IPM KERNALS P, F, and K ############

### SURVIVAL AND GROWTH KERNAL

p.fnc <- function(y, y1, num, growth.par, enviro.par, 
                  atemp.num, aprecip.num, 
                  gp.num=gp.num,
                  mean.size=mean.size, sd.size=sd.size, 
                  surv.par, sv.num=sv.num,
                  itemp.num, iswe.num,
                  small.phi, min.phi){
  return(surv.fnc(y, num, mean.size=mean.size, sd.size=sd.size, 
                  surv.par, enviro.par, sv.num=sv.num,
                  atemp.num, aprecip.num,
                  itemp.num, iswe.num,
                  small.phi, min.phi)*
           growth.fnc(y, y1, num, growth.par, enviro.par, 
                      atemp.num, aprecip.num, 
                      gp.num=gp.num))
}


## FECUNDITY KERNAL
f.fnc <- function(y, y1, num, repro.prob.par, rp.num, 
                  mean.size.repro, sd.size.repro, bio4.scale,
                  mean.eggs,
                  recruit.size.par){
  repro <- repro.prob.fnc(y, num, repro.prob.par, rp.num, 
                          mean.size.repro, sd.size.repro, bio4.scale)*
           fecund.fnc(num, mean.eggs, y)*
           recruit.size.fnc(num, y1, recruit.size.par)*(1/2)
  return(repro)
}

## BUILD THE DISCRETIZED KERNAL
kernal.fnc <- function(m, L, U,
                 num, mean.size, sd.size, 
                 surv.par, enviro.par, sv.num,
                 atemp.num, aprecip.num, itemp.num, 
                 iswe.num,  
                 growth.par, gp.num,
                 repro.prob.par, rp.num, 
                 mean.size.repro, sd.size.repro, bio4.scale,
                 mean.eggs, 
                 recruit.size.par,
                 small.phi, min.phi) {
  # mesh points 
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2)*h
  surv.growth <- h*(outer(meshpts, meshpts, p.fnc, growth.par=growth.par, enviro.par=enviro.par,
                          gp.num=gp.num, num=num,
                          surv.par=surv.par, sv.num=sv.num, mean.size=mean.size, sd.size=sd.size,
                          itemp.num=itemp.num, iswe.num=iswe.num, aprecip.num=aprecip.num, atemp.num=atemp.num,
                          small.phi, min.phi))
  fecundity <- h *(outer(meshpts, meshpts, f.fnc, repro.prob.par=repro.prob.par,
                         rp.num=rp.num, mean.size.repro=mean.size.repro, 
                         sd.size.repro=sd.size.repro, bio4.scale=bio4.scale,
                         mean.eggs=mean.eggs, num=num,
                         recruit.size.par=recruit.size.par))
  kernal <- surv.growth + fecundity
  return(list(kernal=kernal, meshpts=meshpts, 
              surv.growth=surv.growth, fecundity=fecundity))
}



####### POPULATION GROWTH RATE EIGENS, UNNORMALIZED SIZED AND REPRO DISTS ##########

eigens.fnc <- function(num, ipm){
  return(eigen(ipm[[num]]$kernal))
}

lambda.fnc <- function(num, eig.sys){
  return(Re(eig.sys[[num]]$values[1]))
}

wz.fnc <- function(num, eig.sys){
  return(Re(eig.sys[[num]]$vectors[,1]))
}

vz1.fnc <- function(num, ipm){
  return(Re(eigen(t(ipm[[num]]$kernal))$vectors[,1]))
}


################# KERNAL SENSITIVITY AND ELASTICITY #################

sens.fnc <- function(num, v.z1, w.z, h){
  return(outer(v.z1[[num]], w.z[[num]], "*")/
           sum(v.z1[[num]]*w.z[[num]]*h))
}

elas.fnc <- function(num, sens, ipm, h, lambda){
  return(sens[[num]]*(ipm[[num]]$kernal/h)/lambda[[num]])
}


surv.elas.fnc <- function(num, ipm, h, sens, lambda){
  return((ipm[[num]]$surv.growth/h)*sens[[num]]/lambda[[num]])
}



fecund.elas.fnc <- function(num, ipm, h, sens, lambda){
  return((ipm[[num]]$fecundity/h)*sens[[num]]/lambda[[num]])
}
