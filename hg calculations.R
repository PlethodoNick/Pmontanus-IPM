

source('Pmont demog funs_8june2018.R')

load('hg ipm data.Rdata')


###### COMMON MODEL PARAMS ##############

set.seed(8675309)

L <- 0
U <- 90
m <- 100

h <- (U - L)/m
meshpts <- L + ((1:m) - 1/2)*h

reps <- 1000 # number of bootstrap replicates
gp.num <- sample(1:nrow(hg.ipm.dat[['growth.par']]), size=reps, replace=TRUE)
sv.num <- sample(1:nrow(hg.ipm.dat[['surv.par']]), size=reps, replace=TRUE)
rp.num <- sample(1:nrow(hg.ipm.dat[['repro.par']]), size=reps, replace=TRUE)

###### ACTIVE ###########

################## 1-1-2-2 ###############

system.time(hg.ipm1122 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=1, aprecip.num=1, itemp.num=2, iswe.num=2,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num,
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens1122 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm1122))

system.time(hg.lambda1122 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens1122))

system.time(hg.wz1122 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens1122))

system.time(hg.vz1122 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm1122))


## SENSITIVITY
system.time(hg.sens1122 <- lapply(1:reps, sens.fnc, v.z1=hg.vz1122, w.z=hg.wz1122, h=h))

# hg.sum.sens1122 <- Reduce('+', hg.sens1122)
# hg.mean.sens1122 <- hg.sum.sens1122/reps


## ELASTICITY
hg.elas1122 <- lapply(1:reps, elas.fnc, sens=hg.sens1122,
                      ipm=hg.ipm1122, lambda=hg.lambda1122, h=h)

hg.sum.elas1122 <- Reduce('+', hg.elas1122)
hg.mean.elas1122 <- hg.sum.elas1122/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas1122 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm1122, h=h, 
                                       sens=hg.sens1122, lambda=hg.lambda1122))

hg.sum.surv.elas1122 <- Reduce('+', hg.surv.elas1122)
hg.mean.surv.elas1122 <- hg.sum.surv.elas1122/reps

sum(hg.mean.surv.elas1122)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas1122 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm1122, h=h, 
                                         sens=hg.sens1122, lambda=hg.lambda1122))

hg.sum.fecund.elas1122 <- Reduce('+', hg.fecund.elas1122)
hg.mean.fecund.elas1122 <- hg.sum.fecund.elas1122/reps

sum(hg.mean.fecund.elas1122)*h^2


hg1122.list <- list()
hg1122.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg1122.list[['ipm']] <- hg.ipm1122
hg1122.list[['eigens']] <- hg.eigens1122
hg1122.list[['lambda']] <- hg.lambda1122
hg1122.list[['wz']] <- hg.wz1122
hg1122.list[['vz']] <- hg.vz1122
hg1122.list[['sens']] <- hg.sens1122
hg1122.list[['elas']] <- hg.elas1122
hg1122.list[['mean.elas']] <- hg.mean.elas1122
hg1122.list[['sg.elas']] <- hg.surv.elas1122
hg1122.list[['mean.sg.elas']] <- hg.mean.surv.elas1122
hg1122.list[['fec.elas']] <- hg.fecund.elas1122
hg1122.list[['mean.fec.elas']] <- hg.mean.fecund.elas1122

save(hg1122.list, file='Model Results/hg1122.Rdata')


################## 1-2-2-2 ###############

system.time(hg.ipm1222 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=1, aprecip.num=2, itemp.num=2, iswe.num=2,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens1222 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm1222))

system.time(hg.lambda1222 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens1222))

system.time(hg.wz1222 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens1222))

system.time(hg.vz1222 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm1222))


## SENSITIVITY
system.time(hg.sens1222 <- lapply(1:reps, sens.fnc, v.z1=hg.vz1222, w.z=hg.wz1222, h=h))

# hg.sum.sens1222 <- Reduce('+', hg.sens1222)
# hg.mean.sens1222 <- hg.sum.sens1222/reps


## ELASTICITY
hg.elas1222 <- lapply(1:reps, elas.fnc, sens=hg.sens1222,
                      ipm=hg.ipm1222, lambda=hg.lambda1222, h=h)

hg.sum.elas1222 <- Reduce('+', hg.elas1222)
hg.mean.elas1222 <- hg.sum.elas1222/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas1222 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm1222, h=h, 
                                       sens=hg.sens1222, lambda=hg.lambda1222))

hg.sum.surv.elas1222 <- Reduce('+', hg.surv.elas1222)
hg.mean.surv.elas1222 <- hg.sum.surv.elas1222/reps

sum(hg.mean.surv.elas1222)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas1222 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm1222, h=h, 
                                         sens=hg.sens1222, lambda=hg.lambda1222))

hg.sum.fecund.elas1222 <- Reduce('+', hg.fecund.elas1222)
hg.mean.fecund.elas1222 <- hg.sum.fecund.elas1222/reps

sum(hg.mean.fecund.elas1222)*h^2


hg1222.list <- list()
hg1222.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg1222.list[['ipm']] <- hg.ipm1222
hg1222.list[['eigens']] <- hg.eigens1222
hg1222.list[['lambda']] <- hg.lambda1222
hg1222.list[['wz']] <- hg.wz1222
hg1222.list[['vz']] <- hg.vz1222
hg1222.list[['sens']] <- hg.sens1222
hg1222.list[['elas']] <- hg.elas1222
hg1222.list[['mean.elas']] <- hg.mean.elas1222
hg1222.list[['sg.elas']] <- hg.surv.elas1222
hg1222.list[['mean.sg.elas']] <- hg.mean.surv.elas1222
hg1222.list[['fec.elas']] <- hg.fecund.elas1222
hg1222.list[['mean.fec.elas']] <- hg.mean.fecund.elas1222

save(hg1222.list, file='Model Results/hg1222.Rdata')


################## 2-1-2-2 ###############

system.time(hg.ipm2122 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=1, itemp.num=2, iswe.num=2,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens2122 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm2122))

system.time(hg.lambda2122 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens2122))

system.time(hg.wz2122 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens2122))

system.time(hg.vz2122 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm2122))


## SENSITIVITY
system.time(hg.sens2122 <- lapply(1:reps, sens.fnc, v.z1=hg.vz2122, w.z=hg.wz2122, h=h))

# hg.sum.sens2122 <- Reduce('+', hg.sens2122)
# hg.mean.sens2122 <- hg.sum.sens2122/reps


## ELASTICITY
hg.elas2122 <- lapply(1:reps, elas.fnc, sens=hg.sens2122,
                      ipm=hg.ipm2122, lambda=hg.lambda2122, h=h)

hg.sum.elas2122 <- Reduce('+', hg.elas2122)
hg.mean.elas2122 <- hg.sum.elas2122/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas2122 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm2122, h=h, 
                                       sens=hg.sens2122, lambda=hg.lambda2122))

hg.sum.surv.elas2122 <- Reduce('+', hg.surv.elas2122)
hg.mean.surv.elas2122 <- hg.sum.surv.elas2122/reps

sum(hg.mean.surv.elas2122)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas2122 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm2122, h=h, 
                                         sens=hg.sens2122, lambda=hg.lambda2122))

hg.sum.fecund.elas2122 <- Reduce('+', hg.fecund.elas2122)
hg.mean.fecund.elas2122 <- hg.sum.fecund.elas2122/reps

sum(hg.mean.fecund.elas2122)*h^2


hg2122.list <- list()
hg2122.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg2122.list[['ipm']] <- hg.ipm2122
hg2122.list[['eigens']] <- hg.eigens2122
hg2122.list[['lambda']] <- hg.lambda2122
hg2122.list[['wz']] <- hg.wz2122
hg2122.list[['vz']] <- hg.vz2122
hg2122.list[['sens']] <- hg.sens2122
hg2122.list[['elas']] <- hg.elas2122
hg2122.list[['mean.elas']] <- hg.mean.elas2122
hg2122.list[['sg.elas']] <- hg.surv.elas2122
hg2122.list[['mean.sg.elas']] <- hg.mean.surv.elas2122
hg2122.list[['fec.elas']] <- hg.fecund.elas2122
hg2122.list[['mean.fec.elas']] <- hg.mean.fecund.elas2122

save(hg2122.list, file='Model Results/hg2122.Rdata')



################## 2-2-2-2 ###############

system.time(hg.ipm2222 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=2, iswe.num=2,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens2222 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm2222))

system.time(hg.lambda2222 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens2222))

system.time(hg.wz2222 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens2222))

system.time(hg.vz2222 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm2222))


## SENSITIVITY
system.time(hg.sens2222 <- lapply(1:reps, sens.fnc, v.z1=hg.vz2222, w.z=hg.wz2222, h=h))

# hg.sum.sens2222 <- Reduce('+', hg.sens2222)
# hg.mean.sens2222 <- hg.sum.sens2222/reps


## ELASTICITY
hg.elas2222 <- lapply(1:reps, elas.fnc, sens=hg.sens2222,
                      ipm=hg.ipm2222, lambda=hg.lambda2222, h=h)

hg.sum.elas2222 <- Reduce('+', hg.elas2222)
hg.mean.elas2222 <- hg.sum.elas2222/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas2222 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm2222, h=h, 
                                       sens=hg.sens2222, lambda=hg.lambda2222))

hg.sum.surv.elas2222 <- Reduce('+', hg.surv.elas2222)
hg.mean.surv.elas2222 <- hg.sum.surv.elas2222/reps

sum(hg.mean.surv.elas2222)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas2222 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm2222, h=h, 
                                         sens=hg.sens2222, lambda=hg.lambda2222))

hg.sum.fecund.elas2222 <- Reduce('+', hg.fecund.elas2222)
hg.mean.fecund.elas2222 <- hg.sum.fecund.elas2222/reps

sum(hg.mean.fecund.elas2222)*h^2


hg2222.list <- list()
hg2222.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg2222.list[['ipm']] <- hg.ipm2222
hg2222.list[['eigens']] <- hg.eigens2222
hg2222.list[['lambda']] <- hg.lambda2222
hg2222.list[['wz']] <- hg.wz2222
hg2222.list[['vz']] <- hg.vz2222
hg2222.list[['sens']] <- hg.sens2222
hg2222.list[['elas']] <- hg.elas2222
hg2222.list[['mean.elas']] <- hg.mean.elas2222
hg2222.list[['sg.elas']] <- hg.surv.elas2222
hg2222.list[['mean.sg.elas']] <- hg.mean.surv.elas2222
hg2222.list[['fec.elas']] <- hg.fecund.elas2222
hg2222.list[['mean.fec.elas']] <- hg.mean.fecund.elas2222

save(hg2222.list, file='Model Results/hg2222.Rdata')



################## 1-3-2-2 ###############

system.time(hg.ipm1322 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=1, aprecip.num=3, itemp.num=2, iswe.num=2,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens1322 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm1322))

system.time(hg.lambda1322 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens1322))

system.time(hg.wz1322 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens1322))

system.time(hg.vz1322 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm1322))


## SENSITIVITY
system.time(hg.sens1322 <- lapply(1:reps, sens.fnc, v.z1=hg.vz1322, w.z=hg.wz1322, h=h))

# hg.sum.sens1322 <- Reduce('+', hg.sens1322)
# hg.mean.sens1322 <- hg.sum.sens1322/reps


## ELASTICITY
hg.elas1322 <- lapply(1:reps, elas.fnc, sens=hg.sens1322,
                      ipm=hg.ipm1322, lambda=hg.lambda1322, h=h)

hg.sum.elas1322 <- Reduce('+', hg.elas1322)
hg.mean.elas1322 <- hg.sum.elas1322/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas1322 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm1322, h=h, 
                                       sens=hg.sens1322, lambda=hg.lambda1322))

hg.sum.surv.elas1322 <- Reduce('+', hg.surv.elas1322)
hg.mean.surv.elas1322 <- hg.sum.surv.elas1322/reps

sum(hg.mean.surv.elas1322)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas1322 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm1322, h=h, 
                                         sens=hg.sens1322, lambda=hg.lambda1322))

hg.sum.fecund.elas1322 <- Reduce('+', hg.fecund.elas1322)
hg.mean.fecund.elas1322 <- hg.sum.fecund.elas1322/reps

sum(hg.mean.fecund.elas1322)*h^2


hg1322.list <- list()
hg1322.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg1322.list[['ipm']] <- hg.ipm1322
hg1322.list[['eigens']] <- hg.eigens1322
hg1322.list[['lambda']] <- hg.lambda1322
hg1322.list[['wz']] <- hg.wz1322
hg1322.list[['vz']] <- hg.vz1322
hg1322.list[['sens']] <- hg.sens1322
hg1322.list[['elas']] <- hg.elas1322
hg1322.list[['mean.elas']] <- hg.mean.elas1322
hg1322.list[['sg.elas']] <- hg.surv.elas1322
hg1322.list[['mean.sg.elas']] <- hg.mean.surv.elas1322
hg1322.list[['fec.elas']] <- hg.fecund.elas1322
hg1322.list[['mean.fec.elas']] <- hg.mean.fecund.elas1322

save(hg1322.list, file='Model Results/hg1322.Rdata')



################## 3-1-2-2 ###############

system.time(hg.ipm3122 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=3, aprecip.num=1, itemp.num=2, iswe.num=2,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens3122 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm3122))

system.time(hg.lambda3122 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens3122))

system.time(hg.wz3122 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens3122))

system.time(hg.vz3122 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm3122))


## SENSITIVITY
system.time(hg.sens3122 <- lapply(1:reps, sens.fnc, v.z1=hg.vz3122, w.z=hg.wz3122, h=h))

# hg.sum.sens3122 <- Reduce('+', hg.sens3122)
# hg.mean.sens3122 <- hg.sum.sens3122/reps


## ELASTICITY
hg.elas3122 <- lapply(1:reps, elas.fnc, sens=hg.sens3122,
                      ipm=hg.ipm3122, lambda=hg.lambda3122, h=h)

hg.sum.elas3122 <- Reduce('+', hg.elas3122)
hg.mean.elas3122 <- hg.sum.elas3122/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas3122 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm3122, h=h, 
                                       sens=hg.sens3122, lambda=hg.lambda3122))

hg.sum.surv.elas3122 <- Reduce('+', hg.surv.elas3122)
hg.mean.surv.elas3122 <- hg.sum.surv.elas3122/reps

sum(hg.mean.surv.elas3122)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas3122 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm3122, h=h, 
                                         sens=hg.sens3122, lambda=hg.lambda3122))

hg.sum.fecund.elas3122 <- Reduce('+', hg.fecund.elas3122)
hg.mean.fecund.elas3122 <- hg.sum.fecund.elas3122/reps

sum(hg.mean.fecund.elas3122)*h^2


hg3122.list <- list()
hg3122.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg3122.list[['ipm']] <- hg.ipm3122
hg3122.list[['eigens']] <- hg.eigens3122
hg3122.list[['lambda']] <- hg.lambda3122
hg3122.list[['wz']] <- hg.wz3122
hg3122.list[['vz']] <- hg.vz3122
hg3122.list[['sens']] <- hg.sens3122
hg3122.list[['elas']] <- hg.elas3122
hg3122.list[['mean.elas']] <- hg.mean.elas3122
hg3122.list[['sg.elas']] <- hg.surv.elas3122
hg3122.list[['mean.sg.elas']] <- hg.mean.surv.elas3122
hg3122.list[['fec.elas']] <- hg.fecund.elas3122
hg3122.list[['mean.fec.elas']] <- hg.mean.fecund.elas3122

save(hg3122.list, file='Model Results/hg3122.Rdata')


################## 2-3-2-2 ###############

system.time(hg.ipm2322 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=3, itemp.num=2, iswe.num=2,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens2322 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm2322))

system.time(hg.lambda2322 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens2322))

system.time(hg.wz2322 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens2322))

system.time(hg.vz2322 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm2322))


## SENSITIVITY
system.time(hg.sens2322 <- lapply(1:reps, sens.fnc, v.z1=hg.vz2322, w.z=hg.wz2322, h=h))

# hg.sum.sens2322 <- Reduce('+', hg.sens2322)
# hg.mean.sens2322 <- hg.sum.sens2322/reps


## ELASTICITY
hg.elas2322 <- lapply(1:reps, elas.fnc, sens=hg.sens2322,
                      ipm=hg.ipm2322, lambda=hg.lambda2322, h=h)

hg.sum.elas2322 <- Reduce('+', hg.elas2322)
hg.mean.elas2322 <- hg.sum.elas2322/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas2322 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm2322, h=h, 
                                       sens=hg.sens2322, lambda=hg.lambda2322))

hg.sum.surv.elas2322 <- Reduce('+', hg.surv.elas2322)
hg.mean.surv.elas2322 <- hg.sum.surv.elas2322/reps

sum(hg.mean.surv.elas2322)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas2322 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm2322, h=h, 
                                         sens=hg.sens2322, lambda=hg.lambda2322))

hg.sum.fecund.elas2322 <- Reduce('+', hg.fecund.elas2322)
hg.mean.fecund.elas2322 <- hg.sum.fecund.elas2322/reps

sum(hg.mean.fecund.elas2322)*h^2


hg2322.list <- list()
hg2322.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg2322.list[['ipm']] <- hg.ipm2322
hg2322.list[['eigens']] <- hg.eigens2322
hg2322.list[['lambda']] <- hg.lambda2322
hg2322.list[['wz']] <- hg.wz2322
hg2322.list[['vz']] <- hg.vz2322
hg2322.list[['sens']] <- hg.sens2322
hg2322.list[['elas']] <- hg.elas2322
hg2322.list[['mean.elas']] <- hg.mean.elas2322
hg2322.list[['sg.elas']] <- hg.surv.elas2322
hg2322.list[['mean.sg.elas']] <- hg.mean.surv.elas2322
hg2322.list[['fec.elas']] <- hg.fecund.elas2322
hg2322.list[['mean.fec.elas']] <- hg.mean.fecund.elas2322

save(hg2322.list, file='Model Results/hg2322.Rdata')

################## 3-2-2-2 ###############

system.time(hg.ipm3222 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=3, aprecip.num=2, itemp.num=2, iswe.num=2,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens3222 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm3222))

system.time(hg.lambda3222 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens3222))

system.time(hg.wz3222 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens3222))

system.time(hg.vz3222 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm3222))


## SENSITIVITY
system.time(hg.sens3222 <- lapply(1:reps, sens.fnc, v.z1=hg.vz3222, w.z=hg.wz3222, h=h))

# hg.sum.sens3222 <- Reduce('+', hg.sens3222)
# hg.mean.sens3222 <- hg.sum.sens3222/reps


## ELASTICITY
hg.elas3222 <- lapply(1:reps, elas.fnc, sens=hg.sens3222,
                      ipm=hg.ipm3222, lambda=hg.lambda3222, h=h)

hg.sum.elas3222 <- Reduce('+', hg.elas3222)
hg.mean.elas3222 <- hg.sum.elas3222/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas3222 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm3222, h=h, 
                                       sens=hg.sens3222, lambda=hg.lambda3222))

hg.sum.surv.elas3222 <- Reduce('+', hg.surv.elas3222)
hg.mean.surv.elas3222 <- hg.sum.surv.elas3222/reps

sum(hg.mean.surv.elas3222)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas3222 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm3222, h=h, 
                                         sens=hg.sens3222, lambda=hg.lambda3222))

hg.sum.fecund.elas3222 <- Reduce('+', hg.fecund.elas3222)
hg.mean.fecund.elas3222 <- hg.sum.fecund.elas3222/reps

sum(hg.mean.fecund.elas3222)*h^2


hg3222.list <- list()
hg3222.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg3222.list[['ipm']] <- hg.ipm3222
hg3222.list[['eigens']] <- hg.eigens3222
hg3222.list[['lambda']] <- hg.lambda3222
hg3222.list[['wz']] <- hg.wz3222
hg3222.list[['vz']] <- hg.vz3222
hg3222.list[['sens']] <- hg.sens3222
hg3222.list[['elas']] <- hg.elas3222
hg3222.list[['mean.elas']] <- hg.mean.elas3222
hg3222.list[['sg.elas']] <- hg.surv.elas3222
hg3222.list[['mean.sg.elas']] <- hg.mean.surv.elas3222
hg3222.list[['fec.elas']] <- hg.fecund.elas3222
hg3222.list[['mean.fec.elas']] <- hg.mean.fecund.elas3222

save(hg3222.list, file='Model Results/hg3222.Rdata')


################## 3-3-2-2 ###############

system.time(hg.ipm3322 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=3, aprecip.num=3, itemp.num=2, iswe.num=2,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens3322 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm3322))

system.time(hg.lambda3322 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens3322))

system.time(hg.wz3322 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens3322))

system.time(hg.vz3322 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm3322))


## SENSITIVITY
system.time(hg.sens3322 <- lapply(1:reps, sens.fnc, v.z1=hg.vz3322, w.z=hg.wz3322, h=h))

# hg.sum.sens3322 <- Reduce('+', hg.sens3322)
# hg.mean.sens3322 <- hg.sum.sens3322/reps


## ELASTICITY
hg.elas3322 <- lapply(1:reps, elas.fnc, sens=hg.sens3322,
                      ipm=hg.ipm3322, lambda=hg.lambda3322, h=h)

hg.sum.elas3322 <- Reduce('+', hg.elas3322)
hg.mean.elas3322 <- hg.sum.elas3322/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas3322 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm3322, h=h, 
                                       sens=hg.sens3322, lambda=hg.lambda3322))

hg.sum.surv.elas3322 <- Reduce('+', hg.surv.elas3322)
hg.mean.surv.elas3322 <- hg.sum.surv.elas3322/reps

sum(hg.mean.surv.elas3322)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas3322 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm3322, h=h, 
                                         sens=hg.sens3322, lambda=hg.lambda3322))

hg.sum.fecund.elas3322 <- Reduce('+', hg.fecund.elas3322)
hg.mean.fecund.elas3322 <- hg.sum.fecund.elas3322/reps

sum(hg.mean.fecund.elas3322)*h^2


hg3322.list <- list()
hg3322.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg3322.list[['ipm']] <- hg.ipm3322
hg3322.list[['eigens']] <- hg.eigens3322
hg3322.list[['lambda']] <- hg.lambda3322
hg3322.list[['wz']] <- hg.wz3322
hg3322.list[['vz']] <- hg.vz3322
hg3322.list[['sens']] <- hg.sens3322
hg3322.list[['elas']] <- hg.elas3322
hg3322.list[['mean.elas']] <- hg.mean.elas3322
hg3322.list[['sg.elas']] <- hg.surv.elas3322
hg3322.list[['mean.sg.elas']] <- hg.mean.surv.elas3322
hg3322.list[['fec.elas']] <- hg.fecund.elas3322
hg3322.list[['mean.fec.elas']] <- hg.mean.fecund.elas3322

save(hg3322.list, file='Model Results/hg3322.Rdata')

#### INACTIVE ###############

################## 2-2-1-1 ###############

system.time(hg.ipm2211 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=1, iswe.num=1,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens2211 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm2211))

system.time(hg.lambda2211 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens2211))

system.time(hg.wz2211 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens2211))

system.time(hg.vz2211 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm2211))


## SENSITIVITY
system.time(hg.sens2211 <- lapply(1:reps, sens.fnc, v.z1=hg.vz2211, w.z=hg.wz2211, h=h))

# hg.sum.sens2211 <- Reduce('+', hg.sens2211)
# hg.mean.sens2211 <- hg.sum.sens2211/reps


## ELASTICITY
hg.elas2211 <- lapply(1:reps, elas.fnc, sens=hg.sens2211,
                      ipm=hg.ipm2211, lambda=hg.lambda2211, h=h)

hg.sum.elas2211 <- Reduce('+', hg.elas2211)
hg.mean.elas2211 <- hg.sum.elas2211/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas2211 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm2211, h=h, 
                                       sens=hg.sens2211, lambda=hg.lambda2211))

hg.sum.surv.elas2211 <- Reduce('+', hg.surv.elas2211)
hg.mean.surv.elas2211 <- hg.sum.surv.elas2211/reps

sum(hg.mean.surv.elas2211)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas2211 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm2211, h=h, 
                                         sens=hg.sens2211, lambda=hg.lambda2211))

hg.sum.fecund.elas2211 <- Reduce('+', hg.fecund.elas2211)
hg.mean.fecund.elas2211 <- hg.sum.fecund.elas2211/reps

sum(hg.mean.fecund.elas2211)*h^2


hg2211.list <- list()
hg2211.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg2211.list[['ipm']] <- hg.ipm2211
hg2211.list[['eigens']] <- hg.eigens2211
hg2211.list[['lambda']] <- hg.lambda2211
hg2211.list[['wz']] <- hg.wz2211
hg2211.list[['vz']] <- hg.vz2211
hg2211.list[['sens']] <- hg.sens2211
hg2211.list[['elas']] <- hg.elas2211
hg2211.list[['mean.elas']] <- hg.mean.elas2211
hg2211.list[['sg.elas']] <- hg.surv.elas2211
hg2211.list[['mean.sg.elas']] <- hg.mean.surv.elas2211
hg2211.list[['fec.elas']] <- hg.fecund.elas2211
hg2211.list[['mean.fec.elas']] <- hg.mean.fecund.elas2211

save(hg2211.list, file='Model Results/hg2211.Rdata')


################## 2-2-1-2 ###############


system.time(hg.ipm2212 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=1, iswe.num=2,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens2212 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm2212))

system.time(hg.lambda2212 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens2212))

system.time(hg.wz2212 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens2212))

system.time(hg.vz2212 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm2212))


## SENSITIVITY
system.time(hg.sens2212 <- lapply(1:reps, sens.fnc, v.z1=hg.vz2212, w.z=hg.wz2212, h=h))

# hg.sum.sens2212 <- Reduce('+', hg.sens2212)
# hg.mean.sens2212 <- hg.sum.sens2212/reps


## ELASTICITY
hg.elas2212 <- lapply(1:reps, elas.fnc, sens=hg.sens2212,
                      ipm=hg.ipm2212, lambda=hg.lambda2212, h=h)

hg.sum.elas2212 <- Reduce('+', hg.elas2212)
hg.mean.elas2212 <- hg.sum.elas2212/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas2212 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm2212, h=h, 
                                       sens=hg.sens2212, lambda=hg.lambda2212))

hg.sum.surv.elas2212 <- Reduce('+', hg.surv.elas2212)
hg.mean.surv.elas2212 <- hg.sum.surv.elas2212/reps

sum(hg.mean.surv.elas2212)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas2212 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm2212, h=h, 
                                         sens=hg.sens2212, lambda=hg.lambda2212))

hg.sum.fecund.elas2212 <- Reduce('+', hg.fecund.elas2212)
hg.mean.fecund.elas2212 <- hg.sum.fecund.elas2212/reps

sum(hg.mean.fecund.elas2212)*h^2


hg2212.list <- list()
hg2212.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg2212.list[['ipm']] <- hg.ipm2212
hg2212.list[['eigens']] <- hg.eigens2212
hg2212.list[['lambda']] <- hg.lambda2212
hg2212.list[['wz']] <- hg.wz2212
hg2212.list[['vz']] <- hg.vz2212
hg2212.list[['sens']] <- hg.sens2212
hg2212.list[['elas']] <- hg.elas2212
hg2212.list[['mean.elas']] <- hg.mean.elas2212
hg2212.list[['sg.elas']] <- hg.surv.elas2212
hg2212.list[['mean.sg.elas']] <- hg.mean.surv.elas2212
hg2212.list[['fec.elas']] <- hg.fecund.elas2212
hg2212.list[['mean.fec.elas']] <- hg.mean.fecund.elas2212

save(hg2212.list, file='Model Results/hg2212.Rdata')


################## 2-2-2-1 ###############


system.time(hg.ipm2221 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=2, iswe.num=1,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens2221 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm2221))

system.time(hg.lambda2221 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens2221))

system.time(hg.wz2221 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens2221))

system.time(hg.vz2221 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm2221))


## SENSITIVITY
system.time(hg.sens2221 <- lapply(1:reps, sens.fnc, v.z1=hg.vz2221, w.z=hg.wz2221, h=h))

# hg.sum.sens2221 <- Reduce('+', hg.sens2221)
# hg.mean.sens2221 <- hg.sum.sens2221/reps


## ELASTICITY
hg.elas2221 <- lapply(1:reps, elas.fnc, sens=hg.sens2221,
                      ipm=hg.ipm2221, lambda=hg.lambda2221, h=h)

hg.sum.elas2221 <- Reduce('+', hg.elas2221)
hg.mean.elas2221 <- hg.sum.elas2221/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas2221 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm2221, h=h, 
                                       sens=hg.sens2221, lambda=hg.lambda2221))

hg.sum.surv.elas2221 <- Reduce('+', hg.surv.elas2221)
hg.mean.surv.elas2221 <- hg.sum.surv.elas2221/reps

sum(hg.mean.surv.elas2221)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas2221 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm2221, h=h, 
                                         sens=hg.sens2221, lambda=hg.lambda2221))

hg.sum.fecund.elas2221 <- Reduce('+', hg.fecund.elas2221)
hg.mean.fecund.elas2221 <- hg.sum.fecund.elas2221/reps

sum(hg.mean.fecund.elas2221)*h^2


hg2221.list <- list()
hg2221.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg2221.list[['ipm']] <- hg.ipm2221
hg2221.list[['eigens']] <- hg.eigens2221
hg2221.list[['lambda']] <- hg.lambda2221
hg2221.list[['wz']] <- hg.wz2221
hg2221.list[['vz']] <- hg.vz2221
hg2221.list[['sens']] <- hg.sens2221
hg2221.list[['elas']] <- hg.elas2221
hg2221.list[['mean.elas']] <- hg.mean.elas2221
hg2221.list[['sg.elas']] <- hg.surv.elas2221
hg2221.list[['mean.sg.elas']] <- hg.mean.surv.elas2221
hg2221.list[['fec.elas']] <- hg.fecund.elas2221
hg2221.list[['mean.fec.elas']] <- hg.mean.fecund.elas2221

save(hg2221.list, file='Model Results/hg2221.Rdata')


################## 2-2-1-3 ###############


system.time(hg.ipm2213 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=1, iswe.num=3,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens2213 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm2213))

system.time(hg.lambda2213 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens2213))

system.time(hg.wz2213 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens2213))

system.time(hg.vz2213 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm2213))


## SENSITIVITY
system.time(hg.sens2213 <- lapply(1:reps, sens.fnc, v.z1=hg.vz2213, w.z=hg.wz2213, h=h))

# hg.sum.sens2213 <- Reduce('+', hg.sens2213)
# hg.mean.sens2213 <- hg.sum.sens2213/reps


## ELASTICITY
hg.elas2213 <- lapply(1:reps, elas.fnc, sens=hg.sens2213,
                      ipm=hg.ipm2213, lambda=hg.lambda2213, h=h)

hg.sum.elas2213 <- Reduce('+', hg.elas2213)
hg.mean.elas2213 <- hg.sum.elas2213/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas2213 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm2213, h=h, 
                                       sens=hg.sens2213, lambda=hg.lambda2213))

hg.sum.surv.elas2213 <- Reduce('+', hg.surv.elas2213)
hg.mean.surv.elas2213 <- hg.sum.surv.elas2213/reps

sum(hg.mean.surv.elas2213)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas2213 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm2213, h=h, 
                                         sens=hg.sens2213, lambda=hg.lambda2213))

hg.sum.fecund.elas2213 <- Reduce('+', hg.fecund.elas2213)
hg.mean.fecund.elas2213 <- hg.sum.fecund.elas2213/reps

sum(hg.mean.fecund.elas2213)*h^2


hg2213.list <- list()
hg2213.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg2213.list[['ipm']] <- hg.ipm2213
hg2213.list[['eigens']] <- hg.eigens2213
hg2213.list[['lambda']] <- hg.lambda2213
hg2213.list[['wz']] <- hg.wz2213
hg2213.list[['vz']] <- hg.vz2213
hg2213.list[['sens']] <- hg.sens2213
hg2213.list[['elas']] <- hg.elas2213
hg2213.list[['mean.elas']] <- hg.mean.elas2213
hg2213.list[['sg.elas']] <- hg.surv.elas2213
hg2213.list[['mean.sg.elas']] <- hg.mean.surv.elas2213
hg2213.list[['fec.elas']] <- hg.fecund.elas2213
hg2213.list[['mean.fec.elas']] <- hg.mean.fecund.elas2213

save(hg2213.list, file='Model Results/hg2213.Rdata')


################## 2-2-3-1 ###############


system.time(hg.ipm2231 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=3, iswe.num=1,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens2231 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm2231))

system.time(hg.lambda2231 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens2231))

system.time(hg.wz2231 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens2231))

system.time(hg.vz2231 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm2231))


## SENSITIVITY
system.time(hg.sens2231 <- lapply(1:reps, sens.fnc, v.z1=hg.vz2231, w.z=hg.wz2231, h=h))

# hg.sum.sens2231 <- Reduce('+', hg.sens2231)
# hg.mean.sens2231 <- hg.sum.sens2231/reps


## ELASTICITY
hg.elas2231 <- lapply(1:reps, elas.fnc, sens=hg.sens2231,
                      ipm=hg.ipm2231, lambda=hg.lambda2231, h=h)

hg.sum.elas2231 <- Reduce('+', hg.elas2231)
hg.mean.elas2231 <- hg.sum.elas2231/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas2231 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm2231, h=h, 
                                       sens=hg.sens2231, lambda=hg.lambda2231))

hg.sum.surv.elas2231 <- Reduce('+', hg.surv.elas2231)
hg.mean.surv.elas2231 <- hg.sum.surv.elas2231/reps

sum(hg.mean.surv.elas2231)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas2231 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm2231, h=h, 
                                         sens=hg.sens2231, lambda=hg.lambda2231))

hg.sum.fecund.elas2231 <- Reduce('+', hg.fecund.elas2231)
hg.mean.fecund.elas2231 <- hg.sum.fecund.elas2231/reps

sum(hg.mean.fecund.elas2231)*h^2


hg2231.list <- list()
hg2231.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg2231.list[['ipm']] <- hg.ipm2231
hg2231.list[['eigens']] <- hg.eigens2231
hg2231.list[['lambda']] <- hg.lambda2231
hg2231.list[['wz']] <- hg.wz2231
hg2231.list[['vz']] <- hg.vz2231
hg2231.list[['sens']] <- hg.sens2231
hg2231.list[['elas']] <- hg.elas2231
hg2231.list[['mean.elas']] <- hg.mean.elas2231
hg2231.list[['sg.elas']] <- hg.surv.elas2231
hg2231.list[['mean.sg.elas']] <- hg.mean.surv.elas2231
hg2231.list[['fec.elas']] <- hg.fecund.elas2231
hg2231.list[['mean.fec.elas']] <- hg.mean.fecund.elas2231

save(hg2231.list, file='Model Results/hg2231.Rdata')


################## 2-2-2-3 ###############


system.time(hg.ipm2223 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=2, iswe.num=3,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens2223 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm2223))

system.time(hg.lambda2223 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens2223))

system.time(hg.wz2223 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens2223))

system.time(hg.vz2223 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm2223))


## SENSITIVITY
system.time(hg.sens2223 <- lapply(1:reps, sens.fnc, v.z1=hg.vz2223, w.z=hg.wz2223, h=h))

# hg.sum.sens2223 <- Reduce('+', hg.sens2223)
# hg.mean.sens2223 <- hg.sum.sens2223/reps


## ELASTICITY
hg.elas2223 <- lapply(1:reps, elas.fnc, sens=hg.sens2223,
                      ipm=hg.ipm2223, lambda=hg.lambda2223, h=h)

hg.sum.elas2223 <- Reduce('+', hg.elas2223)
hg.mean.elas2223 <- hg.sum.elas2223/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas2223 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm2223, h=h, 
                                       sens=hg.sens2223, lambda=hg.lambda2223))

hg.sum.surv.elas2223 <- Reduce('+', hg.surv.elas2223)
hg.mean.surv.elas2223 <- hg.sum.surv.elas2223/reps

sum(hg.mean.surv.elas2223)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas2223 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm2223, h=h, 
                                         sens=hg.sens2223, lambda=hg.lambda2223))

hg.sum.fecund.elas2223 <- Reduce('+', hg.fecund.elas2223)
hg.mean.fecund.elas2223 <- hg.sum.fecund.elas2223/reps

sum(hg.mean.fecund.elas2223)*h^2


hg2223.list <- list()
hg2223.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg2223.list[['ipm']] <- hg.ipm2223
hg2223.list[['eigens']] <- hg.eigens2223
hg2223.list[['lambda']] <- hg.lambda2223
hg2223.list[['wz']] <- hg.wz2223
hg2223.list[['vz']] <- hg.vz2223
hg2223.list[['sens']] <- hg.sens2223
hg2223.list[['elas']] <- hg.elas2223
hg2223.list[['mean.elas']] <- hg.mean.elas2223
hg2223.list[['sg.elas']] <- hg.surv.elas2223
hg2223.list[['mean.sg.elas']] <- hg.mean.surv.elas2223
hg2223.list[['fec.elas']] <- hg.fecund.elas2223
hg2223.list[['mean.fec.elas']] <- hg.mean.fecund.elas2223

save(hg2223.list, file='Model Results/hg2223.Rdata')



################## 2-2-3-2 ###############


system.time(hg.ipm2232 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=3, iswe.num=2,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens2232 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm2232))

system.time(hg.lambda2232 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens2232))

system.time(hg.wz2232 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens2232))

system.time(hg.vz2232 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm2232))


## SENSITIVITY
system.time(hg.sens2232 <- lapply(1:reps, sens.fnc, v.z1=hg.vz2232, w.z=hg.wz2232, h=h))

# hg.sum.sens2232 <- Reduce('+', hg.sens2232)
# hg.mean.sens2232 <- hg.sum.sens2232/reps


## ELASTICITY
hg.elas2232 <- lapply(1:reps, elas.fnc, sens=hg.sens2232,
                      ipm=hg.ipm2232, lambda=hg.lambda2232, h=h)

hg.sum.elas2232 <- Reduce('+', hg.elas2232)
hg.mean.elas2232 <- hg.sum.elas2232/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas2232 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm2232, h=h, 
                                       sens=hg.sens2232, lambda=hg.lambda2232))

hg.sum.surv.elas2232 <- Reduce('+', hg.surv.elas2232)
hg.mean.surv.elas2232 <- hg.sum.surv.elas2232/reps

sum(hg.mean.surv.elas2232)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas2232 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm2232, h=h, 
                                         sens=hg.sens2232, lambda=hg.lambda2232))

hg.sum.fecund.elas2232 <- Reduce('+', hg.fecund.elas2232)
hg.mean.fecund.elas2232 <- hg.sum.fecund.elas2232/reps

sum(hg.mean.fecund.elas2232)*h^2


hg2232.list <- list()
hg2232.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg2232.list[['ipm']] <- hg.ipm2232
hg2232.list[['eigens']] <- hg.eigens2232
hg2232.list[['lambda']] <- hg.lambda2232
hg2232.list[['wz']] <- hg.wz2232
hg2232.list[['vz']] <- hg.vz2232
hg2232.list[['sens']] <- hg.sens2232
hg2232.list[['elas']] <- hg.elas2232
hg2232.list[['mean.elas']] <- hg.mean.elas2232
hg2232.list[['sg.elas']] <- hg.surv.elas2232
hg2232.list[['mean.sg.elas']] <- hg.mean.surv.elas2232
hg2232.list[['fec.elas']] <- hg.fecund.elas2232
hg2232.list[['mean.fec.elas']] <- hg.mean.fecund.elas2232

save(hg2232.list, file='Model Results/hg2232.Rdata')




################## 2-2-3-3 ###############


system.time(hg.ipm2233 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=hg.ipm.dat[['growth.par']], 
                                 enviro.par=hg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=hg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=hg.ipm.dat[['mean.size']], 
                                 sd.size=hg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=3, iswe.num=3,
                                 repro.prob.par=hg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=hg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=hg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=hg.ipm.dat[['bio4.scale']],
                                 mean.eggs=hg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=hg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(hg.eigens2233 <- lapply(1:reps, eigens.fnc, ipm=hg.ipm2233))

system.time(hg.lambda2233 <- lapply(1:reps, lambda.fnc, eig.sys=hg.eigens2233))

system.time(hg.wz2233 <- lapply(1:reps, wz.fnc, eig.sys=hg.eigens2233))

system.time(hg.vz2233 <- lapply(1:reps, vz1.fnc, ipm=hg.ipm2233))


## SENSITIVITY
system.time(hg.sens2233 <- lapply(1:reps, sens.fnc, v.z1=hg.vz2233, w.z=hg.wz2233, h=h))

# hg.sum.sens2233 <- Reduce('+', hg.sens2233)
# hg.mean.sens2233 <- hg.sum.sens2233/reps


## ELASTICITY
hg.elas2233 <- lapply(1:reps, elas.fnc, sens=hg.sens2233,
                      ipm=hg.ipm2233, lambda=hg.lambda2233, h=h)

hg.sum.elas2233 <- Reduce('+', hg.elas2233)
hg.mean.elas2233 <- hg.sum.elas2233/reps


## SURV/GROWTH ELASTICITY
system.time(hg.surv.elas2233 <- lapply(1:reps, surv.elas.fnc, ipm=hg.ipm2233, h=h, 
                                       sens=hg.sens2233, lambda=hg.lambda2233))

hg.sum.surv.elas2233 <- Reduce('+', hg.surv.elas2233)
hg.mean.surv.elas2233 <- hg.sum.surv.elas2233/reps

sum(hg.mean.surv.elas2233)*h^2


## FECUNDITY ELASTICITY
system.time(hg.fecund.elas2233 <- lapply(1:reps, fecund.elas.fnc, ipm=hg.ipm2233, h=h, 
                                         sens=hg.sens2233, lambda=hg.lambda2233))

hg.sum.fecund.elas2233 <- Reduce('+', hg.fecund.elas2233)
hg.mean.fecund.elas2233 <- hg.sum.fecund.elas2233/reps

sum(hg.mean.fecund.elas2233)*h^2


hg2233.list <- list()
hg2233.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
hg2233.list[['ipm']] <- hg.ipm2233
hg2233.list[['eigens']] <- hg.eigens2233
hg2233.list[['lambda']] <- hg.lambda2233
hg2233.list[['wz']] <- hg.wz2233
hg2233.list[['vz']] <- hg.vz2233
hg2233.list[['sens']] <- hg.sens2233
hg2233.list[['elas']] <- hg.elas2233
hg2233.list[['mean.elas']] <- hg.mean.elas2233
hg2233.list[['sg.elas']] <- hg.surv.elas2233
hg2233.list[['mean.sg.elas']] <- hg.mean.surv.elas2233
hg2233.list[['fec.elas']] <- hg.fecund.elas2233
hg2233.list[['mean.fec.elas']] <- hg.mean.fecund.elas2233

save(hg2233.list, file='Model Results/hg2233.Rdata')
