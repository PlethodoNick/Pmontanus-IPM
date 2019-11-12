

source('Pmont demog funs_8june2018.R')

load('cg ipm data.Rdata')


###### COMMON MODEL PARAMS ##############

set.seed(8675309)

L <- 0
U <- 90
m <- 100

h <- (U - L)/m
meshpts <- L + ((1:m) - 1/2)*h

reps <- 1000 # number of bootstrap replicates
gp.num <- sample(1:nrow(cg.ipm.dat[['growth.par']]), size=reps, replace=TRUE)
sv.num <- sample(1:nrow(cg.ipm.dat[['surv.par']]), size=reps, replace=TRUE)
rp.num <- sample(1:nrow(cg.ipm.dat[['repro.par']]), size=reps, replace=TRUE)


###### ACTIVE ###########

################## 1-1-2-2 ###############

system.time(cg.ipm1122 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=1, aprecip.num=1, itemp.num=2, iswe.num=2,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens1122 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm1122))

system.time(cg.lambda1122 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens1122))

system.time(cg.wz1122 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens1122))

system.time(cg.vz1122 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm1122))


## SENSITIVITY
system.time(cg.sens1122 <- lapply(1:reps, sens.fnc, v.z1=cg.vz1122, w.z=cg.wz1122, h=h))

# cg.sum.sens1122 <- Reduce('+', cg.sens1122)
# cg.mean.sens1122 <- cg.sum.sens1122/reps


## ELASTICITY
cg.elas1122 <- lapply(1:reps, elas.fnc, sens=cg.sens1122,
                      ipm=cg.ipm1122, lambda=cg.lambda1122, h=h)

cg.sum.elas1122 <- Reduce('+', cg.elas1122)
cg.mean.elas1122 <- cg.sum.elas1122/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas1122 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm1122, h=h, 
                                       sens=cg.sens1122, lambda=cg.lambda1122))

cg.sum.surv.elas1122 <- Reduce('+', cg.surv.elas1122)
cg.mean.surv.elas1122 <- cg.sum.surv.elas1122/reps

sum(cg.mean.surv.elas1122)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas1122 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm1122, h=h, 
                                         sens=cg.sens1122, lambda=cg.lambda1122))

cg.sum.fecund.elas1122 <- Reduce('+', cg.fecund.elas1122)
cg.mean.fecund.elas1122 <- cg.sum.fecund.elas1122/reps

sum(cg.mean.fecund.elas1122)*h^2


cg1122.list <- list()
cg1122.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg1122.list[['ipm']] <- cg.ipm1122
cg1122.list[['eigens']] <- cg.eigens1122
cg1122.list[['lambda']] <- cg.lambda1122
cg1122.list[['wz']] <- cg.wz1122
cg1122.list[['vz']] <- cg.vz1122
cg1122.list[['sens']] <- cg.sens1122
cg1122.list[['elas']] <- cg.elas1122
cg1122.list[['mean.elas']] <- cg.mean.elas1122
cg1122.list[['sg.elas']] <- cg.surv.elas1122
cg1122.list[['mean.sg.elas']] <- cg.mean.surv.elas1122
cg1122.list[['fec.elas']] <- cg.fecund.elas1122
cg1122.list[['mean.fec.elas']] <- cg.mean.fecund.elas1122

save(cg1122.list, file='Model Results/cg1122.Rdata')


################## 1-2-2-2 ###############

system.time(cg.ipm1222 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=1, aprecip.num=2, itemp.num=2, iswe.num=2,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens1222 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm1222))

system.time(cg.lambda1222 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens1222))

system.time(cg.wz1222 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens1222))

system.time(cg.vz1222 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm1222))


## SENSITIVITY
system.time(cg.sens1222 <- lapply(1:reps, sens.fnc, v.z1=cg.vz1222, w.z=cg.wz1222, h=h))

# cg.sum.sens1222 <- Reduce('+', cg.sens1222)
# cg.mean.sens1222 <- cg.sum.sens1222/reps


## ELASTICITY
cg.elas1222 <- lapply(1:reps, elas.fnc, sens=cg.sens1222,
                      ipm=cg.ipm1222, lambda=cg.lambda1222, h=h)

cg.sum.elas1222 <- Reduce('+', cg.elas1222)
cg.mean.elas1222 <- cg.sum.elas1222/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas1222 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm1222, h=h, 
                                       sens=cg.sens1222, lambda=cg.lambda1222))

cg.sum.surv.elas1222 <- Reduce('+', cg.surv.elas1222)
cg.mean.surv.elas1222 <- cg.sum.surv.elas1222/reps

sum(cg.mean.surv.elas1222)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas1222 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm1222, h=h, 
                                         sens=cg.sens1222, lambda=cg.lambda1222))

cg.sum.fecund.elas1222 <- Reduce('+', cg.fecund.elas1222)
cg.mean.fecund.elas1222 <- cg.sum.fecund.elas1222/reps

sum(cg.mean.fecund.elas1222)*h^2


cg1222.list <- list()
cg1222.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg1222.list[['ipm']] <- cg.ipm1222
cg1222.list[['eigens']] <- cg.eigens1222
cg1222.list[['lambda']] <- cg.lambda1222
cg1222.list[['wz']] <- cg.wz1222
cg1222.list[['vz']] <- cg.vz1222
cg1222.list[['sens']] <- cg.sens1222
cg1222.list[['elas']] <- cg.elas1222
cg1222.list[['mean.elas']] <- cg.mean.elas1222
cg1222.list[['sg.elas']] <- cg.surv.elas1222
cg1222.list[['mean.sg.elas']] <- cg.mean.surv.elas1222
cg1222.list[['fec.elas']] <- cg.fecund.elas1222
cg1222.list[['mean.fec.elas']] <- cg.mean.fecund.elas1222

save(cg1222.list, file='Model Results/cg1222.Rdata')


################## 2-1-2-2 ###############

system.time(cg.ipm2122 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=1, itemp.num=2, iswe.num=2,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens2122 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm2122))

system.time(cg.lambda2122 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens2122))

system.time(cg.wz2122 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens2122))

system.time(cg.vz2122 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm2122))


## SENSITIVITY
system.time(cg.sens2122 <- lapply(1:reps, sens.fnc, v.z1=cg.vz2122, w.z=cg.wz2122, h=h))

# cg.sum.sens2122 <- Reduce('+', cg.sens2122)
# cg.mean.sens2122 <- cg.sum.sens2122/reps


## ELASTICITY
cg.elas2122 <- lapply(1:reps, elas.fnc, sens=cg.sens2122,
                      ipm=cg.ipm2122, lambda=cg.lambda2122, h=h)

cg.sum.elas2122 <- Reduce('+', cg.elas2122)
cg.mean.elas2122 <- cg.sum.elas2122/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas2122 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm2122, h=h, 
                                       sens=cg.sens2122, lambda=cg.lambda2122))

cg.sum.surv.elas2122 <- Reduce('+', cg.surv.elas2122)
cg.mean.surv.elas2122 <- cg.sum.surv.elas2122/reps

sum(cg.mean.surv.elas2122)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas2122 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm2122, h=h, 
                                         sens=cg.sens2122, lambda=cg.lambda2122))

cg.sum.fecund.elas2122 <- Reduce('+', cg.fecund.elas2122)
cg.mean.fecund.elas2122 <- cg.sum.fecund.elas2122/reps

sum(cg.mean.fecund.elas2122)*h^2


cg2122.list <- list()
cg2122.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg2122.list[['ipm']] <- cg.ipm2122
cg2122.list[['eigens']] <- cg.eigens2122
cg2122.list[['lambda']] <- cg.lambda2122
cg2122.list[['wz']] <- cg.wz2122
cg2122.list[['vz']] <- cg.vz2122
cg2122.list[['sens']] <- cg.sens2122
cg2122.list[['elas']] <- cg.elas2122
cg2122.list[['mean.elas']] <- cg.mean.elas2122
cg2122.list[['sg.elas']] <- cg.surv.elas2122
cg2122.list[['mean.sg.elas']] <- cg.mean.surv.elas2122
cg2122.list[['fec.elas']] <- cg.fecund.elas2122
cg2122.list[['mean.fec.elas']] <- cg.mean.fecund.elas2122

save(cg2122.list, file='Model Results/cg2122.Rdata')



################## 2-2-2-2 ###############

system.time(cg.ipm2222 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=2, iswe.num=2,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens2222 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm2222))

system.time(cg.lambda2222 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens2222))

system.time(cg.wz2222 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens2222))

system.time(cg.vz2222 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm2222))


## SENSITIVITY
system.time(cg.sens2222 <- lapply(1:reps, sens.fnc, v.z1=cg.vz2222, w.z=cg.wz2222, h=h))

# cg.sum.sens2222 <- Reduce('+', cg.sens2222)
# cg.mean.sens2222 <- cg.sum.sens2222/reps


## ELASTICITY
cg.elas2222 <- lapply(1:reps, elas.fnc, sens=cg.sens2222,
                      ipm=cg.ipm2222, lambda=cg.lambda2222, h=h)

cg.sum.elas2222 <- Reduce('+', cg.elas2222)
cg.mean.elas2222 <- cg.sum.elas2222/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas2222 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm2222, h=h, 
                                       sens=cg.sens2222, lambda=cg.lambda2222))

cg.sum.surv.elas2222 <- Reduce('+', cg.surv.elas2222)
cg.mean.surv.elas2222 <- cg.sum.surv.elas2222/reps

sum(cg.mean.surv.elas2222)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas2222 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm2222, h=h, 
                                         sens=cg.sens2222, lambda=cg.lambda2222))

cg.sum.fecund.elas2222 <- Reduce('+', cg.fecund.elas2222)
cg.mean.fecund.elas2222 <- cg.sum.fecund.elas2222/reps

sum(cg.mean.fecund.elas2222)*h^2


cg2222.list <- list()
cg2222.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg2222.list[['ipm']] <- cg.ipm2222
cg2222.list[['eigens']] <- cg.eigens2222
cg2222.list[['lambda']] <- cg.lambda2222
cg2222.list[['wz']] <- cg.wz2222
cg2222.list[['vz']] <- cg.vz2222
cg2222.list[['sens']] <- cg.sens2222
cg2222.list[['elas']] <- cg.elas2222
cg2222.list[['mean.elas']] <- cg.mean.elas2222
cg2222.list[['sg.elas']] <- cg.surv.elas2222
cg2222.list[['mean.sg.elas']] <- cg.mean.surv.elas2222
cg2222.list[['fec.elas']] <- cg.fecund.elas2222
cg2222.list[['mean.fec.elas']] <- cg.mean.fecund.elas2222

save(cg2222.list, file='Model Results/cg2222.Rdata')



################## 1-3-2-2 ###############

system.time(cg.ipm1322 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=1, aprecip.num=3, itemp.num=2, iswe.num=2,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens1322 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm1322))

system.time(cg.lambda1322 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens1322))

system.time(cg.wz1322 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens1322))

system.time(cg.vz1322 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm1322))


## SENSITIVITY
system.time(cg.sens1322 <- lapply(1:reps, sens.fnc, v.z1=cg.vz1322, w.z=cg.wz1322, h=h))

# cg.sum.sens1322 <- Reduce('+', cg.sens1322)
# cg.mean.sens1322 <- cg.sum.sens1322/reps


## ELASTICITY
cg.elas1322 <- lapply(1:reps, elas.fnc, sens=cg.sens1322,
                      ipm=cg.ipm1322, lambda=cg.lambda1322, h=h)

cg.sum.elas1322 <- Reduce('+', cg.elas1322)
cg.mean.elas1322 <- cg.sum.elas1322/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas1322 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm1322, h=h, 
                                       sens=cg.sens1322, lambda=cg.lambda1322))

cg.sum.surv.elas1322 <- Reduce('+', cg.surv.elas1322)
cg.mean.surv.elas1322 <- cg.sum.surv.elas1322/reps

sum(cg.mean.surv.elas1322)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas1322 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm1322, h=h, 
                                         sens=cg.sens1322, lambda=cg.lambda1322))

cg.sum.fecund.elas1322 <- Reduce('+', cg.fecund.elas1322)
cg.mean.fecund.elas1322 <- cg.sum.fecund.elas1322/reps

sum(cg.mean.fecund.elas1322)*h^2


cg1322.list <- list()
cg1322.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg1322.list[['ipm']] <- cg.ipm1322
cg1322.list[['eigens']] <- cg.eigens1322
cg1322.list[['lambda']] <- cg.lambda1322
cg1322.list[['wz']] <- cg.wz1322
cg1322.list[['vz']] <- cg.vz1322
cg1322.list[['sens']] <- cg.sens1322
cg1322.list[['elas']] <- cg.elas1322
cg1322.list[['mean.elas']] <- cg.mean.elas1322
cg1322.list[['sg.elas']] <- cg.surv.elas1322
cg1322.list[['mean.sg.elas']] <- cg.mean.surv.elas1322
cg1322.list[['fec.elas']] <- cg.fecund.elas1322
cg1322.list[['mean.fec.elas']] <- cg.mean.fecund.elas1322

save(cg1322.list, file='Model Results/cg1322.Rdata')


################## 3-1-2-2 ###############

system.time(cg.ipm3122 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=3, aprecip.num=1, itemp.num=2, iswe.num=2,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens3122 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm3122))

system.time(cg.lambda3122 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens3122))

system.time(cg.wz3122 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens3122))

system.time(cg.vz3122 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm3122))


## SENSITIVITY
system.time(cg.sens3122 <- lapply(1:reps, sens.fnc, v.z1=cg.vz3122, w.z=cg.wz3122, h=h))

# cg.sum.sens3122 <- Reduce('+', cg.sens3122)
# cg.mean.sens3122 <- cg.sum.sens3122/reps


## ELASTICITY
cg.elas3122 <- lapply(1:reps, elas.fnc, sens=cg.sens3122,
                      ipm=cg.ipm3122, lambda=cg.lambda3122, h=h)

cg.sum.elas3122 <- Reduce('+', cg.elas3122)
cg.mean.elas3122 <- cg.sum.elas3122/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas3122 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm3122, h=h, 
                                       sens=cg.sens3122, lambda=cg.lambda3122))

cg.sum.surv.elas3122 <- Reduce('+', cg.surv.elas3122)
cg.mean.surv.elas3122 <- cg.sum.surv.elas3122/reps

sum(cg.mean.surv.elas3122)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas3122 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm3122, h=h, 
                                         sens=cg.sens3122, lambda=cg.lambda3122))

cg.sum.fecund.elas3122 <- Reduce('+', cg.fecund.elas3122)
cg.mean.fecund.elas3122 <- cg.sum.fecund.elas3122/reps

sum(cg.mean.fecund.elas3122)*h^2


cg3122.list <- list()
cg3122.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg3122.list[['ipm']] <- cg.ipm3122
cg3122.list[['eigens']] <- cg.eigens3122
cg3122.list[['lambda']] <- cg.lambda3122
cg3122.list[['wz']] <- cg.wz3122
cg3122.list[['vz']] <- cg.vz3122
cg3122.list[['sens']] <- cg.sens3122
cg3122.list[['elas']] <- cg.elas3122
cg3122.list[['mean.elas']] <- cg.mean.elas3122
cg3122.list[['sg.elas']] <- cg.surv.elas3122
cg3122.list[['mean.sg.elas']] <- cg.mean.surv.elas3122
cg3122.list[['fec.elas']] <- cg.fecund.elas3122
cg3122.list[['mean.fec.elas']] <- cg.mean.fecund.elas3122

save(cg3122.list, file='Model Results/cg3122.Rdata')


################## 2-3-2-2 ###############

system.time(cg.ipm2322 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=3, itemp.num=2, iswe.num=2,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens2322 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm2322))

system.time(cg.lambda2322 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens2322))

system.time(cg.wz2322 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens2322))

system.time(cg.vz2322 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm2322))


## SENSITIVITY
system.time(cg.sens2322 <- lapply(1:reps, sens.fnc, v.z1=cg.vz2322, w.z=cg.wz2322, h=h))

# cg.sum.sens2322 <- Reduce('+', cg.sens2322)
# cg.mean.sens2322 <- cg.sum.sens2322/reps


## ELASTICITY
cg.elas2322 <- lapply(1:reps, elas.fnc, sens=cg.sens2322,
                      ipm=cg.ipm2322, lambda=cg.lambda2322, h=h)

cg.sum.elas2322 <- Reduce('+', cg.elas2322)
cg.mean.elas2322 <- cg.sum.elas2322/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas2322 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm2322, h=h, 
                                       sens=cg.sens2322, lambda=cg.lambda2322))

cg.sum.surv.elas2322 <- Reduce('+', cg.surv.elas2322)
cg.mean.surv.elas2322 <- cg.sum.surv.elas2322/reps

sum(cg.mean.surv.elas2322)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas2322 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm2322, h=h, 
                                         sens=cg.sens2322, lambda=cg.lambda2322))

cg.sum.fecund.elas2322 <- Reduce('+', cg.fecund.elas2322)
cg.mean.fecund.elas2322 <- cg.sum.fecund.elas2322/reps

sum(cg.mean.fecund.elas2322)*h^2


cg2322.list <- list()
cg2322.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg2322.list[['ipm']] <- cg.ipm2322
cg2322.list[['eigens']] <- cg.eigens2322
cg2322.list[['lambda']] <- cg.lambda2322
cg2322.list[['wz']] <- cg.wz2322
cg2322.list[['vz']] <- cg.vz2322
cg2322.list[['sens']] <- cg.sens2322
cg2322.list[['elas']] <- cg.elas2322
cg2322.list[['mean.elas']] <- cg.mean.elas2322
cg2322.list[['sg.elas']] <- cg.surv.elas2322
cg2322.list[['mean.sg.elas']] <- cg.mean.surv.elas2322
cg2322.list[['fec.elas']] <- cg.fecund.elas2322
cg2322.list[['mean.fec.elas']] <- cg.mean.fecund.elas2322

save(cg2322.list, file='Model Results/cg2322.Rdata')

################## 3-2-2-2 ###############

system.time(cg.ipm3222 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=3, aprecip.num=2, itemp.num=2, iswe.num=2,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens3222 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm3222))

system.time(cg.lambda3222 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens3222))

system.time(cg.wz3222 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens3222))

system.time(cg.vz3222 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm3222))


## SENSITIVITY
system.time(cg.sens3222 <- lapply(1:reps, sens.fnc, v.z1=cg.vz3222, w.z=cg.wz3222, h=h))

# cg.sum.sens3222 <- Reduce('+', cg.sens3222)
# cg.mean.sens3222 <- cg.sum.sens3222/reps


## ELASTICITY
cg.elas3222 <- lapply(1:reps, elas.fnc, sens=cg.sens3222,
                      ipm=cg.ipm3222, lambda=cg.lambda3222, h=h)

cg.sum.elas3222 <- Reduce('+', cg.elas3222)
cg.mean.elas3222 <- cg.sum.elas3222/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas3222 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm3222, h=h, 
                                       sens=cg.sens3222, lambda=cg.lambda3222))

cg.sum.surv.elas3222 <- Reduce('+', cg.surv.elas3222)
cg.mean.surv.elas3222 <- cg.sum.surv.elas3222/reps

sum(cg.mean.surv.elas3222)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas3222 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm3222, h=h, 
                                         sens=cg.sens3222, lambda=cg.lambda3222))

cg.sum.fecund.elas3222 <- Reduce('+', cg.fecund.elas3222)
cg.mean.fecund.elas3222 <- cg.sum.fecund.elas3222/reps

sum(cg.mean.fecund.elas3222)*h^2


cg3222.list <- list()
cg3222.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg3222.list[['ipm']] <- cg.ipm3222
cg3222.list[['eigens']] <- cg.eigens3222
cg3222.list[['lambda']] <- cg.lambda3222
cg3222.list[['wz']] <- cg.wz3222
cg3222.list[['vz']] <- cg.vz3222
cg3222.list[['sens']] <- cg.sens3222
cg3222.list[['elas']] <- cg.elas3222
cg3222.list[['mean.elas']] <- cg.mean.elas3222
cg3222.list[['sg.elas']] <- cg.surv.elas3222
cg3222.list[['mean.sg.elas']] <- cg.mean.surv.elas3222
cg3222.list[['fec.elas']] <- cg.fecund.elas3222
cg3222.list[['mean.fec.elas']] <- cg.mean.fecund.elas3222

save(cg3222.list, file='Model Results/cg3222.Rdata')


################## 3-3-2-2 ###############

system.time(cg.ipm3322 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=3, aprecip.num=3, itemp.num=2, iswe.num=2,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens3322 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm3322))

system.time(cg.lambda3322 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens3322))

system.time(cg.wz3322 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens3322))

system.time(cg.vz3322 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm3322))


## SENSITIVITY
system.time(cg.sens3322 <- lapply(1:reps, sens.fnc, v.z1=cg.vz3322, w.z=cg.wz3322, h=h))

# cg.sum.sens3322 <- Reduce('+', cg.sens3322)
# cg.mean.sens3322 <- cg.sum.sens3322/reps


## ELASTICITY
cg.elas3322 <- lapply(1:reps, elas.fnc, sens=cg.sens3322,
                      ipm=cg.ipm3322, lambda=cg.lambda3322, h=h)

cg.sum.elas3322 <- Reduce('+', cg.elas3322)
cg.mean.elas3322 <- cg.sum.elas3322/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas3322 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm3322, h=h, 
                                       sens=cg.sens3322, lambda=cg.lambda3322))

cg.sum.surv.elas3322 <- Reduce('+', cg.surv.elas3322)
cg.mean.surv.elas3322 <- cg.sum.surv.elas3322/reps

sum(cg.mean.surv.elas3322)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas3322 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm3322, h=h, 
                                         sens=cg.sens3322, lambda=cg.lambda3322))

cg.sum.fecund.elas3322 <- Reduce('+', cg.fecund.elas3322)
cg.mean.fecund.elas3322 <- cg.sum.fecund.elas3322/reps

sum(cg.mean.fecund.elas3322)*h^2


cg3322.list <- list()
cg3322.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg3322.list[['ipm']] <- cg.ipm3322
cg3322.list[['eigens']] <- cg.eigens3322
cg3322.list[['lambda']] <- cg.lambda3322
cg3322.list[['wz']] <- cg.wz3322
cg3322.list[['vz']] <- cg.vz3322
cg3322.list[['sens']] <- cg.sens3322
cg3322.list[['elas']] <- cg.elas3322
cg3322.list[['mean.elas']] <- cg.mean.elas3322
cg3322.list[['sg.elas']] <- cg.surv.elas3322
cg3322.list[['mean.sg.elas']] <- cg.mean.surv.elas3322
cg3322.list[['fec.elas']] <- cg.fecund.elas3322
cg3322.list[['mean.fec.elas']] <- cg.mean.fecund.elas3322

save(cg3322.list, file='Model Results/cg3322.Rdata')

#### INACTIVE ###############

################## 2-2-1-1 ###############

system.time(cg.ipm2211 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=1, iswe.num=1,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens2211 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm2211))

system.time(cg.lambda2211 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens2211))

system.time(cg.wz2211 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens2211))

system.time(cg.vz2211 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm2211))


## SENSITIVITY
system.time(cg.sens2211 <- lapply(1:reps, sens.fnc, v.z1=cg.vz2211, w.z=cg.wz2211, h=h))

# cg.sum.sens2211 <- Reduce('+', cg.sens2211)
# cg.mean.sens2211 <- cg.sum.sens2211/reps


## ELASTICITY
cg.elas2211 <- lapply(1:reps, elas.fnc, sens=cg.sens2211,
                      ipm=cg.ipm2211, lambda=cg.lambda2211, h=h)

cg.sum.elas2211 <- Reduce('+', cg.elas2211)
cg.mean.elas2211 <- cg.sum.elas2211/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas2211 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm2211, h=h, 
                                       sens=cg.sens2211, lambda=cg.lambda2211))

cg.sum.surv.elas2211 <- Reduce('+', cg.surv.elas2211)
cg.mean.surv.elas2211 <- cg.sum.surv.elas2211/reps

sum(cg.mean.surv.elas2211)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas2211 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm2211, h=h, 
                                         sens=cg.sens2211, lambda=cg.lambda2211))

cg.sum.fecund.elas2211 <- Reduce('+', cg.fecund.elas2211)
cg.mean.fecund.elas2211 <- cg.sum.fecund.elas2211/reps

sum(cg.mean.fecund.elas2211)*h^2


cg2211.list <- list()
cg2211.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg2211.list[['ipm']] <- cg.ipm2211
cg2211.list[['eigens']] <- cg.eigens2211
cg2211.list[['lambda']] <- cg.lambda2211
cg2211.list[['wz']] <- cg.wz2211
cg2211.list[['vz']] <- cg.vz2211
cg2211.list[['sens']] <- cg.sens2211
cg2211.list[['elas']] <- cg.elas2211
cg2211.list[['mean.elas']] <- cg.mean.elas2211
cg2211.list[['sg.elas']] <- cg.surv.elas2211
cg2211.list[['mean.sg.elas']] <- cg.mean.surv.elas2211
cg2211.list[['fec.elas']] <- cg.fecund.elas2211
cg2211.list[['mean.fec.elas']] <- cg.mean.fecund.elas2211

save(cg2211.list, file='Model Results/cg2211.Rdata')


################## 2-2-1-2 ###############


system.time(cg.ipm2212 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=1, iswe.num=2,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens2212 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm2212))

system.time(cg.lambda2212 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens2212))

system.time(cg.wz2212 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens2212))

system.time(cg.vz2212 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm2212))


## SENSITIVITY
system.time(cg.sens2212 <- lapply(1:reps, sens.fnc, v.z1=cg.vz2212, w.z=cg.wz2212, h=h))

# cg.sum.sens2212 <- Reduce('+', cg.sens2212)
# cg.mean.sens2212 <- cg.sum.sens2212/reps


## ELASTICITY
cg.elas2212 <- lapply(1:reps, elas.fnc, sens=cg.sens2212,
                      ipm=cg.ipm2212, lambda=cg.lambda2212, h=h)

cg.sum.elas2212 <- Reduce('+', cg.elas2212)
cg.mean.elas2212 <- cg.sum.elas2212/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas2212 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm2212, h=h, 
                                       sens=cg.sens2212, lambda=cg.lambda2212))

cg.sum.surv.elas2212 <- Reduce('+', cg.surv.elas2212)
cg.mean.surv.elas2212 <- cg.sum.surv.elas2212/reps

sum(cg.mean.surv.elas2212)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas2212 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm2212, h=h, 
                                         sens=cg.sens2212, lambda=cg.lambda2212))

cg.sum.fecund.elas2212 <- Reduce('+', cg.fecund.elas2212)
cg.mean.fecund.elas2212 <- cg.sum.fecund.elas2212/reps

sum(cg.mean.fecund.elas2212)*h^2


cg2212.list <- list()
cg2212.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg2212.list[['ipm']] <- cg.ipm2212
cg2212.list[['eigens']] <- cg.eigens2212
cg2212.list[['lambda']] <- cg.lambda2212
cg2212.list[['wz']] <- cg.wz2212
cg2212.list[['vz']] <- cg.vz2212
cg2212.list[['sens']] <- cg.sens2212
cg2212.list[['elas']] <- cg.elas2212
cg2212.list[['mean.elas']] <- cg.mean.elas2212
cg2212.list[['sg.elas']] <- cg.surv.elas2212
cg2212.list[['mean.sg.elas']] <- cg.mean.surv.elas2212
cg2212.list[['fec.elas']] <- cg.fecund.elas2212
cg2212.list[['mean.fec.elas']] <- cg.mean.fecund.elas2212

save(cg2212.list, file='Model Results/cg2212.Rdata')


################## 2-2-2-1 ###############


system.time(cg.ipm2221 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=2, iswe.num=1,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens2221 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm2221))

system.time(cg.lambda2221 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens2221))

system.time(cg.wz2221 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens2221))

system.time(cg.vz2221 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm2221))


## SENSITIVITY
system.time(cg.sens2221 <- lapply(1:reps, sens.fnc, v.z1=cg.vz2221, w.z=cg.wz2221, h=h))

# cg.sum.sens2221 <- Reduce('+', cg.sens2221)
# cg.mean.sens2221 <- cg.sum.sens2221/reps


## ELASTICITY
cg.elas2221 <- lapply(1:reps, elas.fnc, sens=cg.sens2221,
                      ipm=cg.ipm2221, lambda=cg.lambda2221, h=h)

cg.sum.elas2221 <- Reduce('+', cg.elas2221)
cg.mean.elas2221 <- cg.sum.elas2221/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas2221 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm2221, h=h, 
                                       sens=cg.sens2221, lambda=cg.lambda2221))

cg.sum.surv.elas2221 <- Reduce('+', cg.surv.elas2221)
cg.mean.surv.elas2221 <- cg.sum.surv.elas2221/reps

sum(cg.mean.surv.elas2221)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas2221 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm2221, h=h, 
                                         sens=cg.sens2221, lambda=cg.lambda2221))

cg.sum.fecund.elas2221 <- Reduce('+', cg.fecund.elas2221)
cg.mean.fecund.elas2221 <- cg.sum.fecund.elas2221/reps

sum(cg.mean.fecund.elas2221)*h^2


cg2221.list <- list()
cg2221.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg2221.list[['ipm']] <- cg.ipm2221
cg2221.list[['eigens']] <- cg.eigens2221
cg2221.list[['lambda']] <- cg.lambda2221
cg2221.list[['wz']] <- cg.wz2221
cg2221.list[['vz']] <- cg.vz2221
cg2221.list[['sens']] <- cg.sens2221
cg2221.list[['elas']] <- cg.elas2221
cg2221.list[['mean.elas']] <- cg.mean.elas2221
cg2221.list[['sg.elas']] <- cg.surv.elas2221
cg2221.list[['mean.sg.elas']] <- cg.mean.surv.elas2221
cg2221.list[['fec.elas']] <- cg.fecund.elas2221
cg2221.list[['mean.fec.elas']] <- cg.mean.fecund.elas2221

save(cg2221.list, file='Model Results/cg2221.Rdata')


################## 2-2-1-3 ###############


system.time(cg.ipm2213 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=1, iswe.num=3,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens2213 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm2213))

system.time(cg.lambda2213 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens2213))

system.time(cg.wz2213 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens2213))

system.time(cg.vz2213 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm2213))


## SENSITIVITY
system.time(cg.sens2213 <- lapply(1:reps, sens.fnc, v.z1=cg.vz2213, w.z=cg.wz2213, h=h))

# cg.sum.sens2213 <- Reduce('+', cg.sens2213)
# cg.mean.sens2213 <- cg.sum.sens2213/reps


## ELASTICITY
cg.elas2213 <- lapply(1:reps, elas.fnc, sens=cg.sens2213,
                      ipm=cg.ipm2213, lambda=cg.lambda2213, h=h)

cg.sum.elas2213 <- Reduce('+', cg.elas2213)
cg.mean.elas2213 <- cg.sum.elas2213/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas2213 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm2213, h=h, 
                                       sens=cg.sens2213, lambda=cg.lambda2213))

cg.sum.surv.elas2213 <- Reduce('+', cg.surv.elas2213)
cg.mean.surv.elas2213 <- cg.sum.surv.elas2213/reps

sum(cg.mean.surv.elas2213)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas2213 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm2213, h=h, 
                                         sens=cg.sens2213, lambda=cg.lambda2213))

cg.sum.fecund.elas2213 <- Reduce('+', cg.fecund.elas2213)
cg.mean.fecund.elas2213 <- cg.sum.fecund.elas2213/reps

sum(cg.mean.fecund.elas2213)*h^2


cg2213.list <- list()
cg2213.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg2213.list[['ipm']] <- cg.ipm2213
cg2213.list[['eigens']] <- cg.eigens2213
cg2213.list[['lambda']] <- cg.lambda2213
cg2213.list[['wz']] <- cg.wz2213
cg2213.list[['vz']] <- cg.vz2213
cg2213.list[['sens']] <- cg.sens2213
cg2213.list[['elas']] <- cg.elas2213
cg2213.list[['mean.elas']] <- cg.mean.elas2213
cg2213.list[['sg.elas']] <- cg.surv.elas2213
cg2213.list[['mean.sg.elas']] <- cg.mean.surv.elas2213
cg2213.list[['fec.elas']] <- cg.fecund.elas2213
cg2213.list[['mean.fec.elas']] <- cg.mean.fecund.elas2213

save(cg2213.list, file='Model Results/cg2213.Rdata')


################## 2-2-3-1 ###############


system.time(cg.ipm2231 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=3, iswe.num=1,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens2231 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm2231))

system.time(cg.lambda2231 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens2231))

system.time(cg.wz2231 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens2231))

system.time(cg.vz2231 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm2231))


## SENSITIVITY
system.time(cg.sens2231 <- lapply(1:reps, sens.fnc, v.z1=cg.vz2231, w.z=cg.wz2231, h=h))

# cg.sum.sens2231 <- Reduce('+', cg.sens2231)
# cg.mean.sens2231 <- cg.sum.sens2231/reps


## ELASTICITY
cg.elas2231 <- lapply(1:reps, elas.fnc, sens=cg.sens2231,
                      ipm=cg.ipm2231, lambda=cg.lambda2231, h=h)

cg.sum.elas2231 <- Reduce('+', cg.elas2231)
cg.mean.elas2231 <- cg.sum.elas2231/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas2231 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm2231, h=h, 
                                       sens=cg.sens2231, lambda=cg.lambda2231))

cg.sum.surv.elas2231 <- Reduce('+', cg.surv.elas2231)
cg.mean.surv.elas2231 <- cg.sum.surv.elas2231/reps

sum(cg.mean.surv.elas2231)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas2231 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm2231, h=h, 
                                         sens=cg.sens2231, lambda=cg.lambda2231))

cg.sum.fecund.elas2231 <- Reduce('+', cg.fecund.elas2231)
cg.mean.fecund.elas2231 <- cg.sum.fecund.elas2231/reps

sum(cg.mean.fecund.elas2231)*h^2


cg2231.list <- list()
cg2231.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg2231.list[['ipm']] <- cg.ipm2231
cg2231.list[['eigens']] <- cg.eigens2231
cg2231.list[['lambda']] <- cg.lambda2231
cg2231.list[['wz']] <- cg.wz2231
cg2231.list[['vz']] <- cg.vz2231
cg2231.list[['sens']] <- cg.sens2231
cg2231.list[['elas']] <- cg.elas2231
cg2231.list[['mean.elas']] <- cg.mean.elas2231
cg2231.list[['sg.elas']] <- cg.surv.elas2231
cg2231.list[['mean.sg.elas']] <- cg.mean.surv.elas2231
cg2231.list[['fec.elas']] <- cg.fecund.elas2231
cg2231.list[['mean.fec.elas']] <- cg.mean.fecund.elas2231

save(cg2231.list, file='/Model Results/cg2231.Rdata')


################## 2-2-2-3 ###############


system.time(cg.ipm2223 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=2, iswe.num=3,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens2223 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm2223))

system.time(cg.lambda2223 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens2223))

system.time(cg.wz2223 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens2223))

system.time(cg.vz2223 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm2223))


## SENSITIVITY
system.time(cg.sens2223 <- lapply(1:reps, sens.fnc, v.z1=cg.vz2223, w.z=cg.wz2223, h=h))

# cg.sum.sens2223 <- Reduce('+', cg.sens2223)
# cg.mean.sens2223 <- cg.sum.sens2223/reps


## ELASTICITY
cg.elas2223 <- lapply(1:reps, elas.fnc, sens=cg.sens2223,
                      ipm=cg.ipm2223, lambda=cg.lambda2223, h=h)

cg.sum.elas2223 <- Reduce('+', cg.elas2223)
cg.mean.elas2223 <- cg.sum.elas2223/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas2223 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm2223, h=h, 
                                       sens=cg.sens2223, lambda=cg.lambda2223))

cg.sum.surv.elas2223 <- Reduce('+', cg.surv.elas2223)
cg.mean.surv.elas2223 <- cg.sum.surv.elas2223/reps

sum(cg.mean.surv.elas2223)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas2223 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm2223, h=h, 
                                         sens=cg.sens2223, lambda=cg.lambda2223))

cg.sum.fecund.elas2223 <- Reduce('+', cg.fecund.elas2223)
cg.mean.fecund.elas2223 <- cg.sum.fecund.elas2223/reps

sum(cg.mean.fecund.elas2223)*h^2


cg2223.list <- list()
cg2223.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg2223.list[['ipm']] <- cg.ipm2223
cg2223.list[['eigens']] <- cg.eigens2223
cg2223.list[['lambda']] <- cg.lambda2223
cg2223.list[['wz']] <- cg.wz2223
cg2223.list[['vz']] <- cg.vz2223
cg2223.list[['sens']] <- cg.sens2223
cg2223.list[['elas']] <- cg.elas2223
cg2223.list[['mean.elas']] <- cg.mean.elas2223
cg2223.list[['sg.elas']] <- cg.surv.elas2223
cg2223.list[['mean.sg.elas']] <- cg.mean.surv.elas2223
cg2223.list[['fec.elas']] <- cg.fecund.elas2223
cg2223.list[['mean.fec.elas']] <- cg.mean.fecund.elas2223

save(cg2223.list, file='Model Results/cg2223.Rdata')



################## 2-2-3-2 ###############


system.time(cg.ipm2232 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=3, iswe.num=2,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens2232 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm2232))

system.time(cg.lambda2232 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens2232))

system.time(cg.wz2232 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens2232))

system.time(cg.vz2232 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm2232))


## SENSITIVITY
system.time(cg.sens2232 <- lapply(1:reps, sens.fnc, v.z1=cg.vz2232, w.z=cg.wz2232, h=h))

# cg.sum.sens2232 <- Reduce('+', cg.sens2232)
# cg.mean.sens2232 <- cg.sum.sens2232/reps


## ELASTICITY
cg.elas2232 <- lapply(1:reps, elas.fnc, sens=cg.sens2232,
                      ipm=cg.ipm2232, lambda=cg.lambda2232, h=h)

cg.sum.elas2232 <- Reduce('+', cg.elas2232)
cg.mean.elas2232 <- cg.sum.elas2232/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas2232 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm2232, h=h, 
                                       sens=cg.sens2232, lambda=cg.lambda2232))

cg.sum.surv.elas2232 <- Reduce('+', cg.surv.elas2232)
cg.mean.surv.elas2232 <- cg.sum.surv.elas2232/reps

sum(cg.mean.surv.elas2232)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas2232 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm2232, h=h, 
                                         sens=cg.sens2232, lambda=cg.lambda2232))

cg.sum.fecund.elas2232 <- Reduce('+', cg.fecund.elas2232)
cg.mean.fecund.elas2232 <- cg.sum.fecund.elas2232/reps

sum(cg.mean.fecund.elas2232)*h^2


cg2232.list <- list()
cg2232.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg2232.list[['ipm']] <- cg.ipm2232
cg2232.list[['eigens']] <- cg.eigens2232
cg2232.list[['lambda']] <- cg.lambda2232
cg2232.list[['wz']] <- cg.wz2232
cg2232.list[['vz']] <- cg.vz2232
cg2232.list[['sens']] <- cg.sens2232
cg2232.list[['elas']] <- cg.elas2232
cg2232.list[['mean.elas']] <- cg.mean.elas2232
cg2232.list[['sg.elas']] <- cg.surv.elas2232
cg2232.list[['mean.sg.elas']] <- cg.mean.surv.elas2232
cg2232.list[['fec.elas']] <- cg.fecund.elas2232
cg2232.list[['mean.fec.elas']] <- cg.mean.fecund.elas2232

save(cg2232.list, file='Model Results/cg2232.Rdata')


################## 2-2-3-3 ###############


system.time(cg.ipm2233 <- lapply(1:reps, FUN=kernal.fnc, 
                                 growth.par=cg.ipm.dat[['growth.par']], 
                                 enviro.par=cg.ipm.dat[['enviro.par']],
                                 gp.num=gp.num,
                                 surv.par=cg.ipm.dat[['surv.par']], 
                                 sv.num=sv.num, 
                                 mean.size=cg.ipm.dat[['mean.size']], 
                                 sd.size=cg.ipm.dat[['sd.size']],
                                 U=U, L=L, m=m, 
                                 atemp.num=2, aprecip.num=2, itemp.num=3, iswe.num=3,
                                 repro.prob.par=cg.ipm.dat[['repro.par']],
                                 rp.num=rp.num, 
                                 mean.size.repro=cg.ipm.dat[['mean.size.repro']], 
                                 sd.size.repro=cg.ipm.dat[['sd.size.repro']], 
                                 bio4.scale=cg.ipm.dat[['bio4.scale']],
                                 mean.eggs=cg.ipm.dat[['mean.eggs']], 
                                 recruit.size.par=cg.ipm.dat[['recruit.size.par']],
                                 small.phi=0, min.phi=0))



system.time(cg.eigens2233 <- lapply(1:reps, eigens.fnc, ipm=cg.ipm2233))

system.time(cg.lambda2233 <- lapply(1:reps, lambda.fnc, eig.sys=cg.eigens2233))

system.time(cg.wz2233 <- lapply(1:reps, wz.fnc, eig.sys=cg.eigens2233))

system.time(cg.vz2233 <- lapply(1:reps, vz1.fnc, ipm=cg.ipm2233))


## SENSITIVITY
system.time(cg.sens2233 <- lapply(1:reps, sens.fnc, v.z1=cg.vz2233, w.z=cg.wz2233, h=h))

# cg.sum.sens2233 <- Reduce('+', cg.sens2233)
# cg.mean.sens2233 <- cg.sum.sens2233/reps


## ELASTICITY
cg.elas2233 <- lapply(1:reps, elas.fnc, sens=cg.sens2233,
                      ipm=cg.ipm2233, lambda=cg.lambda2233, h=h)

cg.sum.elas2233 <- Reduce('+', cg.elas2233)
cg.mean.elas2233 <- cg.sum.elas2233/reps


## SURV/GROWTH ELASTICITY
system.time(cg.surv.elas2233 <- lapply(1:reps, surv.elas.fnc, ipm=cg.ipm2233, h=h, 
                                       sens=cg.sens2233, lambda=cg.lambda2233))

cg.sum.surv.elas2233 <- Reduce('+', cg.surv.elas2233)
cg.mean.surv.elas2233 <- cg.sum.surv.elas2233/reps

sum(cg.mean.surv.elas2233)*h^2


## FECUNDITY ELASTICITY
system.time(cg.fecund.elas2233 <- lapply(1:reps, fecund.elas.fnc, ipm=cg.ipm2233, h=h, 
                                         sens=cg.sens2233, lambda=cg.lambda2233))

cg.sum.fecund.elas2233 <- Reduce('+', cg.fecund.elas2233)
cg.mean.fecund.elas2233 <- cg.sum.fecund.elas2233/reps

sum(cg.mean.fecund.elas2233)*h^2


cg2233.list <- list()
cg2233.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
cg2233.list[['ipm']] <- cg.ipm2233
cg2233.list[['eigens']] <- cg.eigens2233
cg2233.list[['lambda']] <- cg.lambda2233
cg2233.list[['wz']] <- cg.wz2233
cg2233.list[['vz']] <- cg.vz2233
cg2233.list[['sens']] <- cg.sens2233
cg2233.list[['elas']] <- cg.elas2233
cg2233.list[['mean.elas']] <- cg.mean.elas2233
cg2233.list[['sg.elas']] <- cg.surv.elas2233
cg2233.list[['mean.sg.elas']] <- cg.mean.surv.elas2233
cg2233.list[['fec.elas']] <- cg.fecund.elas2233
cg2233.list[['mean.fec.elas']] <- cg.mean.fecund.elas2233

save(cg2233.list, file='Model Results/cg2233.Rdata')
