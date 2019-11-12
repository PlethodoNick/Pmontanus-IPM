

source('Pmont demog funs_8june2018.R')

load('spg ipm data.Rdata')


###### COMMON MODEL PARAMS ##############

set.seed(8675309)

L <- 0
U <- 90
m <- 100

h <- (U - L)/m
meshpts <- L + ((1:m) - 1/2)*h

reps <- 1000 # number of bootstrap replicates
gp.num <- sample(1:nrow(spg.ipm.dat[['growth.par']]), size=reps, replace=TRUE)
sv.num <- sample(1:nrow(spg.ipm.dat[['surv.par']]), size=reps, replace=TRUE)
rp.num <- sample(1:nrow(spg.ipm.dat[['repro.par']]), size=reps, replace=TRUE)

###### ACTIVE ###########

################## 1-1-2-2 ###############

system.time(spg.ipm1122 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=1, aprecip.num=1, itemp.num=2, iswe.num=2,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num,
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens1122 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm1122))

system.time(spg.lambda1122 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens1122))

system.time(spg.wz1122 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens1122))

system.time(spg.vz1122 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm1122))


## SENSITIVITY
system.time(spg.sens1122 <- lapply(1:reps, sens.fnc, v.z1=spg.vz1122, w.z=spg.wz1122, h=h))

# spg.sum.sens1122 <- Reduce('+', spg.sens1122)
# spg.mean.sens1122 <- spg.sum.sens1122/reps


## ELASTICITY
spg.elas1122 <- lapply(1:reps, elas.fnc, sens=spg.sens1122,
                       ipm=spg.ipm1122, lambda=spg.lambda1122, h=h)

spg.sum.elas1122 <- Reduce('+', spg.elas1122)
spg.mean.elas1122 <- spg.sum.elas1122/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas1122 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm1122, h=h, 
                                        sens=spg.sens1122, lambda=spg.lambda1122))

spg.sum.surv.elas1122 <- Reduce('+', spg.surv.elas1122)
spg.mean.surv.elas1122 <- spg.sum.surv.elas1122/reps

sum(spg.mean.surv.elas1122)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas1122 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm1122, h=h, 
                                          sens=spg.sens1122, lambda=spg.lambda1122))

spg.sum.fecund.elas1122 <- Reduce('+', spg.fecund.elas1122)
spg.mean.fecund.elas1122 <- spg.sum.fecund.elas1122/reps

sum(spg.mean.fecund.elas1122)*h^2


spg1122.list <- list()
spg1122.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg1122.list[['ipm']] <- spg.ipm1122
spg1122.list[['eigens']] <- spg.eigens1122
spg1122.list[['lambda']] <- spg.lambda1122
spg1122.list[['wz']] <- spg.wz1122
spg1122.list[['vz']] <- spg.vz1122
spg1122.list[['sens']] <- spg.sens1122
spg1122.list[['elas']] <- spg.elas1122
spg1122.list[['mean.elas']] <- spg.mean.elas1122
spg1122.list[['sg.elas']] <- spg.surv.elas1122
spg1122.list[['mean.sg.elas']] <- spg.mean.surv.elas1122
spg1122.list[['fec.elas']] <- spg.fecund.elas1122
spg1122.list[['mean.fec.elas']] <- spg.mean.fecund.elas1122

save(spg1122.list, file='Model Results/spg1122.Rdata')


################## 1-2-2-2 ###############

system.time(spg.ipm1222 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=1, aprecip.num=2, itemp.num=2, iswe.num=2,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens1222 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm1222))

system.time(spg.lambda1222 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens1222))

system.time(spg.wz1222 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens1222))

system.time(spg.vz1222 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm1222))


## SENSITIVITY
system.time(spg.sens1222 <- lapply(1:reps, sens.fnc, v.z1=spg.vz1222, w.z=spg.wz1222, h=h))

# spg.sum.sens1222 <- Reduce('+', spg.sens1222)
# spg.mean.sens1222 <- spg.sum.sens1222/reps


## ELASTICITY
spg.elas1222 <- lapply(1:reps, elas.fnc, sens=spg.sens1222,
                       ipm=spg.ipm1222, lambda=spg.lambda1222, h=h)

spg.sum.elas1222 <- Reduce('+', spg.elas1222)
spg.mean.elas1222 <- spg.sum.elas1222/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas1222 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm1222, h=h, 
                                        sens=spg.sens1222, lambda=spg.lambda1222))

spg.sum.surv.elas1222 <- Reduce('+', spg.surv.elas1222)
spg.mean.surv.elas1222 <- spg.sum.surv.elas1222/reps

sum(spg.mean.surv.elas1222)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas1222 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm1222, h=h, 
                                          sens=spg.sens1222, lambda=spg.lambda1222))

spg.sum.fecund.elas1222 <- Reduce('+', spg.fecund.elas1222)
spg.mean.fecund.elas1222 <- spg.sum.fecund.elas1222/reps

sum(spg.mean.fecund.elas1222)*h^2


spg1222.list <- list()
spg1222.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg1222.list[['ipm']] <- spg.ipm1222
spg1222.list[['eigens']] <- spg.eigens1222
spg1222.list[['lambda']] <- spg.lambda1222
spg1222.list[['wz']] <- spg.wz1222
spg1222.list[['vz']] <- spg.vz1222
spg1222.list[['sens']] <- spg.sens1222
spg1222.list[['elas']] <- spg.elas1222
spg1222.list[['mean.elas']] <- spg.mean.elas1222
spg1222.list[['sg.elas']] <- spg.surv.elas1222
spg1222.list[['mean.sg.elas']] <- spg.mean.surv.elas1222
spg1222.list[['fec.elas']] <- spg.fecund.elas1222
spg1222.list[['mean.fec.elas']] <- spg.mean.fecund.elas1222

save(spg1222.list, file='Model Results/spg1222.Rdata')


################## 2-1-2-2 ###############

system.time(spg.ipm2122 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=1, itemp.num=2, iswe.num=2,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens2122 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm2122))

system.time(spg.lambda2122 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens2122))

system.time(spg.wz2122 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens2122))

system.time(spg.vz2122 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm2122))


## SENSITIVITY
system.time(spg.sens2122 <- lapply(1:reps, sens.fnc, v.z1=spg.vz2122, w.z=spg.wz2122, h=h))

# spg.sum.sens2122 <- Reduce('+', spg.sens2122)
# spg.mean.sens2122 <- spg.sum.sens2122/reps


## ELASTICITY
spg.elas2122 <- lapply(1:reps, elas.fnc, sens=spg.sens2122,
                       ipm=spg.ipm2122, lambda=spg.lambda2122, h=h)

spg.sum.elas2122 <- Reduce('+', spg.elas2122)
spg.mean.elas2122 <- spg.sum.elas2122/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas2122 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm2122, h=h, 
                                        sens=spg.sens2122, lambda=spg.lambda2122))

spg.sum.surv.elas2122 <- Reduce('+', spg.surv.elas2122)
spg.mean.surv.elas2122 <- spg.sum.surv.elas2122/reps

sum(spg.mean.surv.elas2122)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas2122 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm2122, h=h, 
                                          sens=spg.sens2122, lambda=spg.lambda2122))

spg.sum.fecund.elas2122 <- Reduce('+', spg.fecund.elas2122)
spg.mean.fecund.elas2122 <- spg.sum.fecund.elas2122/reps

sum(spg.mean.fecund.elas2122)*h^2


spg2122.list <- list()
spg2122.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg2122.list[['ipm']] <- spg.ipm2122
spg2122.list[['eigens']] <- spg.eigens2122
spg2122.list[['lambda']] <- spg.lambda2122
spg2122.list[['wz']] <- spg.wz2122
spg2122.list[['vz']] <- spg.vz2122
spg2122.list[['sens']] <- spg.sens2122
spg2122.list[['elas']] <- spg.elas2122
spg2122.list[['mean.elas']] <- spg.mean.elas2122
spg2122.list[['sg.elas']] <- spg.surv.elas2122
spg2122.list[['mean.sg.elas']] <- spg.mean.surv.elas2122
spg2122.list[['fec.elas']] <- spg.fecund.elas2122
spg2122.list[['mean.fec.elas']] <- spg.mean.fecund.elas2122

save(spg2122.list, file='Model Results/spg2122.Rdata')


################## 2-2-2-2 ###############

system.time(spg.ipm2222 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=2, iswe.num=2,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens2222 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm2222))

system.time(spg.lambda2222 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens2222))

system.time(spg.wz2222 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens2222))

system.time(spg.vz2222 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm2222))


## SENSITIVITY
system.time(spg.sens2222 <- lapply(1:reps, sens.fnc, v.z1=spg.vz2222, w.z=spg.wz2222, h=h))

# spg.sum.sens2222 <- Reduce('+', spg.sens2222)
# spg.mean.sens2222 <- spg.sum.sens2222/reps


## ELASTICITY
spg.elas2222 <- lapply(1:reps, elas.fnc, sens=spg.sens2222,
                       ipm=spg.ipm2222, lambda=spg.lambda2222, h=h)

spg.sum.elas2222 <- Reduce('+', spg.elas2222)
spg.mean.elas2222 <- spg.sum.elas2222/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas2222 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm2222, h=h, 
                                        sens=spg.sens2222, lambda=spg.lambda2222))

spg.sum.surv.elas2222 <- Reduce('+', spg.surv.elas2222)
spg.mean.surv.elas2222 <- spg.sum.surv.elas2222/reps

sum(spg.mean.surv.elas2222)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas2222 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm2222, h=h, 
                                          sens=spg.sens2222, lambda=spg.lambda2222))

spg.sum.fecund.elas2222 <- Reduce('+', spg.fecund.elas2222)
spg.mean.fecund.elas2222 <- spg.sum.fecund.elas2222/reps

sum(spg.mean.fecund.elas2222)*h^2


spg2222.list <- list()
spg2222.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg2222.list[['ipm']] <- spg.ipm2222
spg2222.list[['eigens']] <- spg.eigens2222
spg2222.list[['lambda']] <- spg.lambda2222
spg2222.list[['wz']] <- spg.wz2222
spg2222.list[['vz']] <- spg.vz2222
spg2222.list[['sens']] <- spg.sens2222
spg2222.list[['elas']] <- spg.elas2222
spg2222.list[['mean.elas']] <- spg.mean.elas2222
spg2222.list[['sg.elas']] <- spg.surv.elas2222
spg2222.list[['mean.sg.elas']] <- spg.mean.surv.elas2222
spg2222.list[['fec.elas']] <- spg.fecund.elas2222
spg2222.list[['mean.fec.elas']] <- spg.mean.fecund.elas2222

save(spg2222.list, file='Model Results/spg2222.Rdata')


################## 1-3-2-2 ###############

system.time(spg.ipm1322 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=1, aprecip.num=3, itemp.num=2, iswe.num=2,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens1322 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm1322))

system.time(spg.lambda1322 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens1322))

system.time(spg.wz1322 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens1322))

system.time(spg.vz1322 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm1322))


## SENSITIVITY
system.time(spg.sens1322 <- lapply(1:reps, sens.fnc, v.z1=spg.vz1322, w.z=spg.wz1322, h=h))

# spg.sum.sens1322 <- Reduce('+', spg.sens1322)
# spg.mean.sens1322 <- spg.sum.sens1322/reps


## ELASTICITY
spg.elas1322 <- lapply(1:reps, elas.fnc, sens=spg.sens1322,
                       ipm=spg.ipm1322, lambda=spg.lambda1322, h=h)

spg.sum.elas1322 <- Reduce('+', spg.elas1322)
spg.mean.elas1322 <- spg.sum.elas1322/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas1322 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm1322, h=h, 
                                        sens=spg.sens1322, lambda=spg.lambda1322))

spg.sum.surv.elas1322 <- Reduce('+', spg.surv.elas1322)
spg.mean.surv.elas1322 <- spg.sum.surv.elas1322/reps

sum(spg.mean.surv.elas1322)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas1322 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm1322, h=h, 
                                          sens=spg.sens1322, lambda=spg.lambda1322))

spg.sum.fecund.elas1322 <- Reduce('+', spg.fecund.elas1322)
spg.mean.fecund.elas1322 <- spg.sum.fecund.elas1322/reps

sum(spg.mean.fecund.elas1322)*h^2


spg1322.list <- list()
spg1322.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg1322.list[['ipm']] <- spg.ipm1322
spg1322.list[['eigens']] <- spg.eigens1322
spg1322.list[['lambda']] <- spg.lambda1322
spg1322.list[['wz']] <- spg.wz1322
spg1322.list[['vz']] <- spg.vz1322
spg1322.list[['sens']] <- spg.sens1322
spg1322.list[['elas']] <- spg.elas1322
spg1322.list[['mean.elas']] <- spg.mean.elas1322
spg1322.list[['sg.elas']] <- spg.surv.elas1322
spg1322.list[['mean.sg.elas']] <- spg.mean.surv.elas1322
spg1322.list[['fec.elas']] <- spg.fecund.elas1322
spg1322.list[['mean.fec.elas']] <- spg.mean.fecund.elas1322

save(spg1322.list, file='Model Results/spg1322.Rdata')



################## 3-1-2-2 ###############

system.time(spg.ipm3122 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=3, aprecip.num=1, itemp.num=2, iswe.num=2,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens3122 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm3122))

system.time(spg.lambda3122 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens3122))

system.time(spg.wz3122 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens3122))

system.time(spg.vz3122 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm3122))


## SENSITIVITY
system.time(spg.sens3122 <- lapply(1:reps, sens.fnc, v.z1=spg.vz3122, w.z=spg.wz3122, h=h))

# spg.sum.sens3122 <- Reduce('+', spg.sens3122)
# spg.mean.sens3122 <- spg.sum.sens3122/reps


## ELASTICITY
spg.elas3122 <- lapply(1:reps, elas.fnc, sens=spg.sens3122,
                       ipm=spg.ipm3122, lambda=spg.lambda3122, h=h)

spg.sum.elas3122 <- Reduce('+', spg.elas3122)
spg.mean.elas3122 <- spg.sum.elas3122/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas3122 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm3122, h=h, 
                                        sens=spg.sens3122, lambda=spg.lambda3122))

spg.sum.surv.elas3122 <- Reduce('+', spg.surv.elas3122)
spg.mean.surv.elas3122 <- spg.sum.surv.elas3122/reps

sum(spg.mean.surv.elas3122)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas3122 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm3122, h=h, 
                                          sens=spg.sens3122, lambda=spg.lambda3122))

spg.sum.fecund.elas3122 <- Reduce('+', spg.fecund.elas3122)
spg.mean.fecund.elas3122 <- spg.sum.fecund.elas3122/reps

sum(spg.mean.fecund.elas3122)*h^2


spg3122.list <- list()
spg3122.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg3122.list[['ipm']] <- spg.ipm3122
spg3122.list[['eigens']] <- spg.eigens3122
spg3122.list[['lambda']] <- spg.lambda3122
spg3122.list[['wz']] <- spg.wz3122
spg3122.list[['vz']] <- spg.vz3122
spg3122.list[['sens']] <- spg.sens3122
spg3122.list[['elas']] <- spg.elas3122
spg3122.list[['mean.elas']] <- spg.mean.elas3122
spg3122.list[['sg.elas']] <- spg.surv.elas3122
spg3122.list[['mean.sg.elas']] <- spg.mean.surv.elas3122
spg3122.list[['fec.elas']] <- spg.fecund.elas3122
spg3122.list[['mean.fec.elas']] <- spg.mean.fecund.elas3122

save(spg3122.list, file='Model Results/spg3122.Rdata')


################## 2-3-2-2 ###############

system.time(spg.ipm2322 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=3, itemp.num=2, iswe.num=2,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens2322 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm2322))

system.time(spg.lambda2322 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens2322))

system.time(spg.wz2322 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens2322))

system.time(spg.vz2322 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm2322))


## SENSITIVITY
system.time(spg.sens2322 <- lapply(1:reps, sens.fnc, v.z1=spg.vz2322, w.z=spg.wz2322, h=h))

# spg.sum.sens2322 <- Reduce('+', spg.sens2322)
# spg.mean.sens2322 <- spg.sum.sens2322/reps


## ELASTICITY
spg.elas2322 <- lapply(1:reps, elas.fnc, sens=spg.sens2322,
                       ipm=spg.ipm2322, lambda=spg.lambda2322, h=h)

spg.sum.elas2322 <- Reduce('+', spg.elas2322)
spg.mean.elas2322 <- spg.sum.elas2322/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas2322 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm2322, h=h, 
                                        sens=spg.sens2322, lambda=spg.lambda2322))

spg.sum.surv.elas2322 <- Reduce('+', spg.surv.elas2322)
spg.mean.surv.elas2322 <- spg.sum.surv.elas2322/reps

sum(spg.mean.surv.elas2322)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas2322 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm2322, h=h, 
                                          sens=spg.sens2322, lambda=spg.lambda2322))

spg.sum.fecund.elas2322 <- Reduce('+', spg.fecund.elas2322)
spg.mean.fecund.elas2322 <- spg.sum.fecund.elas2322/reps

sum(spg.mean.fecund.elas2322)*h^2


spg2322.list <- list()
spg2322.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg2322.list[['ipm']] <- spg.ipm2322
spg2322.list[['eigens']] <- spg.eigens2322
spg2322.list[['lambda']] <- spg.lambda2322
spg2322.list[['wz']] <- spg.wz2322
spg2322.list[['vz']] <- spg.vz2322
spg2322.list[['sens']] <- spg.sens2322
spg2322.list[['elas']] <- spg.elas2322
spg2322.list[['mean.elas']] <- spg.mean.elas2322
spg2322.list[['sg.elas']] <- spg.surv.elas2322
spg2322.list[['mean.sg.elas']] <- spg.mean.surv.elas2322
spg2322.list[['fec.elas']] <- spg.fecund.elas2322
spg2322.list[['mean.fec.elas']] <- spg.mean.fecund.elas2322

save(spg2322.list, file='Model Results/spg2322.Rdata')

################## 3-2-2-2 ###############

system.time(spg.ipm3222 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=3, aprecip.num=2, itemp.num=2, iswe.num=2,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens3222 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm3222))

system.time(spg.lambda3222 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens3222))

system.time(spg.wz3222 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens3222))

system.time(spg.vz3222 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm3222))


## SENSITIVITY
system.time(spg.sens3222 <- lapply(1:reps, sens.fnc, v.z1=spg.vz3222, w.z=spg.wz3222, h=h))

# spg.sum.sens3222 <- Reduce('+', spg.sens3222)
# spg.mean.sens3222 <- spg.sum.sens3222/reps


## ELASTICITY
spg.elas3222 <- lapply(1:reps, elas.fnc, sens=spg.sens3222,
                       ipm=spg.ipm3222, lambda=spg.lambda3222, h=h)

spg.sum.elas3222 <- Reduce('+', spg.elas3222)
spg.mean.elas3222 <- spg.sum.elas3222/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas3222 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm3222, h=h, 
                                        sens=spg.sens3222, lambda=spg.lambda3222))

spg.sum.surv.elas3222 <- Reduce('+', spg.surv.elas3222)
spg.mean.surv.elas3222 <- spg.sum.surv.elas3222/reps

sum(spg.mean.surv.elas3222)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas3222 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm3222, h=h, 
                                          sens=spg.sens3222, lambda=spg.lambda3222))

spg.sum.fecund.elas3222 <- Reduce('+', spg.fecund.elas3222)
spg.mean.fecund.elas3222 <- spg.sum.fecund.elas3222/reps

sum(spg.mean.fecund.elas3222)*h^2


spg3222.list <- list()
spg3222.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg3222.list[['ipm']] <- spg.ipm3222
spg3222.list[['eigens']] <- spg.eigens3222
spg3222.list[['lambda']] <- spg.lambda3222
spg3222.list[['wz']] <- spg.wz3222
spg3222.list[['vz']] <- spg.vz3222
spg3222.list[['sens']] <- spg.sens3222
spg3222.list[['elas']] <- spg.elas3222
spg3222.list[['mean.elas']] <- spg.mean.elas3222
spg3222.list[['sg.elas']] <- spg.surv.elas3222
spg3222.list[['mean.sg.elas']] <- spg.mean.surv.elas3222
spg3222.list[['fec.elas']] <- spg.fecund.elas3222
spg3222.list[['mean.fec.elas']] <- spg.mean.fecund.elas3222

save(spg3222.list, file='Model Results/spg3222.Rdata')


################## 3-3-2-2 ###############

system.time(spg.ipm3322 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=3, aprecip.num=3, itemp.num=2, iswe.num=2,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens3322 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm3322))

system.time(spg.lambda3322 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens3322))

system.time(spg.wz3322 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens3322))

system.time(spg.vz3322 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm3322))


## SENSITIVITY
system.time(spg.sens3322 <- lapply(1:reps, sens.fnc, v.z1=spg.vz3322, w.z=spg.wz3322, h=h))

# spg.sum.sens3322 <- Reduce('+', spg.sens3322)
# spg.mean.sens3322 <- spg.sum.sens3322/reps


## ELASTICITY
spg.elas3322 <- lapply(1:reps, elas.fnc, sens=spg.sens3322,
                       ipm=spg.ipm3322, lambda=spg.lambda3322, h=h)

spg.sum.elas3322 <- Reduce('+', spg.elas3322)
spg.mean.elas3322 <- spg.sum.elas3322/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas3322 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm3322, h=h, 
                                        sens=spg.sens3322, lambda=spg.lambda3322))

spg.sum.surv.elas3322 <- Reduce('+', spg.surv.elas3322)
spg.mean.surv.elas3322 <- spg.sum.surv.elas3322/reps

sum(spg.mean.surv.elas3322)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas3322 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm3322, h=h, 
                                          sens=spg.sens3322, lambda=spg.lambda3322))

spg.sum.fecund.elas3322 <- Reduce('+', spg.fecund.elas3322)
spg.mean.fecund.elas3322 <- spg.sum.fecund.elas3322/reps

sum(spg.mean.fecund.elas3322)*h^2


spg3322.list <- list()
spg3322.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg3322.list[['ipm']] <- spg.ipm3322
spg3322.list[['eigens']] <- spg.eigens3322
spg3322.list[['lambda']] <- spg.lambda3322
spg3322.list[['wz']] <- spg.wz3322
spg3322.list[['vz']] <- spg.vz3322
spg3322.list[['sens']] <- spg.sens3322
spg3322.list[['elas']] <- spg.elas3322
spg3322.list[['mean.elas']] <- spg.mean.elas3322
spg3322.list[['sg.elas']] <- spg.surv.elas3322
spg3322.list[['mean.sg.elas']] <- spg.mean.surv.elas3322
spg3322.list[['fec.elas']] <- spg.fecund.elas3322
spg3322.list[['mean.fec.elas']] <- spg.mean.fecund.elas3322

save(spg3322.list, file='Model Results/spg3322.Rdata')

#### INACTIVE ###############

################## 2-2-1-1 ###############

system.time(spg.ipm2211 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=1, iswe.num=1,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens2211 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm2211))

system.time(spg.lambda2211 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens2211))

system.time(spg.wz2211 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens2211))

system.time(spg.vz2211 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm2211))


## SENSITIVITY
system.time(spg.sens2211 <- lapply(1:reps, sens.fnc, v.z1=spg.vz2211, w.z=spg.wz2211, h=h))

# spg.sum.sens2211 <- Reduce('+', spg.sens2211)
# spg.mean.sens2211 <- spg.sum.sens2211/reps


## ELASTICITY
spg.elas2211 <- lapply(1:reps, elas.fnc, sens=spg.sens2211,
                       ipm=spg.ipm2211, lambda=spg.lambda2211, h=h)

spg.sum.elas2211 <- Reduce('+', spg.elas2211)
spg.mean.elas2211 <- spg.sum.elas2211/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas2211 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm2211, h=h, 
                                        sens=spg.sens2211, lambda=spg.lambda2211))

spg.sum.surv.elas2211 <- Reduce('+', spg.surv.elas2211)
spg.mean.surv.elas2211 <- spg.sum.surv.elas2211/reps

sum(spg.mean.surv.elas2211)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas2211 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm2211, h=h, 
                                          sens=spg.sens2211, lambda=spg.lambda2211))

spg.sum.fecund.elas2211 <- Reduce('+', spg.fecund.elas2211)
spg.mean.fecund.elas2211 <- spg.sum.fecund.elas2211/reps

sum(spg.mean.fecund.elas2211)*h^2


spg2211.list <- list()
spg2211.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg2211.list[['ipm']] <- spg.ipm2211
spg2211.list[['eigens']] <- spg.eigens2211
spg2211.list[['lambda']] <- spg.lambda2211
spg2211.list[['wz']] <- spg.wz2211
spg2211.list[['vz']] <- spg.vz2211
spg2211.list[['sens']] <- spg.sens2211
spg2211.list[['elas']] <- spg.elas2211
spg2211.list[['mean.elas']] <- spg.mean.elas2211
spg2211.list[['sg.elas']] <- spg.surv.elas2211
spg2211.list[['mean.sg.elas']] <- spg.mean.surv.elas2211
spg2211.list[['fec.elas']] <- spg.fecund.elas2211
spg2211.list[['mean.fec.elas']] <- spg.mean.fecund.elas2211

save(spg2211.list, file='Model Results/spg2211.Rdata')


################## 2-2-1-2 ###############


system.time(spg.ipm2212 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=1, iswe.num=2,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens2212 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm2212))

system.time(spg.lambda2212 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens2212))

system.time(spg.wz2212 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens2212))

system.time(spg.vz2212 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm2212))


## SENSITIVITY
system.time(spg.sens2212 <- lapply(1:reps, sens.fnc, v.z1=spg.vz2212, w.z=spg.wz2212, h=h))

# spg.sum.sens2212 <- Reduce('+', spg.sens2212)
# spg.mean.sens2212 <- spg.sum.sens2212/reps


## ELASTICITY
spg.elas2212 <- lapply(1:reps, elas.fnc, sens=spg.sens2212,
                       ipm=spg.ipm2212, lambda=spg.lambda2212, h=h)

spg.sum.elas2212 <- Reduce('+', spg.elas2212)
spg.mean.elas2212 <- spg.sum.elas2212/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas2212 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm2212, h=h, 
                                        sens=spg.sens2212, lambda=spg.lambda2212))

spg.sum.surv.elas2212 <- Reduce('+', spg.surv.elas2212)
spg.mean.surv.elas2212 <- spg.sum.surv.elas2212/reps

sum(spg.mean.surv.elas2212)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas2212 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm2212, h=h, 
                                          sens=spg.sens2212, lambda=spg.lambda2212))

spg.sum.fecund.elas2212 <- Reduce('+', spg.fecund.elas2212)
spg.mean.fecund.elas2212 <- spg.sum.fecund.elas2212/reps

sum(spg.mean.fecund.elas2212)*h^2


spg2212.list <- list()
spg2212.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg2212.list[['ipm']] <- spg.ipm2212
spg2212.list[['eigens']] <- spg.eigens2212
spg2212.list[['lambda']] <- spg.lambda2212
spg2212.list[['wz']] <- spg.wz2212
spg2212.list[['vz']] <- spg.vz2212
spg2212.list[['sens']] <- spg.sens2212
spg2212.list[['elas']] <- spg.elas2212
spg2212.list[['mean.elas']] <- spg.mean.elas2212
spg2212.list[['sg.elas']] <- spg.surv.elas2212
spg2212.list[['mean.sg.elas']] <- spg.mean.surv.elas2212
spg2212.list[['fec.elas']] <- spg.fecund.elas2212
spg2212.list[['mean.fec.elas']] <- spg.mean.fecund.elas2212

save(spg2212.list, file='Model Results/spg2212.Rdata')


################## 2-2-2-1 ###############


system.time(spg.ipm2221 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=2, iswe.num=1,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens2221 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm2221))

system.time(spg.lambda2221 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens2221))

system.time(spg.wz2221 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens2221))

system.time(spg.vz2221 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm2221))


## SENSITIVITY
system.time(spg.sens2221 <- lapply(1:reps, sens.fnc, v.z1=spg.vz2221, w.z=spg.wz2221, h=h))

# spg.sum.sens2221 <- Reduce('+', spg.sens2221)
# spg.mean.sens2221 <- spg.sum.sens2221/reps


## ELASTICITY
spg.elas2221 <- lapply(1:reps, elas.fnc, sens=spg.sens2221,
                       ipm=spg.ipm2221, lambda=spg.lambda2221, h=h)

spg.sum.elas2221 <- Reduce('+', spg.elas2221)
spg.mean.elas2221 <- spg.sum.elas2221/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas2221 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm2221, h=h, 
                                        sens=spg.sens2221, lambda=spg.lambda2221))

spg.sum.surv.elas2221 <- Reduce('+', spg.surv.elas2221)
spg.mean.surv.elas2221 <- spg.sum.surv.elas2221/reps

sum(spg.mean.surv.elas2221)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas2221 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm2221, h=h, 
                                          sens=spg.sens2221, lambda=spg.lambda2221))

spg.sum.fecund.elas2221 <- Reduce('+', spg.fecund.elas2221)
spg.mean.fecund.elas2221 <- spg.sum.fecund.elas2221/reps

sum(spg.mean.fecund.elas2221)*h^2


spg2221.list <- list()
spg2221.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg2221.list[['ipm']] <- spg.ipm2221
spg2221.list[['eigens']] <- spg.eigens2221
spg2221.list[['lambda']] <- spg.lambda2221
spg2221.list[['wz']] <- spg.wz2221
spg2221.list[['vz']] <- spg.vz2221
spg2221.list[['sens']] <- spg.sens2221
spg2221.list[['elas']] <- spg.elas2221
spg2221.list[['mean.elas']] <- spg.mean.elas2221
spg2221.list[['sg.elas']] <- spg.surv.elas2221
spg2221.list[['mean.sg.elas']] <- spg.mean.surv.elas2221
spg2221.list[['fec.elas']] <- spg.fecund.elas2221
spg2221.list[['mean.fec.elas']] <- spg.mean.fecund.elas2221

save(spg2221.list, file='Model Results/spg2221.Rdata')


################## 2-2-1-3 ###############


system.time(spg.ipm2213 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=1, iswe.num=3,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens2213 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm2213))

system.time(spg.lambda2213 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens2213))

system.time(spg.wz2213 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens2213))

system.time(spg.vz2213 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm2213))


## SENSITIVITY
system.time(spg.sens2213 <- lapply(1:reps, sens.fnc, v.z1=spg.vz2213, w.z=spg.wz2213, h=h))

# spg.sum.sens2213 <- Reduce('+', spg.sens2213)
# spg.mean.sens2213 <- spg.sum.sens2213/reps


## ELASTICITY
spg.elas2213 <- lapply(1:reps, elas.fnc, sens=spg.sens2213,
                       ipm=spg.ipm2213, lambda=spg.lambda2213, h=h)

spg.sum.elas2213 <- Reduce('+', spg.elas2213)
spg.mean.elas2213 <- spg.sum.elas2213/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas2213 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm2213, h=h, 
                                        sens=spg.sens2213, lambda=spg.lambda2213))

spg.sum.surv.elas2213 <- Reduce('+', spg.surv.elas2213)
spg.mean.surv.elas2213 <- spg.sum.surv.elas2213/reps

sum(spg.mean.surv.elas2213)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas2213 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm2213, h=h, 
                                          sens=spg.sens2213, lambda=spg.lambda2213))

spg.sum.fecund.elas2213 <- Reduce('+', spg.fecund.elas2213)
spg.mean.fecund.elas2213 <- spg.sum.fecund.elas2213/reps

sum(spg.mean.fecund.elas2213)*h^2


spg2213.list <- list()
spg2213.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg2213.list[['ipm']] <- spg.ipm2213
spg2213.list[['eigens']] <- spg.eigens2213
spg2213.list[['lambda']] <- spg.lambda2213
spg2213.list[['wz']] <- spg.wz2213
spg2213.list[['vz']] <- spg.vz2213
spg2213.list[['sens']] <- spg.sens2213
spg2213.list[['elas']] <- spg.elas2213
spg2213.list[['mean.elas']] <- spg.mean.elas2213
spg2213.list[['sg.elas']] <- spg.surv.elas2213
spg2213.list[['mean.sg.elas']] <- spg.mean.surv.elas2213
spg2213.list[['fec.elas']] <- spg.fecund.elas2213
spg2213.list[['mean.fec.elas']] <- spg.mean.fecund.elas2213

save(spg2213.list, file='Model Results/spg2213.Rdata')


################## 2-2-3-1 ###############


system.time(spg.ipm2231 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=3, iswe.num=1,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens2231 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm2231))

system.time(spg.lambda2231 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens2231))

system.time(spg.wz2231 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens2231))

system.time(spg.vz2231 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm2231))


## SENSITIVITY
system.time(spg.sens2231 <- lapply(1:reps, sens.fnc, v.z1=spg.vz2231, w.z=spg.wz2231, h=h))

# spg.sum.sens2231 <- Reduce('+', spg.sens2231)
# spg.mean.sens2231 <- spg.sum.sens2231/reps


## ELASTICITY
spg.elas2231 <- lapply(1:reps, elas.fnc, sens=spg.sens2231,
                       ipm=spg.ipm2231, lambda=spg.lambda2231, h=h)

spg.sum.elas2231 <- Reduce('+', spg.elas2231)
spg.mean.elas2231 <- spg.sum.elas2231/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas2231 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm2231, h=h, 
                                        sens=spg.sens2231, lambda=spg.lambda2231))

spg.sum.surv.elas2231 <- Reduce('+', spg.surv.elas2231)
spg.mean.surv.elas2231 <- spg.sum.surv.elas2231/reps

sum(spg.mean.surv.elas2231)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas2231 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm2231, h=h, 
                                          sens=spg.sens2231, lambda=spg.lambda2231))

spg.sum.fecund.elas2231 <- Reduce('+', spg.fecund.elas2231)
spg.mean.fecund.elas2231 <- spg.sum.fecund.elas2231/reps

sum(spg.mean.fecund.elas2231)*h^2


spg2231.list <- list()
spg2231.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg2231.list[['ipm']] <- spg.ipm2231
spg2231.list[['eigens']] <- spg.eigens2231
spg2231.list[['lambda']] <- spg.lambda2231
spg2231.list[['wz']] <- spg.wz2231
spg2231.list[['vz']] <- spg.vz2231
spg2231.list[['sens']] <- spg.sens2231
spg2231.list[['elas']] <- spg.elas2231
spg2231.list[['mean.elas']] <- spg.mean.elas2231
spg2231.list[['sg.elas']] <- spg.surv.elas2231
spg2231.list[['mean.sg.elas']] <- spg.mean.surv.elas2231
spg2231.list[['fec.elas']] <- spg.fecund.elas2231
spg2231.list[['mean.fec.elas']] <- spg.mean.fecund.elas2231

save(spg2231.list, file='Model Results/spg2231.Rdata')


################## 2-2-2-3 ###############


system.time(spg.ipm2223 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=2, iswe.num=3,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens2223 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm2223))

system.time(spg.lambda2223 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens2223))

system.time(spg.wz2223 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens2223))

system.time(spg.vz2223 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm2223))


## SENSITIVITY
system.time(spg.sens2223 <- lapply(1:reps, sens.fnc, v.z1=spg.vz2223, w.z=spg.wz2223, h=h))

# spg.sum.sens2223 <- Reduce('+', spg.sens2223)
# spg.mean.sens2223 <- spg.sum.sens2223/reps


## ELASTICITY
spg.elas2223 <- lapply(1:reps, elas.fnc, sens=spg.sens2223,
                       ipm=spg.ipm2223, lambda=spg.lambda2223, h=h)

spg.sum.elas2223 <- Reduce('+', spg.elas2223)
spg.mean.elas2223 <- spg.sum.elas2223/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas2223 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm2223, h=h, 
                                        sens=spg.sens2223, lambda=spg.lambda2223))

spg.sum.surv.elas2223 <- Reduce('+', spg.surv.elas2223)
spg.mean.surv.elas2223 <- spg.sum.surv.elas2223/reps

sum(spg.mean.surv.elas2223)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas2223 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm2223, h=h, 
                                          sens=spg.sens2223, lambda=spg.lambda2223))

spg.sum.fecund.elas2223 <- Reduce('+', spg.fecund.elas2223)
spg.mean.fecund.elas2223 <- spg.sum.fecund.elas2223/reps

sum(spg.mean.fecund.elas2223)*h^2


spg2223.list <- list()
spg2223.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg2223.list[['ipm']] <- spg.ipm2223
spg2223.list[['eigens']] <- spg.eigens2223
spg2223.list[['lambda']] <- spg.lambda2223
spg2223.list[['wz']] <- spg.wz2223
spg2223.list[['vz']] <- spg.vz2223
spg2223.list[['sens']] <- spg.sens2223
spg2223.list[['elas']] <- spg.elas2223
spg2223.list[['mean.elas']] <- spg.mean.elas2223
spg2223.list[['sg.elas']] <- spg.surv.elas2223
spg2223.list[['mean.sg.elas']] <- spg.mean.surv.elas2223
spg2223.list[['fec.elas']] <- spg.fecund.elas2223
spg2223.list[['mean.fec.elas']] <- spg.mean.fecund.elas2223

save(spg2223.list, file='Model Results/spg2223.Rdata')



################## 2-2-3-2 ###############


system.time(spg.ipm2232 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=3, iswe.num=2,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens2232 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm2232))

system.time(spg.lambda2232 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens2232))

system.time(spg.wz2232 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens2232))

system.time(spg.vz2232 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm2232))


## SENSITIVITY
system.time(spg.sens2232 <- lapply(1:reps, sens.fnc, v.z1=spg.vz2232, w.z=spg.wz2232, h=h))

# spg.sum.sens2232 <- Reduce('+', spg.sens2232)
# spg.mean.sens2232 <- spg.sum.sens2232/reps


## ELASTICITY
spg.elas2232 <- lapply(1:reps, elas.fnc, sens=spg.sens2232,
                       ipm=spg.ipm2232, lambda=spg.lambda2232, h=h)

spg.sum.elas2232 <- Reduce('+', spg.elas2232)
spg.mean.elas2232 <- spg.sum.elas2232/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas2232 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm2232, h=h, 
                                        sens=spg.sens2232, lambda=spg.lambda2232))

spg.sum.surv.elas2232 <- Reduce('+', spg.surv.elas2232)
spg.mean.surv.elas2232 <- spg.sum.surv.elas2232/reps

sum(spg.mean.surv.elas2232)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas2232 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm2232, h=h, 
                                          sens=spg.sens2232, lambda=spg.lambda2232))

spg.sum.fecund.elas2232 <- Reduce('+', spg.fecund.elas2232)
spg.mean.fecund.elas2232 <- spg.sum.fecund.elas2232/reps

sum(spg.mean.fecund.elas2232)*h^2


spg2232.list <- list()
spg2232.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg2232.list[['ipm']] <- spg.ipm2232
spg2232.list[['eigens']] <- spg.eigens2232
spg2232.list[['lambda']] <- spg.lambda2232
spg2232.list[['wz']] <- spg.wz2232
spg2232.list[['vz']] <- spg.vz2232
spg2232.list[['sens']] <- spg.sens2232
spg2232.list[['elas']] <- spg.elas2232
spg2232.list[['mean.elas']] <- spg.mean.elas2232
spg2232.list[['sg.elas']] <- spg.surv.elas2232
spg2232.list[['mean.sg.elas']] <- spg.mean.surv.elas2232
spg2232.list[['fec.elas']] <- spg.fecund.elas2232
spg2232.list[['mean.fec.elas']] <- spg.mean.fecund.elas2232

save(spg2232.list, file='Model Results/spg2232.Rdata')




################## 2-2-3-3 ###############


system.time(spg.ipm2233 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=spg.ipm.dat[['growth.par']], 
                                  enviro.par=spg.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=spg.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=spg.ipm.dat[['mean.size']], 
                                  sd.size=spg.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=3, iswe.num=3,
                                  repro.prob.par=spg.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=spg.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=spg.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=spg.ipm.dat[['bio4.scale']],
                                  mean.eggs=spg.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=spg.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(spg.eigens2233 <- lapply(1:reps, eigens.fnc, ipm=spg.ipm2233))

system.time(spg.lambda2233 <- lapply(1:reps, lambda.fnc, eig.sys=spg.eigens2233))

system.time(spg.wz2233 <- lapply(1:reps, wz.fnc, eig.sys=spg.eigens2233))

system.time(spg.vz2233 <- lapply(1:reps, vz1.fnc, ipm=spg.ipm2233))


## SENSITIVITY
system.time(spg.sens2233 <- lapply(1:reps, sens.fnc, v.z1=spg.vz2233, w.z=spg.wz2233, h=h))

# spg.sum.sens2233 <- Reduce('+', spg.sens2233)
# spg.mean.sens2233 <- spg.sum.sens2233/reps


## ELASTICITY
spg.elas2233 <- lapply(1:reps, elas.fnc, sens=spg.sens2233,
                       ipm=spg.ipm2233, lambda=spg.lambda2233, h=h)

spg.sum.elas2233 <- Reduce('+', spg.elas2233)
spg.mean.elas2233 <- spg.sum.elas2233/reps


## SURV/GROWTH ELASTICITY
system.time(spg.surv.elas2233 <- lapply(1:reps, surv.elas.fnc, ipm=spg.ipm2233, h=h, 
                                        sens=spg.sens2233, lambda=spg.lambda2233))

spg.sum.surv.elas2233 <- Reduce('+', spg.surv.elas2233)
spg.mean.surv.elas2233 <- spg.sum.surv.elas2233/reps

sum(spg.mean.surv.elas2233)*h^2


## FECUNDITY ELASTICITY
system.time(spg.fecund.elas2233 <- lapply(1:reps, fecund.elas.fnc, ipm=spg.ipm2233, h=h, 
                                          sens=spg.sens2233, lambda=spg.lambda2233))

spg.sum.fecund.elas2233 <- Reduce('+', spg.fecund.elas2233)
spg.mean.fecund.elas2233 <- spg.sum.fecund.elas2233/reps

sum(spg.mean.fecund.elas2233)*h^2


spg2233.list <- list()
spg2233.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
spg2233.list[['ipm']] <- spg.ipm2233
spg2233.list[['eigens']] <- spg.eigens2233
spg2233.list[['lambda']] <- spg.lambda2233
spg2233.list[['wz']] <- spg.wz2233
spg2233.list[['vz']] <- spg.vz2233
spg2233.list[['sens']] <- spg.sens2233
spg2233.list[['elas']] <- spg.elas2233
spg2233.list[['mean.elas']] <- spg.mean.elas2233
spg2233.list[['sg.elas']] <- spg.surv.elas2233
spg2233.list[['mean.sg.elas']] <- spg.mean.surv.elas2233
spg2233.list[['fec.elas']] <- spg.fecund.elas2233
spg2233.list[['mean.fec.elas']] <- spg.mean.fecund.elas2233

save(spg2233.list, file='Model Results/spg2233.Rdata')
