

source('Pmont demog funs_8june2018.R')

load('bbt ipm data.Rdata')


###### COMMON MODEL PARAMS ##############

set.seed(8675309)

L <- 0
U <- 90
m <- 100

h <- (U - L)/m
meshpts <- L + ((1:m) - 1/2)*h

reps <- 1000 # number of bootstrap replicates
gp.num <- sample(1:nrow(bbt.ipm.dat[['growth.par']]), size=reps, replace=TRUE)
sv.num <- sample(1:nrow(bbt.ipm.dat[['surv.par']]), size=reps, replace=TRUE)
rp.num <- sample(1:nrow(bbt.ipm.dat[['repro.par']]), size=reps, replace=TRUE)


###### ACTIVE ###########

################## 1-1-2-2 ###############

system.time(bbt.ipm1122 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=1, aprecip.num=1, itemp.num=2, iswe.num=2,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens1122 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm1122))

system.time(bbt.lambda1122 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens1122))

system.time(bbt.wz1122 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens1122))

system.time(bbt.vz1122 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm1122))


## SENSITIVITY
system.time(bbt.sens1122 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz1122, w.z=bbt.wz1122, h=h))

# bbt.sum.sens1122 <- Reduce('+', bbt.sens1122)
# bbt.mean.sens1122 <- bbt.sum.sens1122/reps


## ELASTICITY
bbt.elas1122 <- lapply(1:reps, elas.fnc, sens=bbt.sens1122,
                       ipm=bbt.ipm1122, lambda=bbt.lambda1122, h=h)

bbt.sum.elas1122 <- Reduce('+', bbt.elas1122)
bbt.mean.elas1122 <- bbt.sum.elas1122/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas1122 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm1122, h=h, 
                                        sens=bbt.sens1122, lambda=bbt.lambda1122))

bbt.sum.surv.elas1122 <- Reduce('+', bbt.surv.elas1122)
bbt.mean.surv.elas1122 <- bbt.sum.surv.elas1122/reps

sum(bbt.mean.surv.elas1122)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas1122 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm1122, h=h, 
                                          sens=bbt.sens1122, lambda=bbt.lambda1122))

bbt.sum.fecund.elas1122 <- Reduce('+', bbt.fecund.elas1122)
bbt.mean.fecund.elas1122 <- bbt.sum.fecund.elas1122/reps

sum(bbt.mean.fecund.elas1122)*h^2


bbt1122.list <- list()
bbt1122.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt1122.list[['ipm']] <- bbt.ipm1122
bbt1122.list[['eigens']] <- bbt.eigens1122
bbt1122.list[['lambda']] <- bbt.lambda1122
bbt1122.list[['wz']] <- bbt.wz1122
bbt1122.list[['vz']] <- bbt.vz1122
bbt1122.list[['sens']] <- bbt.sens1122
bbt1122.list[['elas']] <- bbt.elas1122
bbt1122.list[['mean.elas']] <- bbt.mean.elas1122
bbt1122.list[['sg.elas']] <- bbt.surv.elas1122
bbt1122.list[['mean.sg.elas']] <- bbt.mean.surv.elas1122
bbt1122.list[['fec.elas']] <- bbt.fecund.elas1122
bbt1122.list[['mean.fec.elas']] <- bbt.mean.fecund.elas1122

save(bbt1122.list, file='Model Results/bbt1122.Rdata')


################## 1-2-2-2 ###############

system.time(bbt.ipm1222 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=1, aprecip.num=2, itemp.num=2, iswe.num=2,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens1222 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm1222))

system.time(bbt.lambda1222 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens1222))

system.time(bbt.wz1222 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens1222))

system.time(bbt.vz1222 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm1222))


## SENSITIVITY
system.time(bbt.sens1222 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz1222, w.z=bbt.wz1222, h=h))

# bbt.sum.sens1222 <- Reduce('+', bbt.sens1222)
# bbt.mean.sens1222 <- bbt.sum.sens1222/reps


## ELASTICITY
bbt.elas1222 <- lapply(1:reps, elas.fnc, sens=bbt.sens1222,
                       ipm=bbt.ipm1222, lambda=bbt.lambda1222, h=h)

bbt.sum.elas1222 <- Reduce('+', bbt.elas1222)
bbt.mean.elas1222 <- bbt.sum.elas1222/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas1222 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm1222, h=h, 
                                        sens=bbt.sens1222, lambda=bbt.lambda1222))

bbt.sum.surv.elas1222 <- Reduce('+', bbt.surv.elas1222)
bbt.mean.surv.elas1222 <- bbt.sum.surv.elas1222/reps

sum(bbt.mean.surv.elas1222)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas1222 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm1222, h=h, 
                                          sens=bbt.sens1222, lambda=bbt.lambda1222))

bbt.sum.fecund.elas1222 <- Reduce('+', bbt.fecund.elas1222)
bbt.mean.fecund.elas1222 <- bbt.sum.fecund.elas1222/reps

sum(bbt.mean.fecund.elas1222)*h^2


bbt1222.list <- list()
bbt1222.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt1222.list[['ipm']] <- bbt.ipm1222
bbt1222.list[['eigens']] <- bbt.eigens1222
bbt1222.list[['lambda']] <- bbt.lambda1222
bbt1222.list[['wz']] <- bbt.wz1222
bbt1222.list[['vz']] <- bbt.vz1222
bbt1222.list[['sens']] <- bbt.sens1222
bbt1222.list[['elas']] <- bbt.elas1222
bbt1222.list[['mean.elas']] <- bbt.mean.elas1222
bbt1222.list[['sg.elas']] <- bbt.surv.elas1222
bbt1222.list[['mean.sg.elas']] <- bbt.mean.surv.elas1222
bbt1222.list[['fec.elas']] <- bbt.fecund.elas1222
bbt1222.list[['mean.fec.elas']] <- bbt.mean.fecund.elas1222

save(bbt1222.list, file='Model Results/bbt1222.Rdata')


################## 2-1-2-2 ###############

system.time(bbt.ipm2122 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=1, itemp.num=2, iswe.num=2,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens2122 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm2122))

system.time(bbt.lambda2122 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens2122))

system.time(bbt.wz2122 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens2122))

system.time(bbt.vz2122 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm2122))


## SENSITIVITY
system.time(bbt.sens2122 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz2122, w.z=bbt.wz2122, h=h))

# bbt.sum.sens2122 <- Reduce('+', bbt.sens2122)
# bbt.mean.sens2122 <- bbt.sum.sens2122/reps


## ELASTICITY
bbt.elas2122 <- lapply(1:reps, elas.fnc, sens=bbt.sens2122,
                       ipm=bbt.ipm2122, lambda=bbt.lambda2122, h=h)

bbt.sum.elas2122 <- Reduce('+', bbt.elas2122)
bbt.mean.elas2122 <- bbt.sum.elas2122/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas2122 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm2122, h=h, 
                                        sens=bbt.sens2122, lambda=bbt.lambda2122))

bbt.sum.surv.elas2122 <- Reduce('+', bbt.surv.elas2122)
bbt.mean.surv.elas2122 <- bbt.sum.surv.elas2122/reps

sum(bbt.mean.surv.elas2122)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas2122 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm2122, h=h, 
                                          sens=bbt.sens2122, lambda=bbt.lambda2122))

bbt.sum.fecund.elas2122 <- Reduce('+', bbt.fecund.elas2122)
bbt.mean.fecund.elas2122 <- bbt.sum.fecund.elas2122/reps

sum(bbt.mean.fecund.elas2122)*h^2


bbt2122.list <- list()
bbt2122.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt2122.list[['ipm']] <- bbt.ipm2122
bbt2122.list[['eigens']] <- bbt.eigens2122
bbt2122.list[['lambda']] <- bbt.lambda2122
bbt2122.list[['wz']] <- bbt.wz2122
bbt2122.list[['vz']] <- bbt.vz2122
bbt2122.list[['sens']] <- bbt.sens2122
bbt2122.list[['elas']] <- bbt.elas2122
bbt2122.list[['mean.elas']] <- bbt.mean.elas2122
bbt2122.list[['sg.elas']] <- bbt.surv.elas2122
bbt2122.list[['mean.sg.elas']] <- bbt.mean.surv.elas2122
bbt2122.list[['fec.elas']] <- bbt.fecund.elas2122
bbt2122.list[['mean.fec.elas']] <- bbt.mean.fecund.elas2122

save(bbt2122.list, file='Model Results/bbt2122.Rdata')



################## 2-2-2-2 ###############

system.time(bbt.ipm2222 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=2, iswe.num=2,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens2222 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm2222))

system.time(bbt.lambda2222 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens2222))

system.time(bbt.wz2222 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens2222))

system.time(bbt.vz2222 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm2222))


## SENSITIVITY
system.time(bbt.sens2222 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz2222, w.z=bbt.wz2222, h=h))

# bbt.sum.sens2222 <- Reduce('+', bbt.sens2222)
# bbt.mean.sens2222 <- bbt.sum.sens2222/reps


## ELASTICITY
bbt.elas2222 <- lapply(1:reps, elas.fnc, sens=bbt.sens2222,
                       ipm=bbt.ipm2222, lambda=bbt.lambda2222, h=h)

bbt.sum.elas2222 <- Reduce('+', bbt.elas2222)
bbt.mean.elas2222 <- bbt.sum.elas2222/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas2222 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm2222, h=h, 
                                        sens=bbt.sens2222, lambda=bbt.lambda2222))

bbt.sum.surv.elas2222 <- Reduce('+', bbt.surv.elas2222)
bbt.mean.surv.elas2222 <- bbt.sum.surv.elas2222/reps

sum(bbt.mean.surv.elas2222)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas2222 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm2222, h=h, 
                                          sens=bbt.sens2222, lambda=bbt.lambda2222))

bbt.sum.fecund.elas2222 <- Reduce('+', bbt.fecund.elas2222)
bbt.mean.fecund.elas2222 <- bbt.sum.fecund.elas2222/reps

sum(bbt.mean.fecund.elas2222)*h^2


bbt2222.list <- list()
bbt2222.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt2222.list[['ipm']] <- bbt.ipm2222
bbt2222.list[['eigens']] <- bbt.eigens2222
bbt2222.list[['lambda']] <- bbt.lambda2222
bbt2222.list[['wz']] <- bbt.wz2222
bbt2222.list[['vz']] <- bbt.vz2222
bbt2222.list[['sens']] <- bbt.sens2222
bbt2222.list[['elas']] <- bbt.elas2222
bbt2222.list[['mean.elas']] <- bbt.mean.elas2222
bbt2222.list[['sg.elas']] <- bbt.surv.elas2222
bbt2222.list[['mean.sg.elas']] <- bbt.mean.surv.elas2222
bbt2222.list[['fec.elas']] <- bbt.fecund.elas2222
bbt2222.list[['mean.fec.elas']] <- bbt.mean.fecund.elas2222

save(bbt2222.list, file='Model Results/bbt2222.Rdata')



################## 1-3-2-2 ###############

system.time(bbt.ipm1322 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=1, aprecip.num=3, itemp.num=2, iswe.num=2,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens1322 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm1322))

system.time(bbt.lambda1322 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens1322))

system.time(bbt.wz1322 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens1322))

system.time(bbt.vz1322 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm1322))


## SENSITIVITY
system.time(bbt.sens1322 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz1322, w.z=bbt.wz1322, h=h))

# bbt.sum.sens1322 <- Reduce('+', bbt.sens1322)
# bbt.mean.sens1322 <- bbt.sum.sens1322/reps


## ELASTICITY
bbt.elas1322 <- lapply(1:reps, elas.fnc, sens=bbt.sens1322,
                       ipm=bbt.ipm1322, lambda=bbt.lambda1322, h=h)

bbt.sum.elas1322 <- Reduce('+', bbt.elas1322)
bbt.mean.elas1322 <- bbt.sum.elas1322/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas1322 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm1322, h=h, 
                                        sens=bbt.sens1322, lambda=bbt.lambda1322))

bbt.sum.surv.elas1322 <- Reduce('+', bbt.surv.elas1322)
bbt.mean.surv.elas1322 <- bbt.sum.surv.elas1322/reps

sum(bbt.mean.surv.elas1322)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas1322 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm1322, h=h, 
                                          sens=bbt.sens1322, lambda=bbt.lambda1322))

bbt.sum.fecund.elas1322 <- Reduce('+', bbt.fecund.elas1322)
bbt.mean.fecund.elas1322 <- bbt.sum.fecund.elas1322/reps

sum(bbt.mean.fecund.elas1322)*h^2


bbt1322.list <- list()
bbt1322.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt1322.list[['ipm']] <- bbt.ipm1322
bbt1322.list[['eigens']] <- bbt.eigens1322
bbt1322.list[['lambda']] <- bbt.lambda1322
bbt1322.list[['wz']] <- bbt.wz1322
bbt1322.list[['vz']] <- bbt.vz1322
bbt1322.list[['sens']] <- bbt.sens1322
bbt1322.list[['elas']] <- bbt.elas1322
bbt1322.list[['mean.elas']] <- bbt.mean.elas1322
bbt1322.list[['sg.elas']] <- bbt.surv.elas1322
bbt1322.list[['mean.sg.elas']] <- bbt.mean.surv.elas1322
bbt1322.list[['fec.elas']] <- bbt.fecund.elas1322
bbt1322.list[['mean.fec.elas']] <- bbt.mean.fecund.elas1322

save(bbt1322.list, file='Model Results/bbt1322.Rdata')



################## 3-1-2-2 ###############

system.time(bbt.ipm3122 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=3, aprecip.num=1, itemp.num=2, iswe.num=2,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens3122 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm3122))

system.time(bbt.lambda3122 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens3122))

system.time(bbt.wz3122 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens3122))

system.time(bbt.vz3122 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm3122))


## SENSITIVITY
system.time(bbt.sens3122 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz3122, w.z=bbt.wz3122, h=h))

# bbt.sum.sens3122 <- Reduce('+', bbt.sens3122)
# bbt.mean.sens3122 <- bbt.sum.sens3122/reps


## ELASTICITY
bbt.elas3122 <- lapply(1:reps, elas.fnc, sens=bbt.sens3122,
                       ipm=bbt.ipm3122, lambda=bbt.lambda3122, h=h)

bbt.sum.elas3122 <- Reduce('+', bbt.elas3122)
bbt.mean.elas3122 <- bbt.sum.elas3122/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas3122 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm3122, h=h, 
                                        sens=bbt.sens3122, lambda=bbt.lambda3122))

bbt.sum.surv.elas3122 <- Reduce('+', bbt.surv.elas3122)
bbt.mean.surv.elas3122 <- bbt.sum.surv.elas3122/reps

sum(bbt.mean.surv.elas3122)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas3122 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm3122, h=h, 
                                          sens=bbt.sens3122, lambda=bbt.lambda3122))

bbt.sum.fecund.elas3122 <- Reduce('+', bbt.fecund.elas3122)
bbt.mean.fecund.elas3122 <- bbt.sum.fecund.elas3122/reps

sum(bbt.mean.fecund.elas3122)*h^2


bbt3122.list <- list()
bbt3122.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt3122.list[['ipm']] <- bbt.ipm3122
bbt3122.list[['eigens']] <- bbt.eigens3122
bbt3122.list[['lambda']] <- bbt.lambda3122
bbt3122.list[['wz']] <- bbt.wz3122
bbt3122.list[['vz']] <- bbt.vz3122
bbt3122.list[['sens']] <- bbt.sens3122
bbt3122.list[['elas']] <- bbt.elas3122
bbt3122.list[['mean.elas']] <- bbt.mean.elas3122
bbt3122.list[['sg.elas']] <- bbt.surv.elas3122
bbt3122.list[['mean.sg.elas']] <- bbt.mean.surv.elas3122
bbt3122.list[['fec.elas']] <- bbt.fecund.elas3122
bbt3122.list[['mean.fec.elas']] <- bbt.mean.fecund.elas3122

save(bbt3122.list, file='Model Results/bbt3122.Rdata')


################## 2-3-2-2 ###############

system.time(bbt.ipm2322 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=3, itemp.num=2, iswe.num=2,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens2322 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm2322))

system.time(bbt.lambda2322 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens2322))

system.time(bbt.wz2322 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens2322))

system.time(bbt.vz2322 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm2322))


## SENSITIVITY
system.time(bbt.sens2322 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz2322, w.z=bbt.wz2322, h=h))

# bbt.sum.sens2322 <- Reduce('+', bbt.sens2322)
# bbt.mean.sens2322 <- bbt.sum.sens2322/reps


## ELASTICITY
bbt.elas2322 <- lapply(1:reps, elas.fnc, sens=bbt.sens2322,
                       ipm=bbt.ipm2322, lambda=bbt.lambda2322, h=h)

bbt.sum.elas2322 <- Reduce('+', bbt.elas2322)
bbt.mean.elas2322 <- bbt.sum.elas2322/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas2322 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm2322, h=h, 
                                        sens=bbt.sens2322, lambda=bbt.lambda2322))

bbt.sum.surv.elas2322 <- Reduce('+', bbt.surv.elas2322)
bbt.mean.surv.elas2322 <- bbt.sum.surv.elas2322/reps

sum(bbt.mean.surv.elas2322)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas2322 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm2322, h=h, 
                                          sens=bbt.sens2322, lambda=bbt.lambda2322))

bbt.sum.fecund.elas2322 <- Reduce('+', bbt.fecund.elas2322)
bbt.mean.fecund.elas2322 <- bbt.sum.fecund.elas2322/reps

sum(bbt.mean.fecund.elas2322)*h^2


bbt2322.list <- list()
bbt2322.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt2322.list[['ipm']] <- bbt.ipm2322
bbt2322.list[['eigens']] <- bbt.eigens2322
bbt2322.list[['lambda']] <- bbt.lambda2322
bbt2322.list[['wz']] <- bbt.wz2322
bbt2322.list[['vz']] <- bbt.vz2322
bbt2322.list[['sens']] <- bbt.sens2322
bbt2322.list[['elas']] <- bbt.elas2322
bbt2322.list[['mean.elas']] <- bbt.mean.elas2322
bbt2322.list[['sg.elas']] <- bbt.surv.elas2322
bbt2322.list[['mean.sg.elas']] <- bbt.mean.surv.elas2322
bbt2322.list[['fec.elas']] <- bbt.fecund.elas2322
bbt2322.list[['mean.fec.elas']] <- bbt.mean.fecund.elas2322

save(bbt2322.list, file='Model Results/bbt2322.Rdata')

################## 3-2-2-2 ###############

system.time(bbt.ipm3222 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=3, aprecip.num=2, itemp.num=2, iswe.num=2,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens3222 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm3222))

system.time(bbt.lambda3222 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens3222))

system.time(bbt.wz3222 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens3222))

system.time(bbt.vz3222 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm3222))


## SENSITIVITY
system.time(bbt.sens3222 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz3222, w.z=bbt.wz3222, h=h))

# bbt.sum.sens3222 <- Reduce('+', bbt.sens3222)
# bbt.mean.sens3222 <- bbt.sum.sens3222/reps


## ELASTICITY
bbt.elas3222 <- lapply(1:reps, elas.fnc, sens=bbt.sens3222,
                       ipm=bbt.ipm3222, lambda=bbt.lambda3222, h=h)

bbt.sum.elas3222 <- Reduce('+', bbt.elas3222)
bbt.mean.elas3222 <- bbt.sum.elas3222/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas3222 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm3222, h=h, 
                                        sens=bbt.sens3222, lambda=bbt.lambda3222))

bbt.sum.surv.elas3222 <- Reduce('+', bbt.surv.elas3222)
bbt.mean.surv.elas3222 <- bbt.sum.surv.elas3222/reps

sum(bbt.mean.surv.elas3222)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas3222 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm3222, h=h, 
                                          sens=bbt.sens3222, lambda=bbt.lambda3222))

bbt.sum.fecund.elas3222 <- Reduce('+', bbt.fecund.elas3222)
bbt.mean.fecund.elas3222 <- bbt.sum.fecund.elas3222/reps

sum(bbt.mean.fecund.elas3222)*h^2


bbt3222.list <- list()
bbt3222.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt3222.list[['ipm']] <- bbt.ipm3222
bbt3222.list[['eigens']] <- bbt.eigens3222
bbt3222.list[['lambda']] <- bbt.lambda3222
bbt3222.list[['wz']] <- bbt.wz3222
bbt3222.list[['vz']] <- bbt.vz3222
bbt3222.list[['sens']] <- bbt.sens3222
bbt3222.list[['elas']] <- bbt.elas3222
bbt3222.list[['mean.elas']] <- bbt.mean.elas3222
bbt3222.list[['sg.elas']] <- bbt.surv.elas3222
bbt3222.list[['mean.sg.elas']] <- bbt.mean.surv.elas3222
bbt3222.list[['fec.elas']] <- bbt.fecund.elas3222
bbt3222.list[['mean.fec.elas']] <- bbt.mean.fecund.elas3222

save(bbt3222.list, file='Model Results/bbt3222.Rdata')


################## 3-3-2-2 ###############

system.time(bbt.ipm3322 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=3, aprecip.num=3, itemp.num=2, iswe.num=2,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens3322 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm3322))

system.time(bbt.lambda3322 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens3322))

system.time(bbt.wz3322 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens3322))

system.time(bbt.vz3322 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm3322))


## SENSITIVITY
system.time(bbt.sens3322 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz3322, w.z=bbt.wz3322, h=h))

# bbt.sum.sens3322 <- Reduce('+', bbt.sens3322)
# bbt.mean.sens3322 <- bbt.sum.sens3322/reps


## ELASTICITY
bbt.elas3322 <- lapply(1:reps, elas.fnc, sens=bbt.sens3322,
                       ipm=bbt.ipm3322, lambda=bbt.lambda3322, h=h)

bbt.sum.elas3322 <- Reduce('+', bbt.elas3322)
bbt.mean.elas3322 <- bbt.sum.elas3322/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas3322 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm3322, h=h, 
                                        sens=bbt.sens3322, lambda=bbt.lambda3322))

bbt.sum.surv.elas3322 <- Reduce('+', bbt.surv.elas3322)
bbt.mean.surv.elas3322 <- bbt.sum.surv.elas3322/reps

sum(bbt.mean.surv.elas3322)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas3322 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm3322, h=h, 
                                          sens=bbt.sens3322, lambda=bbt.lambda3322))

bbt.sum.fecund.elas3322 <- Reduce('+', bbt.fecund.elas3322)
bbt.mean.fecund.elas3322 <- bbt.sum.fecund.elas3322/reps

sum(bbt.mean.fecund.elas3322)*h^2


bbt3322.list <- list()
bbt3322.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt3322.list[['ipm']] <- bbt.ipm3322
bbt3322.list[['eigens']] <- bbt.eigens3322
bbt3322.list[['lambda']] <- bbt.lambda3322
bbt3322.list[['wz']] <- bbt.wz3322
bbt3322.list[['vz']] <- bbt.vz3322
bbt3322.list[['sens']] <- bbt.sens3322
bbt3322.list[['elas']] <- bbt.elas3322
bbt3322.list[['mean.elas']] <- bbt.mean.elas3322
bbt3322.list[['sg.elas']] <- bbt.surv.elas3322
bbt3322.list[['mean.sg.elas']] <- bbt.mean.surv.elas3322
bbt3322.list[['fec.elas']] <- bbt.fecund.elas3322
bbt3322.list[['mean.fec.elas']] <- bbt.mean.fecund.elas3322

save(bbt3322.list, file='Model Results/bbt3322.Rdata')

#### INACTIVE ###############

################## 2-2-1-1 ###############

system.time(bbt.ipm2211 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=1, iswe.num=1,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens2211 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm2211))

system.time(bbt.lambda2211 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens2211))

system.time(bbt.wz2211 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens2211))

system.time(bbt.vz2211 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm2211))


## SENSITIVITY
system.time(bbt.sens2211 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz2211, w.z=bbt.wz2211, h=h))

# bbt.sum.sens2211 <- Reduce('+', bbt.sens2211)
# bbt.mean.sens2211 <- bbt.sum.sens2211/reps


## ELASTICITY
bbt.elas2211 <- lapply(1:reps, elas.fnc, sens=bbt.sens2211,
                       ipm=bbt.ipm2211, lambda=bbt.lambda2211, h=h)

bbt.sum.elas2211 <- Reduce('+', bbt.elas2211)
bbt.mean.elas2211 <- bbt.sum.elas2211/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas2211 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm2211, h=h, 
                                        sens=bbt.sens2211, lambda=bbt.lambda2211))

bbt.sum.surv.elas2211 <- Reduce('+', bbt.surv.elas2211)
bbt.mean.surv.elas2211 <- bbt.sum.surv.elas2211/reps

sum(bbt.mean.surv.elas2211)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas2211 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm2211, h=h, 
                                          sens=bbt.sens2211, lambda=bbt.lambda2211))

bbt.sum.fecund.elas2211 <- Reduce('+', bbt.fecund.elas2211)
bbt.mean.fecund.elas2211 <- bbt.sum.fecund.elas2211/reps

sum(bbt.mean.fecund.elas2211)*h^2


bbt2211.list <- list()
bbt2211.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt2211.list[['ipm']] <- bbt.ipm2211
bbt2211.list[['eigens']] <- bbt.eigens2211
bbt2211.list[['lambda']] <- bbt.lambda2211
bbt2211.list[['wz']] <- bbt.wz2211
bbt2211.list[['vz']] <- bbt.vz2211
bbt2211.list[['sens']] <- bbt.sens2211
bbt2211.list[['elas']] <- bbt.elas2211
bbt2211.list[['mean.elas']] <- bbt.mean.elas2211
bbt2211.list[['sg.elas']] <- bbt.surv.elas2211
bbt2211.list[['mean.sg.elas']] <- bbt.mean.surv.elas2211
bbt2211.list[['fec.elas']] <- bbt.fecund.elas2211
bbt2211.list[['mean.fec.elas']] <- bbt.mean.fecund.elas2211

save(bbt2211.list, file='Model Results/bbt2211.Rdata')


################## 2-2-1-2 ###############


system.time(bbt.ipm2212 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=1, iswe.num=2,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens2212 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm2212))

system.time(bbt.lambda2212 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens2212))

system.time(bbt.wz2212 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens2212))

system.time(bbt.vz2212 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm2212))


## SENSITIVITY
system.time(bbt.sens2212 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz2212, w.z=bbt.wz2212, h=h))

# bbt.sum.sens2212 <- Reduce('+', bbt.sens2212)
# bbt.mean.sens2212 <- bbt.sum.sens2212/reps


## ELASTICITY
bbt.elas2212 <- lapply(1:reps, elas.fnc, sens=bbt.sens2212,
                       ipm=bbt.ipm2212, lambda=bbt.lambda2212, h=h)

bbt.sum.elas2212 <- Reduce('+', bbt.elas2212)
bbt.mean.elas2212 <- bbt.sum.elas2212/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas2212 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm2212, h=h, 
                                        sens=bbt.sens2212, lambda=bbt.lambda2212))

bbt.sum.surv.elas2212 <- Reduce('+', bbt.surv.elas2212)
bbt.mean.surv.elas2212 <- bbt.sum.surv.elas2212/reps

sum(bbt.mean.surv.elas2212)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas2212 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm2212, h=h, 
                                          sens=bbt.sens2212, lambda=bbt.lambda2212))

bbt.sum.fecund.elas2212 <- Reduce('+', bbt.fecund.elas2212)
bbt.mean.fecund.elas2212 <- bbt.sum.fecund.elas2212/reps

sum(bbt.mean.fecund.elas2212)*h^2


bbt2212.list <- list()
bbt2212.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt2212.list[['ipm']] <- bbt.ipm2212
bbt2212.list[['eigens']] <- bbt.eigens2212
bbt2212.list[['lambda']] <- bbt.lambda2212
bbt2212.list[['wz']] <- bbt.wz2212
bbt2212.list[['vz']] <- bbt.vz2212
bbt2212.list[['sens']] <- bbt.sens2212
bbt2212.list[['elas']] <- bbt.elas2212
bbt2212.list[['mean.elas']] <- bbt.mean.elas2212
bbt2212.list[['sg.elas']] <- bbt.surv.elas2212
bbt2212.list[['mean.sg.elas']] <- bbt.mean.surv.elas2212
bbt2212.list[['fec.elas']] <- bbt.fecund.elas2212
bbt2212.list[['mean.fec.elas']] <- bbt.mean.fecund.elas2212

save(bbt2212.list, file='Model Results/bbt2212.Rdata')


################## 2-2-2-1 ###############


system.time(bbt.ipm2221 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=2, iswe.num=1,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens2221 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm2221))

system.time(bbt.lambda2221 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens2221))

system.time(bbt.wz2221 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens2221))

system.time(bbt.vz2221 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm2221))


## SENSITIVITY
system.time(bbt.sens2221 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz2221, w.z=bbt.wz2221, h=h))

# bbt.sum.sens2221 <- Reduce('+', bbt.sens2221)
# bbt.mean.sens2221 <- bbt.sum.sens2221/reps


## ELASTICITY
bbt.elas2221 <- lapply(1:reps, elas.fnc, sens=bbt.sens2221,
                       ipm=bbt.ipm2221, lambda=bbt.lambda2221, h=h)

bbt.sum.elas2221 <- Reduce('+', bbt.elas2221)
bbt.mean.elas2221 <- bbt.sum.elas2221/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas2221 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm2221, h=h, 
                                        sens=bbt.sens2221, lambda=bbt.lambda2221))

bbt.sum.surv.elas2221 <- Reduce('+', bbt.surv.elas2221)
bbt.mean.surv.elas2221 <- bbt.sum.surv.elas2221/reps

sum(bbt.mean.surv.elas2221)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas2221 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm2221, h=h, 
                                          sens=bbt.sens2221, lambda=bbt.lambda2221))

bbt.sum.fecund.elas2221 <- Reduce('+', bbt.fecund.elas2221)
bbt.mean.fecund.elas2221 <- bbt.sum.fecund.elas2221/reps

sum(bbt.mean.fecund.elas2221)*h^2


bbt2221.list <- list()
bbt2221.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt2221.list[['ipm']] <- bbt.ipm2221
bbt2221.list[['eigens']] <- bbt.eigens2221
bbt2221.list[['lambda']] <- bbt.lambda2221
bbt2221.list[['wz']] <- bbt.wz2221
bbt2221.list[['vz']] <- bbt.vz2221
bbt2221.list[['sens']] <- bbt.sens2221
bbt2221.list[['elas']] <- bbt.elas2221
bbt2221.list[['mean.elas']] <- bbt.mean.elas2221
bbt2221.list[['sg.elas']] <- bbt.surv.elas2221
bbt2221.list[['mean.sg.elas']] <- bbt.mean.surv.elas2221
bbt2221.list[['fec.elas']] <- bbt.fecund.elas2221
bbt2221.list[['mean.fec.elas']] <- bbt.mean.fecund.elas2221

save(bbt2221.list, file='Model Results/bbt2221.Rdata')


################## 2-2-1-3 ###############


system.time(bbt.ipm2213 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=1, iswe.num=3,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens2213 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm2213))

system.time(bbt.lambda2213 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens2213))

system.time(bbt.wz2213 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens2213))

system.time(bbt.vz2213 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm2213))


## SENSITIVITY
system.time(bbt.sens2213 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz2213, w.z=bbt.wz2213, h=h))

# bbt.sum.sens2213 <- Reduce('+', bbt.sens2213)
# bbt.mean.sens2213 <- bbt.sum.sens2213/reps


## ELASTICITY
bbt.elas2213 <- lapply(1:reps, elas.fnc, sens=bbt.sens2213,
                       ipm=bbt.ipm2213, lambda=bbt.lambda2213, h=h)

bbt.sum.elas2213 <- Reduce('+', bbt.elas2213)
bbt.mean.elas2213 <- bbt.sum.elas2213/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas2213 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm2213, h=h, 
                                        sens=bbt.sens2213, lambda=bbt.lambda2213))

bbt.sum.surv.elas2213 <- Reduce('+', bbt.surv.elas2213)
bbt.mean.surv.elas2213 <- bbt.sum.surv.elas2213/reps

sum(bbt.mean.surv.elas2213)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas2213 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm2213, h=h, 
                                          sens=bbt.sens2213, lambda=bbt.lambda2213))

bbt.sum.fecund.elas2213 <- Reduce('+', bbt.fecund.elas2213)
bbt.mean.fecund.elas2213 <- bbt.sum.fecund.elas2213/reps

sum(bbt.mean.fecund.elas2213)*h^2


bbt2213.list <- list()
bbt2213.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt2213.list[['ipm']] <- bbt.ipm2213
bbt2213.list[['eigens']] <- bbt.eigens2213
bbt2213.list[['lambda']] <- bbt.lambda2213
bbt2213.list[['wz']] <- bbt.wz2213
bbt2213.list[['vz']] <- bbt.vz2213
bbt2213.list[['sens']] <- bbt.sens2213
bbt2213.list[['elas']] <- bbt.elas2213
bbt2213.list[['mean.elas']] <- bbt.mean.elas2213
bbt2213.list[['sg.elas']] <- bbt.surv.elas2213
bbt2213.list[['mean.sg.elas']] <- bbt.mean.surv.elas2213
bbt2213.list[['fec.elas']] <- bbt.fecund.elas2213
bbt2213.list[['mean.fec.elas']] <- bbt.mean.fecund.elas2213

save(bbt2213.list, file='Model Results/bbt2213.Rdata')


################## 2-2-3-1 ###############


system.time(bbt.ipm2231 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=3, iswe.num=1,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens2231 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm2231))

system.time(bbt.lambda2231 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens2231))

system.time(bbt.wz2231 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens2231))

system.time(bbt.vz2231 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm2231))


## SENSITIVITY
system.time(bbt.sens2231 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz2231, w.z=bbt.wz2231, h=h))

# bbt.sum.sens2231 <- Reduce('+', bbt.sens2231)
# bbt.mean.sens2231 <- bbt.sum.sens2231/reps


## ELASTICITY
bbt.elas2231 <- lapply(1:reps, elas.fnc, sens=bbt.sens2231,
                       ipm=bbt.ipm2231, lambda=bbt.lambda2231, h=h)

bbt.sum.elas2231 <- Reduce('+', bbt.elas2231)
bbt.mean.elas2231 <- bbt.sum.elas2231/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas2231 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm2231, h=h, 
                                        sens=bbt.sens2231, lambda=bbt.lambda2231))

bbt.sum.surv.elas2231 <- Reduce('+', bbt.surv.elas2231)
bbt.mean.surv.elas2231 <- bbt.sum.surv.elas2231/reps

sum(bbt.mean.surv.elas2231)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas2231 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm2231, h=h, 
                                          sens=bbt.sens2231, lambda=bbt.lambda2231))

bbt.sum.fecund.elas2231 <- Reduce('+', bbt.fecund.elas2231)
bbt.mean.fecund.elas2231 <- bbt.sum.fecund.elas2231/reps

sum(bbt.mean.fecund.elas2231)*h^2


bbt2231.list <- list()
bbt2231.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt2231.list[['ipm']] <- bbt.ipm2231
bbt2231.list[['eigens']] <- bbt.eigens2231
bbt2231.list[['lambda']] <- bbt.lambda2231
bbt2231.list[['wz']] <- bbt.wz2231
bbt2231.list[['vz']] <- bbt.vz2231
bbt2231.list[['sens']] <- bbt.sens2231
bbt2231.list[['elas']] <- bbt.elas2231
bbt2231.list[['mean.elas']] <- bbt.mean.elas2231
bbt2231.list[['sg.elas']] <- bbt.surv.elas2231
bbt2231.list[['mean.sg.elas']] <- bbt.mean.surv.elas2231
bbt2231.list[['fec.elas']] <- bbt.fecund.elas2231
bbt2231.list[['mean.fec.elas']] <- bbt.mean.fecund.elas2231

save(bbt2231.list, file='Model Results/bbt2231.Rdata')


################## 2-2-2-3 ###############


system.time(bbt.ipm2223 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=2, iswe.num=3,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens2223 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm2223))

system.time(bbt.lambda2223 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens2223))

system.time(bbt.wz2223 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens2223))

system.time(bbt.vz2223 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm2223))


## SENSITIVITY
system.time(bbt.sens2223 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz2223, w.z=bbt.wz2223, h=h))

# bbt.sum.sens2223 <- Reduce('+', bbt.sens2223)
# bbt.mean.sens2223 <- bbt.sum.sens2223/reps


## ELASTICITY
bbt.elas2223 <- lapply(1:reps, elas.fnc, sens=bbt.sens2223,
                       ipm=bbt.ipm2223, lambda=bbt.lambda2223, h=h)

bbt.sum.elas2223 <- Reduce('+', bbt.elas2223)
bbt.mean.elas2223 <- bbt.sum.elas2223/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas2223 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm2223, h=h, 
                                        sens=bbt.sens2223, lambda=bbt.lambda2223))

bbt.sum.surv.elas2223 <- Reduce('+', bbt.surv.elas2223)
bbt.mean.surv.elas2223 <- bbt.sum.surv.elas2223/reps

sum(bbt.mean.surv.elas2223)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas2223 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm2223, h=h, 
                                          sens=bbt.sens2223, lambda=bbt.lambda2223))

bbt.sum.fecund.elas2223 <- Reduce('+', bbt.fecund.elas2223)
bbt.mean.fecund.elas2223 <- bbt.sum.fecund.elas2223/reps

sum(bbt.mean.fecund.elas2223)*h^2


bbt2223.list <- list()
bbt2223.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt2223.list[['ipm']] <- bbt.ipm2223
bbt2223.list[['eigens']] <- bbt.eigens2223
bbt2223.list[['lambda']] <- bbt.lambda2223
bbt2223.list[['wz']] <- bbt.wz2223
bbt2223.list[['vz']] <- bbt.vz2223
bbt2223.list[['sens']] <- bbt.sens2223
bbt2223.list[['elas']] <- bbt.elas2223
bbt2223.list[['mean.elas']] <- bbt.mean.elas2223
bbt2223.list[['sg.elas']] <- bbt.surv.elas2223
bbt2223.list[['mean.sg.elas']] <- bbt.mean.surv.elas2223
bbt2223.list[['fec.elas']] <- bbt.fecund.elas2223
bbt2223.list[['mean.fec.elas']] <- bbt.mean.fecund.elas2223

save(bbt2223.list, file='Model Results/bbt2223.Rdata')



################## 2-2-3-2 ###############


system.time(bbt.ipm2232 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=3, iswe.num=2,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens2232 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm2232))

system.time(bbt.lambda2232 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens2232))

system.time(bbt.wz2232 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens2232))

system.time(bbt.vz2232 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm2232))


## SENSITIVITY
system.time(bbt.sens2232 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz2232, w.z=bbt.wz2232, h=h))

# bbt.sum.sens2232 <- Reduce('+', bbt.sens2232)
# bbt.mean.sens2232 <- bbt.sum.sens2232/reps


## ELASTICITY
bbt.elas2232 <- lapply(1:reps, elas.fnc, sens=bbt.sens2232,
                       ipm=bbt.ipm2232, lambda=bbt.lambda2232, h=h)

bbt.sum.elas2232 <- Reduce('+', bbt.elas2232)
bbt.mean.elas2232 <- bbt.sum.elas2232/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas2232 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm2232, h=h, 
                                        sens=bbt.sens2232, lambda=bbt.lambda2232))

bbt.sum.surv.elas2232 <- Reduce('+', bbt.surv.elas2232)
bbt.mean.surv.elas2232 <- bbt.sum.surv.elas2232/reps

sum(bbt.mean.surv.elas2232)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas2232 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm2232, h=h, 
                                          sens=bbt.sens2232, lambda=bbt.lambda2232))

bbt.sum.fecund.elas2232 <- Reduce('+', bbt.fecund.elas2232)
bbt.mean.fecund.elas2232 <- bbt.sum.fecund.elas2232/reps

sum(bbt.mean.fecund.elas2232)*h^2


bbt2232.list <- list()
bbt2232.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt2232.list[['ipm']] <- bbt.ipm2232
bbt2232.list[['eigens']] <- bbt.eigens2232
bbt2232.list[['lambda']] <- bbt.lambda2232
bbt2232.list[['wz']] <- bbt.wz2232
bbt2232.list[['vz']] <- bbt.vz2232
bbt2232.list[['sens']] <- bbt.sens2232
bbt2232.list[['elas']] <- bbt.elas2232
bbt2232.list[['mean.elas']] <- bbt.mean.elas2232
bbt2232.list[['sg.elas']] <- bbt.surv.elas2232
bbt2232.list[['mean.sg.elas']] <- bbt.mean.surv.elas2232
bbt2232.list[['fec.elas']] <- bbt.fecund.elas2232
bbt2232.list[['mean.fec.elas']] <- bbt.mean.fecund.elas2232

save(bbt2232.list, file='Model Results/bbt2232.Rdata')


################## 2-2-3-3 ###############


system.time(bbt.ipm2233 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=bbt.ipm.dat[['growth.par']], 
                                  enviro.par=bbt.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=bbt.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=bbt.ipm.dat[['mean.size']], 
                                  sd.size=bbt.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=3, iswe.num=3,
                                  repro.prob.par=bbt.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=bbt.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=bbt.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=bbt.ipm.dat[['bio4.scale']],
                                  mean.eggs=bbt.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=bbt.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(bbt.eigens2233 <- lapply(1:reps, eigens.fnc, ipm=bbt.ipm2233))

system.time(bbt.lambda2233 <- lapply(1:reps, lambda.fnc, eig.sys=bbt.eigens2233))

system.time(bbt.wz2233 <- lapply(1:reps, wz.fnc, eig.sys=bbt.eigens2233))

system.time(bbt.vz2233 <- lapply(1:reps, vz1.fnc, ipm=bbt.ipm2233))


## SENSITIVITY
system.time(bbt.sens2233 <- lapply(1:reps, sens.fnc, v.z1=bbt.vz2233, w.z=bbt.wz2233, h=h))

# bbt.sum.sens2233 <- Reduce('+', bbt.sens2233)
# bbt.mean.sens2233 <- bbt.sum.sens2233/reps


## ELASTICITY
bbt.elas2233 <- lapply(1:reps, elas.fnc, sens=bbt.sens2233,
                       ipm=bbt.ipm2233, lambda=bbt.lambda2233, h=h)

bbt.sum.elas2233 <- Reduce('+', bbt.elas2233)
bbt.mean.elas2233 <- bbt.sum.elas2233/reps


## SURV/GROWTH ELASTICITY
system.time(bbt.surv.elas2233 <- lapply(1:reps, surv.elas.fnc, ipm=bbt.ipm2233, h=h, 
                                        sens=bbt.sens2233, lambda=bbt.lambda2233))

bbt.sum.surv.elas2233 <- Reduce('+', bbt.surv.elas2233)
bbt.mean.surv.elas2233 <- bbt.sum.surv.elas2233/reps

sum(bbt.mean.surv.elas2233)*h^2


## FECUNDITY ELASTICITY
system.time(bbt.fecund.elas2233 <- lapply(1:reps, fecund.elas.fnc, ipm=bbt.ipm2233, h=h, 
                                          sens=bbt.sens2233, lambda=bbt.lambda2233))

bbt.sum.fecund.elas2233 <- Reduce('+', bbt.fecund.elas2233)
bbt.mean.fecund.elas2233 <- bbt.sum.fecund.elas2233/reps

sum(bbt.mean.fecund.elas2233)*h^2


bbt2233.list <- list()
bbt2233.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
bbt2233.list[['ipm']] <- bbt.ipm2233
bbt2233.list[['eigens']] <- bbt.eigens2233
bbt2233.list[['lambda']] <- bbt.lambda2233
bbt2233.list[['wz']] <- bbt.wz2233
bbt2233.list[['vz']] <- bbt.vz2233
bbt2233.list[['sens']] <- bbt.sens2233
bbt2233.list[['elas']] <- bbt.elas2233
bbt2233.list[['mean.elas']] <- bbt.mean.elas2233
bbt2233.list[['sg.elas']] <- bbt.surv.elas2233
bbt2233.list[['mean.sg.elas']] <- bbt.mean.surv.elas2233
bbt2233.list[['fec.elas']] <- bbt.fecund.elas2233
bbt2233.list[['mean.fec.elas']] <- bbt.mean.fecund.elas2233

save(bbt2233.list, file='Model Results/bbt2233.Rdata')
