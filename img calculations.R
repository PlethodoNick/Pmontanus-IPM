

source('Pmont demog funs_8june2018.R')

load('img ipm data.Rdata')


###### COMMON MODEL PARAMS ##############

set.seed(8675309)

L <- 0
U <- 90
m <- 100

h <- (U - L)/m
meshpts <- L + ((1:m) - 1/2)*h

reps <- 1000 # number of bootstrap replicates
gp.num <- sample(1:nrow(img.ipm.dat[['growth.par']]), size=reps, replace=TRUE)
sv.num <- sample(1:nrow(img.ipm.dat[['surv.par']]), size=reps, replace=TRUE)
rp.num <- sample(1:nrow(img.ipm.dat[['repro.par']]), size=reps, replace=TRUE)

###### ACTIVE ###########

################## 1-1-2-2 ###############

system.time(img.ipm1122 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=1, aprecip.num=1, itemp.num=2, iswe.num=2,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens1122 <- lapply(1:reps, eigens.fnc, ipm=img.ipm1122))

system.time(img.lambda1122 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens1122))

system.time(img.wz1122 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens1122))

system.time(img.vz1122 <- lapply(1:reps, vz1.fnc, ipm=img.ipm1122))


## SENSITIVITY
system.time(img.sens1122 <- lapply(1:reps, sens.fnc, v.z1=img.vz1122, w.z=img.wz1122, h=h))

# img.sum.sens1122 <- Reduce('+', img.sens1122)
# img.mean.sens1122 <- img.sum.sens1122/reps


## ELASTICITY
img.elas1122 <- lapply(1:reps, elas.fnc, sens=img.sens1122,
                       ipm=img.ipm1122, lambda=img.lambda1122, h=h)

img.sum.elas1122 <- Reduce('+', img.elas1122)
img.mean.elas1122 <- img.sum.elas1122/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas1122 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm1122, h=h, 
                                        sens=img.sens1122, lambda=img.lambda1122))

img.sum.surv.elas1122 <- Reduce('+', img.surv.elas1122)
img.mean.surv.elas1122 <- img.sum.surv.elas1122/reps

sum(img.mean.surv.elas1122)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas1122 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm1122, h=h, 
                                          sens=img.sens1122, lambda=img.lambda1122))

img.sum.fecund.elas1122 <- Reduce('+', img.fecund.elas1122)
img.mean.fecund.elas1122 <- img.sum.fecund.elas1122/reps

sum(img.mean.fecund.elas1122)*h^2


img1122.list <- list()
img1122.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img1122.list[['ipm']] <- img.ipm1122
img1122.list[['eigens']] <- img.eigens1122
img1122.list[['lambda']] <- img.lambda1122
img1122.list[['wz']] <- img.wz1122
img1122.list[['vz']] <- img.vz1122
img1122.list[['sens']] <- img.sens1122
img1122.list[['elas']] <- img.elas1122
img1122.list[['mean.elas']] <- img.mean.elas1122
img1122.list[['sg.elas']] <- img.surv.elas1122
img1122.list[['mean.sg.elas']] <- img.mean.surv.elas1122
img1122.list[['fec.elas']] <- img.fecund.elas1122
img1122.list[['mean.fec.elas']] <- img.mean.fecund.elas1122

save(img1122.list, file='Model Results/img1122.Rdata')


################## 1-2-2-2 ###############

system.time(img.ipm1222 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=1, aprecip.num=2, itemp.num=2, iswe.num=2,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens1222 <- lapply(1:reps, eigens.fnc, ipm=img.ipm1222))

system.time(img.lambda1222 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens1222))

system.time(img.wz1222 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens1222))

system.time(img.vz1222 <- lapply(1:reps, vz1.fnc, ipm=img.ipm1222))


## SENSITIVITY
system.time(img.sens1222 <- lapply(1:reps, sens.fnc, v.z1=img.vz1222, w.z=img.wz1222, h=h))

# img.sum.sens1222 <- Reduce('+', img.sens1222)
# img.mean.sens1222 <- img.sum.sens1222/reps


## ELASTICITY
img.elas1222 <- lapply(1:reps, elas.fnc, sens=img.sens1222,
                       ipm=img.ipm1222, lambda=img.lambda1222, h=h)

img.sum.elas1222 <- Reduce('+', img.elas1222)
img.mean.elas1222 <- img.sum.elas1222/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas1222 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm1222, h=h, 
                                        sens=img.sens1222, lambda=img.lambda1222))

img.sum.surv.elas1222 <- Reduce('+', img.surv.elas1222)
img.mean.surv.elas1222 <- img.sum.surv.elas1222/reps

sum(img.mean.surv.elas1222)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas1222 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm1222, h=h, 
                                          sens=img.sens1222, lambda=img.lambda1222))

img.sum.fecund.elas1222 <- Reduce('+', img.fecund.elas1222)
img.mean.fecund.elas1222 <- img.sum.fecund.elas1222/reps

sum(img.mean.fecund.elas1222)*h^2


img1222.list <- list()
img1222.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img1222.list[['ipm']] <- img.ipm1222
img1222.list[['eigens']] <- img.eigens1222
img1222.list[['lambda']] <- img.lambda1222
img1222.list[['wz']] <- img.wz1222
img1222.list[['vz']] <- img.vz1222
img1222.list[['sens']] <- img.sens1222
img1222.list[['elas']] <- img.elas1222
img1222.list[['mean.elas']] <- img.mean.elas1222
img1222.list[['sg.elas']] <- img.surv.elas1222
img1222.list[['mean.sg.elas']] <- img.mean.surv.elas1222
img1222.list[['fec.elas']] <- img.fecund.elas1222
img1222.list[['mean.fec.elas']] <- img.mean.fecund.elas1222

save(img1222.list, file='Model Results/img1222.Rdata')


################## 2-1-2-2 ###############

system.time(img.ipm2122 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=1, itemp.num=2, iswe.num=2,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens2122 <- lapply(1:reps, eigens.fnc, ipm=img.ipm2122))

system.time(img.lambda2122 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens2122))

system.time(img.wz2122 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens2122))

system.time(img.vz2122 <- lapply(1:reps, vz1.fnc, ipm=img.ipm2122))


## SENSITIVITY
system.time(img.sens2122 <- lapply(1:reps, sens.fnc, v.z1=img.vz2122, w.z=img.wz2122, h=h))

# img.sum.sens2122 <- Reduce('+', img.sens2122)
# img.mean.sens2122 <- img.sum.sens2122/reps


## ELASTICITY
img.elas2122 <- lapply(1:reps, elas.fnc, sens=img.sens2122,
                       ipm=img.ipm2122, lambda=img.lambda2122, h=h)

img.sum.elas2122 <- Reduce('+', img.elas2122)
img.mean.elas2122 <- img.sum.elas2122/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas2122 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm2122, h=h, 
                                        sens=img.sens2122, lambda=img.lambda2122))

img.sum.surv.elas2122 <- Reduce('+', img.surv.elas2122)
img.mean.surv.elas2122 <- img.sum.surv.elas2122/reps

sum(img.mean.surv.elas2122)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas2122 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm2122, h=h, 
                                          sens=img.sens2122, lambda=img.lambda2122))

img.sum.fecund.elas2122 <- Reduce('+', img.fecund.elas2122)
img.mean.fecund.elas2122 <- img.sum.fecund.elas2122/reps

sum(img.mean.fecund.elas2122)*h^2


img2122.list <- list()
img2122.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img2122.list[['ipm']] <- img.ipm2122
img2122.list[['eigens']] <- img.eigens2122
img2122.list[['lambda']] <- img.lambda2122
img2122.list[['wz']] <- img.wz2122
img2122.list[['vz']] <- img.vz2122
img2122.list[['sens']] <- img.sens2122
img2122.list[['elas']] <- img.elas2122
img2122.list[['mean.elas']] <- img.mean.elas2122
img2122.list[['sg.elas']] <- img.surv.elas2122
img2122.list[['mean.sg.elas']] <- img.mean.surv.elas2122
img2122.list[['fec.elas']] <- img.fecund.elas2122
img2122.list[['mean.fec.elas']] <- img.mean.fecund.elas2122

save(img2122.list, file='Model Results/img2122.Rdata')



################## 2-2-2-2 ###############

system.time(img.ipm2222 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=2, iswe.num=2,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens2222 <- lapply(1:reps, eigens.fnc, ipm=img.ipm2222))

system.time(img.lambda2222 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens2222))

system.time(img.wz2222 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens2222))

system.time(img.vz2222 <- lapply(1:reps, vz1.fnc, ipm=img.ipm2222))


## SENSITIVITY
system.time(img.sens2222 <- lapply(1:reps, sens.fnc, v.z1=img.vz2222, w.z=img.wz2222, h=h))

# img.sum.sens2222 <- Reduce('+', img.sens2222)
# img.mean.sens2222 <- img.sum.sens2222/reps


## ELASTICITY
img.elas2222 <- lapply(1:reps, elas.fnc, sens=img.sens2222,
                       ipm=img.ipm2222, lambda=img.lambda2222, h=h)

img.sum.elas2222 <- Reduce('+', img.elas2222)
img.mean.elas2222 <- img.sum.elas2222/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas2222 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm2222, h=h, 
                                        sens=img.sens2222, lambda=img.lambda2222))

img.sum.surv.elas2222 <- Reduce('+', img.surv.elas2222)
img.mean.surv.elas2222 <- img.sum.surv.elas2222/reps

sum(img.mean.surv.elas2222)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas2222 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm2222, h=h, 
                                          sens=img.sens2222, lambda=img.lambda2222))

img.sum.fecund.elas2222 <- Reduce('+', img.fecund.elas2222)
img.mean.fecund.elas2222 <- img.sum.fecund.elas2222/reps

sum(img.mean.fecund.elas2222)*h^2


img2222.list <- list()
img2222.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img2222.list[['ipm']] <- img.ipm2222
img2222.list[['eigens']] <- img.eigens2222
img2222.list[['lambda']] <- img.lambda2222
img2222.list[['wz']] <- img.wz2222
img2222.list[['vz']] <- img.vz2222
img2222.list[['sens']] <- img.sens2222
img2222.list[['elas']] <- img.elas2222
img2222.list[['mean.elas']] <- img.mean.elas2222
img2222.list[['sg.elas']] <- img.surv.elas2222
img2222.list[['mean.sg.elas']] <- img.mean.surv.elas2222
img2222.list[['fec.elas']] <- img.fecund.elas2222
img2222.list[['mean.fec.elas']] <- img.mean.fecund.elas2222

save(img2222.list, file='Model Results/img2222.Rdata')



################## 1-3-2-2 ###############

system.time(img.ipm1322 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=1, aprecip.num=3, itemp.num=2, iswe.num=2,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens1322 <- lapply(1:reps, eigens.fnc, ipm=img.ipm1322))

system.time(img.lambda1322 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens1322))

system.time(img.wz1322 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens1322))

system.time(img.vz1322 <- lapply(1:reps, vz1.fnc, ipm=img.ipm1322))


## SENSITIVITY
system.time(img.sens1322 <- lapply(1:reps, sens.fnc, v.z1=img.vz1322, w.z=img.wz1322, h=h))

# img.sum.sens1322 <- Reduce('+', img.sens1322)
# img.mean.sens1322 <- img.sum.sens1322/reps


## ELASTICITY
img.elas1322 <- lapply(1:reps, elas.fnc, sens=img.sens1322,
                       ipm=img.ipm1322, lambda=img.lambda1322, h=h)

img.sum.elas1322 <- Reduce('+', img.elas1322)
img.mean.elas1322 <- img.sum.elas1322/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas1322 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm1322, h=h, 
                                        sens=img.sens1322, lambda=img.lambda1322))

img.sum.surv.elas1322 <- Reduce('+', img.surv.elas1322)
img.mean.surv.elas1322 <- img.sum.surv.elas1322/reps

sum(img.mean.surv.elas1322)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas1322 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm1322, h=h, 
                                          sens=img.sens1322, lambda=img.lambda1322))

img.sum.fecund.elas1322 <- Reduce('+', img.fecund.elas1322)
img.mean.fecund.elas1322 <- img.sum.fecund.elas1322/reps

sum(img.mean.fecund.elas1322)*h^2


img1322.list <- list()
img1322.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img1322.list[['ipm']] <- img.ipm1322
img1322.list[['eigens']] <- img.eigens1322
img1322.list[['lambda']] <- img.lambda1322
img1322.list[['wz']] <- img.wz1322
img1322.list[['vz']] <- img.vz1322
img1322.list[['sens']] <- img.sens1322
img1322.list[['elas']] <- img.elas1322
img1322.list[['mean.elas']] <- img.mean.elas1322
img1322.list[['sg.elas']] <- img.surv.elas1322
img1322.list[['mean.sg.elas']] <- img.mean.surv.elas1322
img1322.list[['fec.elas']] <- img.fecund.elas1322
img1322.list[['mean.fec.elas']] <- img.mean.fecund.elas1322

save(img1322.list, file='Model Results/img1322.Rdata')



################## 3-1-2-2 ###############

system.time(img.ipm3122 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=3, aprecip.num=1, itemp.num=2, iswe.num=2,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens3122 <- lapply(1:reps, eigens.fnc, ipm=img.ipm3122))

system.time(img.lambda3122 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens3122))

system.time(img.wz3122 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens3122))

system.time(img.vz3122 <- lapply(1:reps, vz1.fnc, ipm=img.ipm3122))


## SENSITIVITY
system.time(img.sens3122 <- lapply(1:reps, sens.fnc, v.z1=img.vz3122, w.z=img.wz3122, h=h))

# img.sum.sens3122 <- Reduce('+', img.sens3122)
# img.mean.sens3122 <- img.sum.sens3122/reps


## ELASTICITY
img.elas3122 <- lapply(1:reps, elas.fnc, sens=img.sens3122,
                       ipm=img.ipm3122, lambda=img.lambda3122, h=h)

img.sum.elas3122 <- Reduce('+', img.elas3122)
img.mean.elas3122 <- img.sum.elas3122/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas3122 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm3122, h=h, 
                                        sens=img.sens3122, lambda=img.lambda3122))

img.sum.surv.elas3122 <- Reduce('+', img.surv.elas3122)
img.mean.surv.elas3122 <- img.sum.surv.elas3122/reps

sum(img.mean.surv.elas3122)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas3122 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm3122, h=h, 
                                          sens=img.sens3122, lambda=img.lambda3122))

img.sum.fecund.elas3122 <- Reduce('+', img.fecund.elas3122)
img.mean.fecund.elas3122 <- img.sum.fecund.elas3122/reps

sum(img.mean.fecund.elas3122)*h^2


img3122.list <- list()
img3122.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img3122.list[['ipm']] <- img.ipm3122
img3122.list[['eigens']] <- img.eigens3122
img3122.list[['lambda']] <- img.lambda3122
img3122.list[['wz']] <- img.wz3122
img3122.list[['vz']] <- img.vz3122
img3122.list[['sens']] <- img.sens3122
img3122.list[['elas']] <- img.elas3122
img3122.list[['mean.elas']] <- img.mean.elas3122
img3122.list[['sg.elas']] <- img.surv.elas3122
img3122.list[['mean.sg.elas']] <- img.mean.surv.elas3122
img3122.list[['fec.elas']] <- img.fecund.elas3122
img3122.list[['mean.fec.elas']] <- img.mean.fecund.elas3122

save(img3122.list, file='Model Results/img3122.Rdata')


################## 2-3-2-2 ###############

system.time(img.ipm2322 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=3, itemp.num=2, iswe.num=2,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens2322 <- lapply(1:reps, eigens.fnc, ipm=img.ipm2322))

system.time(img.lambda2322 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens2322))

system.time(img.wz2322 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens2322))

system.time(img.vz2322 <- lapply(1:reps, vz1.fnc, ipm=img.ipm2322))


## SENSITIVITY
system.time(img.sens2322 <- lapply(1:reps, sens.fnc, v.z1=img.vz2322, w.z=img.wz2322, h=h))

# img.sum.sens2322 <- Reduce('+', img.sens2322)
# img.mean.sens2322 <- img.sum.sens2322/reps


## ELASTICITY
img.elas2322 <- lapply(1:reps, elas.fnc, sens=img.sens2322,
                       ipm=img.ipm2322, lambda=img.lambda2322, h=h)

img.sum.elas2322 <- Reduce('+', img.elas2322)
img.mean.elas2322 <- img.sum.elas2322/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas2322 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm2322, h=h, 
                                        sens=img.sens2322, lambda=img.lambda2322))

img.sum.surv.elas2322 <- Reduce('+', img.surv.elas2322)
img.mean.surv.elas2322 <- img.sum.surv.elas2322/reps

sum(img.mean.surv.elas2322)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas2322 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm2322, h=h, 
                                          sens=img.sens2322, lambda=img.lambda2322))

img.sum.fecund.elas2322 <- Reduce('+', img.fecund.elas2322)
img.mean.fecund.elas2322 <- img.sum.fecund.elas2322/reps

sum(img.mean.fecund.elas2322)*h^2


img2322.list <- list()
img2322.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img2322.list[['ipm']] <- img.ipm2322
img2322.list[['eigens']] <- img.eigens2322
img2322.list[['lambda']] <- img.lambda2322
img2322.list[['wz']] <- img.wz2322
img2322.list[['vz']] <- img.vz2322
img2322.list[['sens']] <- img.sens2322
img2322.list[['elas']] <- img.elas2322
img2322.list[['mean.elas']] <- img.mean.elas2322
img2322.list[['sg.elas']] <- img.surv.elas2322
img2322.list[['mean.sg.elas']] <- img.mean.surv.elas2322
img2322.list[['fec.elas']] <- img.fecund.elas2322
img2322.list[['mean.fec.elas']] <- img.mean.fecund.elas2322

save(img2322.list, file='Model Results/img2322.Rdata')

################## 3-2-2-2 ###############

system.time(img.ipm3222 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=3, aprecip.num=2, itemp.num=2, iswe.num=2,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens3222 <- lapply(1:reps, eigens.fnc, ipm=img.ipm3222))

system.time(img.lambda3222 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens3222))

system.time(img.wz3222 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens3222))

system.time(img.vz3222 <- lapply(1:reps, vz1.fnc, ipm=img.ipm3222))


## SENSITIVITY
system.time(img.sens3222 <- lapply(1:reps, sens.fnc, v.z1=img.vz3222, w.z=img.wz3222, h=h))

# img.sum.sens3222 <- Reduce('+', img.sens3222)
# img.mean.sens3222 <- img.sum.sens3222/reps


## ELASTICITY
img.elas3222 <- lapply(1:reps, elas.fnc, sens=img.sens3222,
                       ipm=img.ipm3222, lambda=img.lambda3222, h=h)

img.sum.elas3222 <- Reduce('+', img.elas3222)
img.mean.elas3222 <- img.sum.elas3222/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas3222 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm3222, h=h, 
                                        sens=img.sens3222, lambda=img.lambda3222))

img.sum.surv.elas3222 <- Reduce('+', img.surv.elas3222)
img.mean.surv.elas3222 <- img.sum.surv.elas3222/reps

sum(img.mean.surv.elas3222)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas3222 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm3222, h=h, 
                                          sens=img.sens3222, lambda=img.lambda3222))

img.sum.fecund.elas3222 <- Reduce('+', img.fecund.elas3222)
img.mean.fecund.elas3222 <- img.sum.fecund.elas3222/reps

sum(img.mean.fecund.elas3222)*h^2


img3222.list <- list()
img3222.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img3222.list[['ipm']] <- img.ipm3222
img3222.list[['eigens']] <- img.eigens3222
img3222.list[['lambda']] <- img.lambda3222
img3222.list[['wz']] <- img.wz3222
img3222.list[['vz']] <- img.vz3222
img3222.list[['sens']] <- img.sens3222
img3222.list[['elas']] <- img.elas3222
img3222.list[['mean.elas']] <- img.mean.elas3222
img3222.list[['sg.elas']] <- img.surv.elas3222
img3222.list[['mean.sg.elas']] <- img.mean.surv.elas3222
img3222.list[['fec.elas']] <- img.fecund.elas3222
img3222.list[['mean.fec.elas']] <- img.mean.fecund.elas3222

save(img3222.list, file='Model Results/img3222.Rdata')


################## 3-3-2-2 ###############

system.time(img.ipm3322 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=3, aprecip.num=3, itemp.num=2, iswe.num=2,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens3322 <- lapply(1:reps, eigens.fnc, ipm=img.ipm3322))

system.time(img.lambda3322 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens3322))

system.time(img.wz3322 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens3322))

system.time(img.vz3322 <- lapply(1:reps, vz1.fnc, ipm=img.ipm3322))


## SENSITIVITY
system.time(img.sens3322 <- lapply(1:reps, sens.fnc, v.z1=img.vz3322, w.z=img.wz3322, h=h))

# img.sum.sens3322 <- Reduce('+', img.sens3322)
# img.mean.sens3322 <- img.sum.sens3322/reps


## ELASTICITY
img.elas3322 <- lapply(1:reps, elas.fnc, sens=img.sens3322,
                       ipm=img.ipm3322, lambda=img.lambda3322, h=h)

img.sum.elas3322 <- Reduce('+', img.elas3322)
img.mean.elas3322 <- img.sum.elas3322/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas3322 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm3322, h=h, 
                                        sens=img.sens3322, lambda=img.lambda3322))

img.sum.surv.elas3322 <- Reduce('+', img.surv.elas3322)
img.mean.surv.elas3322 <- img.sum.surv.elas3322/reps

sum(img.mean.surv.elas3322)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas3322 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm3322, h=h, 
                                          sens=img.sens3322, lambda=img.lambda3322))

img.sum.fecund.elas3322 <- Reduce('+', img.fecund.elas3322)
img.mean.fecund.elas3322 <- img.sum.fecund.elas3322/reps

sum(img.mean.fecund.elas3322)*h^2


img3322.list <- list()
img3322.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img3322.list[['ipm']] <- img.ipm3322
img3322.list[['eigens']] <- img.eigens3322
img3322.list[['lambda']] <- img.lambda3322
img3322.list[['wz']] <- img.wz3322
img3322.list[['vz']] <- img.vz3322
img3322.list[['sens']] <- img.sens3322
img3322.list[['elas']] <- img.elas3322
img3322.list[['mean.elas']] <- img.mean.elas3322
img3322.list[['sg.elas']] <- img.surv.elas3322
img3322.list[['mean.sg.elas']] <- img.mean.surv.elas3322
img3322.list[['fec.elas']] <- img.fecund.elas3322
img3322.list[['mean.fec.elas']] <- img.mean.fecund.elas3322

save(img3322.list, file='Model Results/img3322.Rdata')

#### INACTIVE ###############

################## 2-2-1-1 ###############

system.time(img.ipm2211 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=1, iswe.num=1,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens2211 <- lapply(1:reps, eigens.fnc, ipm=img.ipm2211))

system.time(img.lambda2211 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens2211))

system.time(img.wz2211 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens2211))

system.time(img.vz2211 <- lapply(1:reps, vz1.fnc, ipm=img.ipm2211))


## SENSITIVITY
system.time(img.sens2211 <- lapply(1:reps, sens.fnc, v.z1=img.vz2211, w.z=img.wz2211, h=h))

# img.sum.sens2211 <- Reduce('+', img.sens2211)
# img.mean.sens2211 <- img.sum.sens2211/reps


## ELASTICITY
img.elas2211 <- lapply(1:reps, elas.fnc, sens=img.sens2211,
                       ipm=img.ipm2211, lambda=img.lambda2211, h=h)

img.sum.elas2211 <- Reduce('+', img.elas2211)
img.mean.elas2211 <- img.sum.elas2211/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas2211 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm2211, h=h, 
                                        sens=img.sens2211, lambda=img.lambda2211))

img.sum.surv.elas2211 <- Reduce('+', img.surv.elas2211)
img.mean.surv.elas2211 <- img.sum.surv.elas2211/reps

sum(img.mean.surv.elas2211)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas2211 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm2211, h=h, 
                                          sens=img.sens2211, lambda=img.lambda2211))

img.sum.fecund.elas2211 <- Reduce('+', img.fecund.elas2211)
img.mean.fecund.elas2211 <- img.sum.fecund.elas2211/reps

sum(img.mean.fecund.elas2211)*h^2


img2211.list <- list()
img2211.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img2211.list[['ipm']] <- img.ipm2211
img2211.list[['eigens']] <- img.eigens2211
img2211.list[['lambda']] <- img.lambda2211
img2211.list[['wz']] <- img.wz2211
img2211.list[['vz']] <- img.vz2211
img2211.list[['sens']] <- img.sens2211
img2211.list[['elas']] <- img.elas2211
img2211.list[['mean.elas']] <- img.mean.elas2211
img2211.list[['sg.elas']] <- img.surv.elas2211
img2211.list[['mean.sg.elas']] <- img.mean.surv.elas2211
img2211.list[['fec.elas']] <- img.fecund.elas2211
img2211.list[['mean.fec.elas']] <- img.mean.fecund.elas2211

save(img2211.list, file='Model Results/img2211.Rdata')


################## 2-2-1-2 ###############


system.time(img.ipm2212 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=1, iswe.num=2,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens2212 <- lapply(1:reps, eigens.fnc, ipm=img.ipm2212))

system.time(img.lambda2212 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens2212))

system.time(img.wz2212 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens2212))

system.time(img.vz2212 <- lapply(1:reps, vz1.fnc, ipm=img.ipm2212))


## SENSITIVITY
system.time(img.sens2212 <- lapply(1:reps, sens.fnc, v.z1=img.vz2212, w.z=img.wz2212, h=h))

# img.sum.sens2212 <- Reduce('+', img.sens2212)
# img.mean.sens2212 <- img.sum.sens2212/reps


## ELASTICITY
img.elas2212 <- lapply(1:reps, elas.fnc, sens=img.sens2212,
                       ipm=img.ipm2212, lambda=img.lambda2212, h=h)

img.sum.elas2212 <- Reduce('+', img.elas2212)
img.mean.elas2212 <- img.sum.elas2212/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas2212 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm2212, h=h, 
                                        sens=img.sens2212, lambda=img.lambda2212))

img.sum.surv.elas2212 <- Reduce('+', img.surv.elas2212)
img.mean.surv.elas2212 <- img.sum.surv.elas2212/reps

sum(img.mean.surv.elas2212)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas2212 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm2212, h=h, 
                                          sens=img.sens2212, lambda=img.lambda2212))

img.sum.fecund.elas2212 <- Reduce('+', img.fecund.elas2212)
img.mean.fecund.elas2212 <- img.sum.fecund.elas2212/reps

sum(img.mean.fecund.elas2212)*h^2


img2212.list <- list()
img2212.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img2212.list[['ipm']] <- img.ipm2212
img2212.list[['eigens']] <- img.eigens2212
img2212.list[['lambda']] <- img.lambda2212
img2212.list[['wz']] <- img.wz2212
img2212.list[['vz']] <- img.vz2212
img2212.list[['sens']] <- img.sens2212
img2212.list[['elas']] <- img.elas2212
img2212.list[['mean.elas']] <- img.mean.elas2212
img2212.list[['sg.elas']] <- img.surv.elas2212
img2212.list[['mean.sg.elas']] <- img.mean.surv.elas2212
img2212.list[['fec.elas']] <- img.fecund.elas2212
img2212.list[['mean.fec.elas']] <- img.mean.fecund.elas2212

save(img2212.list, file='Model Results/img2212.Rdata')


################## 2-2-2-1 ###############


system.time(img.ipm2221 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=2, iswe.num=1,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens2221 <- lapply(1:reps, eigens.fnc, ipm=img.ipm2221))

system.time(img.lambda2221 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens2221))

system.time(img.wz2221 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens2221))

system.time(img.vz2221 <- lapply(1:reps, vz1.fnc, ipm=img.ipm2221))


## SENSITIVITY
system.time(img.sens2221 <- lapply(1:reps, sens.fnc, v.z1=img.vz2221, w.z=img.wz2221, h=h))

# img.sum.sens2221 <- Reduce('+', img.sens2221)
# img.mean.sens2221 <- img.sum.sens2221/reps


## ELASTICITY
img.elas2221 <- lapply(1:reps, elas.fnc, sens=img.sens2221,
                       ipm=img.ipm2221, lambda=img.lambda2221, h=h)

img.sum.elas2221 <- Reduce('+', img.elas2221)
img.mean.elas2221 <- img.sum.elas2221/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas2221 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm2221, h=h, 
                                        sens=img.sens2221, lambda=img.lambda2221))

img.sum.surv.elas2221 <- Reduce('+', img.surv.elas2221)
img.mean.surv.elas2221 <- img.sum.surv.elas2221/reps

sum(img.mean.surv.elas2221)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas2221 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm2221, h=h, 
                                          sens=img.sens2221, lambda=img.lambda2221))

img.sum.fecund.elas2221 <- Reduce('+', img.fecund.elas2221)
img.mean.fecund.elas2221 <- img.sum.fecund.elas2221/reps

sum(img.mean.fecund.elas2221)*h^2


img2221.list <- list()
img2221.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img2221.list[['ipm']] <- img.ipm2221
img2221.list[['eigens']] <- img.eigens2221
img2221.list[['lambda']] <- img.lambda2221
img2221.list[['wz']] <- img.wz2221
img2221.list[['vz']] <- img.vz2221
img2221.list[['sens']] <- img.sens2221
img2221.list[['elas']] <- img.elas2221
img2221.list[['mean.elas']] <- img.mean.elas2221
img2221.list[['sg.elas']] <- img.surv.elas2221
img2221.list[['mean.sg.elas']] <- img.mean.surv.elas2221
img2221.list[['fec.elas']] <- img.fecund.elas2221
img2221.list[['mean.fec.elas']] <- img.mean.fecund.elas2221

save(img2221.list, file='Model Results/img2221.Rdata')


################## 2-2-1-3 ###############


system.time(img.ipm2213 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=1, iswe.num=3,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens2213 <- lapply(1:reps, eigens.fnc, ipm=img.ipm2213))

system.time(img.lambda2213 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens2213))

system.time(img.wz2213 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens2213))

system.time(img.vz2213 <- lapply(1:reps, vz1.fnc, ipm=img.ipm2213))


## SENSITIVITY
system.time(img.sens2213 <- lapply(1:reps, sens.fnc, v.z1=img.vz2213, w.z=img.wz2213, h=h))

# img.sum.sens2213 <- Reduce('+', img.sens2213)
# img.mean.sens2213 <- img.sum.sens2213/reps


## ELASTICITY
img.elas2213 <- lapply(1:reps, elas.fnc, sens=img.sens2213,
                       ipm=img.ipm2213, lambda=img.lambda2213, h=h)

img.sum.elas2213 <- Reduce('+', img.elas2213)
img.mean.elas2213 <- img.sum.elas2213/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas2213 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm2213, h=h, 
                                        sens=img.sens2213, lambda=img.lambda2213))

img.sum.surv.elas2213 <- Reduce('+', img.surv.elas2213)
img.mean.surv.elas2213 <- img.sum.surv.elas2213/reps

sum(img.mean.surv.elas2213)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas2213 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm2213, h=h, 
                                          sens=img.sens2213, lambda=img.lambda2213))

img.sum.fecund.elas2213 <- Reduce('+', img.fecund.elas2213)
img.mean.fecund.elas2213 <- img.sum.fecund.elas2213/reps

sum(img.mean.fecund.elas2213)*h^2


img2213.list <- list()
img2213.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img2213.list[['ipm']] <- img.ipm2213
img2213.list[['eigens']] <- img.eigens2213
img2213.list[['lambda']] <- img.lambda2213
img2213.list[['wz']] <- img.wz2213
img2213.list[['vz']] <- img.vz2213
img2213.list[['sens']] <- img.sens2213
img2213.list[['elas']] <- img.elas2213
img2213.list[['mean.elas']] <- img.mean.elas2213
img2213.list[['sg.elas']] <- img.surv.elas2213
img2213.list[['mean.sg.elas']] <- img.mean.surv.elas2213
img2213.list[['fec.elas']] <- img.fecund.elas2213
img2213.list[['mean.fec.elas']] <- img.mean.fecund.elas2213

save(img2213.list, file='Model Results/img2213.Rdata')


################## 2-2-3-1 ###############


system.time(img.ipm2231 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=3, iswe.num=1,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens2231 <- lapply(1:reps, eigens.fnc, ipm=img.ipm2231))

system.time(img.lambda2231 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens2231))

system.time(img.wz2231 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens2231))

system.time(img.vz2231 <- lapply(1:reps, vz1.fnc, ipm=img.ipm2231))


## SENSITIVITY
system.time(img.sens2231 <- lapply(1:reps, sens.fnc, v.z1=img.vz2231, w.z=img.wz2231, h=h))

# img.sum.sens2231 <- Reduce('+', img.sens2231)
# img.mean.sens2231 <- img.sum.sens2231/reps


## ELASTICITY
img.elas2231 <- lapply(1:reps, elas.fnc, sens=img.sens2231,
                       ipm=img.ipm2231, lambda=img.lambda2231, h=h)

img.sum.elas2231 <- Reduce('+', img.elas2231)
img.mean.elas2231 <- img.sum.elas2231/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas2231 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm2231, h=h, 
                                        sens=img.sens2231, lambda=img.lambda2231))

img.sum.surv.elas2231 <- Reduce('+', img.surv.elas2231)
img.mean.surv.elas2231 <- img.sum.surv.elas2231/reps

sum(img.mean.surv.elas2231)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas2231 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm2231, h=h, 
                                          sens=img.sens2231, lambda=img.lambda2231))

img.sum.fecund.elas2231 <- Reduce('+', img.fecund.elas2231)
img.mean.fecund.elas2231 <- img.sum.fecund.elas2231/reps

sum(img.mean.fecund.elas2231)*h^2


img2231.list <- list()
img2231.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img2231.list[['ipm']] <- img.ipm2231
img2231.list[['eigens']] <- img.eigens2231
img2231.list[['lambda']] <- img.lambda2231
img2231.list[['wz']] <- img.wz2231
img2231.list[['vz']] <- img.vz2231
img2231.list[['sens']] <- img.sens2231
img2231.list[['elas']] <- img.elas2231
img2231.list[['mean.elas']] <- img.mean.elas2231
img2231.list[['sg.elas']] <- img.surv.elas2231
img2231.list[['mean.sg.elas']] <- img.mean.surv.elas2231
img2231.list[['fec.elas']] <- img.fecund.elas2231
img2231.list[['mean.fec.elas']] <- img.mean.fecund.elas2231

save(img2231.list, file='Model Results/img2231.Rdata')


################## 2-2-2-3 ###############


system.time(img.ipm2223 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=2, iswe.num=3,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens2223 <- lapply(1:reps, eigens.fnc, ipm=img.ipm2223))

system.time(img.lambda2223 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens2223))

system.time(img.wz2223 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens2223))

system.time(img.vz2223 <- lapply(1:reps, vz1.fnc, ipm=img.ipm2223))


## SENSITIVITY
system.time(img.sens2223 <- lapply(1:reps, sens.fnc, v.z1=img.vz2223, w.z=img.wz2223, h=h))

# img.sum.sens2223 <- Reduce('+', img.sens2223)
# img.mean.sens2223 <- img.sum.sens2223/reps


## ELASTICITY
img.elas2223 <- lapply(1:reps, elas.fnc, sens=img.sens2223,
                       ipm=img.ipm2223, lambda=img.lambda2223, h=h)

img.sum.elas2223 <- Reduce('+', img.elas2223)
img.mean.elas2223 <- img.sum.elas2223/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas2223 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm2223, h=h, 
                                        sens=img.sens2223, lambda=img.lambda2223))

img.sum.surv.elas2223 <- Reduce('+', img.surv.elas2223)
img.mean.surv.elas2223 <- img.sum.surv.elas2223/reps

sum(img.mean.surv.elas2223)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas2223 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm2223, h=h, 
                                          sens=img.sens2223, lambda=img.lambda2223))

img.sum.fecund.elas2223 <- Reduce('+', img.fecund.elas2223)
img.mean.fecund.elas2223 <- img.sum.fecund.elas2223/reps

sum(img.mean.fecund.elas2223)*h^2


img2223.list <- list()
img2223.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img2223.list[['ipm']] <- img.ipm2223
img2223.list[['eigens']] <- img.eigens2223
img2223.list[['lambda']] <- img.lambda2223
img2223.list[['wz']] <- img.wz2223
img2223.list[['vz']] <- img.vz2223
img2223.list[['sens']] <- img.sens2223
img2223.list[['elas']] <- img.elas2223
img2223.list[['mean.elas']] <- img.mean.elas2223
img2223.list[['sg.elas']] <- img.surv.elas2223
img2223.list[['mean.sg.elas']] <- img.mean.surv.elas2223
img2223.list[['fec.elas']] <- img.fecund.elas2223
img2223.list[['mean.fec.elas']] <- img.mean.fecund.elas2223

save(img2223.list, file='Model Results/img2223.Rdata')



################## 2-2-3-2 ###############


system.time(img.ipm2232 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=3, iswe.num=2,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens2232 <- lapply(1:reps, eigens.fnc, ipm=img.ipm2232))

system.time(img.lambda2232 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens2232))

system.time(img.wz2232 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens2232))

system.time(img.vz2232 <- lapply(1:reps, vz1.fnc, ipm=img.ipm2232))


## SENSITIVITY
system.time(img.sens2232 <- lapply(1:reps, sens.fnc, v.z1=img.vz2232, w.z=img.wz2232, h=h))

# img.sum.sens2232 <- Reduce('+', img.sens2232)
# img.mean.sens2232 <- img.sum.sens2232/reps


## ELASTICITY
img.elas2232 <- lapply(1:reps, elas.fnc, sens=img.sens2232,
                       ipm=img.ipm2232, lambda=img.lambda2232, h=h)

img.sum.elas2232 <- Reduce('+', img.elas2232)
img.mean.elas2232 <- img.sum.elas2232/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas2232 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm2232, h=h, 
                                        sens=img.sens2232, lambda=img.lambda2232))

img.sum.surv.elas2232 <- Reduce('+', img.surv.elas2232)
img.mean.surv.elas2232 <- img.sum.surv.elas2232/reps

sum(img.mean.surv.elas2232)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas2232 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm2232, h=h, 
                                          sens=img.sens2232, lambda=img.lambda2232))

img.sum.fecund.elas2232 <- Reduce('+', img.fecund.elas2232)
img.mean.fecund.elas2232 <- img.sum.fecund.elas2232/reps

sum(img.mean.fecund.elas2232)*h^2


img2232.list <- list()
img2232.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img2232.list[['ipm']] <- img.ipm2232
img2232.list[['eigens']] <- img.eigens2232
img2232.list[['lambda']] <- img.lambda2232
img2232.list[['wz']] <- img.wz2232
img2232.list[['vz']] <- img.vz2232
img2232.list[['sens']] <- img.sens2232
img2232.list[['elas']] <- img.elas2232
img2232.list[['mean.elas']] <- img.mean.elas2232
img2232.list[['sg.elas']] <- img.surv.elas2232
img2232.list[['mean.sg.elas']] <- img.mean.surv.elas2232
img2232.list[['fec.elas']] <- img.fecund.elas2232
img2232.list[['mean.fec.elas']] <- img.mean.fecund.elas2232

save(img2232.list, file='Model Results/img2232.Rdata')


################## 2-2-3-3 ###############


system.time(img.ipm2233 <- lapply(1:reps, FUN=kernal.fnc, 
                                  growth.par=img.ipm.dat[['growth.par']], 
                                  enviro.par=img.ipm.dat[['enviro.par']],
                                  gp.num=gp.num,
                                  surv.par=img.ipm.dat[['surv.par']], 
                                  sv.num=sv.num, 
                                  mean.size=img.ipm.dat[['mean.size']], 
                                  sd.size=img.ipm.dat[['sd.size']],
                                  U=U, L=L, m=m, 
                                  atemp.num=2, aprecip.num=2, itemp.num=3, iswe.num=3,
                                  repro.prob.par=img.ipm.dat[['repro.par']],
                                  rp.num=rp.num, 
                                  mean.size.repro=img.ipm.dat[['mean.size.repro']], 
                                  sd.size.repro=img.ipm.dat[['sd.size.repro']], 
                                  bio4.scale=img.ipm.dat[['bio4.scale']],
                                  mean.eggs=img.ipm.dat[['mean.eggs']], 
                                  recruit.size.par=img.ipm.dat[['recruit.size.par']],
                                  small.phi=0, min.phi=0))



system.time(img.eigens2233 <- lapply(1:reps, eigens.fnc, ipm=img.ipm2233))

system.time(img.lambda2233 <- lapply(1:reps, lambda.fnc, eig.sys=img.eigens2233))

system.time(img.wz2233 <- lapply(1:reps, wz.fnc, eig.sys=img.eigens2233))

system.time(img.vz2233 <- lapply(1:reps, vz1.fnc, ipm=img.ipm2233))


## SENSITIVITY
system.time(img.sens2233 <- lapply(1:reps, sens.fnc, v.z1=img.vz2233, w.z=img.wz2233, h=h))

# img.sum.sens2233 <- Reduce('+', img.sens2233)
# img.mean.sens2233 <- img.sum.sens2233/reps


## ELASTICITY
img.elas2233 <- lapply(1:reps, elas.fnc, sens=img.sens2233,
                       ipm=img.ipm2233, lambda=img.lambda2233, h=h)

img.sum.elas2233 <- Reduce('+', img.elas2233)
img.mean.elas2233 <- img.sum.elas2233/reps


## SURV/GROWTH ELASTICITY
system.time(img.surv.elas2233 <- lapply(1:reps, surv.elas.fnc, ipm=img.ipm2233, h=h, 
                                        sens=img.sens2233, lambda=img.lambda2233))

img.sum.surv.elas2233 <- Reduce('+', img.surv.elas2233)
img.mean.surv.elas2233 <- img.sum.surv.elas2233/reps

sum(img.mean.surv.elas2233)*h^2


## FECUNDITY ELASTICITY
system.time(img.fecund.elas2233 <- lapply(1:reps, fecund.elas.fnc, ipm=img.ipm2233, h=h, 
                                          sens=img.sens2233, lambda=img.lambda2233))

img.sum.fecund.elas2233 <- Reduce('+', img.fecund.elas2233)
img.mean.fecund.elas2233 <- img.sum.fecund.elas2233/reps

sum(img.mean.fecund.elas2233)*h^2


img2233.list <- list()
img2233.list[['common.model.params']] <- c(L=L, U=U, m=m, h=h, meshpts=meshpts)
img2233.list[['ipm']] <- img.ipm2233
img2233.list[['eigens']] <- img.eigens2233
img2233.list[['lambda']] <- img.lambda2233
img2233.list[['wz']] <- img.wz2233
img2233.list[['vz']] <- img.vz2233
img2233.list[['sens']] <- img.sens2233
img2233.list[['elas']] <- img.elas2233
img2233.list[['mean.elas']] <- img.mean.elas2233
img2233.list[['sg.elas']] <- img.surv.elas2233
img2233.list[['mean.sg.elas']] <- img.mean.surv.elas2233
img2233.list[['fec.elas']] <- img.fecund.elas2233
img2233.list[['mean.fec.elas']] <- img.mean.fecund.elas2233

save(img2233.list, file='Model Results/img2233.Rdata')
