rm(list = ls())

library(limma)
library(parallel)
library(scales)
####T1 cov 0.8####
#polygenic
dat_flist = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.8_prop1.0",recursive = T,
                       pattern = "phenotype_p1_p*",full.names = T)
dat_list = mclapply(dat_flist,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix = simplify2array(dat_list)
gen = rep(seq(2,102,10),each = 300)

avg = apply(dat_matrix,2,function(x) tapply(x,gen,median))
sigma2 = apply(dat_matrix,2,function(x) tapply(x,gen,var))

avg_std = apply(avg,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std = apply(sqrt(sigma2),2,function(x) (x/x[1]))

# #moderate
# dat_flist2 = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.8_prop0.5/",recursive = T,
#                         pattern = "phenotype_p1_p*",full.names = T)
# dat_list2 = mclapply(dat_flist2,function(x) read.table(x)[,1],mc.cores = 20)
# 
# dat_matrix2 = simplify2array(dat_list2)
# 
# avg2 = apply(dat_matrix2,2,function(x) tapply(x,gen,median))
# sigma2_2 = apply(dat_matrix2,2,function(x) tapply(x,gen,var))
# 
# avg_std2 = apply(avg2,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
# sigma_std2 = apply(sqrt(sigma2_2),2,function(x) (x/x[1]))

#oligonenic
dat_flist3 = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.8_prop0.05/",recursive = T,
                        pattern = "phenotype_p1_p*",full.names = T)
dat_list3 = mclapply(dat_flist3,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix3 = simplify2array(dat_list3)

avg3 = apply(dat_matrix3,2,function(x) tapply(x,gen,median))
sigma2_3 = apply(dat_matrix3,2,function(x) tapply(x,gen,var))

avg_std3 = apply(avg3,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std3 = apply(sqrt(sigma2_3),2,function(x) (x/x[1]))

# avg_sigma_std1 = rbind(apply(sigma_std,1,mean),apply(sigma_std2,1,mean),apply(sigma_std3,1,mean))
avg_sigma_std4 = rbind(apply(sigma_std,1,mean),apply(sigma_std3,1,mean))


####T1 cov 0.5####
#polygenic
dat_flist = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.5_prop1.0",recursive = T,
                       pattern = "phenotype_p1_p*",full.names = T)
dat_list = mclapply(dat_flist,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix = simplify2array(dat_list)
gen = rep(seq(2,102,10),each = 300)

avg = apply(dat_matrix,2,function(x) tapply(x,gen,median))
sigma2 = apply(dat_matrix,2,function(x) tapply(x,gen,var))

avg_std = apply(avg,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std = apply(sqrt(sigma2),2,function(x) (x/x[1]))

# #moderate
# dat_flist2 = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.5_prop0.5/",recursive = T,
#                         pattern = "phenotype_p1_p*",full.names = T)
# dat_list2 = mclapply(dat_flist2,function(x) read.table(x)[,1],mc.cores = 20)
# 
# dat_matrix2 = simplify2array(dat_list2)
# 
# avg2 = apply(dat_matrix2,2,function(x) tapply(x,gen,median))
# sigma2_2 = apply(dat_matrix2,2,function(x) tapply(x,gen,var))
# 
# avg_std2 = apply(avg2,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
# sigma_std2 = apply(sqrt(sigma2_2),2,function(x) (x/x[1]))

#oligonenic
dat_flist3 = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.5_prop0.05/",recursive = T,
                        pattern = "phenotype_p1_p*",full.names = T)
dat_list3 = mclapply(dat_flist3,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix3 = simplify2array(dat_list3)

avg3 = apply(dat_matrix3,2,function(x) tapply(x,gen,median))
sigma2_3 = apply(dat_matrix3,2,function(x) tapply(x,gen,var))

avg_std3 = apply(avg3,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std3 = apply(sqrt(sigma2_3),2,function(x) (x/x[1]))

# avg_sigma_std1 = rbind(apply(sigma_std,1,mean),apply(sigma_std2,1,mean),apply(sigma_std3,1,mean))
avg_sigma_std3 = rbind(apply(sigma_std,1,mean),apply(sigma_std3,1,mean))


####T1 cov 0.1####
#polygenic
dat_flist = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.1_prop1.0",recursive = T,
                       pattern = "phenotype_p1_p*",full.names = T)
dat_list = mclapply(dat_flist,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix = simplify2array(dat_list)
gen = rep(seq(2,102,10),each = 300)

avg = apply(dat_matrix,2,function(x) tapply(x,gen,median))
sigma2 = apply(dat_matrix,2,function(x) tapply(x,gen,var))

avg_std = apply(avg,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std = apply(sqrt(sigma2),2,function(x) (x/x[1]))

#moderate
dat_flist2 = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.1_prop0.5/",recursive = T,
                       pattern = "phenotype_p1_p*",full.names = T)
dat_list2 = mclapply(dat_flist2,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix2 = simplify2array(dat_list2)

avg2 = apply(dat_matrix2,2,function(x) tapply(x,gen,median))
sigma2_2 = apply(dat_matrix2,2,function(x) tapply(x,gen,var))

avg_std2 = apply(avg2,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std2 = apply(sqrt(sigma2_2),2,function(x) (x/x[1]))

#oligonenic
dat_flist3 = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.1_prop0.05/",recursive = T,
                        pattern = "phenotype_p1_p*",full.names = T)
dat_list3 = mclapply(dat_flist3,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix3 = simplify2array(dat_list3)

avg3 = apply(dat_matrix3,2,function(x) tapply(x,gen,median))
sigma2_3 = apply(dat_matrix3,2,function(x) tapply(x,gen,var))

avg_std3 = apply(avg3,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std3 = apply(sqrt(sigma2_3),2,function(x) (x/x[1]))

avg_sigma_std1 = rbind(apply(sigma_std,1,mean),apply(sigma_std2,1,mean),apply(sigma_std3,1,mean))

####T1 cov 0.0####
#polygenic
dat_flist = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.0_prop1.0",recursive = T,
                       pattern = "phenotype_p1_p*",full.names = T)
dat_list = mclapply(dat_flist,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix = simplify2array(dat_list)
gen = rep(seq(2,102,10),each = 300)

avg = apply(dat_matrix,2,function(x) tapply(x,gen,median))
sigma2 = apply(dat_matrix,2,function(x) tapply(x,gen,var))

avg_std = apply(avg,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std = apply(sqrt(sigma2),2,function(x) (x/x[1]))

#moderate
dat_flist2 = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.0_prop0.5/",recursive = T,
                        pattern = "phenotype_p1_p*",full.names = T)
dat_list2 = mclapply(dat_flist2,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix2 = simplify2array(dat_list2)

avg2 = apply(dat_matrix2,2,function(x) tapply(x,gen,median))
sigma2_2 = apply(dat_matrix2,2,function(x) tapply(x,gen,var))

avg_std2 = apply(avg2,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std2 = apply(sqrt(sigma2_2),2,function(x) (x/x[1]))

#oligonenic
dat_flist3 = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.0_prop0.05/",recursive = T,
                        pattern = "phenotype_p1_p*",full.names = T)
dat_list3 = mclapply(dat_flist3,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix3 = simplify2array(dat_list3)

avg3 = apply(dat_matrix3,2,function(x) tapply(x,gen,median))
sigma2_3 = apply(dat_matrix3,2,function(x) tapply(x,gen,var))

avg_std3 = apply(avg3,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std3 = apply(sqrt(sigma2_3),2,function(x) (x/x[1]))

avg_sigma_std2 = rbind(apply(sigma_std,1,mean),apply(sigma_std2,1,mean),apply(sigma_std3,1,mean))


####T2 cov 0.1####
#polygenic
dat_flist = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/T2_cov0.1_prop1.0",recursive = T,
                       pattern = "phenotype_p1_p*",full.names = T)
dat_list = mclapply(dat_flist,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix = simplify2array(dat_list)
gen = rep(seq(2,102,10),each = 300)

avg = apply(dat_matrix,2,function(x) tapply(x,gen,median))
sigma2 = apply(dat_matrix,2,function(x) tapply(x,gen,var))

avg_std = apply(avg,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std = apply(sqrt(sigma2),2,function(x) (x/x[1]))

#moderate
dat_flist2 = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/T2_cov0.1_prop0.5/",recursive = T,
                        pattern = "phenotype_p1_p*",full.names = T)
dat_list2 = mclapply(dat_flist2,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix2 = simplify2array(dat_list2)

avg2 = apply(dat_matrix2,2,function(x) tapply(x,gen,median))
sigma2_2 = apply(dat_matrix2,2,function(x) tapply(x,gen,var))

avg_std2 = apply(avg2,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std2 = apply(sqrt(sigma2_2),2,function(x) (x/x[1]))

#oligonenic
dat_flist3 = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/T2_cov0.1_prop0.05/",recursive = T,
                        pattern = "phenotype_p1_p*",full.names = T)
dat_list3 = mclapply(dat_flist3,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix3 = simplify2array(dat_list3)

avg3 = apply(dat_matrix3,2,function(x) tapply(x,gen,median))
sigma2_3 = apply(dat_matrix3,2,function(x) tapply(x,gen,var))

avg_std3 = apply(avg3,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std3 = apply(sqrt(sigma2_3),2,function(x) (x/x[1]))

avg_sigma_std3 = rbind(apply(sigma_std,1,mean),apply(sigma_std2,1,mean),apply(sigma_std3,1,mean))

####T2 cov 0.0####
#polygenic
dat_flist = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/T2_cov0.0_prop1.0",recursive = T,
                       pattern = "phenotype_p1_p*",full.names = T)
dat_list = mclapply(dat_flist,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix = simplify2array(dat_list)
gen = rep(seq(2,102,10),each = 300)

avg = apply(dat_matrix,2,function(x) tapply(x,gen,median))
sigma2 = apply(dat_matrix,2,function(x) tapply(x,gen,var))

avg_std = apply(avg,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std = apply(sqrt(sigma2),2,function(x) (x/x[1]))

#moderate
dat_flist2 = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/T2_cov0.0_prop0.5/",recursive = T,
                        pattern = "phenotype_p1_p*",full.names = T)
dat_list2 = mclapply(dat_flist2,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix2 = simplify2array(dat_list2)

avg2 = apply(dat_matrix2,2,function(x) tapply(x,gen,median))
sigma2_2 = apply(dat_matrix2,2,function(x) tapply(x,gen,var))

avg_std2 = apply(avg2,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std2 = apply(sqrt(sigma2_2),2,function(x) (x/x[1]))

#oligonenic
dat_flist3 = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/T2_cov0.0_prop0.05/",recursive = T,
                        pattern = "phenotype_p1_p*",full.names = T)
dat_list3 = mclapply(dat_flist3,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix3 = simplify2array(dat_list3)

avg3 = apply(dat_matrix3,2,function(x) tapply(x,gen,median))
sigma2_3 = apply(dat_matrix3,2,function(x) tapply(x,gen,var))

avg_std3 = apply(avg3,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std3 = apply(sqrt(sigma2_3),2,function(x) (x/x[1]))

avg_sigma_std4 = rbind(apply(sigma_std,1,mean),apply(sigma_std2,1,mean),apply(sigma_std3,1,mean))




#### plot ####
png("~/Dropbox (PopGen)/application_for_popgen_vienna/figureSX.png",width = 17.4,height = 8.7,units = "cm",
    pointsize = 10,res = 600)
par(mfrow =c(2,2))
plot(seq(0,100,10),avg_sigma_std4[1,],type ='l',ylim = c(0.5,1),xlim = c(0, 100),lty = 3,lwd = 2,
     xlab = "Generation",ylab = "Phenotypic variance change (F)",main = "Strong genetic correlation")
points(seq(0,100,10),avg_sigma_std4[2,],type ='l',lty = 1, lwd =2)
legend("bottomleft", lty = c(1,3),lwd =2, legend = c("oligogenic", "polygenic"),bty = "n")
plot(seq(0,100,10),avg_sigma_std3[1,],type ='l',ylim = c(0.5,1),xlim = c(0, 100),lty = 3,lwd = 2,
     xlab = "Generation",ylab = "Phenotypic variance change (F)",main = "Moderate genetic correlation")
points(seq(0,100,10),avg_sigma_std3[2,],type ='l',lty = 1, lwd =2)
plot(seq(0,100,10),avg_sigma_std1[1,],type ='l',ylim = c(0.5,1),xlim = c(0, 100),lty = 3,lwd = 2,
     xlab = "Generation",ylab = "Phenotypic variance change (F)",main = "Weak genetic correlation")
#points(seq(0,100,10),avg_sigma_std1[2,],type ='l',lty = 2, lwd = 2)
points(seq(0,100,10),avg_sigma_std1[3,],type ='l',lty = 1, lwd =2)
plot(seq(0,100,10),avg_sigma_std2[1,],type ='l',ylim = c(0.5,1),xlim = c(0, 100),lty = 3,lwd = 2,
     xlab = "Generation",ylab = "Phenotypic variance change (F)",main = "No genetic correlation")
#points(seq(0,100,10),avg_sigma_std2[2,],type ='l',lty = 2, lwd = 2)
points(seq(0,100,10),avg_sigma_std2[3,],type ='l',lty = 1, lwd =2)
dev.off()

png("~/Dropbox (PopGen)/application_for_popgen_vienna/figureSX_update.png",width = 8.7*1.5,height = 8.7,units = "cm",
    pointsize = 8,res = 600)
plot(seq(0,100,10),avg_sigma_std4[1,],type ='l',ylim = c(0.4,1),xlim = c(0, 100),
     lty = 3,lwd = 2,col = "brown"
     ,xlab = "Generation",ylab = "Phenotypic variance change (F)")
points(seq(0,100,10),avg_sigma_std4[2,],type ='l',lty = 1, lwd =2,col = "brown")
legend("bottomleft", lty = c(1,3,1,1,1,1),lwd =2, 
       col = c("black","black","brown","red","orange","grey"), 
       legend = c("oligogenic","polygenic","strong (0.8)", "moderate (0.5)","weak (0.1)","no (0)"),
       bty = "n")

points(seq(0,100,10),avg_sigma_std3[1,],type ='l',ylim = c(0.5,1),xlim = c(0, 100),
       lty = 3,lwd = 2,col = "red")
points(seq(0,100,10),avg_sigma_std3[2,],type ='l',lty = 1, lwd =2,col = "red")

points(seq(0,100,10),avg_sigma_std1[1,],type ='l',ylim = c(0.5,1),xlim = c(0, 100),
       lty = 3,lwd = 2,col = "orange")
points(seq(0,100,10),avg_sigma_std1[3,],type ='l',lty = 1, lwd =2, col ="orange")

points(seq(0,100,10),avg_sigma_std2[1,],type ='l',ylim = c(0.5,1),xlim = c(0, 100),
       lty = 3,lwd = 2,col = "grey")
points(seq(0,100,10),avg_sigma_std2[3,],type ='l',lty = 1, lwd =2, col = "grey")
dev.off()


png("~/Dropbox (PopGen)/application_for_popgen_vienna/figureSX_v2.png",width = 17.4,height = 17.4,units = "cm",
    pointsize = 10,res = 600)
par(mfrow =c(2,2))
plot(seq(0,100,10),avg_sigma_std1[1,],type ='l',ylim = c(0.5,1),xlim = c(0, 100),lty = 3,lwd = 2,
     xlab = "Generation",ylab = "Phenotypic variance change (F)",main = "uncorrelated effects; uncorrelated fitness")
# points(seq(0,100,10),avg_sigma_std1[2,],type ='l',lty = 2, lwd = 2)
points(seq(0,100,10),avg_sigma_std1[3,],type ='l',lty = 1, lwd =2)
# legend("bottomleft", lty = 1:3,lwd =2, legend = c("oligogenic", "polygenic","omnigenic"),bty = "n")
plot(seq(0,100,10),avg_sigma_std2[1,],type ='l',ylim = c(0.5,1),xlim = c(0, 100),lty = 3,lwd = 2,
     xlab = "Generation",ylab = "Phenotypic variance change (F)",main = "correlated effects; uncorrelated fitness")
# points(seq(0,100,10),avg_sigma_std2[2,],type ='l',lty = 2, lwd = 2)
points(seq(0,100,10),avg_sigma_std2[3,],type ='l',lty = 1, lwd =2)

plot(seq(0,100,10),avg_sigma_std3[1,],type ='l',ylim = c(0.5,1),xlim = c(0, 100),lty = 3,lwd = 2,
     xlab = "Generation",ylab = "Phenotypic variance change (F)",main = "uncorrelated effects; correlated fitness")
# points(seq(0,100,10),avg_sigma_std3[2,],type ='l',lty = 2, lwd = 2)
points(seq(0,100,10),avg_sigma_std3[3,],type ='l',lty = 1, lwd =2)
legend("bottomleft", lty = c(1,3),lwd =2, legend = c("oligogenic", "polygenic"),bty = "n")
plot(seq(0,100,10),avg_sigma_std4[1,],type ='l',ylim = c(0.5,1),xlim = c(0, 100),lty = 3,lwd = 2,
     xlab = "Generation",ylab = "Phenotypic variance change (F)",main = "correlated effects; correlated fitness")
# points(seq(0,100,10),avg_sigma_std4[2,],type ='l',lty = 2, lwd = 2)
points(seq(0,100,10),avg_sigma_std4[3,],type ='l',lty = 1, lwd =2)
dev.off()


####between-within module 0.1 0.5####
dat_flist = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.1-0.5_prop0.25",recursive = T,
                       pattern = "phenotype_p1_p*",full.names = T)
dat_list = mclapply(dat_flist,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix = simplify2array(dat_list)
gen = rep(seq(2,102,10),each = 300)

avg = apply(dat_matrix,2,function(x) tapply(x,gen,median))
sigma2 = apply(dat_matrix,2,function(x) tapply(x,gen,var))

avg_std = apply(avg,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std = apply(sqrt(sigma2),2,function(x) (x/x[1]))


dat_flist2 = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.1-0.5_prop0.05",recursive = T,
                       pattern = "phenotype_p1_p*",full.names = T)
dat_list2 = mclapply(dat_flist2,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix2 = simplify2array(dat_list2)

avg = apply(dat_matrix2,2,function(x) tapply(x,gen,median))
sigma2 = apply(dat_matrix2,2,function(x) tapply(x,gen,var))

avg_std2 = apply(avg,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std2 = apply(sqrt(sigma2),2,function(x) (x/x[1]))

avg_sigma_std5 = rbind(apply(sigma_std,1,mean),apply(sigma_std2,1,mean))

####between-within module 0.0 0.5####
dat_flist = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.0-0.5_prop0.25",recursive = T,
                       pattern = "phenotype_p1_p*",full.names = T)
dat_list = mclapply(dat_flist,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix = simplify2array(dat_list)
gen = rep(seq(2,102,10),each = 300)

avg = apply(dat_matrix,2,function(x) tapply(x,gen,median))
sigma2 = apply(dat_matrix,2,function(x) tapply(x,gen,var))

avg_std = apply(avg,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std = apply(sqrt(sigma2),2,function(x) (x/x[1]))


dat_flist2 = list.files("/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov0.0-0.5_prop0.05",recursive = T,
                        pattern = "phenotype_p1_p*",full.names = T)
dat_list2 = mclapply(dat_flist2,function(x) read.table(x)[,1],mc.cores = 20)

dat_matrix2 = simplify2array(dat_list2)

avg = apply(dat_matrix2,2,function(x) tapply(x,gen,median))
sigma2 = apply(dat_matrix2,2,function(x) tapply(x,gen,var))

avg_std2 = apply(avg,2,function(x) (x - x[1]))/sqrt(sigma2[1,])
sigma_std2 = apply(sqrt(sigma2),2,function(x) (x/x[1]))

avg_sigma_std6 = rbind(apply(sigma_std,1,mean),apply(sigma_std2,1,mean))


####visualization####
png("~/Dropbox (PopGen)/application_for_popgen_vienna/figureSY.png",width = 17.4,height = 8.7,units = "cm",
    pointsize = 10,res = 600)
par(mfrow = c(1,2))
plot(seq(0,100,10),avg_sigma_std5[1,],type ='l',ylim = c(0.2,1),xlim = c(0, 100),lty = 3,lwd = 2,
     xlab = "Generation",ylab = "Phenotypic variance change (F)",main = "between 0.1; within 0.5")
points(seq(0,100,10),avg_sigma_std5[2,],type ='l',lty = 1, lwd =2)
# legend("bottomleft", lty = c(1,3),lwd =2, legend = c("oligogenic", "polygenic"),bty = "n")

plot(seq(0,100,10),avg_sigma_std6[1,],type ='l',ylim = c(0.2,1),xlim = c(0, 100),lty = 3,lwd = 2,
     xlab = "Generation",ylab = "Phenotypic variance change (F)",main = "between 0.0; within 0.5")
points(seq(0,100,10),avg_sigma_std6[2,],type ='l',lty = 1, lwd =2)
legend("topright", lty = c(1,3),lwd =2, legend = c("oligogenic", "polygenic"),bty = "n")
dev.off()



plot(seq(2,102,10),apply(avg_std,1,mean),type = "l",ylim = c(-1,3),col = "red")
for (i in 1:1000){
  points(seq(2,102,10),avg_std[,i],type = "l",col =alpha("black",.1))
}
points(seq(2,102,10),apply(avg_std,1,mean),type = "l",ylim = c(0,3),col = "red",lwd =2)


plot(seq(2,102,10),apply(avg_std2,1,mean),type = "l",ylim = c(-1,3),col = "red")
for (i in 1:1000){
  points(seq(2,102,10),avg_std2[,i],type = "l",col =alpha("black",.1))
}
points(seq(2,102,10),apply(avg_std2,1,mean),type = "l",ylim = c(0,3),col = "red",lwd =2)


plot(seq(2,102,10),apply(avg_std3,1,mean),type = "l",ylim = c(-1,3),col = "red")
for (i in 1:1000){
  points(seq(2,102,10),avg_std3[,i],type = "l",col =alpha("black",.1))
}
points(seq(2,102,10),apply(avg_std3,1,mean),type = "l",ylim = c(0,3),col = "red",lwd =2)

#### other ####
boxplot(t(sigma_std))
boxplot(t(sigma_std2))
boxplot(t(sigma_std3))


plot(avg_std[,1],type = "b",ylim =c (-0.5,1.5),col = alpha("black",.1))
for ( i in 2:100){
  points(avg_std[,i],type = "b",col = alpha("black",.1))
}

plot(sigma_std[,1],type = "b",ylim =c (.5,1.2),col = alpha("black",.1))
for ( i in 2:100){
  points(sigma_std[,i],type = "b",col = alpha("black",.1))
}


dat = read.table("/Volumes/Temp1/shengkai/evolVar_pleiotropy/phenotype_p1_p0.txt")
tapply(dat[,1],gen,mean)
