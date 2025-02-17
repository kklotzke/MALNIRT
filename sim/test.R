dat4.1 <-  simdataLNIRT(N = 500, K = 10, delta = c(0.2,0), tau = c(0.2,0), nu = c(seq(0.2, -0.05, length.out = 3),
                                                                                   seq(-0.1, -0.2, length.out = 4), seq(-0.3, -0.4, length.out = 3)))
dat4.2 <- simdataLNIRT(N = 500, K = 10, delta = c(0.3,0), tau = c(0.3,0),
                       nu = c(seq(-0.05, 0.2, length.out = 3), c(0, -0.2, -0.4, 0.05, -0.1), seq(-0.15, 0.2, length.out = 2)),
                       lambda = dat4.1$lambda, beta = dat4.1$beta, theta.offset = 0, zeta.offset = 0)
print(mean(dat4.2$Z) + mean(dat4.1$beta))

#dat4.1 <-  simdataLNIRT(N = 500, K = 10, delta = c(0.2,0), tau = c(0.2,0), nu = rep(0,10))
#dat4.2 <- simdataLNIRT(N = 500, K = 10, delta = c(0.3,0), tau = c(0.3,0),
#                       nu = rep(0,10), lambda = dat4.1$lambda, beta = dat4.1$beta, theta.offset = -0.2, zeta.offset = 0.4)


#dat4.3 <- simdataLNIRT(N = 350, K = 10, delta = c(0.15,0), tau = c(0.15,0), nu = rep(0,10), lambda = dat4.1$lambda, beta = dat4.1$beta, theta.offset = -0.2, zeta.offset = -0.45)
group <- c(rep(1, 500), rep(2, 500))#, rep(3, 350))
y.all <- rbind(dat4.1$Y, dat4.2$Y)#, dat4.3$Y)
rt.all <- rbind(dat4.1$RT, dat4.2$RT)#, dat4.3$RT)



out4 <- MALNIRT(Y = y.all, RT = rt.all, group = group, XG = 600, est.person = FALSE)
summary((out4$post.means[[1]]$beta + out4$post.means[[2]]$beta)/2 - dat4.1$beta)
summary((out4$post.means[[1]]$lambda + out4$post.means[[2]]$lambda)/2 - dat4.1$lambda)
c(out4$post.means[[1]]$theta, out4$post.means[[2]]$theta)#, out4$post.means[[3]]$theta)
c(out4$post.means[[1]]$zeta, out4$post.means[[2]]$zeta)#, out4$post.means[[3]]$zeta)
c(out4$post.means[[1]]$tau, out4$post.means[[2]]$tau)#, out4$post.means[[3]]$tau)
c(out4$post.means[[1]]$delta, out4$post.means[[2]]$delta)#, out4$post.means[[3]]$delta)
rbind(out4$post.means[[1]]$nu, out4$post.means[[2]]$nu)#, out4$post.means[[3]]$nu)
rowMeans(rbind(out4$post.means[[1]]$nu, out4$post.means[[2]]$nu))#, out4$post.means[[3]]$nu))
BIC.zeta <- out4$BIC.zeta
BIC.zeta$BIC2 - BIC.zeta$BIC1
cor(out4$post.means[[1]]$sig2k, dat4.1$sig2k)
cor(out4$post.means[[2]]$sig2k, dat4.2$sig2k)




library(LNIRT)
#dat.lnirt <- simLNIRT123(N = 500, K = 10, rho = 0.6, WL = FALSE, td = TRUE, id = TRUE)
#dat.lnirt <- simLNIRT(N = 500, K = 10, rho = 0.6, WL = FALSE, td = TRUE, id = TRUE)

dat5 <- simdataLNIRT(N = 500, K = 10, delta = c(0.6,0), tau = c(0.5,0), nu = c(rep(-0.2, 5), rep(0.4, 5)))
dat6 <- simdataLNIRTP(N = 500, K = 10, delta = 1, tau = 1, rho=0, nu = c(rep(-0.6, 5), rep(0.6, 5)))
round(cov(dat6$th, dat6$ze), digits=5)
round(cov(dat6$th, dat6$th), digits=1)
round(cov(dat6$ze, dat6$ze), digits=1)


round(cov(dat6$Z, dat6$RT), digits=1)
round(cov(dat6$Z), digits=1)
round(cov(dat6$RT), digits=1)
round(cov(dat6$theta, dat6$zeta), digits=1)
round(cov(dat6$theta, dat6$theta), digits=1)
round(cov(dat6$zeta, dat6$zeta), digits=1)

round(cov(dat6$theta[,1], dat6$zeta[,1]), digits=2)
var(dat6$theta[,1])
var(dat6$zeta[,1])


cov(rowMeans(dat6$theta), rowMeans(dat6$zeta))
cov(rowMeans(dat6$theta), rowMeans(dat6$theta))
cov(rowMeans(dat6$zeta), rowMeans(dat6$zeta))

cov(rowMeans(dat6$theta[,1:2]), rowMeans(dat6$theta[, 3:4]))
cov(rowMeans(dat6$zeta[,1:2]), rowMeans(dat6$zeta[, 3:4]))
cov(dat6$beta, dat6$lambda)


#out.lnirt <- LNIRT(Y = Y1, RT = RT1, data = dat.lnirt, XG = 600, WL = FALSE, td = TRUE)
out.lnirt2 <- LNIRT(Y = Y, RT = RT, data = dat6, XG = 600, WL = FALSE, td = TRUE)
#s.lnirt2 <- summary(out.lnirt2)

#out5 <- MALNIRT3Steps(Y = dat.lnirt$Y1, RT = dat.lnirt$RT1, XG = 600, est.person = TRUE)
out52 <- MALNIRT3Steps(Y = Y, RT = RT, data = dat6, XG = 800, est.person = TRUE)
#out53 <- MAIRT(Y = dat.lnirt$Y1, XG = 600, est.person = FALSE)

mean(out52$post.means[[1]]$nu)
out52$post.means[[1]]$tau
out52$post.means[[1]]$delta
summary(out52$post.means[[1]]$nu + dat6$nu )#diag(cov(dat6$Z, dat6$RT)))#(dat6$nu))
summary(out52$post.means[[1]]$beta - dat6$beta)
summary(out52$post.means[[1]]$lambda - dat6$lambda)
summary(out52$post.means[[1]]$sig2k - dat6$sig2k)
#cor(out52$post.means[[1]]$theta_i, rowMeans(dat6$theta))
#cor(out52$post.means[[1]]$zeta_i, rowMeans(dat6$zeta))
mean(out52$post.means[[1]]$theta_ik[6,1:10]) - mean(dat6$theta[6,1:10])

cor(rowMeans(out52$post.means[[1]]$theta_ik[,1:10]), rowMeans(dat6$theta[,1:10]))
cor(rowMeans(out52$post.means[[1]]$zeta_ik[,1:10]), rowMeans(dat6$zeta[,1:10]))

theta_all <- rowMeans(out52$post.means[[1]]$theta_ik[,1:10])
theta_1 <- rowMeans(out52$post.means[[1]]$theta_ik[,1:5])
theta_2 <- rowMeans(out52$post.means[[1]]$theta_ik[,6:10])

zeta_all <- rowMeans(out52$post.means[[1]]$zeta_ik[,1:10])
zeta_1 <- rowMeans(out52$post.means[[1]]$zeta_ik[,1:5])
zeta_2 <- rowMeans(out52$post.means[[1]]$zeta_ik[,6:10])

diag(cor(out52$post.means[[1]]$theta_ik, dat6$theta))
diag(cor(out52$post.means[[1]]$zeta_ik, dat6$zeta))
cov(out52$post.means[[1]]$theta_ik, out52$post.means[[1]]$theta_ik)
cov(out52$post.means[[1]]$zeta_ik, out52$post.means[[1]]$zeta_ik)
cov(out52$post.means[[1]]$theta_ik, out52$post.means[[1]]$zeta_ik)

summary(out52$post.means[[1]]$theta_ik[,] - dat6$theta[,])
summary(out52$post.means[[1]]$zeta_ik[,] - dat6$zeta[,])

hist(rowMeans(out52$post.means[[1]]$theta_ik[,1:5]) - rowMeans(dat6$theta[,1:5]), breaks=20)
summary(out.lnirt2$Mtheta[,1] - rowMeans(dat6$theta[,1:15]))
summary(rowMeans(out52$post.means[[1]]$zeta_ik[,1:15]) - rowMeans(dat6$zeta[,1:15]))
summary(out.lnirt2$Mtheta[,2] - rowMeans(dat6$zeta[,1:15]))
cor(out.lnirt2$Mtheta[,1], rowMeans(dat6$theta[,1:5]))
cor(out.lnirt2$Mtheta[,2], rowMeans(dat6$zeta[,1:15]))

round(cov(out52$post.means[[1]]$theta_ik, out52$post.means[[1]]$zeta_ik), digits=2)

cor(out52$post.means[[1]]$theta_i, dat.lnirt$theta[, 1])
cor(out52$post.means[[1]]$zeta_i, dat.lnirt$theta[, 2])
cov(out52$post.means[[1]]$theta_i, out5$post.means[[1]]$zeta_i)
cor(out52$post.means[[1]]$theta_i, s.lnirt2$Mtheta[,1])
cor(out52$post.means[[1]]$zeta_i, s.lnirt2$Mtheta[,2])


out5i <- MAIRT(Y = Y, data = dat6, XG = 500, est.person = FALSE)
summary(out5i$beta - dat4$beta + 1)
print(mean(out5i$beta))
out5i$tau


#nu1_1 <- c(seq(0.2, -0.05, length.out = 3), seq(-0.1, -0.2, length.out = 4), seq(-0.3, -0.4, length.out = 3))
#nu1_2 <- c(seq(-0.05, 0.2, length.out = 3), c(0, -0.2, -0.4, 0.05, -0.1), seq(-0.15, 0.2, length.out = 2))
nu1_1 <- seq(-0.3, 0.3, length.out = 20)
nu1_2 <- seq(-0.35, 0.35, length.out = 20)
dat1 <- simdataLNIRT(N = 500, K = 20, delta = c(0.3,0), tau = c(0.3,0), nu = nu1_1)
dat2 <- simdataLNIRT(N = 500, K = 20, delta = c(0.4,0), tau = c(0.4,0), nu = nu1_2, beta = dat1$beta, lambda = dat1$lambda, theta.offset = 0.1, zeta.offset = 0.1)

group <- c(rep(1, 500), rep(2, 500))#, rep(3, 350))
y.all <- rbind(dat1$Y, dat2$Y)#, dat4.3$Y)
rt.all <- rbind(dat1$RT, dat2$RT)#, dat4.3$RT)
out4 <- MALNIRT3Steps(Y = y.all, RT = rt.all, group = group, XG = 600, est.person = FALSE, doBIC.zeta = TRUE, doBIC.theta = TRUE)
BIC.zeta <- out4$BIC.zeta
BIC.zeta$BIC1 - BIC.zeta$BIC2
BIC.theta <- out4$BIC.theta
BIC.theta$BIC1 - BIC.theta$BIC2


summary((out4$post.means[[1]]$beta + out4$post.means[[2]]$beta)/2 - dat1$beta)
summary((out4$post.means[[1]]$lambda + out4$post.means[[2]]$lambda)/2 - dat1$lambda)
c(out4$post.means[[1]]$theta, out4$post.means[[2]]$theta)#, out4$post.means[[3]]$theta)
c(out4$post.means[[1]]$zeta, out4$post.means[[2]]$zeta)#, out4$post.means[[3]]$zeta)
c(out4$post.means[[1]]$tau, out4$post.means[[2]]$tau)#, out4$post.means[[3]]$tau)
c(out4$post.means[[1]]$delta, out4$post.means[[2]]$delta)#, out4$post.means[[3]]$delta)
summary(out4$post.means[[1]]$nu - nu1_1)
summary(out4$post.means[[2]]$nu - nu1_2)
summary(out4$post.means[[1]]$sig2k - dat1$sig2k)
summary(out4$post.means[[2]]$sig2k - dat2$sig2k)
summary(out4$post.means[[1]]$nu.sd)
summary(out4$post.means[[2]]$nu.sd)

rbind(out4$post.means[[1]]$nu, out4$post.means[[2]]$nu) - rbind(nu1_1, nu1_2)#, out4$post.means[[3]]$nu)
rowMeans(rbind(out4$post.means[[1]]$nu, out4$post.means[[2]]$nu)) #, out4$post.means[[3]]$nu))
cor(out4$post.means[[1]]$sig2k, dat1$sig2k)
cor(out4$post.means[[2]]$sig2k, dat2$sig2k)


plot(1:10000, out4$samples[[1]]$nu.5, type="l")



dattmp <- simdataLNIRT(N = 800, K = 10, delta = c(0.3,0), tau = c(0.3,0),
                     nu = rep(-0.1,10), theta.offset = 0)

nu <- seq(-0.15, 0.15, length.out = 20)
dat4 <- simdataLNIRT(N = 800, K = 20, delta = c(0.3,0), tau = c(0.5,0),
                       nu = nu, theta.offset = 0, zeta.offset = 0)

out5 <- MALNIRT3Steps(Y = Y, RT = RT, data = dat4, XG = 600, burnin = 0.1, est.person = FALSE, doBIC.nu = TRUE)
BIC.nu <- out5$BIC.nu[[1]]

out52 <- MALNIRT3Steps(Y = Y, RT = RT, data = dat4, XG = 600, burnin = 0.1, est.person = FALSE, doBIC.nu = FALSE, doBIC.full = TRUE)
BIC.full <- out52$BIC.full[[1]]

BIC.full$d - BIC.nu$d
BIC.full$BIC - BIC.nu$BIC

summary(out5$post.means[[1]]$beta - dat4$beta)
summary(out5$post.means[[1]]$lambda - dat4$lambda)
summary(out5$post.means[[1]]$nu - dat4$nu)
summary(out52$post.means[[1]]$nu - dat4$nu)

summary(out5$post.means[[1]]$sig2k - dat4$sig2k)
out5$post.means[[1]]$tau
out5$post.means[[1]]$delta
cor(out5$post.means[[1]]$sig2k, dat4$sig2k)
print(mean(out5$post.means[[1]]$beta))
print(mean(out5$post.means[[1]]$lambda))


out5$post.means[[1]]$beta.sd
out5$post.means[[1]]$lambda.sd
out5$post.means[[1]]$nu.sd
out5$post.means[[1]]$sig2k.sd
out5$post.means[[1]]$tau.sd
out5$post.means[[1]]$delta.sd





plot(1:10, out5$post.means[[1]]$beta + 1, ylim=c(-3,3))
points(dat4$beta, col="blue")
abline(0,0)

out5l <- LNIRT(RT = RT, Y = Y, data = dat4, XG = 1000, td = TRUE)
s.out5l <- summary(out5l)
summary(s.out5l$idiff - dat4$beta)
summary(s.out5l$tintens - dat5$lambda)
summary(s.out5l$estsigma2 - dat5$sig2k)
s.out5l$SigmaPcor
s.out5l$SigmaIcor
cor(s.out5l$estsigma2, dat5$sig2k)


out5i <- MAIRT(Y = Y, data = dat2, XG = 500, est.person = FALSE)
summary(out5i$beta - dat4$beta + 1)
print(mean(out5i$beta))
out5i$tau

out5li <- MALNRT(RT = RT, data = dat4.2, XG = 1000, est.person = FALSE)
summary(out5li$lambda - dat4.2$lambda)
summary(out5li$sig2k[,1] - dat4.2$sig2k)
out5li$delta




plot(1:10, dat4$beta, col="black")
points(out5$post.means[[1]]$beta, col="blue", pch=5)
points(s.out5l$idiff, col="green", pch=2)

plot(1:10, dat4$lambda, col="black")
points(out5$post.means[[1]]$lambda, col="blue", pch=5)
points(s.out5l$tintens, col="green", pch=2)

plot(1:10, dat4$sig2k, col="black")
points(out5$post.means[[1]]$sig2k, col="blue", pch=5)
points(s.out5l$estsigma2, col="green", pch=2)

hist(out5$samples[[1]]$tau, breaks=60)
plot(1:1200, out5$samples[[1]]$tau, type = "l")
plot(1:1200, out5$samples[[1]]$delta, type = "l")
plot(1:1200, out5$samples[[1]]$beta.1, type = "l")
plot(1:1200, out5$samples[[1]]$nu.1, type = "l")
plot(c(60:600,661:1200), out5$samples[[1]]$nu.1[c(60:600,661:1200)], type = "l")
plot(1:1200, out5$samples[[1]]$sig2k.2, type = "l")
plot(1:1200, out5$samples[[1]]$beta.8, type = "l")



#dat <- simMALNIRT(1000, 20, 1, 0, 0)
#out <- MALNRT(RT1, Group = NULL, dat, XG = 1000, burnin = 0.15)
#out$post.lambda; dat$ab[,4]
#out$post.sig2k; dat$sigma2
#out$post.delta
#out$post.sig2/20; mean(dat$sigma2); mean(out$post.sig2k)
#out$post.delta/(out$post.delta + out$post.sig2)
# summary(LNRT(RT = RT1, data = dat, XG = 1000, Discrimination = FALSE, WL = FALSE))
#cor(out$post.lambda, dat$ab[,4])
#cor(out$post.sig2k, dat$sigma2)
#out$post.sig2-sum(dat$sigma2)
#out2 <- LNRT(RT = RT1, data = dat, XG = 1000, Discrimination = FALSE, WL = FALSE)
#cor(colMeans(out2$MAB[,,2]), dat$ab[,4])
#cor(colMeans(out2$Msigma2[101:1000,]), dat$sigma2)
#abs(out$post.sig2k - dat$sigma2)
#abs(out$post.lambda - dat$ab[,4])


dat2 <- simdata(100, 5, c(-0.1, 0), zeta.offset=0)
#out1 <- LNRT(dat2$RT, XG = 1000, Discrimination = FALSE, WL = FALSE)
out2 <- MALNRT(RT, Group = NULL, dat2, XG = 1000, burnin = 0.15, est.person = F)
#cor(out2$zeta_i, dat2$zeta_i); cor(out2$sig2k, dat2$sig2k); cor(out2$lambda, dat2$lambda)

out2$FBF
exp(out2$FBF)
log(out2$FBF)

out2$m0
out2$m1


out2$sig2 - dat2$sig2 / 20
out2$delta[1]
summary(out2$lambda - dat2$lambda)
cor(out2$lambda, dat2$lambda)
cor(out2$sig2k, dat2$sig2k)
#summary(out2$sig2k - dat2$sig2k)
summary(out2$zeta_i - dat2$zeta_i)
plot(out2$zeta_i, dat2$zeta_i)
cor(out2$zeta_i, dat2$zeta_i)
cov(dat2$RT)
pred.RT <- matrix(out2$lambda, nrow = 200, ncol = 20, byrow = T) - (matrix(out2$zeta_i, nrow = 200, ncol = 20, byrow = F))
pred2.RT <- matrix(dat2$lambda, nrow = 200, ncol = 20, byrow = T) - (matrix(dat2$zeta_i, nrow = 200, ncol = 20, byrow = F))
cov(pred.RT[,1], pred.RT[,2])
cov(pred2.RT)


#ICC out2$delta[1]/(out2$delta[1] + out2$sig2[1])
#ICC out2$delta[1]/(out2$delta[1] + out2$sig2[1]*200)
#ICC dat2$delta/(dat2$delta + dat2$sig2/20)
#x <- dat2$RT - matrix(out2$lambda, nrow = 1000, ncol = 20, byrow = TRUE) + out2$zeta_i[,1]
#save(x, file = "errors.Rdata")


## Groups
dat3 <- simdata(200, 20, c(0, 1.0), zeta.offset=0); dat4 <- simdata(200, 20, c(0,0.5), lambda=dat3$lambda, zeta.offset=-2); dat5 <- simdata(200, 20, c(0,0.35), lambda=dat3$lambda, zeta.offset=4)
df <- rbind(dat3$RT, dat4$RT, dat5$RT)
Group <- c(rep(1,200), rep(2,200), rep(3,200))
out3 <- MALNRT(df, Group = Group, data = NULL, XG = 500, burnin = 0.20)
out3$zeta
cor(out3$zeta_i[,1], dat3$zeta_i)
cor(out3$zeta_i[,2], dat4$zeta_i)
cor(out3$zeta_i[,3], dat5$zeta_i)

#ICC out3$post.delta[2]/(out3$post.delta[2] + out3$post.sig2[2]); out3$post.delta[3]/(out3$post.delta[3] + out3$post.sig2[3])
out3$delta[2]; out3$delta[3]; out3$delta[4]
out3$lambda - dat3$lambda
#abs(out3$sig2k[,2] - dat3$sig2k); abs(out3$sig2k[,3] - dat4$sig2k); abs(out3$sig2k[,4] - dat5$sig2k)
#out3$zeta
#mean(dat5$RT) - mean(dat5$lambda)


data.lnirt <- simdataLNIRT(N = 1000, K = 5, delta = c(0.1,0), tau = c(0.2,0), nu = rep(-0.25, 5))
out.lnrt <- MALNRT(data.lnirt$RTY[,1:5], Group = NULL, data = NULL, XG = 1000, burnin = 0.15, est.person = F)
out.lnrt$delta
data.lnirt$lambda
out.lnrt$lambda
data.lnirt$sig2k.lnrt
out.lnrt$sig2k[,1]

out.irt <- MAIRT(data.lnirt$RTY[,11:20], Group = NULL, data = NULL, XG = 1000, burnin = 0.15, est.person = F)

sim1 <- simdataLNRT(1000, 5, c(0.1, 0), zeta.offset=zeta1)
out <- MALNRT(sim1$RT, Group = NULL, data = NULL, XG = XG, burnin = burnin, est.person = FALSE, silent = FALSE)
out$delta
sim1$lambda
out$lambda
sim1$sig2k
out$sig2k[,1]

data.irt <- simdataIRT(1000, 20, c(0, 0.55))
out.irt <- MAIRT(data.irt$Y, XG = 500, est.person = TRUE)
#out.irt <- MAIRT(Y, data = dat, XG = 1000, est.person = FALSE)

out.irt$tau
out.irt$beta
data.irt$beta
cor(data.irt$theta_i, out.irt$theta_i)
data.irt$theta_i[1:5]
out.irt$theta_i[1:5]


#data.lnirt <- simdataLNIRT(N = 1000, K = 10, delta = c(0,0.35), tau = c(0,0.7), nu = rep(0.6, 10))
#data.lnirt <- simdataLNIRT(N = 1000, K = 10, delta = c(0,0), tau = c(0,0), nu = rep(0.5, 10))

out.irt <- MAIRT(data.lnirt$RTY[,21:40], XG = 1500, est.person = TRUE)

data.lnirt <- simdataLNIRT(N = 1000, K = 20, delta = c(0.7,0), tau = c(1,0), nu = rep(-0.25, 20))
out.lnirt <- MALNIRT(data.lnirt$RTY[,1:20], data.lnirt$RTY[,21:40], XG = 1000, est.person = FALSE)

mean(out.lnirt$nu)
out.lnirt$tau
out.lnirt$delta

mean(data.lnirt$sig2k.lnrt) - out.lnirt$sig2
summary(data.lnirt$sig2k.lnrt - out.lnirt$sig2k[,1])

summary(out.lnirt$beta - data.lnirt$beta)
summary(out.lnirt$lambda - data.lnirt$lambda)


mean(diag(cor(out.lnirt$Z, data.lnirt$RTZ[,1:20])))
mean(diag(cor(out.lnirt$Z, data.lnirt$RTZ[,21:40])))

mean(diag(cov(data.lnirt$ZRT, data.lnirt$RTY[,21:40])))

mean(diag(cor(data.lnirt$RTZ[,1:20], data.lnirt$RTZ[,21:40])))

summary(out.lnirt$Z - data.lnirt$RTZ[,1:20])




data.lnirt$beta - out.lnirt$beta

out.irt$tau

data.lnirt2 <- simdataLNIRT(N = 1000, K = 20, delta = c(0,0.25), tau = c(0,0.35), nu = rep(0, 20))

data.lnrt <- simdataLNRT(1000, 10, c(0,0.6))
out.lnrt <- MALNRT(data.lnrt$RT, XG = 500, est.person = FALSE)
summary(data.lnrt$sig2k - out.lnrt$sig2k[,1])
summary(data.lnrt$sig2k - out.lnrt$sig2k[,2])
data.lnrt$sig2/10 - out.lnrt$sig2
out.lnrt$delta

cor(data.lnrt$sig2k, out.lnrt$sig2k[,1])
cor(data.lnrt$sig2k, out.lnrt$sig2k[,2])



hist(out.lnirt$nu.obs, breaks = 25, main = "N-fold observations", xlab = "E[(RT2 - E[RT2])*(Z2 - E[Z2])]")
abline(v=mean(out.lnirt$nu.obs), col="red", lwd=2)
mean(out.lnirt$nu.obs)

hist(out.lnirt$nu.prop, breaks = 25)
mean(out.lnirt$nu.prop)

ind.prop <- sample(1:length(out.lnirt$nu.prop), size = 1000, replace = FALSE)
nu.prop.N <- out.lnirt$nu.prop[ind.prop]

hist(nu.prop.N, breaks = 25, main = "Helmert proposal", xlab = "E[(t2 - E[t2])*(z2 - E[z2])]")
abline(v=mean(nu.prop.N), col="red", lwd=2)
mean(nu.prop.N)

ww1 <- out.lnirt$nu.obs/nu.prop.N
qq1 <- ww1 / sum(ww1)
ind1 <- sample(1:1000, prob = qq1, replace = TRUE)
draw.nu <- out.lnirt$nu.prop[ind1]
hist(draw.nu, breaks = 25, main = "SIR with Helmert proposal", xlab = "nu")
abline(v=mean(draw.nu), col="red", lwd=2)
mean(draw.nu)

mean(diag(cov(data.lnirt$RTZ[,1:20], data.lnirt$RTZ[,21:40])))
cor(out.lnirt$Z[,2], data.lnirt$RTZ[,22])



dat2 <- simdataLNIRT2(N = 1000, K = 20, delta = 0.2, tau = 0.35, nu = runif(20, -0.5, 0.5))
out2 <- MALNIRT(Y = Y, RT = RT, data = dat2, XG = 800, est.person = FALSE)


dat3 <- simdataIRT(1000, 10, c(0.25,0))
out3 <- MAIRT(Y = Y, data = dat4, XG = 500, est.person = FALSE)
out3$tau
summary(out3$beta - dat4$beta)

out2$nu - dat2$nu
summary(out2$nu - dat2$nu)

-colMeans(out2$ZT2) - out2$beta
summary(out2$beta - dat2$beta)
summary(out2$lambda - dat2$lambda)


#out4 <- MALNIRT3(Y = Y, RT = RT, data = dat4, XG = 1000, est.person = FALSE)
#dat4 <- simdataLNIRT2(N = 1000, K = 5, delta = 0.2, tau = 0.3, nu = seq(-0.2, 0.2, length.out = 5))
#dat4 <- simdataLNIRT(N = 1000, K = 10, delta = c(0.2,0), tau = c(0.15,0), nu = seq(-0.25, 0.25, length.out = 10))



out5 <- MALNIRT(Y = Y, RT = RT, data = dat4, XG = 200, est.person = FALSE)
summary(out5$nu - dat4$nu)
plot(out5$nu - dat4$nu)
out5$tau
out5$delta
out5$nu
out5$sd.nu
out5$mce.nu
summary(out5$beta - dat4$beta)
summary(out5$lambda - dat4$lambda)
summary(out5$sig2k - dat4$sig2k)

#out52 <- MALNIRT.noMH(Y = Y, RT = RT, data = dat4, XG = 1000, est.person = FALSE)
out52 <- MALNIRT.1StepZT2(Y = Y, RT = RT, data = dat4, XG = 4000, est.person = FALSE)
summary(out52$nu - dat4$nu)
plot(out52$nu - dat4$nu)
out52$tau
out52$delta
out52$nu
out52$sd.nu
out52$mce.nu
summary(out52$beta - dat4$beta)
summary(out52$lambda - dat4$lambda)
summary(out52$sig2k[,1] - dat4$sig2k)

par(mfrow=c(2,2))
plot(200:4000, out52$data.chain1$nu.1[200:4000], type = "l", col = "red", main = "nu (item 1)",
     xlab = "", ylab = "", frame.plot=F, cex.axis=1.1)
lines(200:4000, out52$data.chain2$nu.1[200:4000], col = "blue")

plot(200:4000, out52$data.chain1$tau.1[200:4000], type = "l", col = "red", main = "tau",
     xlab = "", ylab = "", frame.plot=F, cex.axis=1.1)
lines(200:4000, out52$data.chain2$tau.1[200:4000], col = "blue")

plot(200:4000, out52$data.chain1$delta.1[200:4000], type = "l", col = "red", main = "delta",
     xlab = "", ylab = "", frame.plot=F, cex.axis=1.1)
lines(200:4000, out52$data.chain2$delta.1[200:4000], col = "blue")

plot(200:4000, out52$data.chain1$sig2k.1[200:4000], type = "l", col = "red", main = "sig2 (item 1)",
     xlab = "", ylab = "", frame.plot=F, cex.axis=1.1)
lines(200:4000, out52$data.chain2$sig2k.1[200:4000], col = "blue")

par(mfrow=c(1,1))

plot(200:4000, out52$data.chain1$sig2.1[200:4000], type = "l", col = "red", main = "sig2",
     xlab = "", ylab = "", frame.plot=F, cex.axis=1.1)
lines(200:4000, out52$data.chain2$sig2.1[200:4000], col = "blue")


#dat6 <- simdataLNRT(1000, 10, c(0.2,0))
out6 <- MALNRT(RT = RT, data = dat4, est.person = FALSE)
out6$delta
summary(out6$sig2k[,1] - dat4$sig2k)

out7 <- MAIRT(Y = Y, data = dat4, est.person = FALSE)
out7$tau

diag(cov(dat4$RT, out7$Z))

summary(out4$nu - dat4$nu)
summary(out4$beta - dat4$beta)

#mean(out4$nu)

summary(out4$o - dat4$nu)
#mean(out4$o)
out4$o.accept

