dat4.1 <- simdataLNIRT(N = 500, K = 20, delta = c(0.0,0), tau = c(0.05,0), nu = rep(-0.25,20))
dat4.2 <- simdataLNIRT(N = 500, K = 20, delta = c(0.1,0), tau = c(0.15,0), nu = rep(-0.15,20), lambda = dat4.1$lambda, beta = dat4.1$beta, theta.offset = 0.5, zeta.offset = 0.5)
dat4.3 <- simdataLNIRT(N = 500, K = 20, delta = c(0.15,0), tau = c(0.10,0), nu = rep(0,20), lambda = dat4.1$lambda, beta = dat4.1$beta, theta.offset = -0.5, zeta.offset = -0.5)
group <- c(rep(1, 500), rep(2, 500), rep(3, 500))
y.all <- rbind(dat4.1$Y, dat4.2$Y, dat4.3$Y)
rt.all <- rbind(dat4.1$RT, dat4.2$RT, dat4.3$RT)

out4 <- MALNIRT(Y = y.all, RT = rt.all, group = group, XG = 800, est.person = FALSE)
summary(out4$post.means[[1]]$beta - dat4.1$beta)
summary(out4$post.means[[1]]$lambda - dat4.1$lambda)
c(out4$post.means[[1]]$theta, out4$post.means[[2]]$theta, out4$post.means[[3]]$theta)
c(out4$post.means[[1]]$zeta, out4$post.means[[2]]$zeta, out4$post.means[[3]]$zeta)
c(out4$post.means[[1]]$tau, out4$post.means[[2]]$tau, out4$post.means[[3]]$tau)
c(out4$post.means[[1]]$delta, out4$post.means[[2]]$delta, out4$post.means[[3]]$delta)
rbind(out4$post.means[[1]]$nu, out4$post.means[[2]]$nu, out4$post.means[[3]]$nu)
rowMeans(rbind(out4$post.means[[1]]$nu, out4$post.means[[2]]$nu, out4$post.means[[3]]$nu))


out5 <- MALNIRT(Y = Y, RT = RT, data = dat4.3, XG = 500, est.person = FALSE)
summary(out5$post.means[[1]]$beta - dat4.1$beta)
summary(out5$post.means[[1]]$lambda - dat4.1$lambda)



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

