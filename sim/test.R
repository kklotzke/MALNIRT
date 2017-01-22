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



