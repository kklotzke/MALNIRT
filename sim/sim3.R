
nu.pred <- nu.obs <- de <- ta <- sig <- numeric(10)
for (ii in 1:10)
{
  print(ii)
  data.lnirt <- simdataLNIRT(N = 1000, K = 20, delta = c(0.35,0), tau = c(0.45,0), nu = rep(-0.25, 20))
  out <- MALNIRT(data.lnirt$RTY[,1:20], data.lnirt$RTY[,21:40], XG = 500, est.person = FALSE, silent = TRUE)

  nu.pred[ii] <- mean(out$nu[complete.cases(out$nu)])
  nu.obs[ii] <- mean(diag(cov(data.lnirt$RTZ[,1:20], data.lnirt$RTZ[,21:40])))
  de[ii] <- out$delta[1]
  ta[ii] <- out$tau[1]
  sig[ii] <- out$sig2[1]
}

mean(nu.pred)
mean(nu.obs)
sd(nu.pred)
sd(nu.obs)

mean(de)
sd(de)
mean(ta)
sd(ta)
mean(sig)
sd(sig)
