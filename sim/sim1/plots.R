setwd("~/Desktop/github_kklotzke/MALNIRT/sim/sim1")

load(file = "simulation1_020417_1.RData")
load(file = "simulation1_020417_2.RData")
load(file = "simulation1_020417_3.RData")
load(file = "simulation1_020417_4.RData")

out.nu1 <- rbind(sim.nu1_1, sim.nu2_1, sim.nu3_1, sim.nu4_1)
out.nu2 <- rbind(sim.nu1_2, sim.nu2_2, sim.nu3_2, sim.nu4_2)


x <- 1:10
y1 <- colMeans(out.nu1)
sd1 <- apply(out.nu1, FUN = sd, MARGIN = 2)

plot(x, nu1_1, col="blue", ylim=c(-0.6, 0.4))
points(y, col="black", pch=5)
segments(x, y-sd,x, y+sd, col="black")
epsilon = 0.02
segments(x-epsilon,y-sd,x+epsilon,y-sd, col="black")
segments(x-epsilon,y+sd,x+epsilon,y+sd, col="black")


x <- 1:10
y2 <- colMeans(out.nu2)
sd2 <- apply(out.nu2, FUN = sd, MARGIN = 2)

plot(x, nu1_2, col="blue", ylim=c(-0.6, 0.4))
points(y, col="black", pch=5)
segments(x, y-sd,x, y+sd, col="black")
epsilon = 0.02
segments(x-epsilon,y-sd,x+epsilon,y-sd, col="black")
segments(x-epsilon,y+sd,x+epsilon,y+sd, col="black")

df <- data.frame("item" = c(1:10, 1:10), "nu.s" = c(nu1_1, nu1_2), "nu" = c(y1, y2), "sd" = c(sd1, sd2), "group" = c(rep("1", 10), rep("2", 10)))

pd <- position_dodge(0.1) # move them .05 to the left and right
ggplot(df, aes(x=item, y=nu, colour=group, group=group)) +
  geom_errorbar(aes(ymin=nu-sd, ymax=nu+sd), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(data=df, aes(x=item, y=nu.s, shape=group), size=3, colour="black", fill="black") + # 21 is filled circle
  geom_point(position=pd, size=3, shape=1, fill="white") + # 21 is filled circle
  scale_x_discrete("Item", limits=1:10) +
  ylab(expression(nu)) +
  scale_colour_hue(name="",    # Legend label, use darker colors
                   breaks=c("1", "2"),
                   labels=c("Group 1", "Group2"),
                   l=40) +                    # Use darker colors, lightness=40
  ggtitle("Posterior means and SD's of cross-covariance parameters") +
  expand_limits(y=0) +                        # Expand y range
  scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0))               # Position legend in bottom right



ggplot(df[df$group==1, ], aes(x=item, y=nu, colour="blue")) +
  geom_errorbar(aes(ymin=nu-sd, ymax=nu+sd), colour="black", width=.1, size=0.3, position=pd) +
  #geom_line(position=pd, colour="blue") +
  geom_point(data=df[df$group==1, ], aes(x=item, y=nu.s), size=1, shape=21, colour="blue", fill="blue") + # 21 is filled circle
  geom_point(position=pd, size=3, shape=2, colour="black", fill="white") + # 21 is filled circle
  scale_x_discrete("Item", limits=1:10) +
  ylab(expression(nu)) +
  # scale_colour_hue(name="",    # Legend label, use darker colors
  #                  breaks=c("1", "2"),
  #                  labels=c("Group 1", "Group2"),
  #                  l=40) +                    # Use darker colors, lightness=40
  ggtitle("Posterior means and SD's of cross-covariance parameters") +
  # expand_limits(y=0) +                        # Expand y range
  scale_y_continuous(limits=c(-0.6, 0.3), breaks=c(-0.50, -0.25, 0, 0.25)) +      # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0))               # Position legend in bottom right


ggplot(df[df$group==2, ], aes(x=item, y=nu, colour="red")) +
  geom_errorbar(aes(ymin=nu-sd, ymax=nu+sd), colour="black", width=.1, size=0.3, position=pd) +
  #geom_line(position=pd, colour="blue") +
  geom_point(data=df[df$group==2, ], aes(x=item, y=nu.s), size=1, shape=21, colour="red", fill="red") + # 21 is filled circle
  geom_point(position=pd, size=3, shape=2, colour="black", fill="white") + # 21 is filled circle
  scale_x_discrete("Item", limits=1:10) +
  ylab(expression(nu)) +
  scale_colour_hue(name="Group 2",    # Legend label, use darker colors
                    breaks=c("1", "2"),
                    labels=c("Group 1", "Group2"),
                    l=40) +                    # Use darker colors, lightness=40
  #ggtitle("Posterior means and SD's of cross-covariance parameters") +
  # expand_limits(y=0) +                        # Expand y range
  scale_y_continuous(limits=c(-0.6, 0.3), breaks=c(-0.50, -0.25, 0, 0.25)) +      # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position="bottom")               # Position legend in bottom right


