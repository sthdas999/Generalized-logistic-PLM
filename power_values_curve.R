power.r5 = c(0.050 , 0.054 , 0.059 , 0.078 , 0.103 , 0.129 , 0.164 , 0.212 , 0.265 , 0.323 , 0.385 , 0.461 , 0.520 , 0.611 , 0.692 , 0.758 , 0.803 , 0.849 , 0.892 , 0.923 , 0.946 , 0.969 , 0.980 , 0.991 , 0.998 , 1.000 , 1.000 , 1.000 , 1.000 , 1.000 , 1.000)
power.r4 = c(0.050 , 0.054 , 0.065 , 0.078 , 0.101 , 0.126 , 0.168 , 0.211 , 0.259 , 0.325 , 0.385 , 0.441 , 0.525 , 0.588 , 0.649 , 0.712 , 0.771 , 0.829 , 0.872 , 0.898 , 0.921 , 0.942 , 0.959 , 0.968 , 0.976 , 0.986 , 0.988 , 0.988 , 0.994 , 0.997 , 0.997)
power.r3 = c(0.050 , 0.052 , 0.059 , 0.060 , 0.070 , 0.079 , 0.092 , 0.114 , 0.132 , 0.157 , 0.194 , 0.209 , 0.239 , 0.280 , 0.322 , 0.365 , 0.410 , 0.462 , 0.504 , 0.549 , 0.589 , 0.620 , 0.663 , 0.693 , 0.731 , 0.770 , 0.805 , 0.832 , 0.864 , 0.887 , 0.908)
power.r2 = c(0.050 , 0.053 , 0.056 , 0.069 , 0.075 , 0.088 , 0.103 , 0.120 , 0.150 , 0.173 , 0.203 , 0.238 , 0.264 , 0.293 , 0.325 , 0.362 , 0.388 , 0.427 , 0.470 , 0.493 , 0.529 , 0.568 , 0.611 , 0.644 , 0.689 , 0.718 , 0.747 , 0.771 , 0.801 , 0.826 , 0.847)

theta = 0:30

par(mar = c(5.4, 4, 2, 2.2), xpd = TRUE)
plot(theta, power.r2, xlim=c(theta[1],theta[length(theta)]), ylim=c(0,1), type="l", pch=3, col="red", xlab=expression(theta), ylab=expression(Power))
# Add a line
lines(theta, power.r3, xlim=c(theta[1],theta[length(theta)]), ylim=c(0,1), pch=3, col="black", type="l")
# Add a line
lines(theta, power.r4, xlim=c(theta[1],theta[length(theta)]), ylim=c(0,1), pch=3, col="blue", type="l")
# Add a line
lines(theta, power.r5, xlim=c(theta[1],theta[length(theta)]), ylim=c(0,1), pch=3, col="green", type="l")
# Add a legend
legend("bottomright",
       c(expression(paste('r=2'),paste('r=3'),paste('r=4'),paste('r=5'))),
       fill=c("red","black","blue","green"))