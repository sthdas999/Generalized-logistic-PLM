power.r2 = c(0.050 , 0.051 , 0.065 , 0.078 , 0.101 , 0.126 , 0.168 , 0.211 , 0.259 , 0.325 , 0.385 , 0.441 , 0.525 , 0.588 , 0.649 , 0.712 , 0.771 , 0.829 , 0.872 , 0.898 , 0.921 , 0.942 , 0.959 , 0.968 , 0.976 , 0.986 , 0.988 , 0.988 , 0.994 , 0.996 , 0.997)

power.r3 = c(0.050 , 0.054 , 0.059 , 0.078 , 0.103 , 0.129 , 0.164 , 0.212 , 0.265 , 0.323 , 0.385 , 0.461 , 0.520 , 0.611 , 0.692 , 0.758 , 0.803 , 0.849 , 0.892 , 0.923 , 0.946 , 0.964 , 0.974 , 0.986 , 0.990 , 0.991 , 0.995 , 0.997 , 0.997 , 0.999, 1.000)

power.r4 = c(0.050 , 0.054 , 0.066 , 0.085 , 0.113 , 0.147 , 0.187 , 0.248 , 0.304 , 0.369 , 0.456 , 0.531 , 0.611 , 0.666 , 0.732 , 0.777 , 0.822 , 0.869 , 0.900 , 0.926 , 0.9607 , 0.969 , 0.980 , 0.986 , 0.992 , 0.995 , 0.997 , 0.997 , 0.998, 1.000, 1.000)

power.r5 = c(0.050 , 0.059 , 0.067 , 0.088 , 0.116 , 0.155 , 0.207 , 0.270 , 0.342 , 0.403 , 0.480 , 0.540 , 0.625 , 0.677 , 0.749 , 0.798 , 0.847 , 0.881 , 0.914 , 0.938 , 0.960 , 0.973 , 0.982 , 0.991 , 0.998, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000)
theta = 0:30

par(mar = c(5.4, 4, 2, 2.2), xpd = TRUE)
plot(theta, power.r2, xlim=c(theta[1],theta[length(theta)]), ylim=c(0,1), type="l", pch=3, col="red", xlab=expression(theta), ylab=expression(Power))
# Add a line
lines(theta, power.r3, xlim=c(theta[1],theta[length(theta)]), ylim=c(0,1), pch=3, col="blue", type="l")
# Add a line
lines(theta, power.r4, xlim=c(theta[1],theta[length(theta)]), ylim=c(0,1), pch=3, col="green", type="l")
# Add a line
lines(theta, power.r5, xlim=c(theta[1],theta[length(theta)]), ylim=c(0,1), pch=3, col="black", type="l")
# Add a legend
legend("bottomright",
       c(expression(paste('r==2'),paste('r==3'),paste('r==4'),paste('r==5'))),
       fill=c("red","blue","green","black"))