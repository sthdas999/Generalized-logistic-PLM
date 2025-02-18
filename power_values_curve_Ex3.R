power.r2 = c(0.050 , 0.051 , 0.057 , 0.068 , 0.089 , 0.093 , 0.118 , 0.140 , 0.172 , 0.197 , 0.246 , 0.293 , 0.330 , 0.377 , 0.419 , 0.459 , 0.500 , 0.546 , 0.584 , 0.626 , 0.666 , 0.700 , 0.737 , 0.771 , 0.793 , 0.819 , 0.843 , 0.864 , 0.894 , 0.909 , 0.925)

power.r3 = c(0.050 , 0.052 , 0.060 , 0.075 , 0.096 , 0.131 , 0.162 , 0.200 , 0.249 , 0.313 , 0.378 , 0.443 , 0.510 , 0.597 , 0.652 , 0.712 , 0.764 , 0.810 , 0.847 , 0.881 , 0.911 , 0.931 , 0.948 , 0.966 , 0.973 , 0.979 , 0.987 , 0.991 , 0.993 , 0.993 , 0.995)

power.r4 = c(0.050 , 0.053 , 0.063 , 0.077 , 0.099 , 0.134 , 0.172 , 0.212 , 0.259 , 0.324 , 0.389 , 0.465 , 0.531 , 0.598 , 0.661 , 0.714 , 0.770 , 0.811 , 0.863 , 0.890 , 0.923 , 0.942 , 0.951 , 0.970 , 0.980 , 0.985 , 0.989 , 0.994 , 0.995 , 0.996 , 0.999)

power.r5 = c(0.050 , 0.056 , 0.069 , 0.079 , 0.106 , 0.139 , 0.184 , 0.254 , 0.301 , 0.354 , 0.406 , 0.487 , 0.549 , 0.602 , 0.670 , 0.734 , 0.788 , 0.849 , 0.886 , 0.925 , 0.948 , 0.962 , 0.972 , 0.986 , 0.995 , 0.997 , 0.998 , 0.999 , 0.999 , 0.999 , 1.000)

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