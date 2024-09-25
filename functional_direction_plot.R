library(ggplot2)
library(splines)

# Set parameters
theta_0 <- c(0, 1.741539, 0, 1.741539, -1.741539, -1.741539)
a <- 0
b <- 1
nknot.theta <- 3
order.Bspline <- 3

t <- seq(a,b,length.out=100)
Knot.theta <- seq(a, b, length = nknot.theta + 2)[ - c(1, nknot.theta + 2)]
delta.theta <- sort(c(rep(c(a, b),order.Bspline), Knot.theta))
Bspline.theta <- splineDesign(delta.theta,t,order.Bspline)
theta.rec <- Bspline.theta%*%theta_0 
theta_df <- data.frame(t, theta.rec)

# Plot functional direction
ggplot(theta_df, aes(x = t, y = theta.rec)) +
  geom_line(linewidth = 1.5, color = "purple") +
  labs(x = "t", y = expression(theta[0](t)), title = "Functional Direction of Projection") +
  theme_minimal()
