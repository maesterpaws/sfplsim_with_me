library(ggplot2)

# Set seed for reproducibility
set.seed(378)

# Sample size
n <- 6

# Functional data discretized domain
t <- seq(0,1,length.out=100)

# Generate functional data
a <- ifelse(runif(n) < 0.5, runif(n, 5, 10), runif(n, 20, 20.5))
b <- ifelse(runif(n) < 0.5, runif(n, 5, 10), runif(n, 20, 20.5))
c <- ifelse(runif(n) < 0.5, runif(n, 5, 10), runif(n, 20, 20.5))
x <- t(sapply(1:n, function(j) {
  a[j]*cos(2*pi*t) + b[j]*sin(3*pi*t) + 4*c[j]*(t-0.2)^2
}))

# Reshape the data for ggplot2
df <- data.frame(t = rep(t, n), Curve = rep(1:n, each = length(t)), Value = as.vector(x))

# Plot samples
ggplot(df, aes(x = t, y = Value, group = Curve, color = as.factor(Curve))) +
  geom_line() +
  labs(x = "t", y = expression(chi(t)), title = "Functional Data Curves") +
  theme_minimal() +
  theme(legend.position = "none")
