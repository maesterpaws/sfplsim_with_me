library(ggplot2)
library(reshape2)

# Import data
coef_bias <- read.csv("simulation_results_bias.csv", stringsAsFactors = FALSE)

# Filter dataframe for sigma value
coef_bias_filtered <- subset(coef_bias, sigma == 0.01)

# Select columns of one coefficient
coef_bias_filtered <- coef_bias_filtered[, c("fit", "n", "sigma", 
                                             "naive.oracle_1", 
                                             "naive.ME_1", 
                                             "prop.ME_1")]

# Rename columns
colnames(coef_bias_filtered) <- c("fit", "n", "sigma", 
                                  "ORACLE", "NAIVE", "PROP")

# Reshape from wide to long format
coef_bias_long <- melt(coef_bias_filtered, id.vars = c("fit", "n", "sigma"),
                       variable.name = "Group", value.name = "Bias")

# Create plot
ggplot(coef_bias_long, aes(x = interaction(Group, n), y = Bias, fill = fit)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  labs(title = expression("Bias for" ~ beta[0][1] ~ "where" ~ sigma^2 ~ "= 0.1"^2),
       x = "", y = "Bias") +
  scale_x_discrete(labels = rep(c("n=50", "n=100", "n=200"), 3)) +
  facet_wrap(~ Group, scales = "free_x", strip.position = "bottom") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")
