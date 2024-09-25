# Designate file path
file_path_coef_bias <- "simulation_results_bias.csv"

# Open dataframe
if (file.exists(file_path_coef_bias)) {
  coef_bias <- read.csv(file_path_coef_bias, stringsAsFactors = FALSE)
} else {
  coef_bias <- data.frame(
    fit = character(),
    n = numeric(),
    sigma = numeric(),
    naive.oracle_1 = numeric(),
    naive.oracle_2 = numeric(),
    naive.ME_1 = numeric(),
    naive.ME_2 = numeric(),
    prop.ME_1 = numeric(),
    prop.ME_2 = numeric(),
    stringsAsFactors = FALSE
  )
}

# Create row of data with parameters
new_row <- data.frame(
  fit = "kNN",
  n = 200,
  sigma = 0.04,
  naive.oracle_1 = b.naive.oracle[,1],
  naive.oracle_2 = b.naive.oracle[,2],
  naive.ME_1 = b.naive.ME[,1],
  naive.ME_2 = b.naive.ME[,2],
  prop.ME_1 = b.prop.ME[,1],
  prop.ME_2 = b.prop.ME[,2]
)

# Append the results to the existing data frame
coef_bias <- rbind(coef_bias, new_row)

# Update csv file
write.csv(coef_bias, file = file_path_coef_bias, row.names = FALSE)

#parameters <- expand.grid(
#  fit = c("kernel","kNN"),
#  n = c(50,100,150),
#  sigma = c(0.1,0.2^2,0.1^2)
#)