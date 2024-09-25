# Designate file path
file_path_response <- "simulation_results_response_v2.csv"
file_path_coefficients <- "simulation_results_coefficients_v2.csv"

# Open dataframe
if (file.exists(file_path_response)) {
  MSE_response <- read.csv(file_path_response, stringsAsFactors = FALSE)
} else {
  MSE_response <- data.frame(
    fit = character(),
    n = numeric(),
    sigma = numeric(),
    naive.oracle = numeric(),
    naive.ME = numeric(),
    prop.ME = numeric(),
    stringsAsFactors = FALSE
  )
}

if (file.exists(file_path_coefficients)) {
  MSE_coefficients <- read.csv(file_path_coefficients, stringsAsFactors = FALSE)
} else {
  MSE_coefficients <- data.frame(
    fit = character(),
    n = numeric(),
    sigma = numeric(),
    naive.oracle = numeric(),
    naive.ME = numeric(),
    prop.ME = numeric(),
    stringsAsFactors = FALSE
  )
}

# Create row of data with parameters
response_row <- data.frame(
  fit = "kNN",
  n = 200,
  sigma = 0.1^2,
  naive.oracle = avg_MSE[1],
  naive.ME = avg_MSE[2],
  prop.ME = avg_MSE[3]
)

coefficients_row <- data.frame(
  fit = "kNN",
  n = 200,
  sigma = 0.1^2,
  naive.oracle = avg_MSE_b[1],
  naive.ME = avg_MSE_b[2],
  prop.ME = avg_MSE_b[3]
)

# Append the results to the existing data frame
MSE_response <- rbind(MSE_response, response_row)
MSE_coefficients <- rbind(MSE_coefficients, coefficients_row)

# Update csv file
write.csv(MSE_response, file = file_path_response, row.names = FALSE)
write.csv(MSE_coefficients, file = file_path_coefficients, row.names = FALSE)

#parameters <- expand.grid(
#  fit = c("kernel","kNN"),
#  n = c(50,100,150),
#  sigma = c(0.1,0.2^2,0.1^2)
#)