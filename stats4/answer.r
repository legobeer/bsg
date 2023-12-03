library(MASS)
# Question 1
data <- read.table("Chr21.dat", header = FALSE, sep = " ")
num_individuals <- nrow(data) - 1
cat("Number of individuals:", num_individuals, "\n")
variant_cols_indices <- which(grepl("^rs", data[1, ], ignore.case = TRUE))
num_variant <- length(variant_cols_indices)
cat("Number of variants:", num_variant, "\n")
missing_percentage <- mean(is.na(data)) * 100
cat("Percentage of missing data:", missing_percentage, "%\n")

# Question 2
genotype_data <- data[-1, variant_cols_indices]
genotype_data[is.na(genotype_data)] <- 0
manhattan_dist_matrix <- dist(genotype_data, method = "manhattan")
submatrix <- as.matrix(manhattan_dist_matrix)[1:5, 1:5]
print(submatrix)

# Question 3

# Question 4

mds_result <- cmdscale(manhattan_dist_matrix, k = 2)
colors <- ifelse(mds_result[, 1] > 0, "red", "green")
plot(mds_result[, 1], mds_result[, 2], col = colors, pch = 16, main = "MDS Plot", xlab = "Dimension 1", ylab = "Dimension 2")
# Maybe two subpopulations
population1 <- sum(mds_result[, 1] > 0)
population2 <- sum(mds_result[, 1] < 0)
cat("Population 1 : ", population1, "%\n")
cat("Population 2 : ", population2, "%\n")


# Question 5
cat("Stress value:", isoMDS(manhattan_dist_matrix, k = 2)$stress
, "\n")
# Quite big normal we have two subpopulations.

# Question 6
observed_distances <- as.vector(manhattan_dist_matrix)
estimated_distances <- as.vector(dist(mds_result))

plot(observed_distances, estimated_distances,
     main = "Estimated vs Observed Distances",
     xlab = "Observed Distances", ylab = "Estimated Distances")

regression_model <- lm(estimated_distances ~ observed_distances)
abline(regression_model, col = "red")

cat("Coefficient of Determination (R-squared):", summary(regression_model)$r.squared, "\n")

# Question 7

manhattan_matrix <- as.matrix(manhattan_dist_matrix)
n <- nrow(manhattan_matrix)
m <- ncol(manhattan_matrix)

set.seed(12345)
random_matrix <- matrix(runif(n * m), nrow = n, ncol = m)
nonmetric_mds_result <- isoMDS(random_matrix, k = 2)
plot(nonmetric_mds_result$points[, 1], nonmetric_mds_result$points[, 2],
     main = "Non-Metric MDS Plot", xlab = "Dimension 1", ylab = "Dimension 2")
text(nonmetric_mds_result$points[, 1], nonmetric_mds_result$points[, 2], labels = rownames(nonmetric_mds_result$points))
nonmetric_mds_result$stress

# Question 8

num_runs <- 150
stress_result <- c()
for (i in 1:num_runs) {
  set.seed(i)  
  random_matrix <- matrix(runif(n * m), nrow = n, ncol = m)
  mds_result <- isoMDS(random_matrix, k = 2, trace=FALSE)
  
  stress_result <- c(stress_result, mds_result$stress)
}

plot(stress_result)

# Question 9

set.seed(12345)
stress_result <- c()
random_matrix <- matrix(runif(n * m), nrow = n, ncol = m)
for (i in 1:50) {
  set.seed(i)  
  mds_result <- isoMDS(random_matrix, k = i, trace=FALSE)
  
  stress_result <- c(stress_result, mds_result$stress)
}

plot(stress_result)

