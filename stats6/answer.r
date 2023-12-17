library(data.table)

# Question 1
data <- fread("YRI6.raw", header = FALSE, sep = " ")

num_individuals <- nrow(data)
num_snps <- ncol(data) - 6

percentage_missing <- mean(is.na(data)) * 100

cat("Number of individuals:", num_individuals, "\n")
cat("Number of SNPs:", num_snps, "\n")
cat("Percentage of missing data:", percentage_missing, "%\n")

# Question 2

