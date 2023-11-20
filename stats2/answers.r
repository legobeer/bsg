# Question 1
library(HardyWeinberg)
library(data.table)

file <- "TSIChr22v4.raw"
data <- fread(file, header = FALSE, sep = " ")
genetic_data <- data[, 7:ncol(data)]
first_elements <- sapply(genetic_data, function(x) x[1])
rs_columns <- grepl("^rs", first_elements)
genetic_data_rs <- genetic_data[, ..rs_columns, with = FALSE]
num_variants <- ncol(genetic_data_rs)
missing_percentage <- mean(is.na(genetic_data_rs)) * 100

cat("Number of variants in the database:", num_variants, "\n")
cat("Percentage of missing data:", missing_percentage, "%\n")

# Question 2

monomorphic_variants <- sapply(genetic_data_rs, function(col) {
  cleaned_col <- col[-1]  # Remove the first element (variant name)
  length(unique(cleaned_col[!is.na(cleaned_col)])) == 1
})
percentage_monomorphic <- (sum(monomorphic_variants) / length(monomorphic_variants)) * 100
non_monomorphic_data <- genetic_data_rs[, !monomorphic_variants, with = FALSE]
num_remaining_variants <- ncol(non_monomorphic_data)
cat("Percentage of monomorphic variants:", percentage_monomorphic, "%\n")
cat("Number of remaining variants:", num_remaining_variants, "\n")

# Question 3

variant_column_index <- which(genetic_data_rs[1, ] == "rs587756191_T")
variant_column <- genetic_data_rs[[variant_column_index]]

variant_column <- factor(variant_column, levels = c(0, 1, 2), labels = c("AA", "AB", "BB"))
print(table(variant_column))

# Create a table with counts
chi_square_result <- HWChisq(table(variant_column))
print(chi_square_result)

# Chi-square test with continuity correction
chi_square_correction_result <- HWChisq(table(variant_column))
print(chi_square_correction_result)

# Exact test
exact_test_result <- HWExact(table(variant_column))
print(exact_test_result)

# Permutation test
permutation_test_result <- HWPerm(table(variant_column))
print(permutation_test_result)

# Question 4 
genetic_data_rs <- non_monomorphic_data
genotype_counts_matrix <- matrix(NA, ncol = 3, nrow = ncol(genetic_data_rs))
for (i in 1:ncol(genetic_data_rs)) {
  variant_column <- genetic_data_rs[[i]]
  genotype_counts <- table(variant_column)
  
  genotype_counts <- genotype_counts[c("0", "1", "2")]
  
  genotype_counts_matrix[i, ] <- as.numeric(genotype_counts)
}

colnames(genotype_counts_matrix) <- c("AA", "AB", "BB")

# Question 5
genetic_data_no_header <- genetic_data_rs[-1, ]

hwe_results <- c()

# Loop through each SNP column
for (col_name in names(genetic_data_no_header)) {
  # Extract the SNP column
  variant_column <- as.integer(genetic_data_no_header[[col_name]])
  variant_column <- factor(variant_column, levels = c(0, 1, 2), labels = c("AA", "AB", "BB"))
  
  # Create a matrix X with genotype counts for the SNP
  X <- matrix(table(variant_column), ncol = 3, byrow = TRUE)
  
  # Apply HWExactStats to the matrix X
  hwe_result <- HWExactStats(X, x.linked = FALSE, plinkcode = TRUE, midp = FALSE)
  
  # Collect the p-value for each SNP
  hwe_results <- c(hwe_results, hwe_result)
}

# Calculate the percentage of significant SNPs (p-value < 0.05)
significant_snps_percentage <- mean(hwe_results < 0.05) * 100

cat("Percentage of significant SNPs at alpha = 0.05:", significant_snps_percentage, "%\n")

# Question 6

most_significant_index <- which.min(hwe_results)

most_significant_snp_value <- genetic_data_rs[[most_significant_index]][1]
most_significant_p_value <- hwe_results[most_significant_index]

cat("Most Significant SNP:", most_significant_snp_value, "\n")
cat("P-value for Most Significant SNP:", most_significant_p_value, "\n")

# Question 7

inbreeding_coefficients <- c()

# Loop through each SNP column
for (col_name in names(genetic_data_no_header)) {
  # Extract the SNP column
  variant_column <- as.integer(genetic_data_no_header[[col_name]])
  variant_column <- factor(variant_column, levels = c(0, 1, 2), labels = c("AA", "AB", "BB"))
  
  # Create a matrix X with genotype counts for the SNP
  X <- matrix(table(variant_column), ncol = 3, byrow = TRUE)
  
  # Apply HWf to calculate the inbreeding coefficient
  f <- HWf(X)
  
  # Store the inbreeding coefficient
  inbreeding_coefficients <- c(inbreeding_coefficients, f)
}

# Make a histogram of f
hist(inbreeding_coefficients, main = "Histogram of Inbreeding Coefficients", xlab = "Inbreeding Coefficient (f)", col = "lightblue")

# Calculate and print descriptive statistics
mean_f <- mean(inbreeding_coefficients)
sd_f <- sd(inbreeding_coefficients)

cat("Descriptive Statistics for Inbreeding Coefficients:\n")
cat("Mean:", mean_f, "\n")
cat("Standard Deviation:", sd_f, "\n")

# Probability plot
qqnorm(inbreeding_coefficients, main = "Probability Plot of Inbreeding Coefficients")
qqline(inbreeding_coefficients)


# Question 8

