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
print(monomorphic_variants)
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
print(genotype_counts_matrix)
