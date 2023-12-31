# First data set

# Question 1

file <-"TSICHR22RAW.raw"
data <- read.table(file, header=FALSE, sep=" ")
genetic_data <- data[, 7:ncol(data)]
first_elements <- sapply(genetic_data, function(x) x[1])
rs_columns <- grepl("^rs", first_elements)
genetic_data_rs <- genetic_data[, rs_columns]
print(dim(genetic_data_rs))
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
non_monomorphic_data <- genetic_data_rs[, !monomorphic_variants]
num_remaining_variants <- ncol(non_monomorphic_data)

cat("Percentage of monomorphic variants:", percentage_monomorphic, "%\n")
cat("Number of remaining variants:", num_remaining_variants, "\n")

# Question 3

column_index <- which(sapply(genetic_data_rs, function(col) col[1] == "rs8138488_C"))
rs8138488_column <- genetic_data_rs[, column_index]
genotype_counts <- table(rs8138488_column[-1])
genotype_0 <- genotype_counts["0"]
genotype_1 <- genotype_counts["1"]
genotype_2 <- genotype_counts["2"]
minor_allele_count <- min(genotype_0, genotype_1, genotype_2)
total_individuals <- length(rs8138488_column) - 1  
MAF <- minor_allele_count / (2 * total_individuals)
cat("Genotype counts for rs8138488_C:\n")
print(genotype_counts)
cat("Minor Allele Count:", minor_allele_count, "\n")
cat("Minor Allele Frequency (MAF):", MAF, "\n")

# Question 4

MAF_values <- numeric(length(genetic_data_rs)) 

for (i in 1:length(genetic_data_rs)) {
  column <- genetic_data_rs[[i]]
  allele_counts <- table(column[-1])  
  minor_allele_count <- min(allele_counts)
  total_individuals <- length(column) - 1 
  MAF <- minor_allele_count / (2 * total_individuals)  
  MAF_values[i] <- MAF
}

hist(MAF_values, breaks = 20, xlab = "Minor Allele Frequency (MAF)", main = "Histogram of Minor Allele Frequencies")

percentage_below_005 <- sum(MAF_values < 0.05) / length(MAF_values) * 100
percentage_below_001 <- sum(MAF_values < 0.01) / length(MAF_values) * 100

cat("Percentage of markers with MAF below 0.05:", percentage_below_005, "%\n")
cat("Percentage of markers with MAF below 0.01:", percentage_below_001, "%\n")

# Question 5

H0_values <- numeric(length(genetic_data_rs))  

for (i in 1:length(genetic_data_rs)) {
  column <- genetic_data_rs[[i]]
  genotypes <- column[-1]  
  num_heterozygous <- sum(genotypes == 1)  
  total_individuals <- length(genotypes)
  H0 <- num_heterozygous / total_individuals
  H0_values[i] <- H0
}

hist(H0_values, breaks = 20, xlab = "Observed Heterozygosity (H0)", main = "Histogram of Observed Heterozygosity (H0)")

# Theoretical range of variation for H0 is [0, 0.5]

# Question 6

He_values <- numeric(length(genetic_data_rs))  
for (i in 1:length(genetic_data_rs)) {
  column <- genetic_data_rs[[i]]
  genotypes <- column[-1]  # Exclude the first element (variant name)
  allele_freq <- table(genotypes) / length(genotypes)
  He <- 1 - sum(allele_freq^2)
  He_values[i] <- He
}

hist(He_values, breaks = 20, xlab = "Expected Heterozygosity (He)", main = "Histogram of Expected Heterozygosity (He)")

average_He <- mean(He_values)
cat("Average Expected Heterozygosity (He) for this database:", average_He, "\n")

# Second data set

# Question 1

library(HardyWeinberg)
data(NistSTRs)
dimensions <- dim(NistSTRs)
num_individuals <- dimensions[1]  
num_STRs <- dimensions[2]
cat("Number of individuals in the database:", num_individuals, "\n")
cat("Number of STRs in the database:", num_STRs, "\n")

# Question 2

count_alleles <- function(STR_locus) {
  unique_alleles <- unique(STR_locus)
  num_alleles <- length(unique_alleles)
  return(num_alleles)
}
num_alleles_list <- sapply(NistSTRs, count_alleles)

mean_num_alleles <- mean(num_alleles_list)
sd_num_alleles <- sd(num_alleles_list)
median_num_alleles <- median(num_alleles_list)
min_num_alleles <- min(num_alleles_list)
max_num_alleles <- max(num_alleles_list)

cat("Descriptive statistics of the number of alleles:\n")
cat("Mean:", mean_num_alleles, "\n")
cat("Standard Deviation:", sd_num_alleles, "\n")
cat("Median:", median_num_alleles, "\n")
cat("Minimum:", min_num_alleles, "\n")
cat("Maximum:", max_num_alleles, "\n")

# Question 3

table_num_alleles <- table(num_alleles_list)

barplot(table_num_alleles, xlab = "Number of Alleles", ylab = "Number of STRs", main = "Number of STRs for Each Number of Alleles")

most_common_alleles <- names(table_num_alleles)[which.max(table_num_alleles)]
cat("The most common number of alleles for an STR is:", most_common_alleles, "\n")

# Question 4

calculate_He <- function(STR_locus) {
  unique_alleles <- unique(STR_locus)
  allele_freq <- table(STR_locus) / length(STR_locus)
  He <- 1 - sum((allele_freq / sum(allele_freq))^2)
  return(He)
}
He_values <- sapply(NistSTRs, calculate_He)

hist(He_values, breaks = 20, xlab = "Expected Heterozygosity (He)", main = "Histogram of Expected Heterozygosity (He) across STRs")

average_He <- mean(He_values)
cat("Average Expected Heterozygosity over all STRs:", average_He, "\n")

# Question 5

calculate_Ho <- function(STR_locus) {
  num_unique_alleles <- length(unique(STR_locus))
  total_individuals <- length(STR_locus)
  Ho <- 1 - (num_unique_alleles - 1) / (2 * total_individuals)
  return(Ho)
}
Ho_values <- sapply(NistSTRs, calculate_Ho)

heterozygosity_data <- data.frame(He = He_values, Ho = Ho_values)
print(heterozygosity_data)

plot(heterozygosity_data$He, heterozygosity_data$Ho, xlab = "Expected Heterozygosity (He)", ylab = "Observed Heterozygosity (Ho)", 
     main = "Observed vs. Expected Heterozygosity for All STRs")

