# Linkage Disequilibrium
library(genetics)
library(HardyWeinberg)


# Question 1
data <- read.table("FOXP2/FOXP2.dat", header = TRUE, sep = " ")
num_individuals <- nrow(data)
num_snps <- ncol(data)

cat("Number of individuals:", num_individuals, "\n")
cat("Number of SNPs:", num_snps, "\n")

missing_percentage <- mean(is.na(data)) * 100
cat("Percentage of missing data:", missing_percentage, "%\n")

# Question 2
g1 <- genotype(data$rs34684677)
g2 <- genotype(data$rs2894715)
LD(g1,g2)

# Question 4

bim_file <- "FOXP2/FOXP2.bim"
bim_data <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)

genotype_matrix <- data[,-1]
for (col in colnames(genotype_matrix)) {
  genotype_matrix[[col]] <- genotype(genotype_matrix[[col]])
}
calculate_MAF <- function(genotype_counts) {
  allele1 <- 2 * genotype_counts[[1]] + genotype_counts[[2]]
  allele2 <- 2 * genotype_counts[[3]] + genotype_counts[[2]]
  return (min(allele1, allele2) / (allele2 + allele1))
}

filter_index_MAF <- c()
res <- 0
for (i in 1:ncol(genotype_matrix)) {
  alleles <- c(paste(bim_data$V5[[i]], "/", bim_data$V6[[i]], sep = ""))
  count_gentoype <- MakeCounts(genotype_matrix[[i]], alleles = alleles, sep="/")
  MAF <- calculate_MAF(count_gentoype)
  if (MAF < 0.35) {
    filter_index_MAF <- c(filter_index_MAF, i)
  }
  p_value <- HWE.chisq(genotype_matrix[[i]])
  p_value <- p_value$p.value
  if (p_value < 0.05) {
    res <- res + 1
  }
}
print(res)

# Question 5

ld_matrix <- LD(genotype_matrix)
r2_matrix <- ld_matrix$`R^2`
r2_matrix[is.na(r2_matrix)] <- 1
heatmap(r2_matrix, main = "LD Heatmap", xlab = "Marker Index", ylab = "Marker Index")

# Question 6

# Filter all the variants with a MAF below 0.35
r2_matrix <- r2_matrix[-filter_index_MAF, -filter_index_MAF]
heatmap(r2_matrix, main = "LD Heatmap only with MAF >= 0.35", xlab = "Marker Index", ylab = "Marker Index")

# Question 7

# Haplotype estimation

# Question 1

data <- read.table("APOE/APOE.dat", header = TRUE, sep = " ")
cat("Number of individuals:", num_individuals, "\n")
cat("Number of SNPs:", num_snps, "\n")

missing_percentage <- mean(is.na(data)) * 100
cat("Percentage of missing data:", missing_percentage, "%\n")

# Question 2

