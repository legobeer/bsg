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
genomic_data <- data[2:nrow(data), 7:ncol(data)]
genomic_data[, (1:ncol(genomic_data)) := lapply(.SD, as.numeric), .SDcols = 1:ncol(genomic_data)]
shared_mean <- matrix(c(NA), nrow = num_individuals, ncol = num_individuals)
shared_sd <- matrix(c(NA), nrow = num_individuals, ncol = num_individuals)
for (i in 1:num_individuals) {
  for (j in 1:num_individuals) {
    shared <- numeric(0)
    for (k in 1:nrow(genomic_data)){
    shared <- c(shared, 2 - abs(as.matrix(genomic_data[k,i, with=FALSE]) - as.matrix(genomic_data[k,j, with=FALSE])))
    }
    mean <- mean(shared)
    sd <- sd(shared)
    shared_mean[i,j] <- mean
    shared_sd[i,j] <- sd
  }
}

print(shared_mean[1:5, 1:5])
print(shared_sd[1:5, 1:5])


# Question 3

p0 <- matrix(c(NA), nrow = num_individuals, ncol = num_individuals)
p2 <- matrix(c(NA), nrow = num_individuals, ncol = num_individuals)
m <- matrix(c(NA), nrow = num_individuals, ncol = num_individuals)
for (i in 1:num_individuals) {
  for (j in 1:num_individuals) {
    shared <- numeric(0)
    for (k in 1:nrow(genomic_data)){
      shared <- c(shared, 2 - abs(as.matrix(genomic_data[k,i, with=FALSE]) - as.matrix(genomic_data[k,j, with=FALSE])))
    }
    p0[i,j] <- sum(shared == 0)/length(shared)
    p2[i,j] <- sum(shared == 2)/length(shared)
    m[i,j] = 1 - p0[i,j] + p2[i,j]
  }
}
print(p0[1:5, 1:5])
print(p2[1:5, 1:5])
print(m[1:5, 1:5])
print(shared_mean[1:5, 1:5])
print(all.equal(shared_mean, m))

# m = 1 - p0 + p2 holds.

# Question 4

plot(shared_mean, shared_sd, main = "m vs s", xlab = "m", ylab = "s", col = "blue", pch = 16)
plot(p0, p2, main = "p0 vs p2", xlab = "p0", ylab = "p2", col = "green", pch = 16)

# Question 5


m <- shared_mean
s <- shared_sd
family_relationship <- data[, c(3, 4)]
colors <- rep("blue", nrow(data) - 1)  
parent_offspring_indices <- which(family_relationship[, 1] > 0 | family_relationship[, 2] > 0)
colors[parent_offspring_indices - 1] <- "red"
parent <- c()
for (i in 2:nrow(data)) {
  print(data[i, 2] %in% family_relationship[[1]])
  if (data[i, 2] %in% family_relationship[[1]] | data[i, 2] %in% family_relationship[[2]]) {
    parent <- c(parent, i - 1)
  }
}
print(parent)
colors[parent] <- "red"

plot(m, s, main = "Shared Mean vs Shared SD", xlab = "Shared Mean", ylab = "Shared SD", col = colors, pch = 16)
legend("bottomleft", legend = c("Unrelated", "Parent-Offspring"), col = c("blue", "red"), pch = 16, cex = 0.8)


# Question 6

library(SNPRelate)

raw_data <- read.table("YRI6.raw", header = TRUE, sep = " ")
sample_ids <- raw_data[, 2] 
snp_ids <- names(raw_data[, 7:ncol(raw_data)])
genotypes <- as.matrix(raw_data[, 7:ncol(raw_data), drop = FALSE]) 
print(dim(genotypes))
print(length(sample_ids))
print(length(snp_ids))
snpgdsCreateGeno("test.gds", genmat = t(genotypes), sample.id = sample_ids, snp.id = snp_ids)
                 
 
genofile <- snpgdsOpen("test.gds")
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
snpset.id <- unlist(unname(snpset))
ibd <- snpgdsIBDMLE(genofile,maf=0.05, missing.rate=0.05,snp.id=snp_ids, num.thread=2)                
ibd.coeff <- snpgdsIBDSelection(ibd)
family_info <- raw_data[, c(3, 4)]


colors <- ifelse(family_info[, 1] > 0 & family_info[, 2] > 0, "red", "blue")

plot(ibd.coeff$k0, ibd.coeff$k1, col = colors, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "k0", ylab = "k1", main = "YRI samples (MLE)")

legend("topright", legend = c("Unrelated", "Parent-Offspring"), col = c("blue", "red"), pch = 16, cex = 0.8)
lines(c(0,1), c(1,0), col="black", lty=2)

# Question 7

# In IBD probability plot there are pairs with very low k1 values marked as Parent-Offspring. Parent-Offspring relationships are usually associated with a k0 of 0 and a k1 > 0.5. This suggests that there are incorrectly specified family relationships in the data.

