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
colors <- rep("blue", nrow(data))  
print(data[6,2] %in% family_relationship[[1]])
print(data[6,2])
print(family_relationship[, 1])
parent_offspring_indices <- which(family_relationship[, 1] > 0 | family_relationship[, 2] > 0)
print(parent_offspring_indices)
colors[parent_offspring_indices] <- "red"

plot(m, s, main = "Shared Mean vs Shared SD", xlab = "Shared Mean", ylab = "Shared SD", col = colors, pch = 16)
legend("bottomleft", legend = c("Unrelated", "Parent-Offspring"), col = c("blue", "red"), pch = 16, cex = 0.8)


# Question 6

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SNPRelate")
