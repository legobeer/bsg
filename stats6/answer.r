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