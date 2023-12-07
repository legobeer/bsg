# Question 1

X <- matrix(c(112, 278, 150, 206, 348, 150), nrow = 2, byrow = TRUE)
colnames(X) <- c("AA","Aa","aa")
rownames(X) <- c("Cases","Controls")

Y <- cbind(2*X[,1]+X[,2],2*X[,3]+X[,2])
colnames(Y) <- c("A","a")
chi_square_result <- chisq.test(Y,correct=FALSE)
p_value <- chi_square_result$p.value
odds_ratio <- Y[1, 1] * Y[2, 2] / (Y[1, 2] * Y[2, 1])
print(odds_ratio)


# Print the results
cat("Chi-square test results:\n")
cat("P-value:", p_value, "\n")
cat("Odds Ratio:", odds_ratio, "\n")

# Question 2 and 3

# codominant test
test <- chisq.test(X)
p_value <- test$p.value
cat("Codomminant test results:\n")
cat("P-value:", p_value, "\n")
odds_ratios <- c((X[1, 1] / (X[1, 2] + X[1, 3])) / (X[2, 1] / (X[2, 2] + X[2, 3])))
odds_ratios <- c(odds_ratios, (X[1, 2] / (X[1, 1] + X[1, 3])) / (X[2, 2] / (X[2, 1] + X[2, 3])))
odds_ratios <- c(odds_ratios, (X[1, 3] / (X[1, 1] + X[1, 2])) / (X[2, 3] / (X[2, 1] + X[2, 2])))
risk <- odds_ratios / (1 + odds_ratios)

allele_counts <- c(0, 1, 2)
plot(allele_counts, risk, type = "b", pch = 19, col = "blue", main = "Risk of Disease by Number of M Alleles", xlab = "Genotype", ylab = "Risk")

# Dominant test
Y <- cbind(X[,1],X[,2]+X[,3])
colnames(Y) <- c("AA","Aa or aa")
rownames(Y) <- c("Cases","Control")
test <- chisq.test(Y)
p_value <- test$p.value
cat("Dominant test results:\n")
cat("P-value:", p_value, "\n")

odds_ratios <- c((X[1, 1] / (X[1, 2] + X[1, 3])) / (X[2, 1] / (X[2, 2] + X[2, 3])))
odds_ratios <- c(odds_ratios, (Y[1, 2] / Y[1, 1]) / (Y[2, 2] / Y[2, 1]))
odds_ratios <- c(odds_ratios, (Y[1, 2] / Y[1, 1]) / (Y[2, 2] / Y[2, 1]))

risk <- odds_ratios / (1 + odds_ratios)

allele_counts <- c(0, 1, 2)
plot(allele_counts, risk, type = "b", pch = 19, col = "blue", main = "Risk of Disease by Number of M Alleles", xlab = "Genotype", ylab = "Risk")

# Recessive test
Y <- cbind(X[,1]+X[,2],X[,3])
colnames(Y) <- c("AA or Aa","aa")
rownames(Y) <- c("Cases","Control")
test <- chisq.test(Y)
p_value <- test$p.value
cat("Recessive test results:\n")
cat("P-value:", p_value, "\n")

odds_ratios <- c((Y[1, 1] / Y[1, 2]) / (Y[2, 1] / Y[2, 2]))
odds_ratios <- c(odds_ratios, (Y[1, 1] / Y[1, 2]) / (Y[2, 1] / Y[2, 2]))
odds_ratios <- c(odds_ratios, (Y[1, 2] / Y[1, 1]) / (Y[2, 2] / Y[2, 1]))

risk <- odds_ratios / (1 + odds_ratios)

allele_counts <- c(0, 1, 2)
plot(allele_counts, risk, type = "b", pch = 19, col = "blue", main = "Risk of Disease by Number of M Alleles", xlab = "Genotype", ylab = "Risk")

# Question 4
Cases <- c(112, 278, 150)
Controls <- c(206, 348, 150)
X <- rbind(Cases,Controls)
rownames(X) <- c("Cases","Controls")
colnames(X) <- c("AA","Aa","aa")
n <- sum(X)
cas <- rep(c(0,1,2),Cases)
con <- rep(c(0,1,2),Controls)
y <- c(rep(1,sum(Cases)), rep(0,sum(Controls)))
x <- c(cas,con)
r <- cor(x,y)
A <- n*(r^2)
pvalue <- pchisq(A,df=1,lower.tail=FALSE)

cat("Armitage Trend Tes:\n")
cat("P-value:", pvalue, "\n")
