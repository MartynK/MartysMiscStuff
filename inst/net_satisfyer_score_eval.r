# Sample data
set.seed(123)
data <- sample(0:4, size = 20, replace = TRUE, prob = c(0.05, 0.1, 0.15, 0.2, 0.5))
n <- length(data)

# Proportions
x_top <- sum(data == 4)
p_top <- x_top / n

x_bottom <- sum(data == 0)
p_bottom <- x_bottom / n

# Net Score
net_score <- p_top - p_bottom

# Standard Errors
se_top <- sqrt(p_top * (1 - p_top) / n)
se_bottom <- sqrt(p_bottom * (1 - p_bottom) / n)
se_diff <- sqrt(se_top^2 + se_bottom^2)

# Confidence Interval
z <- 1.96  # For 95% confidence
ci_lower <- net_score - z * se_diff
ci_upper <- net_score + z * se_diff

# Results
print(paste("Net Score:", round(net_score * 100, 2), "%"))
print(paste("95% Confidence Interval for Net Score: [",
            round(ci_lower * 100, 2), "%, ",
            round(ci_upper * 100, 2), "% ]", sep = ""))
