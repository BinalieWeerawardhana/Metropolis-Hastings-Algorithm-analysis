#PART 1#

#(a)#

#Defining the probability density function.
pdf_ <- function(x) {
  0.5 * exp(-abs(x))
}

# Implementing the "Random-walk-Metropolis algorithm":
r_w_metropolis <- function(N, s, x0) {  # N is the no. of generated samples, 's' is the standard-deviation, 'x0' is the initial value.
  
  x <- numeric(N) # Creating a numeric vector of length N, filled with NA values.
  
  x[1] <- x0  # assigning x0 as the first element of vector x.
  
  for (i in 2:N) {  # generating samples starting from the second sample.
    
    x_star_ <- rnorm(1, mean = x[i-1], sd = s) # Simulating a random num 'x*' from a Normal-Distribution with mean 'xi-1' and standard-deviation 's'.
    
    acc_ratio <- f(x_star_) / f(x[i-1]) # Calculating the 'acceptance ratio', where x* is the proposed sample and x[i-1] is the current sample.
    
    u_ <- runif(1) # Producing a single random number from a Uniform-Distribution.
    
    if (log(u_) < log(acc_ratio)) { # Deciding whether the proposed sample x* should be accepted or rejected.
      x[i] <- x_star_
    } else {
      x[i] <- x[i-1]
    }
  }
  return(x)
}

# Parameters:
N <- 10000
s <- 1
x0 <- 0

#Generating samples using the 'Random-walk-Metropolis algorithm'.
samples_ <- r_w_metropolis(N, s, x0)

# Plotting a Histogram and Kernel Density plot.

hist(samples_, breaks = 40, probability = TRUE, col = "lightblue", border = "black",  
     xlab = "x", ylab = "f(x)", main = "")

curve(f, from = -10, to = 10, col = "forestgreen", add = TRUE, lwd = 3)
lines(density(samples), col = "red", lty = 2, lwd = 3) 

legend("topright", legend = c("Histogram", "f(x)", "Kernel Density Estimation"), 
       col = c("lightblue", "forestgreen", "red"), lwd = c(NA, 2, 2.5), lty = c(NA, 1, 2),
       inset = 0.05, cex = 0.65)


# Calculating and displaying the sample mean and standard-deviation.
sample_mean <- mean(samples_)
cat("Sample Mean:", sample_mean, "\n")

sample_std <- sd(samples_)
cat("Sample Standard Deviation:", sample_std, "\n")


#(b)#

# The R^ value is denoted as "Rb" in the codes below.

# J represents the number of chains to be generated(for j = 1,2,...J).
# N represents the length of each chain/number of iterations in the Markov Chain.

# Defining a function to calculate the R^ value:
calculating_Rb <- function(N, s, J, x0) {
  
  # Replicating the random-walk-metropolis function N times, and storing them in a list "chains".   
  chains_ <- replicate(J, r_w_metropolis(N, s, x0), simplify = FALSE) 
  
  Mj <- sapply(chains_, mean) # Calculating the mean of the 'jth' chain.
  
  Vj <- sapply(chains_, function(chain_) mean((chain_ - mean(chain_))^2))  # Calculating the 'within sample variance' of the jth chain.
  
  W <- mean(Vj) # Calculating the 'overall within sample variance' across all chains.
  
  M <- mean(Mj) # Calculating the 'overall sample mean' across all chains.
  
  B <- mean((Mj - M)^2) # Calculating the 'between sample variance'.
  
  Rb <- sqrt((B + W) / W)  # Calculating R^ value.
  return(Rb)
}

# Running multiple chains using the random walk Metropolis algorithm and collecting R^ values into chains array.
run_chains_ <- function(N, s, J, x0) {
  chains_ <- array(dim = c(N, J))
  for (j in 1:J) {
    chains_[, j] <- r_w_metropolis(N, s, x0)
  }
  return(chains_)
}

# Parameters for the R^ value:
N <- 2000   # No. of iterations
s <- 0.001  # Standard-deviation
J <- 4      # No. of chains
x0 <- 0     # Initial value of each chain

# Executing multiple chains using the above parameters and storing them in the chains array.
chains_ <- run_chains_(N, s, J, x0)

# Calculating and displaying the required R^ value
Rb <- calculating_Rb(N, s, J, x0)
cat("Rb value:", Rb, "\n")

# Defining a function to run multiple chains and calculate R^ values for a given s-value.
Rb_for_given_s <- function(N, s, J, x0) {
  chains <- run_chains(N, s, J, x0)  # Generate chains using run_chains function
  return(calculating_Rb(N, s, J, x0))   # Calculate R^ using calculating_Rb function
}

#'Creating an array of 100 evenly spaced s-values, between 0.001 and 1 inclusive'.
s_values_ <- seq(0.001, 1, length.out = 100)

Rb_values_ <- numeric(length(s_values_))
for (i in seq_along(s_values_)) {
  Rb_values_[i] <- Rb_for_given_s(N, s_values_[i], J, x0)
}

# Plotting 'R^' values over a grid of s-values in the range between 0.001 and 1'.

library(ggplot2)
options(repr.plot.width=12, repr.plot.height=6)
plot(s_values_, Rb_values_, type = "l", 
     xlab = "s", ylab = "R^",
     main = "R^ values over a grid of s values",
     col = "blue", lwd = 2)
grid()
legend("topright", legend = "Rb values", col = "blue", lty = 1, lwd = 2)
