####################
# Example: Simulating a single path from a 3-state HCTMC up to time T=10, starting at X(0)=1 (Last example in CTMC.pdf)
#####################

# Set time frame
T<-10
# Define the generator
G <- matrix(c(-4, 2, 2, 1, -2, 1, 0.5, 1, -1.5),ncol=3, byrow=T)
G
# Calculate the (conditional) transition probability matrix
P_tilde <- -G / diag(G)
diag(P_tilde) <- 0
P_tilde

# Initialize
set.seed(9878) # Fix random number generator (for reproducibility)
n <- 0 # Counter for total number of jumps
t <- 0 # Current time
times <- 0 # Vector of jump times
x_t <- 1 # Current state: initializing at X(0)=1
x <- x_t # Vector of visited states

# Simulate the first holding time:
t <- rexp(1, -G[x_t, x_t])

# Simulate sequentially until sopping rule:
while (t <= T){
  n <- n+1
  times <- c(times,t)
  x_t <- sample(1:3, 1, prob = P_tilde[x_t,])
  x <- c(x, x_t)
  t <- t + rexp(1, -G[x_t, x_t])
}

# The last jump (after time T) was never added, but we need to add the final time (T=10) and corresponding process value, to complete the plot
times <- c(times,T)
x <- c(x,x[n+1])

# Plot the path:
plot(times, x, type="s",col="blue", yaxt="n")
axis(side=2, at=1:3)

# Alternatively, if we want to see the jump times in the plot:
plot(times, x, type="s",col="blue", yaxt="n", xaxt="n")
axis(side=2, at=1:3)
axis(side=1, at=round(times,2))
# However, depending on the number of jumps and the scale of the plot, this may not be convenient 

################
# Beyond the simulation:
# We can also find the stationary distribution of the process
################

# Define the matrix A by substituting the last column of G with ones and then transposing:
A <- G
A[,3] <- rep(1,3)
A <- t(A)
A

#Define the vector b = e_3
b <- c(rep(0,3-1),1)
b

# Find the stationary distribution by solving the system A?? = b:
p <- solve(A, b)
p
# WARNING: R uses the name "pi" to denote the famous irrational number, so avoid using this name or you will overwrite the constant!
pi

# Notice that, from the vector of jump times, we can recover the n+1 sojourn times:
sojourn <- times[2:(n+2)]-times[1:(n+1)]

# We can then look at the observed avgerage ammount of time spent at each state i (remember to ignor the last value of x, since we don't have the complete sojourn time, since we stopped observing the process at time T=10)
prop_time=rep(0,3)
for(i in 1:3){
  prop_time[i]=sum(sojourn[x[1:(n+1)]==i])/T
}
prop_time

# Estimate expected length of a sojourn at state i
expected_sojourn=rep(0,3)
for(i in 1:3){
  expected_sojourn[i]=mean(sojourn[x[1:(n+1)]==i])
}

# And the rates for the sojourn times (negatives of the elements of the diagonal in G)
estimated_diagG = 1/expected_sojourn

# We could measure our estimation error, for example, through the square root of the sum of squared errors
sqrt(mean((estimated_diagG-expected_sojourn)^2))

# What happens if T is larger? For example T=50 or T=100?

