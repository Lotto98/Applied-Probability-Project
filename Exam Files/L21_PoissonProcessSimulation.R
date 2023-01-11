# Simulating a Poisson process path:
# (The most inefficient but intuitive way)
lambda= 0.8
t=0
n=0
jump_time=rexp(1,lambda)
t=c(t,jump_time)
n=c(n,1)

# Slightly more efficient
t=rep(0,6)
n=rep(0,6)
n_jumps=5
for(i in (1:n_jumps)){
  t[i+1]=rexp(1,lambda)+t[i]
  n[i+1]=1+n[i]
}

# More efficient
set.seed(876)
n = (0:n_jumps)
t = cumsum(c(0,rexp(n_jumps,lambda)))

# ASSIGNMENT 1: How to plot this path? 
# (How do you plot a step function?)

  # Simplest solution:
  plot(t, n, type ="s")

  # Check package ggplot for "nicer" plotting options
  
# ASSIGNMENT 2: Simulate and plot a path from this 
# Poisson process until time T=10
  
  # Solution:
  set.seed(876)
  T=10
  n_jumps=lambda*T
  n = (0:n_jumps)
  t = cumsum(c(0,rexp(n_jumps,lambda)))
  while (max(t) <= T){
    n = c(n, max(n)+1)
    t = c(t, max(t) + rexp(1,lambda))
  }
  if (max(t) > T){
    t[length(t)]=T
    n[length(n)]=n[length(t)]-1
  }
  plot(t, n, type="s")