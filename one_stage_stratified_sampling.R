library(stratallo)

inverse_num <- function(x){
  return(1/x)
}

target <- function(x,A){
  return(sum(1/x * A^2))
}

gen_prop_box <- function(x){
  index_move_1 <- ceiling(runif(1, min=0, max=length(x)))
  index_move_2 <- ceiling(runif(1, min=0, max=length(x)))
  while(index_move_1 == index_move_2){
    index_move_1 <- ceiling(runif(1, min=0, max=length(x)))
    index_move_2 <- ceiling(runif(1, min=0, max=length(x)))
  }
  if(index_move_1 != index_move_2){
    if(x[index_move_1] >= m[index_move_1] & x[index_move_2] > m[index_move_2]){
      if(x[index_move_1] < M[index_move_1] & x[index_move_2] <= M[index_move_2]){
        x_prop <- x
        x_prop[index_move_1] <- x_prop[index_move_1] + 1
        x_prop[index_move_2] <- x_prop[index_move_2] - 1
      }else{
        x_prop <- x
      }
    }else{
      x_prop <- x
    }
  }else{
    x_prop <- x
  }
  return(x_prop)
}

# Generate example:
N <- c(3000, 4000, 5000, 2000) # Strata sizes.
S <- c(48, 79, 76, 16) # Standard deviations of a study variable in strata.
A <- N * S
m <- c(100, 90, 500, 50) # Lower bounds imposed on sample sizes in strata.
M <- c(300, 400, 800, 90) # Upper bounds imposed on sample sizes in strata.
n <- 1284  # Established sample size
T_schedule <- 10:1/1000 # Annealing schedule (temperature)
N_schedule <- (1:length(T_schedule))*1000 # Annealing schedule (time)
x <- c(100, 300, 800, 84) # Seed

# Tests: (all should be TRUE)
n >= sum(m) && n <= sum(M)
all(N-M>0)
all(m>0)
all(sum(x)==n)
all(x-m>=0)
all(M-x>=0)

sim_anal_box <- function(x, N, S, m, M, n, T_schedule, N_schedule){
  A <- N * S
  for(i in 1:length(N_schedule)){
    T_temp <- T_schedule[i]
    for(j in 1:N_schedule[i]){
      x_prop <- gen_prop_box(x)
      print(target(x, A))
      transit_func <- exp(-(target(x, A)-target(x_prop, A))/T_temp)
      prob_stay <- min(1, transit_func)
      prob_move <- 1 - prob_stay
      coin_flip <- runif(1)
      if (coin_flip <= prob_stay){
        x <- x
      }else{
        x <- x_prop
      }
    }
  }
  return(x)
}
x_sim <- sim_anal_box(x, N, S, m, M, n, T_schedule, N_schedule)

# Optimal solution: 
x_opt <- opt(n = n, A = A, m = m, M = M)

# Relative error: (x_opt - dimensional, rounded)

round(abs(1-(x_sim/x_opt)), 5)

# Test:
df <- pop10_mM
N <- df[, 1]
S <- df[,2]
m <- df[,3]
M <- df[,4]
n <- 3000
x <- c(60, 90, 75, 65, 400, 400, 500, 605, 205, 600)
A <- N * S
