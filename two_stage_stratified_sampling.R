set.seed(42069)

# target function:

inverse_num <- function(x){
  return(1/x)
}

target <- function(V_sample, M_PSU, strat_var, df_picked_students_number, df, school_var, df_picked_schools){
    temp1 <- sum((inverse_num(V_sample) - inverse_num(M_PSU)) * M_PSU^2 * strat_var)
    temp2 <- 0
    for(i in 1:length(M_PSU)){
      temp3 <- 0
      for(j in 1:max(M_PSU)){
        if(df_picked_schools[i, j] != 0){
        temp3 <- temp3 + (inverse_num(df_picked_students_number[i, j]) - 
                            inverse_num(df[i, j]) * df[i, j]^2 * school_var[i, j])
        }else{temp3 <- temp3}
      }
      temp2 <- temp2 + temp3 * M_PSU[i]/V_sample[i]
    }
    return(temp1 + temp2)
}

# generate propositions:

gen_prop_box_stratified_PSU <- function(V_sample, m_PSU, M_PSU){
  index_move_PSU_1 <- ceiling(runif(1, min=0, max=length(V_sample)))
  index_move_PSU_2 <- ceiling(runif(1, min=0, max=length(V_sample)))
  if(index_move_PSU_1 == index_move_PSU_2){
    return(V_sample)
  }
  if(V_sample[index_move_PSU_1] >= m_PSU[index_move_PSU_1] &
     V_sample[index_move_PSU_1] < M_PSU[index_move_PSU_1] &
     V_sample[index_move_PSU_2] > m_PSU[index_move_PSU_2] &
     V_sample[index_move_PSU_2] <= M_PSU[index_move_PSU_2]){
    
    V_prop <- V_sample
    V_prop[index_move_PSU_1] <- V_prop[index_move_PSU_1] + 1
    V_prop[index_move_PSU_2] <- V_prop[index_move_PSU_2] - 1
    
  }else{
    return(V_sample)
  }
  return(V_prop)
}

gen_prop_box_stratified_SSU <- function(m_school_bound, df, df_picked_students_number, df_picked_schools){
  index_move_SSU_1 <- ceiling(runif(1, min=0, max=m_school_bound))
  index_move_SSU_2 <- ceiling(runif(1, min=0, max=m_school_bound))
  if(index_move_SSU_1 == index_move_SSU_2){
    return(df_picked_students_number)
  }
  indexes_of_picked_PSU <- unname(which(df_picked_schools == 1, arr.ind = T))
  index_move_SSU_1 <- indexes_of_picked_PSU[index_move_SSU_1,]
  index_move_SSU_2 <- indexes_of_picked_PSU[index_move_SSU_2,]
  if(df_picked_students_number[index_move_SSU_1[1], index_move_SSU_1[2]] < 
     df[index_move_SSU_1[1], index_move_SSU_1[2]] &
     df_picked_students_number[index_move_SSU_1[1], index_move_SSU_1[2]] >= 0 & 
     df_picked_students_number[index_move_SSU_2[1], index_move_SSU_2[2]] <=
     df[index_move_SSU_2[1], index_move_SSU_2[2]] & 
     df_picked_students_number[index_move_SSU_2[1], index_move_SSU_2[2]] > 1){
    df_picked_students_number_prop <- df_picked_students_number
    df_picked_students_number_prop[index_move_SSU_1[1], index_move_SSU_1[2]] <-
      df_picked_students_number_prop[index_move_SSU_1[1], index_move_SSU_1[2]] + 1
    df_picked_students_number_prop[index_move_SSU_2[1], index_move_SSU_2[2]] <- 
      df_picked_students_number_prop[index_move_SSU_2[1], index_move_SSU_2[2]] - 1
  }else{
    return(df_picked_students_number)
  }
  return(df_picked_students_number_prop)
}

################################################################################

# Box constrains:
M_PSU <- c(35, 40, 50, 70, 150, 50, 30, 55, 200, 450, 50, 70, 130, 150, 200) # (V) PSU upper bounds
m_PSU <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) # PSU lower bounds
M_SSU <- 250 # Max number of students in a school
stud_param <- 250 # Max value of a student parameter (e. g. height, used to estimate pop. variance)
m_school_bound <- 500 # num of schools

# Generate school matrix: 

df <- matrix(0, nrow = length(M_PSU), ncol = max(M_PSU))
for(i in 1:length(M_PSU)){
  v <- sample.int(M_SSU, M_PSU[i], replace = T) + 1 # (!!!) schools need to have at least 2 students 
  df[i, 1:M_PSU[i]] <- v
} # df[i, j] - number of 'students' in j'SSU from i'PSU

# Generate students parameters: (for the sake of D^2 & S^2) 

df_t <- 0
n <- length(M_PSU) * max(M_PSU)
df_t <- as.list(numeric(n))
dim(df_t) <- c(length(M_PSU), max(M_PSU))

for(i in 1:length(M_PSU)){
  for(j in 1:max(M_PSU)){
    l <- df[i, j]
    if(l == 0){
      df_t[[i, j]] <- 0
    }else{
      u <- sample.int(stud_param, l, replace = T)
      df_t[[i, j]] <- u
    }
  }
} # df_t[[i, j]][u] - value of the interest parameter form the u-th student form j-th school from i-th discrict

# school totals: t_(h, j):
df_t_totals <- matrix(0, nrow = length(M_PSU), ncol = max(M_PSU))

for(i in 1:length(M_PSU)){
  for(j in 1:max(M_PSU)){
    df_t_totals[i, j] <- sum(df_t[[i, j]])
  }
}

# stratum average: t_dash:
strat_avg <- numeric(length(M_PSU))
for(i in 1:length(M_PSU)){
  strat_avg[i] <- sum(df_t_totals[i, ])/M_PSU[i]
}

# school average:
df_t_school_avg <- matrix(0, nrow = length(M_PSU), ncol = max(M_PSU))
for(i in 1:length(M_PSU)){
  for(j in 1:max(M_PSU)){
    if(length(df_t[[i, j]]) != 0){
      df_t_school_avg[i, j] <- df_t_totals[i, j] / length(df_t[[i, j]])
    }else{df_t_school_avg[i, j] <- 0}
  }
}

# stratum variance (D^2)_h:
strat_var <- numeric(length(M_PSU))
for(i in 1:length(M_PSU)){
  strat_var[i] <- sum((df_t_totals[i,]-strat_avg[i])^2)/(M_PSU[i]-1)
}

# school variance (S^2)_(h, j):
school_var <- matrix(0, nrow = length(M_PSU), ncol = max(M_PSU))
for(i in 1:length(M_PSU)){
  for(j in 1:max(M_PSU)){
    if(length(df_t[[i, j]]) != 1){
      school_var[i, j] <- sum((df_t[[i, j]]-df_t_school_avg[i, j])^2)/(length(df_t[[i, j]])-1)
    }else{school_var[i, j] <- 0}
  }
}

################################################################################
set.seed(42069)
# generate school sample with size bounds:
V_sample <- numeric(length(M_PSU))
while(sum(V_sample) != m_school_bound){
  V_sample <- numeric(length(M_PSU))
  for(i in 1:length(M_PSU)){
    V_sample[i] <- sample.int(M_PSU[i], 1)
  }
  if(any(V_sample < m_PSU)){
    V_sample <- numeric(length(M_PSU))
  }
}
# Test the sample:
sum(V_sample) == m_school_bound
all(V_sample >= m_PSU)
all(V_sample <= M_PSU)

# pick speficic schools (indexes from df)
df_picked_schools <- matrix(0, nrow = length(M_PSU), ncol = max(M_PSU))
for(i in 1:length(M_PSU)){
  s <- sample.int(M_PSU[i], V_sample[i])
  for(j in 1:length(s)){
    df_picked_schools[i, s[j]] <- 1
  }
}

# generate student sample size (from df with respect to V_sample))
stud_sample_size <- 0
for(i in 1:length(M_PSU)){
  temp <- V_sample[i]/M_PSU[i]
  stud_sample_size <- stud_sample_size + temp * sum(df[i, ] * df_picked_schools[i, ]) 
}

# all students:
stud_all <- sum(df * df_picked_schools)

# generate student sample 
df_picked_students_number <- df_picked_schools
for(i in 1:length(M_PSU)){
  for(j in 1:max(M_PSU)){
    df_picked_students_number[i, j] <- df_picked_students_number[i, j] * 
      df[i, j] * stud_sample_size / stud_all
  }
}
# generate student matrix (small error):
for(i in 1:length(M_PSU)){
  for(j in 1:max(M_PSU)){
    if(df_picked_students_number[i, j] - floor(df_picked_students_number[i, j]) > 0.5){
      df_picked_students_number[i, j] <- ceiling(df_picked_students_number[i, j])
    }else{df_picked_students_number[i, j] <- floor(df_picked_students_number[i, j])}
  }
}

################################################################################

# annealing schedule:

eps <- 1e-10 #for the sake of numeric stability (avoiding zero-division)
T_schedule <- 1:20*1000000000 # Annealing schedule (temperature) #100000
N_schedule <- numeric(length(T_schedule)) + 1500
N_schedule[1] <- N_schedule[1] + 5000
N_schedule[2] <- N_schedule[2] + 1000

# for plotting and assessment:

temp_target_PSU_change <- list()
temp_target_SSU_change <- list()
temp_step_PSU <- list()
temp_step_SSU <- list()

trajectories <- list()

# THE annealing

for(i in 1:length(N_schedule)){
  T_temp <- T_schedule[i]
  for(j in 1:N_schedule[i]){
    trajectories <- append(trajectories, V_sample)
    V_prop <- gen_prop_box_stratified_PSU(V_sample, m_PSU, M_PSU)
    df_picked_schools_V_prop <- matrix(0, nrow = length(M_PSU), ncol = max(M_PSU))
    for(i in 1:length(M_PSU)){
      s <- sample.int(M_PSU[i], V_prop[i])
      for(j in 1:length(s)){
        df_picked_schools_V_prop[i, s[j]] <- 1
      }
    }
    stud_sample_size_V_prop <- 0
    for(i in 1:length(M_PSU)){
      temp <- V_prop[i]/M_PSU[i]
      stud_sample_size_V_prop <- stud_sample_size_V_prop + temp * sum(df[i, ] * df_picked_schools_V_prop[i, ]) 
    }
    stud_all_V_prop <- sum(df * df_picked_schools_V_prop)
    df_picked_students_number_V_prop <- df_picked_schools_V_prop
    for(i in 1:length(M_PSU)){
      for(j in 1:max(M_PSU)){
        df_picked_students_number_V_prop[i, j] <- df_picked_students_number_V_prop[i, j] * 
          df[i, j] * stud_sample_size_V_prop / stud_all_V_prop
      }
    }
    for(i in 1:length(M_PSU)){
      for(j in 1:max(M_PSU)){
        if(df_picked_students_number_V_prop[i, j] - floor(df_picked_students_number_V_prop[i, j]) > 0.5){
          df_picked_students_number_V_prop[i, j] <- ceiling(df_picked_students_number_V_prop[i, j])
        }else{df_picked_students_number_V_prop[i, j] <- floor(df_picked_students_number_V_prop[i, j])}
      }
    }
    df_picked_students_number_prop <- gen_prop_box_stratified_SSU(m_school_bound, df, 
                                                                  df_picked_students_number, 
                                                                  df_picked_schools)
    transit_value_PSU <- exp((target(V_sample, M_PSU, strat_var, 
                                     df_picked_students_number, 
                                     df, school_var, df_picked_schools)-
                                target(V_prop, M_PSU, strat_var, 
                                       df_picked_students_number_V_prop, 
                                       df, school_var, df_picked_schools_V_prop))/T_temp)
    transit_value_SSU <- exp((target(V_sample, M_PSU, strat_var, 
                                     df_picked_students_number, 
                                     df, school_var, df_picked_schools)-
                                target(V_sample, M_PSU, strat_var, 
                                       df_picked_students_number_prop, 
                                       df, school_var, df_picked_schools))/T_temp)
    temp_target_PSU_change <- append(temp_target_PSU_change, target(V_sample, M_PSU, strat_var,
                                                                    df_picked_students_number,
                                                                    df, school_var, df_picked_schools))
    temp_target_SSU_change <- append(temp_target_SSU_change, target(V_sample, M_PSU, strat_var,
                                                                    df_picked_students_number,
                                                                    df, school_var, df_picked_schools))
    prob_move_PSU <- min(1, transit_value_PSU)
    prob_move_SSU <- min(1, transit_value_SSU)
    mass <- target(V_prop, M_PSU, strat_var, 
                   df_picked_students_number_V_prop, 
                   df, school_var, df_picked_schools_V_prop) + 
      target(V_sample, M_PSU, strat_var, 
             df_picked_students_number_prop, 
             df, school_var, df_picked_schools)
    print(V_sample)
    if(prob_move_PSU == 1 & prob_move_SSU == 1){
      coin_flip <- runif(1)
      if(target(V_prop, M_PSU, strat_var, 
                df_picked_students_number_V_prop, 
                df, school_var, df_picked_schools_V_prop)/mass > coin_flip){
        V_sample <- V_prop
        temp_step_PSU <- append(temp_step_PSU, length(which(unlist(temp_step_PSU) > 0)) + 1)
        temp_step_SSU <- append(temp_step_SSU, 0)
      }else{
        df_picked_students_number <- df_picked_students_number_prop
        temp_step_SSU <- append(temp_step_SSU, length(which(unlist(temp_step_SSU) > 0)) + 1)
        temp_step_PSU <- append(temp_step_PSU, 0)
      }
    }
    if(prob_move_PSU == 1 & prob_move_SSU < 1){
      V_sample <- V_prop
      temp_step_PSU <- append(temp_step_PSU, length(which(unlist(temp_step_PSU) > 0)) + 1)
      temp_step_SSU <- append(temp_step_SSU, 0)
    }
    if(prob_move_PSU < 1 & prob_move_SSU == 1){
      df_picked_students_number <- df_picked_students_number_prop
      temp_step_SSU <- append(temp_step_SSU, length(which(unlist(temp_step_SSU) > 0)) + 1)
      temp_step_PSU <- append(temp_step_PSU, 0)
    }
    if(prob_move_PSU < 1 & prob_move_SSU < 1){
      if(prob_move_PSU < eps & prob_move_SSU < eps){
        V_sample <- V_sample
        df_picked_students_number <- df_picked_students_number
        temp_step_SSU <- append(temp_step_SSU, -50)
        temp_step_PSU <- append(temp_step_PSU, -50)
      }
      if(prob_move_PSU > prob_move_SSU){
        coin_flip <- runif(1)
        if(coin_flip < prob_move_PSU){
          V_sample <- V_prop
          temp_step_PSU <- append(temp_step_PSU, length(which(unlist(temp_step_PSU) > 0)) + 1)
          temp_step_SSU <- append(temp_step_SSU, 0)
        }else{
          V_sample <- V_sample
          df_picked_students_number <- df_picked_students_number
          temp_step_SSU <- append(temp_step_SSU, -50)
          temp_step_PSU <- append(temp_step_PSU, -50)
        }
      }
      if(prob_move_SSU > prob_move_PSU){
        coin_flip <- runif(1)
        if(coin_flip < prob_move_SSU){
          df_picked_students_number <- df_picked_students_number_prop
          temp_step_SSU <- append(temp_step_SSU, length(which(unlist(temp_step_SSU) > 0)) + 1)
          temp_step_PSU <- append(temp_step_PSU, 0)
        }else{
          V_sample <- V_sample
          df_picked_students_number <- df_picked_students_number
          temp_step_SSU <- append(temp_step_SSU, -50)
          temp_step_PSU <- append(temp_step_PSU, -50)
        }
      }
      if(prob_move_PSU == prob_move_SSU){
        V_sample <- V_sample
        df_picked_students_number <- df_picked_students_number
        temp_step_SSU <- append(temp_step_SSU, -50)
        temp_step_PSU <- append(temp_step_PSU, -50)
      }
    }
  }
}

################################################################################

plot(unlist(temp_target_PSU_change),type = "l", pch = 23, col = "black", lty = 1, 
     xlab = "Liczba wygenerowanych kroków", 
     ylab = "Wartość funkcji celu")
#plot(unlist(temp_target_SSU_change), type = "l")
plot(unlist(temp_step_PSU) + unlist(temp_step_SSU),
     xlab = "Liczba wygenerowanych kroków",
     ylab = "")
V_sample

plot(unlist(temp_step_SSU))
final_var <- unlist(temp_target_PSU_change)[length(unlist(temp_target_PSU_change))]

# Calculate optimal (unbounded) results: 

# gamma: 

gamma <- M_PSU * strat_var
for(i in 1:length(gamma)){
  gamma[i] <- gamma[i] - sum(df[i,]*school_var[i,])
}
all(gamma > 0)

# optimal alloc:

opt_PSU <- m_school_bound * sqrt(M_PSU * gamma)/sum(sqrt(M_PSU * gamma))

opt_SSU <- matrix(0, nrow = length(M_PSU), ncol = max(M_PSU))
temp <- 0
for(i in 1:length(M_PSU)){
  temp <- temp + sum(df[i, ] * sqrt(school_var[i, ]))
}

for(i in 1:length(M_PSU)){
  for(j in 1:max(M_PSU)){
    opt_SSU[i, j] <- stud_sample_size * (M_PSU[i]/opt_PSU[i]) * df[i, j] * sqrt(school_var[i, j])/temp 
  }
}

# optimal variance:
opt_var <- sum(sqrt(M_PSU * gamma))^2/m_school_bound + (temp)^2/stud_sample_size - sum(M_PSU * strat_var)

round(abs(1-(final_var/opt_var)), 5)*100
opt_var
final_var

df_picked_schools
length(which(df_picked_students_number == round(opt_SSU)))
length(M_PSU) * max(M_PSU)
length(which(abs(df_picked_students_number - opt_SSU) < 2))
temp <- abs(df_picked_students_number - opt_SSU)
trac_matrix <- matrix(0, nrow = sum(N_schedule), ncol = length(M_PSU))

for(i in 1:(length(trajectories)/15)){s
  
  trac_matrix[i,] <- unlist(trajectories[(1 + 15*(i-1)):(15 + 15*(i-1))])
}

matplot(1:(length(trajectories)/15), trac_matrix[1:(length(trajectories)/15),], type = "l", 
        xlab = "Iteracja",
        ylab = "Rozmiar próbki")

