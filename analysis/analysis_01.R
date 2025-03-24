rm(list=ls())
hablar::set_wd_to_script_path()
library(tidyverse)

d <- read_csv("../task/analysis/analysis_all_data.csv")
str(d)

# # remove this until you manage to solve issue from matlab:
# d <- d %>%
#   filter(str_detect(participant_id, "^[A-Za-z]"))

# general prep
d <- d %>%
  mutate(PID = str_sub(edf_file, 1, 4)) %>%
  rename(session_id = participant_id) %>%
  mutate(win=ifelse(win==-1, NA, win),
         P=ifelse(win==-1, NA, P)) %>%
  filter(!is.na(win))

# function that put data of a single participant into stan format
#d_i <- d[d$PID==unique(d$PID)[1],]
#str(d_i)

# int<lower=1> J;                    // n block
# int<lower=1> N[J];                 // n trials x block
# int<lower=1> maxN;                 // max number of trials
# int<lower=0,upper=1> C[J,maxN];    // choice (of option with high reward prob)
# int<lower=0,upper=1> R[J,maxN,];   // feedback (positive or negative)

prep_stan_data <- function(pid, data=d){
  
  d_i <- d[d$PID==pid,]
  d_i$block_id <- paste(d_i$session_id,d_i$block_n,sep="_")
  
  J <- length(unique(d_i$block_id))
  N <- tapply(d_i$block_id,d_i$block_id, length)
  
  if(any(N)<10){
    xbid <- names(N)[which(N<10)] 
    d_i <- d_i %>%
      filter(!is.element(block_id, xbid))
    J <- length(unique(d_i$block_id))
    N <- tapply(d_i$block_id,d_i$block_id, length)
  }
  
  maxN <- max(N)
  
  d_i$C <- ifelse(d_i$P>0.5,1,0)
  C <- array(dim=c(J,maxN))
  R <- array(dim=c(J,maxN))
  for(j in 1:J){
    c_j <- unique(d_i$block_id)[j]
    C[j,1:N[j]] <- d_i$C[d_i$block_id==c_j]
    R[j,1:N[j]] <- d_i$win[d_i$block_id==c_j]
  }
  
  # check NA
  C[is.na(C)] <- -99
  R[is.na(R)] <- -99
  
  d_stan <- list(J=J, N=N, maxN=maxN, C=C, R=R)

  str(d_stan)
  return(d_stan)
}


d_stan <- prep_stan_data(unique(d$PID)[1], data=d)

library(rstan)

options(mc.cores = parallel::detectCores()) # indicate stan to use multiple cores if available
# run sampling: 4 chains in parallel on separate cores

test_fit <- stan(file = "Qlearning_M1.stan", data = d_stan, iter = 2000, chains = 4)


