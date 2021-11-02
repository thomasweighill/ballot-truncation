#PL Truncation

#only use this if you are running on a cluster
library(DirichletReg, lib='~/local/R_libs/')
library(STV, lib='~/local/R_libs/')
library(matrixStats, lib='~/local/R_libs/')
args = commandArgs(trailingOnly=TRUE)
both <- as.logical(args[1])
alpha <- as.numeric(args[2])
minority <- as.numeric(args[3])
#-------------------------------
# 
# library(DirichletReg)
# library(STV)
# library(matrixStats)

set.seed(2021)
num_simulations <- 100

LuceWithS <- function(candidate_names, support_vector, k){
  ballot <- sample(candidate_names, size = k, prob=support_vector/sum(support_vector), replace = FALSE)
  while (ballot[1] == 'S'){
    ballot <- sample(candidate_names, size = k, prob=support_vector/sum(support_vector), replace = FALSE)
  }
  ballot
}

print(paste('Running with minority = ', minority, ' and alpha = ', alpha))
#initialize
group_1_candidate_names <- c("A", "B", "C", "D", "E", "F", "S")
group_2_candidate_names <- c("G", "H", "I", "J", "K", "L", "S")
real_candidate_names <- c(group_1_candidate_names[-7], group_2_candidate_names[-7])
representation_results <- data.frame(matrix(nrow = 0, ncol = 4))
names(representation_results) <- c('S', 'group1wins', 'average_ballot_length_all', 'average_ballot_length_group1')

#run simulations
for (run in seq(num_simulations)) {
  print(paste('Running election number', run))
  #create support vectors
  vector_1 <- sort(rdirichlet(1, rep(alpha, 6)), decreasing = TRUE)
  vector_2 <- sort(rdirichlet(1, rep(alpha, 6)), decreasing = TRUE)
  
  #vary truncation levels
  for (S in seq(-0.2,-2.0,by=-0.2)){
    
    #generate group 1 ballots
    k = 7
    n = as.integer(minority*1000)
    new_vector_1 <- c(vector_1, 10^(S))
    ballots_1 <- replicate(n, LuceWithS(group_1_candidate_names, new_vector_1/sum(new_vector_1), k))
    ballots_1 <- t(matrix(ballots_1, nrow=k, byrow = F))
    ranksframe_1 = data.frame(matrix(nrow= n,ncol= 12))
    names(ranksframe_1) <- real_candidate_names
    for (b in c(1:nrow(ballots_1))){
      for (rank in which(!is.na(ballots_1[b,]))){
        if (ballots_1[b,rank] == "S") {
          break
        }
        ranksframe_1[b,ballots_1[b,rank]] = rank
      }
    }
    
    #generate group 2 ballots
    n = as.integer((1-minority)*1000)
    if (both){
      new_vector_2 <- c(vector_2, 10^(S))
      k = 7
      ballots_2 <- replicate(n, LuceWithS(group_2_candidate_names, new_vector_2/sum(new_vector_2), k))
    }
    else {
      new_vector_2 <- vector_2
      k = 6
      ballots_2 <- replicate(n, LuceWithS(group_2_candidate_names[-7], new_vector_2/sum(new_vector_2), k))
    }
    ballots_2 <- t(matrix(ballots_2, nrow = k, byrow = F))
    ranksframe_2 = data.frame(matrix(nrow= n,ncol= 12))
    names(ranksframe_2) <- real_candidate_names
    for (b in c(1:nrow(ballots_2))){
      for (rank in which(!is.na(ballots_2[b,]))){
        if (ballots_2[b,rank] == "S") {
          break
        }
        ranksframe_2[b,ballots_2[b,rank]] = rank
      }
    }
    
    #combine and compute winners
    full_ranksframe <- rbind(ranksframe_1, ranksframe_2)
    Bal <- cleanBallots(full_ranksframe, cand.names = real_candidate_names)
    out <-STV::stv(Bal, seats = 6, surplusMethod = "Cambridge")
    
    #store outcome
    is_elected <- as.integer(real_candidate_names%in%(out$elected))
    group1_is_elected <- sum(as.integer(group_1_candidate_names%in%(out$elected)))
    average_ballot_length_all <- mean(rowMaxs(as.matrix(full_ranksframe), na.rm=T))
    average_ballot_length_group1 <- mean(rowMaxs(as.matrix(ranksframe_1), na.rm=T))
    representation_results[nrow(representation_results)+1,] <- c(S, group1_is_elected, average_ballot_length_all, average_ballot_length_group1)
  }
}

representation_results$S <- round(representation_results$S, digits = 2)
representation_results$S <- as.factor(x=representation_results$S)
if (both){
  write.csv(x = representation_results, file = paste( "PL_both_a",alpha,"_",minority, "_minority_scatterplot.csv"))
} else{
  write.csv(x = representation_results, file = paste( "PL_g1_a",alpha,"_",minority, "_minority_scatterplot.csv"))
}