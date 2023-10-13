#Fixed Truncation

#only use this if you are running on a cluster
library(DirichletReg, lib='~/local/R_libs/')
library(STV, lib='~/local/R_libs/')
library(matrixStats, lib='~/local/R_libs/')
library(Rsolnp, lib='~/local/R_libs/')
library(DiscreteWeibull, lib='~/local/R_libs/')
args = commandArgs(trailingOnly=TRUE)
both <- as.logical(args[1])
alpha <- as.numeric(args[2])
minority <- as.numeric(args[3])
#-------------------------------
# 
# library(DirichletReg)
# library(STV)
# library(matrixStats)
# library(DiscreteWeibull)

set.seed(2021)
num_simulations <- 100
scales <- 1-10^(-seq(1,10)/5)
global_shape <- 1.4

print(paste('Running with minority = ', minority, ' and alpha = ', alpha))
#initialize
group_1_candidate_names <- c("A", "B", "C", "D", "E", "F")
group_2_candidate_names <- c("G", "H", "I", "J", "K", "L")
real_candidate_names <- c(group_1_candidate_names, group_2_candidate_names)
representation_results <- data.frame(matrix(nrow = 0, ncol = 4))
names(representation_results) <- c('scale', 'group1wins', 'average_ballot_length_all', 'average_ballot_length_group1')

#run simulations
for (run in seq(num_simulations)) {
  print(paste('Running election number', run))
  #create support vectors
  vector_1 <- sort(rdirichlet(1, rep(alpha, 6)), decreasing = TRUE)
  vector_2 <- sort(rdirichlet(1, rep(alpha, 6)), decreasing = TRUE)
  
  #vary truncation levels
  for (scale in scales){
    
    #generate group 1 ballots
    k = length(group_1_candidate_names)
    n = as.integer(minority*1000)
    ballots_1 <- replicate(n, sample(group_1_candidate_names, size = k, prob=vector_1/sum(vector_1), replace = FALSE))
    ballots_1 <- t(matrix(ballots_1, nrow=k, byrow = F))
    ranksframe_1 = data.frame(matrix(nrow= n,ncol= 12))
    names(ranksframe_1) <- real_candidate_names
    for (b in c(1:nrow(ballots_1))){
      length <- rdweibull(1, scale, global_shape, zero=F)
      for (rank in which(!is.na(ballots_1[b,]))){
        if (rank > length) {break}
        ranksframe_1[b,ballots_1[b,rank]] = rank
      }
    }
    
    #generate group 2 ballots
    n = as.integer((1-minority)*1000)
    ballots_2 <-replicate(n, sample(group_2_candidate_names, size = k, prob=vector_2/sum(vector_2),replace = FALSE))
    ballots_2 <- t(matrix(ballots_2, nrow = k, byrow = F))
    ranksframe_2 = data.frame(matrix(nrow= n,ncol= 12))
    names(ranksframe_2) <- real_candidate_names
    for (b in c(1:nrow(ballots_2))){
      if (both){
        length <- rdweibull(1, scale, global_shape, zero=F)
      }
      else {
        length <- length(group_2_candidate_names)
      }
      for (rank in which(!is.na(ballots_2[b,]))){
        if (rank > length) {break}
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
    representation_results[nrow(representation_results)+1,] <- c(scale, group1_is_elected, average_ballot_length_all, average_ballot_length_group1)
  }
}

representation_results$scale <- round(representation_results$scale, digits = 2)
representation_results$scale <- as.factor(x=representation_results$scale)
if (both){
  write.csv(x = representation_results, file = paste( "IT_both_a",alpha,"_",minority, "_minority_scatterplot.csv"))
} else{
  write.csv(x = representation_results, file = paste( "IT_g1_a",alpha,"_",minority, "_minority_scatterplot.csv"))
}