#simulation of a prospective trial with a rare event
library(stats)
library(dplyr)

#Hat tip to 
#https://github.com/aalthous/RCT-Simulation-v1/blob/main/RCT_Binary_Outcome.R 


# Trial Design Parameters
nPatients <- 5000 # here is where you specify the number of patients you want included in each trial
prospective.nc <- 0.0006401762 # here is where you specify the prospective nc event rate. Non carriers
prospective.c <- 0 #As they won't get carbamazepine
allele.freq <- 0.9606
histinc <- 0.0023 #This is the historical incidence of the outcome, from Chen et al 2011 (same historical incidence)

# Simulation Parameters
nSims <- 10000 # here is where you specify the number of trials that you want to simulate
#10000 is actual value

set.seed(123) # set for reproducibility
results=rep(NA,nSims*5)
results=matrix(results,nSims,5)
for(i in 1:nSims){

  treatment=rep(1, nPatients) # this creates a vector of "treatment allocations", which is actually just 1
  carriers = rbinom(nPatients, 1, allele.freq) # vector for carrier probability
  
  adr=rep(NA, nPatients)
  l = length(carriers)
    for(j in 1:l){
      if (
        carriers[j]==0){adr[j]=0}
      else 
        {
          adr[j] = rbinom(1, 1, prospective.nc)}
    }
  trialdata=data.frame(cbind(treatment, adr)) # this creates a data frame with pid, treatment allocation, and adr outcome
  a <- nPatients*mean(trialdata$adr)
  b <- nPatients-(nPatients*mean(trialdata$adr))
  c <- nPatients*histinc
  d <- nPatients-(nPatients*histinc)

  or <- 
    (a*d)/
    (b*c)
  
  results_temp <- c(a,b,c,d)
  
  results[i,]=c(results_temp, or)
}

total_cases <- mean(results[,1])
total_controls <- mean(results[,2])
mean_or <- mean(results[,5])
quantiles <- quantile(results[,5], prob=c(0.025, 0.975))
print(data.frame(nPatients, total_cases, total_controls, mean_or, quantiles, 1/mean_or, 1/quantiles))

# This section calculates power using Fisher's exact:
#results_noor <- results[,1:4]
#pvalues <- apply(results_noor,1, function(x) fisher.test(matrix(x,nr=2))$p.value) 
#warnings show that rounding is done automatically
#confints <- apply(results_noor,1,function(x) fisher.test(matrix(x,nr=2))$conf.int)
#confints <- as.data.frame(t(confints))
#colnames(confints) <- c("ll", "ul")
#estimate <- apply(results_noor,1,function(x) fisher.test(matrix(x,nr=2))$estimate)
#results_fisher <- data.frame(confints, estimate)
#power <- (nrow(results_fisher[results_fisher$pvalues < 0.05,]))/nrow(results_fisher)


