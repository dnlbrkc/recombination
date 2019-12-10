packages <- c('tictoc',
              'profvis',
              'compositions',
              'here',
              'gtools',
              'fGarch',
              'tidyverse'
)

lapply(packages, require, character.only = TRUE)

source(here("setup/NK_function.R"))

na <- 10          # number of agents
# N <- 16            # string length
# K <- c(0, 2, 6, 8) # task complexities
# 1 <- 1         # number of replications
nK <- length(K)
local=1

N <- 10           # string length
K <- c(0, 4, 8) # task complexities
repmax <- 1         # number of replications
nK <- length(K)
local=1
reps=30
# load landscapes
#name <- paste0(here("landscapes/"),"LS_",N,'.Rds')
landscape <- permutations(2,N,v=c(0,1),repeats.allowed=TRUE)
nsol <- nrow(landscape)

# all combinations of groups (col1=prop mutators, col2=prop recombinators)
groups <- seq(0,na,na/10)
groups <- matrix(groups, nrow=length(groups), ncol=2)
groups[,2] <- na - groups[,1]
ng <- nrow(groups)

ch_extents <- c(0.05, 0.25, 0.5, 0.75, 1)  #delete 0.95 in case of n=10
#simchoices <- c(0,1, 11, 21, 51, 81, 91)
simchoices <- c(0,2,3,4,5,6,7,8,9)
# tmax may depend on K 
timesteps <- rep(30,nK)

# initialize payoffs
payoffs <- array(0,dim=c(2^N,nK,reps))
# load payoffs

c=0
for( ik in K){
  c=c+1
  for(frep in 1:reps){
    LS<-permutations(2,N,v=c(0,1),repeats.allowed=TRUE)
    LS<-as.data.frame(LS)
    if (ik==0){
      depends <- as.vector(1:N)
      values <- replicate(N,round(runif(2,0,1),1))
      fitness <- values
    } else {
      # INSTEAD OF depends <- rbind(1:N,replicate(N,sample(c(1:N),K,replace=F)))
      depends <- matrix(nrow = ik+1, ncol = N)
      for(i in 1:N) {
        string <- c(1:N)
        string <- string[-i]
        depends[,i] <- c(i,sample(string,ik,replace=F))
      }
      #
      combinations <- permutations(2,ik+1,v=c(0,1),repeats.allowed=TRUE)
      values <- replicate(N,round(runif(nrow(combinations),0,1),1))
      fitness <- cbind(combinations,values)
    }
    
    land <- generate_landscape(N,ik,LS,fitness,depends)
    payoffs[,c,frep] <- (land[,N+1]/max(land[,N+1]))^8
  }
}



for(f in 1:reps){ #repetitions  
  # loop over different ways neighbors selected (0=rand, 1=10% most, 11=20% most etc )
  # REVIEW 2----------------------------------------------------------------
  
  
  ro=0
  for(simchoice in simchoices) {  #with 10: 0,1,3,4,6,8,9 ? 
    
    ro=ro+1
    # loop over different Ks
    for( ik in 1:nK ) {
      
      # tmax depends on K
      tmax <- timesteps[ik]
      
      for( cex in 1:length(ch_extents)){
        tic()
        
        ch_extent <- ceiling(ch_extents[cex]*N) 
        
        # initialize storage variables
        pay_g <- matrix(NA, nrow=ng, ncol=tmax+1)
        truedif_g <- matrix(NA, nrow=ng,ncol=tmax+1)
        truecex_g <- matrix(NA, nrow=ng,ncol=tmax+1)
        unique_g <- matrix(NA, nrow=ng,ncol=tmax+1)
        # loop over all possible proportions of groups composed of individuals using different strategies
        for( g in 1:ng ){
          a_groups <- c( unlist(replicate(groups[g,1], 1)),
                         unlist(replicate(groups[g,2], 2)) )
          a_groups <- sample(a_groups) # order of groups randomized
          
          pay_t <- matrix(NA, 1, tmax)
          dif_true_avg <- matrix(NA, 1, tmax)
          cex_true_avg <- matrix(NA, 1, tmax)
          uniqueSols_avg <- matrix(NA,1,tmax)
          # starting solutions of agents
          sol_no <- ceiling(runif(na)*nsol)
          a_sols <- landscape[sol_no, ]
          
          # storage for this replication
          dif_true <- matrix(NA, nrow=na, ncol=tmax)
          cex_true <- matrix(NA, nrow=na, ncol=tmax)
          uniqueSols <- vector()
          # loop over tmax time steps
          for( t in 1:tmax ){
            # temporary matrices to store new solutions - needed for simultaneous updating
            sol_no_new <- replicate(length(sol_no), NA)
            a_sols_new <- matrix(NA, nrow=nrow(a_sols), ncol=ncol(a_sols))
            
            
            if(simchoice >= 1){
              distances <- matrix(NA, nrow=na, ncol=na)
              for( iia in 1:na ){
                for( ja in 1:na ){
                  distances[iia, ja] <- sum(a_sols[iia, ] != a_sols[ja, ])
                }
              }
            }
            
            # order of agents random
            ord <- sample(1:na)
            # loop over agents - learning depends on the group - some update some not
            for( ia in ord ){
              if( simchoice == 0 ) { # random choice
                neighbor <- sample(1:na,1);
              } else { # most similar neighbor
                x <- sort.int(distances[ia,], index.return=T)
                range <- sample(simchoices[-1][ro-1]:simchoice,1)
                neighbor <- x$ix[range]
              }
              
              test <- a_sols[ia,] # my current solution
              neighbor_sol <- a_sols[neighbor, ] # find neighbor's solution
              dif_true[ia, t] <- sum( test != neighbor_sol)/N
              
              disbits <- which(neighbor_sol != test)
              cex_true[ia, t] <- length(disbits)/N 
              
              if( length(disbits) >= ch_extent ){ 
                disbits_take <- disbits[1:ch_extent]
                cex_true[ia, t] <- ch_extent/N 
                
              }
              if( length(disbits) < ch_extent ){
                disbits_take <- disbits
                cex_true[ia, t] <- length(disbits_take)/N  #record
              } # end of finding dissimilar bits
              
              if( a_groups[ia] == 1){ # implement mutation
                #change ch_extent number of bits
                test <- unlist(a_sols[ia,]) # my current solution
                position <- sample(1:N,ch_extent)#
                test[position] <- abs(test[position]-1) # flip the bits
                
                indx <- 1 + unbinary(paste(test, collapse=''))
                if( payoffs[indx, ik,f] > payoffs[sol_no[ia], ik,f] ) { # find payoff of the new string
                  a_sols_new[ia, ] <- as.numeric(test) # take it if better
                  sol_no_new[ia] <- indx
                }
              }
              
              else if ( a_groups[ia] == 2 ){ # can it be improved by recombination
                test[disbits_take] <- neighbor_sol[disbits_take]
                
                #indx <- which(apply(landscape, 1, function(row) all( row == test)))
                indx <- 1 + unbinary(paste(test, collapse=''))
                if( payoffs[indx, ik,f] > payoffs[sol_no[ia], ik,f] ) { # find payoff of the new string
                  a_sols_new[ia, ] <- as.numeric(test) # take it if better
                  sol_no_new[ia] <- indx
                }
                
              } # end of social learning in groups
            } # end of looping over agents
            
            # update a_sols after all agents made decisions
            a_sols[complete.cases(a_sols_new),] <- a_sols_new[complete.cases(a_sols_new),]
            sol_no[complete.cases(sol_no_new)] <- sol_no_new[complete.cases(sol_no_new)]
            
            # save current payoff - sum of payoffs of all agents
            pay_t[1, t] <- sum( payoffs[sol_no, ik,f])
            
            #track unique solutions over time
            uniqueSols[t] <-length(unique(sol_no))
          } # end looping through time
          
          # save average difference between models
          dif_true_avg[1, ] <- colMeans(dif_true) #: check mean
          cex_true_avg[1, ] <- colMeans(cex_true)
          uniqueSols_avg[1,] <- uniqueSols
          
          pay_g[g, ] <- c(groups[g,2], pay_t)
          truedif_g[g, ] <- c(groups[g,2], dif_true_avg)
          truecex_g[g, ] <- c(groups[g,2], cex_true_avg)
          unique_g[g,] <-c(groups[g,2], uniqueSols_avg)
        } # end of looping through all possible groups
        
        #setwd("results")
        savefile_prefix <- paste(here("results/"),'results_local', local, '_sch', simchoice, '_K', K[ik], '_cex', cex,'_rep=',f ,sep='')
        savefile_pay <- paste(savefile_prefix, '_simult_ds_pay.RData', sep='')
        savefile_truedif <- paste(savefile_prefix, '_simult_ds_truedif.RData', sep='')
        savefile_truecex <- paste(savefile_prefix, '_simult_ds_truecex.RData', sep='')
        savefile_unique <- paste(savefile_prefix, '_simult_ds_unique.RData', sep='')
        save(pay_g, file=savefile_pay)
        save(truedif_g, file=savefile_truedif)
        save(truecex_g, file=savefile_truecex)
        save(unique_g, file=savefile_unique)
        toc()
        #setwd("..")
      } # end of looping over change extents
    } # end of looping over different Ks (complexities)
  } # end of looping over different neighbor selection strategies
}