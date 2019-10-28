#load required packages
packages <- c('gtools',
              'fGarch',
              'tidyverse',
              'here')

lapply(packages, require, character.only = TRUE)


source("setup/NK_function.R")



N=10

for(K in seq(from = 0, to = N-2, by = 2)) {

#GENERATE LANDSCAPE
LS<-permutations(2,N,v=c(0,1),repeats.allowed=TRUE)
LS<-as.data.frame(LS)
if (K==0){
  depends <- as.vector(1:N)
  values <- replicate(N,round(runif(2,0,1),1))
  fitness <- values
} else {
# INSTEAD OF depends <- rbind(1:N,replicate(N,sample(c(1:N),K,replace=F)))
depends <- matrix(nrow = K+1, ncol = N)
for(i in 1:N) {
  string <- c(1:N)
  string <- string[-i]
  depends[,i] <- c(i,sample(string,K,replace=F))
  }
#
  combinations <- permutations(2,K+1,v=c(0,1),repeats.allowed=TRUE)
  values <- replicate(N,round(runif(nrow(combinations),0,1),1))
  fitness <- cbind(combinations,values)
}

landscape <- generate_landscape(N,K,LS,fitness,depends)
landscape[,N+1] <- (landscape[,N+1]/max(landscape[,N+1]))^8
landscape <- cbind(0:(nrow(landscape)-1),landscape[,N+1])


filename <- paste0(here("landscapes"),"/landscape_",N,'_',K,'.Rds')
write_rds(landscape,filename)
#write.table(landscape, filename, col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

filename <- paste0(here("landscapes"),"/LS_",N,'.Rds')
write_rds(LS,filename)
}