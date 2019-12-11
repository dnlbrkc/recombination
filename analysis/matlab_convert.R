library(R.matlab)
library(here)
setwd(here("results"))

for(l in c(1)){
  for(s in c(0,2,3,4,5,6,7,8,9)){ # 
    for(k in c(0,4,8)){
      for(c in c(1,2,3,4,5)){
        
            #pay_g
            mat <- matrix(0,nrow=11,ncol=31)
            files <- list.files(pattern=c(paste0("results_local",l,"_sch",s,"_K",k,
                                     "_cex",c,"_rep_db=\\d{1,3}","_simult_ds_pay")))
            mat <- array(0,dim=c(11,31,length(files))) #nrow=11,ncol=51)
            c=0
            for(f in files){
             c=c+1
               load(f)
              mat[,,c] <- mat[,,c] + pay_g
            }
            mat <- mat / length(files)
            #pay_g <- apply(mat,c(1,2),sd)
            pay_g <- mat
            path <- paste0("~/Desktop/matlab/","results_local",l,"_sch",s,"_K",k,
                             "_cex",c,"_simult_ds_pay.mat")
            
            writeMat(con=path,pay_g=pay_g)
            
            
            
            #truedif
            mat <- matrix(0,nrow=11,ncol=31)
            files <- list.files(pattern=c(paste0("results_local",l,"_sch",s,"_K",k,
                                                 "_cex",c,"_rep_db=\\d{1,3}","_simult_ds_truedif")))
            for(f in files){
              load(f)
              mat <- mat + truedif_g
            }
            mat <- mat / length(files)
            truedif_g <- mat
            path <- paste0("~/Desktop/matlab/","results_local",l,"_sch",s,"_K",k,
                           "_cex",c,"_simult_ds_truedif.mat")
            
            writeMat(con=path,truedif_g=truedif_g)
            
            #truecex
            mat <- matrix(0,nrow=11,ncol=31)
            files <- list.files(pattern=c(paste0("results_local",l,"_sch",s,"_K",k,
                                                 "_cex",c,"_rep_db=\\d{1,3}","_simult_ds_truecex")))
            for(f in files){
              load(f)
              mat <- mat + truecex_g
            }
            mat <- mat / length(files)
            truecex_g <- mat
            path <- paste0("~/Desktop/matlab/","results_local",l,"_sch",s,"_K",k,
                           "_cex",c,"_simult_ds_truecex.mat")
            
            writeMat(con=path,truecex_g=truecex_g)
            
            
            #truecex
            mat <- matrix(0,nrow=11,ncol=31)
            files <- list.files(pattern=c(paste0("results_local",l,"_sch",s,"_K",k,
                                                 "_cex",c,"_rep_dbmatl=\\d{1,3}","_simult_ds_unique")))
            for(f in files){
              load(f)
              mat <- mat + unique_g
            }
            mat <- mat / length(files)
            unique_g <- mat
            path <- paste0("~/Desktop/matlab/","results_local",l,"_sch",s,"_K",k,
                           "_cex",c,"_simult_ds_unique.mat")
            
            writeMat(con=path,unique_g=unique_g)
            
      }
    }
  }
}




