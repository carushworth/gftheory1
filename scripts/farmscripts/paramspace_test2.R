library(readr)

source("/home/crusher/projects/GFs/scripts/gf_sims.R")

params <- commandArgs(trailingOnly = TRUE) %>% as.numeric
str(params)

n.gen <- 10000
#myGFruns <- function(s, prop.replaced0, prop.replaced1, r12, r23, init.M1.freq, n.gen){
  output  <- runGFsim(n.gen = n.gen,  
                      discrim = 1,
                      s   =  params[1],
                      r12 =  params[3],
                      r23 =  params[4],
                      init.freqs =  c(fA_0 = 0, fM_0 = 0, fF_0 = 0,fA_1 = 1, fM_1 = params[5], fF_1 = 0.01),
                      prop.replaced0 = params[2],
                      prop.replaced1 = params[2], 
                      n.unlinked     =  as.numeric(params[6]),
                      delta_hap_components = TRUE)
  save(output, file=paste("CHANGEME",paste(paste("s",params[1],"m",params[2],"r",params[3],"M",params[5]),collapse="_"),".Robj",sep=""))
  write.table(c(output$summary.stats,output$params),file=paste("CHANGEME",paste(paste("s",params[1],"m",params[2],"r",params[3],"M",params[5]),"n.unl",params[6],collapse="_"),".txt",sep="")) 
#}
