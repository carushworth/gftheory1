library(parallel)
library(readr)
# source("gf_sims.R")
source("https://raw.githubusercontent.com/carushworth/gftheory1/master/scripts/gf_sims.R")


n.gen <- 500

# ran below separately to make a txt file
# jeff's param combos 8/1
param.combos <- expand.grid(
  s = c(0.01, 0.1, 0.9) ,
  m = c(0.01,0.1) ,
  r12 = c(0.001,0.5),
  r23 = c(0.001,0.5),
  init.M1.freq = c(0.05,1),
  n_unlinked = c(1,4)
)

write.table(param.combos,file="param_combos_nunlinked.txt",sep="\t",row.names=FALSE,col.names=FALSE)

# yaniv's param combos 8/12
param.combos <- expand.grid(
  s = seq(0,1,.1),
  m12 = seq(0,.2,.01),
  m21 = seq(0,.2,.01), 
  r = 0.0001,
  M0  = 1, 
  n.unlinked = 0)

write.table(param.combos,file="param_combos_YB.txt",sep="\t",row.names=FALSE,col.names=FALSE)


### for R script ###

summaries <- as_tibble(do.call(rbind,mclapply(X = seq_along(param.combos[,1] ), function(Z){
  PARAMS <- param.combos[Z,]
  output  <- runGFsim(n.gen = n.gen, 
                      discrim = 1,
                      s   =     PARAMS[["s"]],
                      r12 =     PARAMS[["r12"]],
                      r23 =     PARAMS[["r23"]],
                      init.freqs =  c(fA_0 = 0, fM_0 = 0, fF_0 = 0,
                                      fA_1 = 1, fM_1 = 1, fF_1 = PARAMS[["init_freq"]]),
                      prop.replaced0 = PARAMS[["m"]],
                      prop.replaced1 = PARAMS[["m"]], 
                      n.unlinked     = PARAMS[["n_unlinked"]] 
  )
  save(output, file =  paste("CHANGEME",paste(paste(names(PARAMS), PARAMS, sep = "="),collapse = "_"),".csv",sep=""))
  c(output$summary.stats,output$params)
})))  %>% 
  mutate_at(setdiff(colnames(summaries),"this.order"), as.numeric  )

write_csv(x = summaries, "CHANGEME.csv")


