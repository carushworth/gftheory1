# Ali
library(dplyr)
# Conflict

# we start with 3 loci that are in various configurations
# in AMF, A is closer to M than F
# in MFA, A is closer to F than M
# in MAF, M and F are far apart and A is equidistant between the two

# we then get a table of meiosis haps for each, where r12 is recombination between the first and second loci
# and r23 is recombination between the second and third

meiosis3loc <- function(this.dip, haplo.names, r12, r23, this.order = "AMF") {
  # Figure out how to change order here
  # order: "A.mat" "M.mat" "F.mat" "A.pat" "M.pat" "F.pat" "pop" 
  # if it was all 0.5 they'd all be on different chromosomes
  if( this.order == "AMF"){
    meiosis.haps <- .5 * (1 - r12)* (1 - r23) *(as.numeric(haplo.names == paste(this.dip[c(1, 2, 3)], collapse = ""))) + # no_rec 
      .5 * (1 - r12)* (1 - r23) *(as.numeric(haplo.names == paste(this.dip[c(4, 5, 6)], collapse = ""))) + # no_rec   
      .5 * (r12)    * (1 - r23) *(as.numeric(haplo.names == paste(this.dip[c(1, 5, 6)], collapse = ""))) + # 12 rec 
      .5 * (r12)    * (1 - r23) *(as.numeric(haplo.names == paste(this.dip[c(4, 2, 3)], collapse = ""))) + # 12 rec   
      .5 * (1 - r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[c(1, 2, 6)], collapse = ""))) + # 23_rec 
      .5 * (1 - r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[c(4, 5, 3)], collapse = ""))) + # 23_rec   
      .5 * (    r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[c(1, 5, 3)], collapse = ""))) + # both_rec 
      .5 * (    r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[c(4, 2, 6)], collapse = "")))   # both_rec  
    return(meiosis.haps)  
  }
  if( this.order == "MFA"){
    meiosis.haps <- .5 * (1 - r12)* (1 - r23) *(as.numeric(haplo.names == paste(this.dip[c(1, 2, 3)], collapse = ""))) + # no_rec 
      .5 * (1 - r12)* (1 - r23) *(as.numeric(haplo.names == paste(this.dip[c(4, 5, 6)], collapse = ""))) + # no_rec   
      .5 * (r12)    * (1 - r23) *(as.numeric(haplo.names == paste(this.dip[c(4, 2, 6)], collapse = ""))) + # 12 rec   
      .5 * (r12)    * (1 - r23) *(as.numeric(haplo.names == paste(this.dip[c(1, 5, 3)], collapse = ""))) + # 12 rec []
      .5 * (1 - r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[c(4, 2, 3)], collapse = ""))) + # 23_rec 
      .5 * (1 - r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[c(1, 5, 6)], collapse = ""))) + # 23_rec   
      .5 * (    r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[c(1, 2, 6)], collapse = ""))) + # both_rec 
      .5 * (    r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[c(4, 5, 3)], collapse = "")))   # both_rec  
    return(meiosis.haps)  
  }
  if( this.order == "FAM"){
    meiosis.haps <- .5 * (1 - r12)* (1 - r23) *(as.numeric(haplo.names == paste(this.dip[c(1, 2, 3)], collapse = ""))) + # no_rec 
      .5 * (1 - r12)* (1 - r23) *(as.numeric(haplo.names == paste(this.dip[c(4, 5, 6)], collapse = ""))) + # no_rec   
      .5 * (r12)    * (1 - r23) *(as.numeric(haplo.names == paste(this.dip[c(4, 5, 3)], collapse = ""))) + # 12 rec   
      .5 * (r12)    * (1 - r23) *(as.numeric(haplo.names == paste(this.dip[c(1, 2, 6)], collapse = ""))) + # 12 rec 
      .5 * (1 - r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[c(1, 5, 3)], collapse = ""))) + # 23_rec 
      .5 * (1 - r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[c(4, 2, 6)], collapse = ""))) + # 23_rec   
      .5 * (    r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[c(4, 2, 3)], collapse = ""))) + # both_rec 
      .5 * (    r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[c(1, 5, 6)], collapse = "")))   # both_rec  
    return(meiosis.haps)  
  }
}

# here we're making a function giving us the pollen pool 
# diplos is a df with freqs and a pop 
# we make a thing called home and a thing called away and in each we select the freqs with As that either match (==) or don't match (!=)
# the focal.pop and multiply by this thing prop.replaced (or not being replaced if match the focal pop)
# then sum them and multiply by meiotic.prod (this is the product of meiosis? what is this)


pollenPool <- function(diplos, focal.pop, prop.replaced, meiotic.prod){ # we should probably do meiosis before we migrate, YANIV
  home <- diplos$freqs[diplos$pop == focal.pop] * (1-prop.replaced) 
  away <- diplos$freqs[diplos$pop != focal.pop] * (prop.replaced)  
  colSums((home + away) *  meiotic.prod)
}

# then we make the mating combos
# my.dips is the females and my.haps is the males
# we have a vector(?) of something that is ... stuck on how.discrim
# 
matingCombos <- function(my.haps, my.dips, discrim, meiotic.prod){
  how.discrim  <- 1 - ((rowSums(my.dips[,c("F.mat","F.pat")]) != 0) * discrim) 
  mating.freqs <- matrix(1,ncol = length(my.haps), nrow = length(how.discrim))
  colnames(mating.freqs) <- names(my.haps) # dad names as column names
  rownames(mating.freqs) <- rownames(meiotic.prod) # mom names as row names
  mating.freqs           <- t(t(mating.freqs) * my.haps) # transposing 2x makes values in columns the rows
  mating.freqs[,grep(pattern = ".0.", names(my.haps))]  <- mating.freqs[,grep(pattern = ".0.", names(my.haps))]  * how.discrim
  mating.freqs / rowSums(mating.freqs) # standardizing so that discriminating moms aren't producing fewer kids
}

# local adaptation step
diploidSel <- function(tmp.freq, geno.labels, s, focal.pop){
  n.maladapt <- rowSums(do.call(rbind,strsplit(geno.labels,""))[,c(1,4)] != focal.pop) # doesn't generalize to different marker order
  w          <- (1 - s)^ n.maladapt
  wbar       <- sum(tmp.freq *  w)
  tmp.freq *  w / wbar
}

# WHY DID YOU HARDCODE discrim AND prop.replaced YANIV
migrateMateReproduceSelect <- function(diplos, focal.pop, prop.replaced, meiotic.prod, discrim, s){
  recover()
  my.haps    <- pollenPool(diplos = diplos, focal.pop = focal.pop, prop.replaced = .3, meiotic.prod = meiotic.prod ) # pollen pool
  my.dips    <- diplos[diplos$pop == focal.pop,]                                                                     # females
  pat.haps   <- matingCombos(my.haps = my.haps , my.dips =  my.dips, discrim =1, meiotic.prod = meiotic.prod)        # mating rules
  mat.haps   <- meiotic.prod * my.dips$freqs                                                                         # meiosis in feamles
  new.genos  <- t(pat.haps) %*% (mat.haps) # rows are paternal, colums are maternal ## YANIV SAYS HE'S NOT SURE THIS IS RIGHT
  before.sel <- c(t(new.genos)) # flip rows and columns and then flatten to make  new diploid freqsQ
  after.sel  <- diploidSel(tmp.freq = before.sel, geno.labels = rownames(mat.haps), s = s , focal.pop = focal.pop)
  return(after.sel)
}  

# WHY DID YOU HARDCODE r12 & r23 YANIV (although that actually seems ok)
runThreLocusSim <-function(n.gen = 1000, init.freqs = c(fA_0 = 0, fM_0 = 0, fF_0 = 0,fA_1 = 1, fM_1 = 1, fF_1 = .001, discrim = 1)){
  # SETUP
  # Put in a flip for diff orders later.. the easiest way to do this is only in rec
  diplos <- expand.grid(data.frame(A.mat = c(0:1), M.mat = c(0:1), F.mat = c(0:1),
                                 A.pat = c(0,1), M.pat = c(0:1), F.pat = c(0,1),
                                 pop = c(0,1)))
  diplos$freqs <-  apply(diplos, 1, function(X){
    this.pop <- X[["pop"]]
    freqs    <- grep( this.pop, names(init.freqs))
    f.big    <- rep(init.freqs[grep( this.pop, names(init.freqs))], 2) 
    prod(as.numeric((1 - X[1:6]) == 0) * f.big  + as.numeric((X[1:6]) == 0) * (1-f.big))
  })
  haplos <- expand.grid(data.frame(A = c(0:1), M = c(0:1),  F = c(0,1), pop = c(0,1)))
  haplo.names  <- apply(unique(haplos[,-4]), 1, paste, collapse = "")
  meiotic.prod <- t(apply(unique(diplos[,1:6]),1, meiosis3loc, haplo.names = haplo.names, r12 = .5, r23 = .5) )
  colnames(meiotic.prod) <- haplo.names
  rownames(meiotic.prod) <- apply(unique(diplos[,1:6]), 1, paste, collapse = "")
  geno.time <- matrix(ncol =  nrow(diplos) + 1, nrow = n.gen)
    colnames(geno.time) <- c("gen", apply(unique(diplos[,1:7]), 1, paste, collapse = ""))

  for(g in 1:n.gen){
    pop0 <- migrateMateReproduceSelect(diplos = diplos, 
                                       focal.pop = 0, 
                                       prop.replaced = .1, 
                                       meiotic.prod = meiotic.prod , 
                                       discrim = discrim , 
                                       s = .4)
    pop1 <- migrateMateReproduceSelect(diplos = diplos, 
                                       focal.pop = 1, 
                                       prop.replaced = .1, 
                                       meiotic.prod = meiotic.prod , 
                                       discrim = discrim, 
                                       s = .4)
    diplos$freqs  <- c(pop0, pop1)
    geno.time[g,] <- c(g,diplos$freqs )
  }
  return(data.frame(geno.time))
}

  
runThreLocusSim()
