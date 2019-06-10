# Best script in the world, for champions
library(tidyverse)
library(gridExtra)
library(ggthemes)
# Conflict




meiosis4loc <- function(this.dip, haplo.names, r12,r23,r34, this.order = "AMF") {
  recover()
  this.order <- c(strsplit(this.order,"")[[1]],"N")  # this is the order we came in
  reorder    <- order(this.order)[c(1,3,2,4)]
  # order: "A.mat" "M.mat" "F.mat" "N.mat" "A.pat" "M.pat" "F.pat"  "N.pat" "pop" 
  # cathy tomorrow and tonight plan
  meiosis.haps <- 
    # no rec between 3 & 4 
    .5 * (1 - r12)* (1 - r23) * (1 - r34) * (as.numeric(haplo.names == paste(this.dip[paste(this.order,c("mat","mat","mat","mat"),sep=".")][reorder], collapse = ""))) + # no_rec 
    .5 * (1 - r12)* (1 - r23) * (1 - r34) * (as.numeric(haplo.names == paste(this.dip[paste(this.order,c("pat","pat","pat","pat"),sep=".")][reorder], collapse = ""))) + # no_rec   
    #
    .5 * (r12)    * (1 - r23) * (1 - r34) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("mat","pat","pat","pat"),sep=".")][reorder], collapse = ""))) + # 12 rec 
    .5 * (r12)    * (1 - r23) * (1 - r34) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("pat","mat","mat","mat"),sep=".")][reorder], collapse = ""))) + # 12 rec   
    #
    .5 * (1 - r12)* (    r23) * (1 - r34) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("mat","mat","pat","pat"),sep=".")][reorder],collapse = "")))  + # 23_rec 
    .5 * (1 - r12)* (    r23) * (1 - r34) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("pat","pat","mat","mat"),sep=".")][reorder], collapse = ""))) + # 23_rec 
    #
    .5 * (    r12)* (    r23) * (1 - r34) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("mat","pat","mat","mat"),sep=".")][reorder],collapse = ""))) + # both_rec 
    .5 * (    r12)* (    r23) * (1 - r34) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("pat","mat","pat","pat"),sep=".")][reorder],collapse = ""))) + # both_rec  
    #
    #
    # rec between 3 & 4 
    .5 * (1 - r12)* (1 - r23) * (    r34) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("mat","mat","mat","pat"),sep=".")][reorder], collapse = ""))) + # no_rec 
    .5 * (1 - r12)* (1 - r23) * (    r34) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("pat","pat","pat","mat"),sep=".")][reorder], collapse = ""))) + # no_rec   
    #
    .5 * (r12)    * (1 - r23) * (    r34) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("mat","pat","pat","mat"),sep=".")][reorder], collapse = ""))) + # 12 rec 
    .5 * (r12)    * (1 - r23) * (    r34) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("pat","mat","mat","pat"),sep=".")][reorder], collapse = ""))) + # 12 rec   
    #
    .5 * (1 - r12)* (    r23) * (    r34) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("mat","mat","pat","mat"),sep=".")][reorder],collapse = "")))  + # 23_rec 
    .5 * (1 - r12)* (    r23) * (    r34) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("pat","pat","mat","pat"),sep=".")][reorder], collapse = ""))) + # 23_rec 
    #
    .5 * (    r12)* (    r23) * (    r34) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("mat","pat","mat","pat"),sep=".")][reorder],collapse = ""))) + # both_rec 
    .5 * (    r12)* (    r23) * (    r34) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("pat","mat","pat","mat"),sep=".")][reorder],collapse = "")))  # both_rec  
    
  
  return(meiosis.haps)  
}


# yb moved pollenPool into main function for accounting purposes
# pollenPool <- function(diplos, focal.pop, prop.replaced, meiotic.prod){ 
  # we should probably do meiosis before we migrate, YANIV
  # we did this, CATHY
#  home.haps <-  colSums(meiotic.prod * diplos$freqs[diplos$pop == focal.pop] ) * (1-prop.replaced) 
#  away.haps <-  colSums(meiotic.prod * diplos$freqs[diplos$pop != focal.pop] ) * (prop.replaced) # prop.replaced is migration rate
#  home.haps + away.haps
#}

# then we make the mating combos
# we set up who's allowed to mate with who
# my.dips is the females and my.haps is the males
# what are we ending up with here?
matingCombos <- function(my.haps, my.dips, discrim, meiotic.prod){
  how.random  <- 1 - ((rowSums(my.dips[,c("F.mat","F.pat")]) != 0) * discrim) 
  # first take F alleles for mat and pat, if there's a F present in either the sum will be >0
  # for those that have a sum >0, multiply by preference intensity ("discrim") and subtract that amt from 1
  # so this says if you have a 1 at either F locus you'll be discriminated against
  mating.freqs <- matrix(1,ncol = length(my.haps), nrow = length(how.random))
  # making a mat with 8 columns and 64 rows--one col for each dad hap and one row for how discrim each mom geno is
  # these are the frequencies of matings if you're a certain mom/dad combo
  # but right now it's just a matrix of all 1s
  colnames(mating.freqs) <- names(my.haps) # dad names as column names
  rownames(mating.freqs) <- rownames(meiotic.prod) # mom names as row names, taken from meiotic.prod
  mating.freqs           <- t(t(mating.freqs) * my.haps) 
  # transposing 2x makes values in columns the rows, wait this just gives us back 8 cols and 64 rows
  mating.freqs[,grep(pattern = ".0..", names(my.haps))]  <- mating.freqs[,grep(pattern = ".0..", names(my.haps))]  * how.random
  # these are grabbing all the little m's cause we're AMF and multiplying by how permissive each mom is/how random the mating is
  denom <- rowSums(mating.freqs)
  mating.freqs / ifelse( denom  == 0,1,denom)# standardizing so that discriminating moms aren't producing fewer kids
}
# This is just the pollen hap frequency on each mom, not taking into account 

diploidSel <- function(tmp.freq, my.dips, s, focal.pop){
  #n.maladapt <- rowSums(my.dips[,c("A.mat","A.pat")] != focal.pop)
  n.maladapt <- rowSums(my.dips[,c("A.mat","A.pat","N.mat","N.pat")] != focal.pop)
  w          <- (1 - s)^ n.maladapt # fitness = (1-s)^n.maladapt
  wbar       <- sum(tmp.freq *  w) # mean fitness = the sum of t
  wrel       <- w / wbar  #relative fitness
  newfreq    <-  wrel  * tmp.freq
  return(newfreq)
}

migrateMateReproduceSelect <- function(diplos, focal.pop, prop.replaced, meiotic.prod, discrim, s){
#  recover()
  home.haps  <- colSums(meiotic.prod * diplos$freqs[diplos$pop == focal.pop] ) * (1-prop.replaced) 
  away.haps  <- colSums(meiotic.prod * diplos$freqs[diplos$pop != focal.pop] ) * (prop.replaced) # prop.replaced is migration rate
  my.haps    <- home.haps + away.haps# pollen pool
  my.dips    <- diplos[diplos$pop == focal.pop,]                                                                     # females
  pat.haps   <- matingCombos(my.haps = my.haps , my.dips =  my.dips, discrim =discrim, meiotic.prod = meiotic.prod)  # mating rules
  mat.haps   <- meiotic.prod * my.dips$freqs                                                                         # meiosis in feamles
  new.genos  <- t(pat.haps) %*% (mat.haps) # rows are paternal, colums are maternal
  #
  p.migrant  <- away.haps  /  my.haps # proportion of each pollen haplotype that is a migrant.. this is useful for downstream accounting. 
  pop.reinforce <- 1 - sum(rowSums(new.genos) * p.migrant ,na.rm=TRUE) /prop.replaced   # quantifying reinforcement as 1 - prob mating with migrant
  #
  before.sel <- c(t(new.genos)) # flip rows and columns and then flatten to make  new diploid freqsQ
  after.sel  <- diploidSel(tmp.freq = before.sel, my.dips = my.dips, s = s , focal.pop = focal.pop)
  return(c(pop.reinforce,after.sel))
  
  tmp.dips <- data.frame(my.dips, before.sel= before.sel,after.sel=after.sel)
}  



# WHY DID YOU HARDCODE r12 & r23 YANIV # I did not. these are default values. not hardcoded values
runFourLocusSim <-function(n.gen = 1000, r12 = .1, r23 = .3,r34 = .5, 
                           init.freqs = c(fA_0 = 0, fM_0 = 0, fF_0 = 0,fN_0=0,fA_1 = 1, fM_1 = 1, fF_1 = .01, fN_1 = 1), 
                           discrim = 1, s = .4, 
                           prop.replaced0 = .1, 
                           prop.replaced1 = .1,
                           this.order = "AMF",
                           get.blank = FALSE){
  # SETUP
  print(r12)
  # Put in a flip for diff orders later.. the easiest way to do this is only in rec
  diplos <- expand.grid(data.frame(A.mat = c(0:1), M.mat = c(0:1), F.mat = c(0:1),N.mat = c(0,1),
                                   A.pat = c(0,1), M.pat = c(0:1), F.pat = c(0,1),N.pat = c(0,1),
                                   pop = c(0,1)))
  diplos$freqs <-  apply(diplos, 1, function(X){
    this.pop <- X[["pop"]]
    freqs    <- grep( this.pop, names(init.freqs))
    f.big    <- rep(init.freqs[grep( this.pop, names(init.freqs))], 2) 
    prod(as.numeric((1 - X[1:8]) == 0) * f.big  + as.numeric((X[1:8]) == 0) * (1-f.big))
  })
  haplos <- expand.grid(data.frame(A = c(0:1), M = c(0:1),  F = c(0,1), N = c(0,1),pop = c(0,1)))
  haplo.names  <- apply(unique(haplos[,-5]), 1, paste, collapse = "")
  meiotic.prod <- t(apply(unique(diplos[,1:8]),1, meiosis4loc, haplo.names = haplo.names, r12 = r12, r23 = r23, r34 = r34, this.order = this.order) )
  colnames(meiotic.prod) <- haplo.names
  rownames(meiotic.prod) <- apply(unique(diplos[,1:8]), 1, paste, collapse = "")
  geno.time <- matrix(ncol =  nrow(diplos) + 3, nrow = n.gen)
  colnames(geno.time) <- c("gen", apply(unique(diplos[,1:9]), 1, paste, collapse = ""),"reinf_p_zero","reinf_p_one")
  if(get.blank){
    geno.time <- data.frame(geno.time)
    geno.time[,"gen"] <- 1:n.gen
    return(geno.time)
  }
  for(g in 1:n.gen){
#    if(g == 5000){recover()}
    pop0 <- migrateMateReproduceSelect(diplos = diplos, 
                                       focal.pop = 0, 
                                       prop.replaced = prop.replaced0, 
                                       meiotic.prod = meiotic.prod , 
                                       discrim = discrim , 
                                       s = s)
    pop1 <- migrateMateReproduceSelect(diplos = diplos, 
                                       focal.pop = 1, 
                                       prop.replaced = prop.replaced1, 
                                       meiotic.prod = meiotic.prod , 
                                       discrim = discrim, 
                                       s = s)
    diplos$freqs  <- c(pop0[-1], pop1[-1])
    geno.time[g,] <- c(g,diplos$freqs,pop0[1], pop1[1] )
  }
  return(data.frame(geno.time))
}





