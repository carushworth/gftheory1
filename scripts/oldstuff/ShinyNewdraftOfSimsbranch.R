# Best script in the world, for champions

library(tidyverse)
library(gridExtra)
library(ggthemes)
# Conflict


meiosis3loc <- function(this.dip, haplo.names, r12,r23, this.order = "AMF") {
  # Figure out how to change order here
  this.order <- strsplit(this.order,"")[[1]]
  reorder <- order(this.order)[c(1,3,2)]
  # order: "A.mat" "M.mat" "F.mat" "A.pat" "M.pat" "F.pat" "pop" 
  # cathy tomorrow and tonight plan
  meiosis.haps <- 
    .5 * (1 - r12)* (1 - r23) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("mat","mat","mat"),sep=".")][reorder], collapse = ""))) + # no_rec 
    .5 * (1 - r12)* (1 - r23) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("pat","pat","pat"),sep=".")][reorder], collapse = ""))) + # no_rec   
    #
    .5 * (r12)    * (1 - r23) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("mat","pat","pat"),sep=".")][reorder], collapse = ""))) + # 12 rec 
    .5 * (r12)    * (1 - r23) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("pat","mat","mat"),sep=".")][reorder], collapse = ""))) + # 12 rec   
    #
    .5 * (1 - r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("mat","mat","pat"),sep=".")][reorder],collapse = ""))) + # 23_rec 
    .5 * (1 - r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("pat","pat","mat"),sep=".")][reorder], collapse = ""))) + # 23_rec 
    #
    .5 * (    r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("mat","pat","mat"),sep=".")][reorder],collapse = ""))) + # both_rec 
    .5 * (    r12)* (    r23) *(as.numeric(haplo.names == paste(this.dip[paste(this.order,c("pat","mat","pat"),sep=".")][reorder],collapse = "")))   # both_rec  
  return(meiosis.haps)  
}

pollenPool <- function(diplos, focal.pop, prop.replaced, meiotic.prod){ # we should probably do meiosis before we migrate, YANIV
  home.haps <-  colSums(meiotic.prod * diplos$freqs[diplos$pop == focal.pop] ) * (1-prop.replaced) 
  away.haps <-  colSums(meiotic.prod * diplos$freqs[diplos$pop != focal.pop] ) * (prop.replaced) # prop.replaced is migration rate
  home.haps + away.haps
}

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
  mating.freqs[,grep(pattern = ".0.", names(my.haps))]  <- mating.freqs[,grep(pattern = ".0.", names(my.haps))]  * how.random
  # these are grabbing all the little m's cause we're AMF and multiplying by how permissive each mom is/how random the mating is
  mating.freqs / rowSums(mating.freqs) # standardizing so that discriminating moms aren't producing fewer kids
}
# This is just the pollen hap frequency on each mom, not taking into account 

diploidSel <- function(tmp.freq, geno.labels, s, focal.pop){
  n.maladapt <- rowSums(do.call(rbind,strsplit(geno.labels,""))[,c(1,4)] != focal.pop)
  w          <- (1 - s)^ n.maladapt # fitness = (1-s)^n.maladapt
  wbar       <- sum(tmp.freq *  w) # mean fitness = the sum of t
  tmp.freq *  w / wbar
}

migrateMateReproduceSelect <- function(diplos, focal.pop, prop.replaced, meiotic.prod, discrim, s){
  my.haps    <- pollenPool(diplos = diplos, focal.pop = focal.pop, prop.replaced = prop.replaced, meiotic.prod = meiotic.prod ) # pollen pool
  my.dips    <- diplos[diplos$pop == focal.pop,]                                                                     # females
  pat.haps   <- matingCombos(my.haps = my.haps , my.dips =  my.dips, discrim =discrim, meiotic.prod = meiotic.prod)        # mating rules
  mat.haps   <- meiotic.prod * my.dips$freqs                                                                         # meiosis in feamles
  new.genos  <- t(pat.haps) %*% (mat.haps) # rows are paternal, colums are maternal
  before.sel <- c(t(new.genos)) # flip rows and columns and then flatten to make  new diploid freqsQ
  after.sel  <- diploidSel(tmp.freq = before.sel, geno.labels = rownames(mat.haps), s = s , focal.pop = focal.pop)
  return(after.sel)
}  



# WHY DID YOU HARDCODE r12 & r23 YANIV # I did not. these are default values. not hardcoded values
runThreLocusSim <-function(n.gen = 1000, r12 = .1, r23 = .3,
                           init.freqs = c(fA_0 = 0, fM_0 = 0, fF_0 = 0,fA_1 = 1, fM_1 = 1, fF_1 = .01), 
                           discrim = 1, s = .4, prop.replaced0 = .1, prop.replaced1 = .1){
  # SETUP
  print(r12)
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
  meiotic.prod <- t(apply(unique(diplos[,1:6]),1, meiosis3loc, haplo.names = haplo.names, r12 = r12, r23 = r23) )
  colnames(meiotic.prod) <- haplo.names
  rownames(meiotic.prod) <- apply(unique(diplos[,1:6]), 1, paste, collapse = "")
  geno.time <- matrix(ncol =  nrow(diplos) + 1, nrow = n.gen)
  colnames(geno.time) <- c("gen", apply(unique(diplos[,1:7]), 1, paste, collapse = ""))
  
  for(g in 1:n.gen){
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
    diplos$freqs  <- c(pop0, pop1)
    geno.time[g,] <- c(g,diplos$freqs )
  }
  return(data.frame(geno.time))
}


#AMF   
output        <- runThreLocusSim(n.gen = 10000, discrim = 1, s=.5, r12=.01, r23 = 0) %>% tidying()
output_freq   <- findFreqs(output)
output_LD     <- findLD(output)

geno_freq_plot <- output_freq %>% 
  select(- am, - mf, - af)    %>%  
  gather(key = locus, value = freq, - pop, -gen) %>%
  ggplot(aes(x=gen,y=freq,color=locus )) + 
  geom_line() + 
  facet_wrap(~pop, labeller = "label_both")+
  scale_x_continuous(trans = "log10") +
  ggtitle("Allele frequencies over time")+
  theme_tufte()
  
  
  
  ggplot(output, aes(x=gen,y=freq,color=pollen_style, group = haplo,  linetype = local_adapt )) + 
  geom_line() + 
  facet_wrap(~pop, labeller = "label_both")+
  scale_x_continuous(trans = "log10") +
  ggtitle("Haplotype frequencies over time")+
  theme_tufte()

ld_plot  <- output_LD %>% 
  select(pop,gen,D_am, D_mf, D_af) %>%
  gather(key = pair, value = D, - pop, - gen) %>%
  ggplot(aes(x = gen, y = D, color = pair)) + 
  geom_line()+
  facet_wrap(~pop, labeller = "label_both")+ 
  scale_x_continuous(trans = "log10")+
  ggtitle("Linkage disequilibrium over time")+
  theme_tufte()

haplo_freq_plot <- ggplot(output, aes(x=gen,y=freq,color=pollen_style, group = haplo,  linetype = local_adapt )) + 
  geom_line() + 
  facet_wrap(~pop, labeller = "label_both")+
  scale_x_continuous(trans = "log10") +
  ggtitle("Haplotype frequencies over time")+
  theme_tufte()


grid.arrange(geno_freq_plot, ld_plot,haplo_freq_plot, ncol = 1)
