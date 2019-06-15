library(tidyverse)

# Meiosis

# Meiosis for the gf complex
meiosisGFcomplex <- function(this.dip, tmp.haplo.names, r12,r23, this.order = "AMF") {
  this.order <- strsplit(this.order,"")[[1]]
  reorder <- order(this.order)[c(1,3,2)]
  # oder in dataframe order: "A.mat" "M.mat" "F.mat" "A.pat" "M.pat" "F.pat" "pop" ... 
  # physical order AMF 
  # genet
  meiosis.haps <- 
    .5 * (1 - r12)* (1 - r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("mat","mat","mat"),sep=".")][reorder], collapse = ""))) + # no_rec 
    .5 * (1 - r12)* (1 - r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("pat","pat","pat"),sep=".")][reorder], collapse = ""))) + # no_rec   
    #
    .5 * (r12)    * (1 - r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("mat","pat","pat"),sep=".")][reorder], collapse = ""))) + # 12 rec 
    .5 * (r12)    * (1 - r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("pat","mat","mat"),sep=".")][reorder], collapse = ""))) + # 12 rec   
    #
    .5 * (1 - r12)* (    r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("mat","mat","pat"),sep=".")][reorder],collapse = ""))) + # 23_rec 
    .5 * (1 - r12)* (    r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("pat","pat","mat"),sep=".")][reorder], collapse = ""))) + # 23_rec 
    #
    .5 * (    r12)* (    r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("mat","pat","mat"),sep=".")][reorder],collapse = ""))) + # both_rec 
    .5 * (    r12)* (    r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("pat","mat","pat"),sep=".")][reorder],collapse = "")))   # both_rec  
  return(meiosis.haps)  
}

# meiosis for unlinked loci
unlinkedMeiosis <- function(tmp.diplo, tmp.haplo.names, haplo.vals){
  tmp.geno <-  abs(1-(tmp.diplo[grep("mat",names(tmp.diplo))] + tmp.diplo[grep("pat",names(tmp.diplo))] )/2)
  apply(abs(tmp.geno - t(haplo.vals)) ,2,prod)
}



# whom mates with whome
matingCombos <- function(my.haps, my.dips, discrim, meiotic.prod){
  how.random  <- 1 - ((rowSums(my.dips[,c("F.mat","F.pat")]) != 0) * discrim) 
  # first take F alleles for mat and pat, if there's a F present in either the sum will be >0
  # for those that have a sum >0, multiply by preference intensity ("discrim") and subtract that amt from 1
  # so this says if you have a 1 at either F locus you'll be discriminated against
  mating.freqs <- matrix(1,ncol = length(my.haps), nrow = length(how.random))
  # these are the frequencies of matings if you're a certain mom/dad combo
  # but right now it's just a matrix of all 1s
  colnames(mating.freqs) <- names(my.haps) # dad names as column names
  rownames(mating.freqs) <- rownames(meiotic.prod) # mom names as row names, taken from meiotic.prod
  mating.freqs           <- t(t(mating.freqs) * my.haps) 
  # transposing 2x makes values in columns the rows, wait this just gives us back 8 cols and 64 rows
  mating.freqs[,grep(pattern = ".0.", substring(names(my.haps),first = 1, last=3))]  <- mating.freqs[,grep(pattern = ".0.", substring(names(my.haps),first = 1, last=3))]  * how.random
  # these are grabbing all the little m's cause we're AMF and multiplying by how permissive each mom is/how random the mating is
  denom <- rowSums(mating.freqs)
  mating.freqs / ifelse( denom  == 0,1,denom)# standardizing so that discriminating moms aren't producing fewer kids
}


diploidSel <- function(tmp.freq, my.dips, s, focal.pop){
  s_linked <- s_unlinked <- s # Currently selection on linked and unlinked loci is eaual... this can be changed
  n.maladapt_unlinked <- rowSums(my.dips[,grep("U",colnames(my.dips))] != focal.pop)
  n.maladapt_linked <- rowSums(my.dips[,grep("A",colnames(my.dips))] != focal.pop)
  w          <- (1 - s_unlinked)^ n.maladapt_unlinked * (1-s_linked)^n.maladapt_linked 
  wbar       <- sum(tmp.freq *  w) # mean fitness = the sum of t
  wrel       <- w / wbar  #relative fitness
  newfreq    <-  wrel  * tmp.freq
  return(newfreq)
}

migrateMateReproduceSelect <- function(diplos, focal.pop, prop.replaced, meiotic.prod, discrim, s, delta_hap_components, hap.ids){
  dhap <- rep(NA, 4*length(hap.ids))
  home.haps  <- colSums(meiotic.prod * diplos$freqs[diplos$pop == focal.pop] ) * (1-prop.replaced) 
  away.haps  <- colSums(meiotic.prod * diplos$freqs[diplos$pop != focal.pop] ) * (prop.replaced) # prop.replaced is migration rate
  my.haps    <- home.haps + away.haps# pollen pool
  if(delta_hap_components){
    init.hap <- colSums(meiotic.prod * diplos$freqs[diplos$pop == focal.pop])
    # change in hap freq by migrations 
    # change in pollen       no change in ovules   average
    dp.migration.pollen <- my.haps - init.hap
  }
  my.dips    <- diplos[diplos$pop == focal.pop,]                                                                     # females
  pat.haps   <- matingCombos(my.haps = my.haps , my.dips =  my.dips, discrim =discrim, meiotic.prod = meiotic.prod)  # mating rules
  mat.haps   <- meiotic.prod * my.dips$freqs                                                                         # meiosis in feamles
  new.genos  <- t(pat.haps) %*% (mat.haps) # rows are paternal, colums are maternal
  #
  if(delta_hap_components){
    # change in hap freq by migrations 
    # change in pollen       no change in ovules   average
    dp.mating.pollen <- rowSums(new.genos) - my.haps
  }
  p.migrant  <- away.haps  /  my.haps # proportion of each pollen haplotype that is a migrant.. this is useful for downstream accounting. 
  pop.reinforce <- 1 - sum(rowSums(new.genos) * p.migrant ,na.rm=TRUE) /prop.replaced   # quantifying reinforcement as 1 - prob mating with migrant
  #
  before.sel <- c(t(new.genos)) # flip rows and columns and then flatten to make  new diploid freqsQ
  after.sel  <- diploidSel(tmp.freq = before.sel, my.dips = my.dips, s = s , focal.pop = focal.pop)
  if(delta_hap_components){
    dhap <- full_join(tibble(
      migration_pollen = dp.migration.pollen, 
      mating_pollen    = dp.mating.pollen,
      hap              = names(dp.mating.pollen),
      gf_hap           = str_sub(names( dp.mating.pollen ), 1 ,3)) %>%
      group_by(gf_hap) %>%
      summarise(migration_pollen = sum(migration_pollen), 
                mating_pollen    = sum(mating_pollen)),
    # change in hap freq by migrations 
    # change in pollen       no change in ovules   average
    my.dips %>% 
      mutate(gf_hap  = paste(A.pat, M.pat, F.pat, sep = ""),
             before.sel = before.sel, 
             after.sel = after.sel) %>% 
      group_by(gf_hap) %>% 
      summarise(dpsel_pat = sum(after.sel) - sum(before.sel)), 
    by = "gf_hap")
    dhap <- full_join(dhap,my.dips %>% 
      mutate(gf_hap  = paste(A.mat, M.mat, F.mat, sep = ""),
             before.sel = before.sel, 
             after.sel  = after.sel) %>% 
      group_by(gf_hap) %>% 
      summarise(dpsel_mat = sum(after.sel) - sum(before.sel)) ,
      by = "gf_hap")
    rownames_dhap <- dhap$gf_hap
    dhap <- as.matrix(dhap[,-1])
    rownames(dhap) <- rownames_dhap
    dhap <- dhap[hap.ids,]
  }
  return(c(pop.reinforce, after.sel, c(dhap)))
}  



# WHY DID YOU HARDCODE r12 & r23 YANIV # I did not. these are default values. not hardcoded values
runGFsim <-function(n.gen = 1000, r12 = .1, r23 = .3,r34 = .5, 
                           init.freqs = c(fA_0 = 0, fM_0 = 0, fF_0 = 0,fA_1 = 1, fM_1 = 1, fF_1 = .01), 
                           discrim = 1, s = .4, 
                           prop.replaced0 = .1, 
                           prop.replaced1 = .1,
                           this.order = "AMF",
                           n.unlinked =0,
                           get.blank = FALSE,
                           delta_hap_components = FALSE){
  # SETUP
  print(sprintf("s = %s, m = %s, r12 = %s, r23 = %s, init_freq = %s", s, prop.replaced0, r12, r23, init.freqs["fF_1"]))
  diplos <- expand.grid(data.frame(rbind(numeric(length = 2*(3+n.unlinked)+1),1)))
  tmp.names <- c("A","M","F",paste("U", 0:n.unlinked,sep ="")[-1])
  colnames(diplos) <- c(paste(rep(tmp.names,times = 2), rep(c("mat","pat"), each = length(tmp.names)), sep="."),"pop")
  diplos$freqs <-  apply(diplos, 1, function(X){
    this.pop        <- X[["pop"]]
    freqs           <- grep( this.pop, names(init.freqs))
    tmp.big         <- c(init.freqs[grep( this.pop, names(init.freqs))],rep(this.pop,n.unlinked))
    names(tmp.big)  <- c(names(tmp.big)[1:3],paste(paste(rep("U",n.unlinked), (0:n.unlinked) , sep=""),this.pop,sep="_"  )[-1]) # assumes we begin as locally adapted at all loci
    f.big           <- rep(tmp.big, 2) 
    prod(as.numeric((1 - X[-length(X)]) == 0) * f.big  + as.numeric((X[-length(X)]) == 0) * (1-f.big))
  })
  haplos <- expand.grid(data.frame(rbind(numeric(length = (3+n.unlinked)+1),1)))
  names(haplos) <- c("A","M","F",paste("U", 0:n.unlinked,sep ="")[-1],"pop")
  haplo.names  <-  apply(unique(haplos[,-ncol(haplos)]), 1, paste, collapse = "")
  #
  ### MEIOSIS
  # First worry about our linked complex
  gf.complex <- unique(diplos[,grep("A|M|F", colnames(diplos))])
  tmp.gf     <- apply(diplos[diplos$pop == 0,grep("A|M|F", colnames(diplos))],1, paste,collapse="") 
  three.loc.meiosis <- t(apply(gf.complex,1, meiosisGFcomplex, tmp.haplo.names = unique(str_sub(haplo.names,1,3)), r12 = r12, r23 = r23,  this.order = this.order) )
  rownames(three.loc.meiosis) <- apply(gf.complex,1, paste,collapse="")
  colnames(three.loc.meiosis) <- unique(str_sub(haplo.names,1,3))
  three.loc.meiosis           <- three.loc.meiosis[ tmp.gf,str_sub(haplo.names,1,3)]
  colnames(three.loc.meiosis) <- haplo.names
  rm(gf.complex ,tmp.gf)
  if(n.unlinked == 0){
    meiotic.prod <- three.loc.meiosis
    rm(three.loc.meiosis)
  }
  # Now meioisis at unlinked loci
  if(n.unlinked > 0){
    unlinked.complex <- unique(diplos[,grep("U", colnames(diplos))])
    tmp.complex      <- apply(diplos[diplos$pop == 0,grep("U", colnames(diplos))],1, paste,collapse="") 
    unlinked.haps    <- unique(str_sub(haplo.names,start = 4))
    unlinked.meiosis <- t(apply(X = unlinked.complex,
                                MARGIN = 1, 
                                FUN = unlinkedMeiosis, 
                                tmp.haplo.names = unlinked.haps, 
                                haplo.vals = apply(do.call(rbind,strsplit(unlinked.haps,"")), 2,as.numeric)
                                ))
    rownames(unlinked.meiosis) <- apply(unlinked.complex,1, paste,collapse="")
    colnames(unlinked.meiosis) <- unlinked.haps
    unlinked.meiosis           <- unlinked.meiosis[tmp.complex ,str_sub(haplo.names,4)]
    meiotic.prod               <- three.loc.meiosis * unlinked.meiosis 
    rownames( meiotic.prod)    <- apply(diplos[diplos$pop == 0,1:(ncol(diplos)-2)], 1, paste, collapse = "")
    rm(unlinked.meiosis, three.loc.meiosis,tmp.complex,unlinked.haps)
  }
  ###
  names.geno.time       <- c(apply(unique(diplos[,grep("U|q",names(diplos),invert = T)]), 1, paste, collapse = ""))
  geno.time             <- matrix(ncol = 3 + length(names.geno.time), nrow = n.gen)
  colnames(geno.time)   <- c("gen","reinf_0","reinf_1", names.geno.time)
  meanUs <- geno.time
  colnames(meanUs)[c(2:3)] <- c("U_0","U_1")
  ##
  hap.ids <- unique(str_sub(haplo.names,1,3))
  dhap_components <- matrix(ncol = 1 + 2 * length(hap.ids) * 4, nrow = n.gen)
  colnames(dhap_components) <- c("gen",paste(rep(c("p0","p1"), each=32),rep(paste(hap.ids,rep(c("migrationPollen", "matingPollen", "selPat", "selMat"),each = 8),sep="_"),2),sep="_"))
  #
  if(get.blank){
    geno.time <- data.frame(geno.time)
    meanUs    <- data.frame(meanUs)
    geno.time[,"gen"] <- 1:n.gen
    meanUs[,"gen"] <- 1:n.gen
    return(list(geno.time = geno.time, meanUs = meanUs))
  }
  for(g in 1:n.gen){
    #    if(g == 5000){recover()}
    pop0 <- migrateMateReproduceSelect(diplos = diplos, 
                                       focal.pop = 0, 
                                       prop.replaced = prop.replaced0, 
                                       meiotic.prod = meiotic.prod , 
                                       discrim = discrim , 
                                       s = s,
                                       delta_hap_components = delta_hap_components,
                                       hap.ids = hap.ids)
    pop1 <- migrateMateReproduceSelect(diplos = diplos, 
                                       focal.pop = 1, 
                                       prop.replaced = prop.replaced1, 
                                       meiotic.prod = meiotic.prod , 
                                       discrim = discrim, 
                                       s = s,
                                       delta_hap_components = delta_hap_components,
                                       hap.ids = hap.ids)
    dhap_components[g,]<- c(g,pop0[(length(pop0)-31):length(pop0)],pop1[(length(pop1)-31):length(pop1)])
    pop0 <- pop0[-((length(pop0)-31):length(pop0))] 
    pop1 <- pop1[-((length(pop1)-31):length(pop1))] 
    print(g)
    diplos$freqs  <- c(pop0[-1], pop1[-1])
    if(n.unlinked == 0){ geno.time[g,]  <- c(g, pop0[1], pop1[1],  diplos$freqs)  }
    if(n.unlinked > 0){
    tmp.diplos <- diplos %>% 
      group_by(A.mat, M.mat, F.mat, A.pat ,M.pat, F.pat, pop) %>%
      summarise(freqs = sum(freqs)) %>%
      ungroup() %>%
      mutate(geno = paste(A.mat, M.mat, F.mat, A.pat, M.pat, F.pat,   pop, sep = "")) 
    
    # FIX THIS
    geno.time[g,]  <- c(g, pop0[1], pop1[1], 
                        unlist(tmp.diplos %>% select(freqs, geno ) %>% 
                                 spread(key = geno, value = freqs))[names.geno.time]   )
    
    tmp.diplos2 <- diplos %>% 
      mutate(mean_U = diplos %>% select(starts_with("U")) %>% rowMeans() )      %>% 
      group_by(A.mat, M.mat, F.mat, A.pat ,M.pat, F.pat, pop) %>%
      summarise(U = sum(freqs * mean_U), freqs = sum(freqs))       %>% 
      ungroup() 
    
    meanUs[g,] <- c(g, 
      tmp.diplos2 %>% group_by(pop) %>% summarise(x = sum(U)) %>%pull(),
      unlist(tmp.diplos2 %>%
               mutate(geno = paste(A.mat, M.mat, F.mat, A.pat, M.pat, F.pat,   pop, sep = ""),
                      meanGenoU    = U / freqs) %>% 
               select(meanGenoU, geno ) %>% 
               spread(key = geno, value = meanGenoU))[names.geno.time])
    }
  }
  return(list(geno.time = data.frame(geno.time), 
              meanUs = data.frame(meanUs),
              dhaps  = data.frame(dhap_components)))
}



