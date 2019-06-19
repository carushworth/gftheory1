# fully worked example
source("gf_sims.R")
source("functions_forsummaries.R")
source("functions_forplots.R")
example.sim <- runGFsim(n.gen = 100, discrim = 1, s   =     .2,
         r12 =     .001, r23 =     .001, 
         init.freqs =  c(fA_0 = 0, fM_0 = 0, fF_0 = 0, fA_1 = 1, fM_1 = 1, fF_1 = 0.01),
         prop.replaced0 = .01,   prop.replaced1 = .01,  n.unlinked     = 2,delta_hap_components = TRUE ) 

#############
# Frequencies and reinfrocement
# haplotypes
tidy.haps           <- tidyingHaps(example.sim$geno.time) # local adapt strictly refers to the a locus, not unlinked loci

# alleles and ld
allele.freqs.and.ld <- findFreqs(tidy.haps)
tidy.allele.freqs   <- tidyingAlleleFreqs(allele.freqs.and.ld)
tidy.ld             <- tidyingLD(allele.freqs.and.ld)

# reinforcement  
tidy.reinforce  <- tidyingReinforcement(example.sim$geno.time)

# freq change components
tidy.freq.change   <- tidyingFreqChange(example.sim$dhaps)
tidy.allele.change <- alleleChange(tidy.freq.change)


# unlinked adaptation 
tidy.MeanU <- tidyingMeanU(runs = example.sim)


## Plots  
# Haplotype Frequency Plot
  hapFreqPlot(tidy.haps, title = "s=.3, r12 = r23 = 1e-04, m = .01, n.unlinked = 2")
# Allele Frequency Plot
  alleleFreqPlot(tidy.allele.freqs, title = "s=.3, r12 = r23 = 1e-04, m = .01, n.unlinked = 2")
# LD Plot
  ldPlot(tidy.ld, title = "s=.3, r12 = r23 = 1e-04, m = .01, n.unlinked = 2")
# reinforcement plot
  reinforcePlot(tidy.reinforce, title = "s=.3, r12 = r23 = 1e-04, m = .01, n.unlinked = 2")
# Haplotype Frequency Change Plot
  hapFreqChangePlot(tidy.freq.change , title = "s=.3, r12 = r23 = 1e-04, m = .01, n.unlinked = 2")
# Allele Frequency Change Plot
  alleleFreqChangePlot(tidy.allele.change, title = "s=.3, r12 = r23 = 1e-04, m = .01, n.unlinked = 2")
  
# Mean U plots
hapUnlinkedPlot(tidy.MeanU$haps)
alleleUnlinkedPlot(tidy.MeanU$alleles)


