# fully worked example
source("gf_sims.R")
source("functions_forsummaries.R")
source("functions_forplots.R")
example.sim <- runGFsim(n.gen = 10000, discrim = 1, s   =     .3,
         r12 =     1e-04, r23 =     1e-04,
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
tidy.reinforce      <- tidyReinforcement(example.sim$geno.time)

# Haplotype Frequency Plot
hapFreqPlot(tidy.haps, title = "s=.3, r12 = r23 = 1e-04, m = .01, n.unlinked = 2")
# Allele Frequency Plot
alleleFreqPlot(tidy.allele.freqs, title = "s=.3, r12 = r23 = 1e-04, m = .01, n.unlinked = 2")
# LD Plot
ldPlot(tidy.ld, title = "s=.3, r12 = r23 = 1e-04, m = .01, n.unlinked = 2")
# reinforcement plot
reinforcePlot(tidy.reinforce, title = "s=.3, r12 = r23 = 1e-04, m = .01, n.unlinked = 2")

tidy.MeanU <- tidyingMeanU(example.sim$meanUs, thing = "meanUs")
