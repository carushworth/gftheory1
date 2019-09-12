# fully worked example
setwd("/Users/catherinerushworth1/projects/gftheory1/scripts")
source("gf_sims.R")
source("functions_forsummaries.R")
source("functions_forplots.R")
# example.sim <- runGFsim(n.gen = 10000, discrim = 1, s   =     .75,
#          r12 =     .0001, r23 =     .0001, 
#          init.freqs =  c(fA_0 = 0, fM_0 = 0, fF_0 = 0, fA_1 = 1, fM_1 = 1, fF_1 = 0.01),
#          prop.replaced0 = .1,   prop.replaced1 = .1,  n.unlinked     = 0,delta_hap_components = TRUE ) 



# load R object and every time it would be called output, or you can name it example.sim
# run once or twice with 0 unlinked vs >0 unlinked

this.file <- "/Users/catherinerushworth1/projects/gftheory1/results/outs 0.01 m 0.01 r12 0.5 r23 0.001 M 0.05 n.unl 4.Robj"
this.title <- str_remove(string=this.file, pattern= "prelimSimResults/out") %>% 
  str_remove(".Robj") %>% 
  str_split(" ")      %>%
  unlist()            %>%
  matrix(ncol=2,byrow=TRUE)%>% 
  data.frame()%>% 
  mutate(vals = paste(X1,X2, sep = " = ")) %>% 
  select(vals) %>% pull %>% paste(collapse = ", ")

load(this.file)
# this is sim #1 from from this cluster run (number 1 in column "sim" in summarystats.csv)
#example.sim <- output

#############
# Frequencies and reinfrocement
# haplotypes, only looking at AMF not at unlinked stuff
tidy.haps           <- tidyingHaps(output$geno.time) # local adapt strictly refers to the a locus, not unlinked loci

# alleles and ld
allele.freqs.and.ld <- findFreqs(tidy.haps)
tidy.allele.freqs   <- tidyingAlleleFreqs(allele.freqs.and.ld)
tidy.ld             <- tidyingLD(allele.freqs.and.ld) # each pairwise D

# reinforcement  
tidy.reinforce  <- tidyingReinforcement(output$geno.time)

# freq change components
tidy.freq.change   <- tidyingFreqChange(output$dhaps)
# dhap is hap change freq in that step for each gen
# 
tidy.allele.change <- alleleChange(tidy.freq.change)


# unlinked adaptation 
tidy.MeanU <- tidyingMeanU(runs = output)


## Plots  
# Haplotype Frequency Plot
hapFreqPlot(tidy.haps, title = this.title)
ggsave(file="hapfreq_test1.pdf")
# Allele Frequency Plot
alleleFreqPlot(tidy.allele.freqs, title = this.title)
ggsave(file="allelefreq_test1.pdf")
# LD Plot
ldPlot(tidy.ld, title = this.title)
ggsave(file="LD_test1.pdf")
# reinforcement plot
reinforcePlot(tidy.reinforce, title = this.title)
ggsave(file="reinforce_test1.pdf")
# Haplotype Frequency Change Plot
hapFreqChangePlot(tidy.freq.change , title = this.title)
ggsave(file="hapFreqChange_test1.pdf")
# Allele Frequency Change Plot
alleleFreqChangePlot(tidy.allele.change, title = this.title)
ggsave(file="alleleFreqChange_test1.pdf")

# Mean U plots
hapUnlinkedPlot(tidy.MeanU$haps)
ggsave(file="hapUnlinked_test1.pdf")
alleleUnlinkedPlot(tidy.MeanU$alleles)
ggsave(file="alleleUnlinked_test1.pdf")