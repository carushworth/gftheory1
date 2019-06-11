# run sims / tidy output 
# run once
source("gf_sims.R")
source("functions_forsummaries.R")
n.gen          <- 10000
blank          <- runGFsim(n.gen = n.gen, get.blank = TRUE)
blank.template <- blankTemplate(blank)

param.combos <- expand.grid(
  s = c(0.01, 0.1, .25, 0.5, 0.75) ,
  m = c(0.01,0.05,0.1) , 
  r12 = c(0.0001),
  r23 = c(0.0001),
  init_freq = c(0.01),
  n_unlinked = c(2,3,4,5)
)



# run params
apply(param.combos, 1, function(PARAMS){
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
  write_csv(output[["geno.time"]], path =  paste("prelimSimResults/outputs/",paste(paste(names(PARAMS), PARAMS, sep = "="),collapse = "_"),".csv",sep=""))
  write_csv(output[["meanUs"]],    path =  paste("prelimSimResults/outputs/",paste("MeanUs",paste(names(PARAMS), PARAMS, sep = "="),collapse = "_"),".csv",sep=""))
})

#make plots
apply(param.combos, 1, function(PARAMS){
  output               <- read_csv(paste("prelimSimResults/",paste(paste(names(PARAMS), PARAMS, sep = "="),collapse = "_"),".csv",sep=""))
  tidy.output          <- tidying(output,blank.template)
  tidy.reinforcement   <- tidyReinforcement(output)
  
  fig1 <- ggplot(ungroup(tidy.output) %>% 
                   filter(pollen_style !="mF")%>%
                   mutate(pop = factor(case_when(pop == 0 ~ "maize",pop == 1~"mex"),levels = c("mex","maize")))%>%
                   mutate(local_adapt = factor(case_when(local_adapt == 0 ~ "maize",local_adapt == 1~"mex"),levels = c("mex","maize"))) , 
                 aes(x = gen, y = freq, 
                          color = pollen_style, 
                          lty = local_adapt)) +
    geom_line() +
    facet_wrap(~pop)+
    scale_x_continuous(trans = "log10")+
    theme_light()  +
    ggtitle(paste(paste(names(PARAMS), PARAMS, sep = "="),collapse = ", "))
  
  fig2 <- ggplot(tidy.reinforcement %>% 
           mutate(pop = factor(case_when(pop == "p0" ~ "maize",pop == "p1"~"mex"),levels = c("mex","maize"))),
         aes(x = gen, y = reinforce_strength, color = pop)) +
    geom_line() +
    scale_x_continuous(trans = "log10") +
    theme_light() 
  freqs.time           <- findFreqs(tidy.output) 
  
  fig3 <- freqs.time %>% 
    select(-starts_with("ld")) %>% 
    rename(A = freq_A, M = freq_M, F = freq_F )    %>%
    gather(key = allele, value = freq, -gen,-pop)  %>%
    mutate(pop = factor(case_when(pop == 0 ~ "maize",pop == 1 ~"mex"),levels = c("mex","maize"))) %>%
    ggplot(aes(x = gen, y = freq, color = allele)) + 
    geom_line() +
    facet_wrap(~pop)+
    scale_x_continuous(trans = "log10")+
    theme_light() 
  
  fig4 <- freqs.time %>% 
    select(-starts_with("freq")) %>% 
    rename(MF = ld_MF, AM = ld_AM, AF = ld_AF )    %>%
    gather(key = pair, value = ld, -gen,-pop)  %>%
    mutate(pop = factor(case_when(pop == 0 ~ "maize",pop == 1 ~"mex"),levels = c("mex","maize"))) %>%
    ggplot(aes(x = gen, y = ld, color = pair)) + 
    geom_line() +
    facet_wrap(~pop)+
    scale_x_continuous(trans = "log10")+
    scale_y_continuous(trans = "log1p")+
    theme_light() 
  my.plot <- grid.arrange(fig1,fig2,fig3,fig4, ncol =1)
  ggsave(filename = paste("prelimSimResults/figures/",paste(paste(names(PARAMS), PARAMS, sep = "="),collapse = "_"),".pdf",sep=""), 
         plot = my.plot,device = "pdf")
})





### more loc
n.gen          <- 10000
blank          <- runTmpSim(n.gen = n.gen, get.blank = TRUE)
blank.template <- blankTemplate(blank)

param.combos <- expand.grid(
  s = c(0.01, 0.1, .25, 0.5, 0.75) ,
  m = c(0.01,0.05,0.1) , 
  r12 = c(0.0001),
  r23 = c(0.0001),
  init_freq = c(0.01),
  n.loc = c(2:5)
)
rm(output)
# run params
apply(param.combos, 1, function(PARAMS){
  output  <- runTmpSim(n.gen = n.gen, 
                       discrim = 1,
                       s   =     PARAMS[["s"]],
                       r12 =     PARAMS[["r12"]],
                       r23 =     PARAMS[["r23"]],
                       init.freqs =  c(fA_0 = 0, fM_0 = 0, fF_0 = 0,
                                       fA_1 = 1, fM_1 = 1, fF_1 = PARAMS[["init_freq"]]),
                       prop.replaced0 = PARAMS[["m"]],
                       prop.replaced1 = PARAMS[["m"]],
                       n.unlinked     = PARAMS[["n.loc"]]
  )
  write_csv(output, path =  paste("prelimSimResults/",paste(paste(names(PARAMS), PARAMS, sep = "="),collapse = "_"),".csv",sep=""))
})











# run in loop of output
output               <- runTmpSim(n.gen = n.gen, 
                                        discrim = 1,
                                        s   =     .9, 
                                        r12 =     0.0001, 
                                        r23 =     0.0,
                                        init.freqs =  c(fA_0 = 0, fM_0 = 0, fF_0 = 0,
                                                        fA_1 = 1, fM_1 = 1, fF_1 = .01),
                                        prop.replaced0 = .3,
                                        prop.replaced1 = .3,
                                        this.order = "AFM")
tidy.output          <- tidying(output,blank.template)
ggplot(tidy.output%>% ungroup(), aes(x = gen, y = freq, color = haplo)) + geom_line() + facet_wrap(~ pop)



tidy.reinforcement   <- tidyReinforcement(output)
ggplot(tidy.reinforcement, aes(x = gen, y = reinforce_strength,color = pop)) + 
  geom_line()

freqs.time %>%
  #select(-am, -af, -mf) %>% 
  select(pop  ,   gen ,a,m,f) %>% 
  gather(key = pair, value = freq, -gen, - pop)%>% 
  ungroup()%>%
  mutate(pop = factor(case_when(pop == 0 ~ "maize",pop == 1~"mex"), 
                      levels = c("mex","maize"))) %>%
  ggplot(aes(x = gen, y = freq, color = pair)) + geom_line()+facet_wrap(~pop)

freqs.time %>%
  select(gen,pop, starts_with("freq")) %>% 
  gather(key = pair, value = freq, -gen, - pop)%>% 
  ungroup()%>%
  mutate(pop=case_when(pop == 0 ~ "maize",pop == 1~"mex")) %>%
  ggplot(aes(x = gen, y = freq, color = pair)) + 
  geom_line()+ 
  scale_x_log10() + 
  facet_wrap(~pop)
  
tidy.reinforcement %>% 
  mutate(pop=case_when(pop == "p0" ~ "maize",pop == "p1"~"mex")) %>%
  ggplot(aes(x = gen, y = reinforce_strength)) +
  geom_hline(yintercept = 0, color = "pink")+
  geom_line() +
  scale_x_log10() + 
  facet_wrap(~pop)+ 
  theme_light()
  





#AMF   
#blank         <- runFourLocusSim(n.gen = 10000, get.blank = TRUE)
#output        <- runFourLocusSim(n.gen = 10000, discrim = 1, s=.5, r12=.01, r23 = 0)


#output %>% tidying()
output_freq   <- findFreqs(output)
output_LD     <- findLD(output)

#geno_freq_plot <- output_freq %>% 
#  select(- am, - mf, - af)    %>%  
#  gather(key = locus, value = freq, - pop, -gen) %>%
#  ggplot(aes(x=gen,y=freq,color=locus )) + 
#  geom_line() + 
#  facet_wrap(~pop, labeller = "label_both")+
#  scale_x_continuous(trans = "log10") +
#  ggtitle("Allele frequencies over time")+
#  theme_tufte()



#ggplot(output, aes(x=gen,y=freq,color=pollen_style, group = haplo,  linetype = local_adapt )) + 
#  geom_line() + 
#  facet_wrap(~pop, labeller = "label_both")+
#  scale_x_continuous(trans = "log10") +
#  ggtitle("Haplotype frequencies over time")+
#  theme_tufte()
#
#ld_plot  <- output_LD %>% 
#  select(pop,gen,D_am, D_mf, D_af) %>%
#  gather(key = pair, value = D, - pop, - gen) %>%
#  ggplot(aes(x = gen, y = D, color = pair)) + 
#  geom_line()+
#  facet_wrap(~pop, labeller = "label_both")+ 
#  scale_x_continuous(trans = "log10")+
#  ggtitle("Linkage disequilibrium over time")+
#  theme_tufte()

#haplo_freq_plot <- ggplot(output, aes(x=gen,y=freq,color=pollen_style, group = haplo,  linetype = local_adapt )) + 
#  geom_line() + 
#  facet_wrap(~pop, labeller = "label_both")+
#  scale_x_continuous(trans = "log10") +
#  ggtitle("Haplotype frequencies over time")+
#  theme_tufte()


#grid.arrange(geno_freq_plot, ld_plot,haplo_freq_plot, ncol = 1)



# so , why doesn't F spread with bideirectional migration... it's because F is indirectly selected against in the other pop because it mates with assholes
# ideas. bad fitness is dom. 
# or unlinked locus under sel