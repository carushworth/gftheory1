# making plots! 

# Haplotype frequencies
hapFreqPlot <- function(tmp.haps, title = NULL){
  my.plot <- tmp.haps %>% ungroup() %>%
    mutate(local_adapt = factor(case_when(local_adapt == 0 ~ "maize",local_adapt == 1~"mex"),levels = c("mex","maize")),
         pop = factor(case_when(pop == 0 ~ "maize",pop == 1~"mex"),levels = c("mex","maize"))) %>%
    ggplot(aes(x = gen, 
             y = freq,    
             group = haplo, 
             color = pollen_style,
             size  = local_adapt))+
    geom_line() +
    facet_wrap(~pop,ncol = 1,labeller = "label_both") +
    scale_size_manual(values = c(.7,1.6))+
    scale_x_continuous(trans = "log10",limits = c(100,10000))+
    ggtitle(label = title)+
    theme_light() 
    return(my.plot)
}



# Allele frequencies
alleleFreqPlot <- function(tmp.alleles, title = NULL){
  my.plot <- tmp.alleles %>% ungroup() %>%
    mutate(pop = factor(case_when(pop == 0 ~ "maize",pop == 1~"mex"),levels = c("mex","maize"))) %>%
    ggplot(aes(x = gen, 
               y = freq,    
               group = allele, 
               color = allele))+
    geom_line() +
    facet_wrap(~pop,ncol = 1,labeller = "label_both") +
    scale_x_continuous(trans = "log10",limits = c(100,10000))+
    ggtitle(label = title)+
    theme_light() 
  return(my.plot)
}


# LD
ldPlot <- function(tmp.ld, title = NULL){
  my.plot <- tmp.ld %>% ungroup() %>%
    mutate(pop = factor(case_when(pop == 0 ~ "maize",pop == 1~"mex"),levels = c("mex","maize"))) %>%
    ggplot(aes(x = gen, 
               y = ld,    
               group = pair, 
               color = pair))+
    geom_line() +
    facet_wrap(~pop,ncol = 1,labeller = "label_both") +
    scale_x_continuous(trans = "log10",limits = c(100,10000))+
    ggtitle(label = title)+
    theme_light() 
  return(my.plot)
}

# Reinforcement frequencies
reinforcePlot <- function(tmp.reinforce, title = NULL){
  my.plot <- tmp.reinforce %>% ungroup() %>%
    mutate(pop = factor(case_when(pop == 0 ~ "maize",pop == 1~"mex"),
                        levels = c("mex","maize"))) %>%
    ggplot(aes(x = gen, 
               y = reinforce_strength,    
               color = pop))                        +
    geom_line()                                     +
    geom_hline(yintercept = 0)                      +
    scale_x_continuous(trans = "log10",
                       limits = c(100,10000))       +
    ggtitle(label = title)                          +
    theme_light() 
  return(my.plot)
}