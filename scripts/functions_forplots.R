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
    scale_x_continuous(trans = "log10")+#,limits = c(100,10000))+
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
    scale_x_continuous(trans = "log10")+#,limits = c(100,10000))+
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
    scale_x_continuous(trans = "log10")+#,limits = c(100,10000))+
    ggtitle(label = title)+
    theme_light() 
  return(my.plot)
}

# Reinforcement 
reinforcePlot <- function(tmp.reinforce, title = NULL){
  my.plot <- tmp.reinforce %>% ungroup() %>%
    mutate(pop = factor(case_when(pop == 0 ~ "maize",pop == 1~"mex"),
                        levels = c("mex","maize"))) %>%
    ggplot(aes(x = gen, 
               y = reinforce_strength,    
               color = pop))                        +
    geom_line()                                     +
    geom_hline(yintercept = 0)                      +
    scale_x_continuous(trans = "log10")+ #,
                      # limits = c(100,10000))       +
    ggtitle(label = title)                          +
    theme_light() 
  return(my.plot)
}

# haplotype change
hapFreqChangePlot <- function(tmp.hap.change, title = NULL){
  my.plot <- tmp.hap.change %>% ungroup() %>%
    mutate(local_adapt = factor(case_when(local_adapt == 0 ~ "maize",local_adapt == 1~"mex"),levels = c("mex","maize")),
           pop = factor(case_when(pop == 0 ~ "maize",pop == 1~"mex"),levels = c("mex","maize"))) %>%
    ggplot(aes(x = gen, 
               y = dhap,    
               group = haplo, 
               color = pollen_style,
               size  = local_adapt))+
    geom_line() +
    facet_grid(pop~time) +
    scale_size_manual(values = c(.7,1.6))+
    scale_x_continuous(trans = "log10")+#,limits = c(100,10000))+
    ggtitle(label = title)+
    geom_hline(yintercept = 0,lty = 2)+
    theme_light() 
  return(my.plot)
}

hapFreqChangePlotAlt <- function(tmp.hap.change, title = NULL){
  my.plot <- tmp.hap.change %>% ungroup() %>% filter(pollen_style != "mF" & !haplo %in% c("aMF","Amf") ) %>%
    mutate(local_adapt = factor(case_when(local_adapt == 0 ~ "maize",local_adapt == 1~"mex"),levels = c("mex","maize")),
           pop = factor(case_when(pop == 0 ~ "maize",pop == 1~"mex"),levels = c("mex","maize"))) %>%
    ggplot(aes(x = gen, 
               y = dhap,    
               group = paste(haplo, time), 
               color = time))+
    geom_line() +
    facet_grid(pop~haplo) +
    scale_x_continuous(trans = "log10")+#,limits = c(100,10000))+
    ggtitle(label = title)+
    geom_hline(yintercept = 0,lty = 2)+
    theme_light() 
  return(my.plot)
}

# allele change 
alleleFreqChangePlot <- function(tmp.allele.change, title = NULL){
  my.plot <- tmp.allele.change %>% ungroup() %>%
    mutate(pop = factor(case_when(pop == 0 ~ "maize",pop == 1~"mex"),levels = c("mex","maize"))) %>%
    ggplot(aes(x = gen, 
               y = delta_p,    
               color = allele))+
    geom_line() +
    facet_grid(pop~time) +
    scale_x_continuous(trans = "log10") + # ,limits = c(100,10000))+
    ggtitle(label = title)+
    geom_hline(yintercept = 0,lty = 2)+
    theme_light()
  return(my.plot)
}


alleleFreqChangePlotAlt <- function(tmp.allele.change, title = NULL){
  my.plot <- tmp.allele.change %>% ungroup() %>%
    mutate(pop = factor(case_when(pop == 0 ~ "maize",pop == 1~"mex"),levels = c("mex","maize"))) %>%
    ggplot(aes(x = gen, 
               y = delta_p,    
               color = time))+
    geom_line() +
    facet_grid(pop~allele) +
    scale_x_continuous(trans = "log10",limits = c(100,10000))+
    ggtitle(label = title)+
    geom_hline(yintercept = 0,lty = 2)+
    theme_light() 
  return(my.plot)
}


hapUnlinkedPlot <- function(tmp.U.hap, title = NULL){
  my.plot <- tmp.U.hap %>% 
    mutate(pollen_style = str_remove(haplo, "A|a"),
           local_adapt  = as.numeric(str_sub(haplo, 1,1) == "A"))%>%
    mutate(local_adapt = factor(case_when(local_adapt == 0 ~ "maize",local_adapt == 1~"mex"),levels = c("mex","maize")),
           pop = factor(case_when(pop == 0 ~ "maize",pop == 1~"mex"),levels = c("mex","maize"))) %>%
    ggplot(aes(x = gen, 
               y = U,    
               group = paste(pop, haplo), 
               color = pollen_style,
               lwd   = local_adapt))+
    geom_line() +
    facet_grid(~pop) +
    scale_x_continuous(trans = "log10")+#,limits = c(100,10000))+
    ggtitle(label = title)+
    scale_size_manual(values = c(.7,1.6))+ 
    theme_light() 
  return(my.plot)
}


alleleUnlinkedPlot <- function(tmp.U.allele, title = NULL){
  my.plot <- tmp.U.allele %>% 
    mutate(pop = factor(case_when(pop == 0 ~ "maize",pop == 1~"mex"),levels = c("mex","maize"))) %>%
    ggplot(aes(x = gen, 
               y = U,    
               group = paste(pop, locus, allele), 
               color = locus))+
    geom_line() +
    facet_grid(allele~pop, labeller = "label_both") +
    scale_x_continuous(trans = "log10")+#,limits = c(100,10000))+
    ggtitle(label = title)+
    scale_size_manual(values = c(.7,1.6))+ 
    theme_light() 
  return(my.plot)
}