# helper functions to summarize haplo freqs, allele freqs, and ld 
# this is to get allele freqs out
blankTemplate <- function(run){
  tmp.output <- gather(data = run[,-grep("reinf",names(run))], key = genopop, value = freq, - gen)
  tmp <- do.call(rbind,strsplit(tmp.output$genopop,""))[, -1]
  tmp[,1] <- ifelse( tmp[,1] == "0", "a","A")
  tmp[,2] <- ifelse( tmp[,2] == "0", "m","M")
  tmp[,3] <- ifelse( tmp[,3] == "0", "f","F")
  tmp[,4] <- ifelse( tmp[,4] == "0", "a","A")
  tmp[,5] <- ifelse( tmp[,5] == "0", "m","M")
  tmp[,6] <- ifelse( tmp[,6] == "0", "f","F")
  tmp.output <- tmp.output %>% 
    mutate(pop = tmp[,7],
           mat = paste( tmp[,1], tmp[,2], tmp[,3], sep =  ""), 
           pat = paste( tmp[,4], tmp[,5], tmp[,6], sep =  "")) 
  return(tmp.output)
}




tidyingHaps   <- function(output, blank.template = NULL){
  if(is.null(blank.template )){
  blank.template <- gather(data = output[,-grep("reinf",names(output))], key = genopop, value = freq, - gen)
  tmp <- do.call(rbind,strsplit(blank.template$genopop,""))[, -1]
  tmp[,1] <- ifelse( tmp[,1] == "0", "a","A")
  tmp[,2] <- ifelse( tmp[,2] == "0", "m","M")
  tmp[,3] <- ifelse( tmp[,3] == "0", "f","F")
  tmp[,4] <- ifelse( tmp[,4] == "0", "a","A")
  tmp[,5] <- ifelse( tmp[,5] == "0", "m","M")
  tmp[,6] <- ifelse( tmp[,6] == "0", "f","F")
  blank.template <- blank.template %>% 
    mutate(pop = tmp[,7],
           mat = paste( tmp[,1], tmp[,2], tmp[,3], sep =  ""), 
           pat = paste( tmp[,4], tmp[,5], tmp[,6], sep =  "")) 
  }
  my.freqs <- gather(data = output[,-grep("reinf",names(output))], key = genopop, value = freq, - gen) %>%select(freq)
  tmp.output <- as_tibble(blank.template) %>%
    mutate(freq = pull(my.freqs))         %>%
    gather(key = parent, value =  haplo, -gen, -genopop, - freq, - pop) %>%
    mutate(freq = freq/2) %>% 
    group_by(gen , pop, haplo) %>% 
    summarise(freq = sum(freq))
  #this bit is slow and non-param result specific so could be resturctured for speed
  tmp.output <- tmp.output %>%
    mutate(pollen_style = str_remove_all(string = haplo ,pattern = "A|a"),
           pollen_adapt = str_remove_all(string = haplo ,pattern = "F|f"),
           style_adapt = str_remove_all(string = haplo ,pattern = "M|m"),
           local_adapt  = as.numeric(grepl(pattern = "A", x = haplo)),
           style_pref   = as.numeric(grepl(pattern = "F", x = haplo)),
           pollen_sig   = as.numeric(grepl(pattern = "M", x = haplo))) 
  return(tmp.output %>% ungroup())
}


tidyingReinforcement <- function(output){
  reinf <- output[,c(1,grep("reinf",names(output)))]
  colnames(reinf) <- c("gen","0","1")
  gather(data = reinf, key = pop, value = reinforce_strength, -gen) 
}

findFreqs <- function(tidied_data){
  tidied_data %>% 
    group_by(gen, pop) %>%
    summarise(freq_A = sum(local_adapt*freq),
              freq_M = sum(pollen_sig*freq),
              freq_F = sum(style_pref*freq),
              ld_MF  = sum(freq * as.numeric(pollen_style == "MF")) - freq_M * freq_F,
              ld_AM  = sum(freq * as.numeric(pollen_adapt == "AM")) - freq_A * freq_M,
              ld_AF  = sum(freq * as.numeric(style_adapt  == "AF")) - freq_A * freq_F) %>% 
    ungroup()
}

tidyingAlleleFreqs <- function(output){
  # require allele freqs and ld format
  output %>% 
    select(-starts_with("ld"))               %>%
    rename(A = freq_A, M = freq_M, F = freq_F) %>%
    gather(key = allele, value = freq, - gen, -pop)
}


tidyingLD <- function(output){
  output %>% 
    select(-starts_with("freq"))               %>%
    rename(AM = ld_AM, MF = ld_MF, AF = ld_AF) %>%
    gather(key = pair, value = ld, - gen, -pop)
}


tidyingFreqChange <- function(output){ 
  tmp.output <- output  %>% 
    gather(key = thing, value = dhap, -gen) %>%
    mutate( pop = ifelse(grepl(pattern= "p0",thing), 0, 1),
            thing = str_remove(thing, pattern = "p0_|p1_"),
            time  = str_remove(thing, pattern = "..._"),
            hap   = str_remove(thing, pattern = "_.*"))     %>%
    select(-thing)                                          %>%
    mutate(haplo = paste(ifelse(str_sub(hap,1,1)==0,"a","A"),
                         ifelse(str_sub(hap,2,2)==0,"m","M"),
                         ifelse(str_sub(hap,3,3)==0,"f","F"), sep =  ""),
           A = str_sub(hap,1,1), 
           M = str_sub(hap,2,2),
           F = str_sub(hap,3,3))
  tmp.output <- tmp.output %>%
    mutate(time = case_when(time == "migrationPollen" ~ "migrationPollen",
                            time == "matingPollen"    ~ "matingPollen",
                            grepl("sel", time)        ~ "sel")) %>%
    group_by(gen, pop, haplo, A, M, F, time) %>%
    summarise(dhap = mean(dhap))   %>% ungroup() %>%
    mutate(time = factor(time,levels = c("migrationPollen","matingPollen","sel")))
  tmp.output <- tmp.output %>%
    mutate(pollen_style = str_remove_all(string = haplo ,pattern = "A|a"),
           pollen_adapt = str_remove_all(string = haplo ,pattern = "F|f"),
           style_adapt = str_remove_all(string = haplo ,pattern = "M|m"),
           local_adapt  = as.numeric(grepl(pattern = "A", x = haplo)),
           style_pref   = as.numeric(grepl(pattern = "F", x = haplo)),
           pollen_sig   = as.numeric(grepl(pattern = "M", x = haplo))) 
  return(tmp.output %>% ungroup())
}


alleleChange <- function(output){
  output %>% 
    group_by(pop, gen, time, A) %>%  mutate(dA = sum(as.numeric(A) * dhap)) %>% ungroup()  %>%
    group_by(pop, gen, time, M) %>%  mutate(dM = sum(as.numeric(M) * dhap)) %>% ungroup()  %>%
    group_by(pop, gen, time, F) %>%  mutate(dF = sum(as.numeric(F) * dhap)) %>% ungroup()  %>%
    select(gen, pop, time, A = dA, M = dM, F = dF)                                    %>%
    gather(key = allele, value = delta_p, -gen, -pop, -time)                          %>% 
    group_by(gen, pop, time, allele) %>% summarise(delta_p = mean(delta_p)) %>%ungroup()
}





#findFreqs <- function(tidied_data){
#  tidied_data %>% 
#    ungroup() %>%
#    mutate(am = freq * as.numeric(str_remove(haplo, "F|f") == "am"),
#           mf = freq * as.numeric(str_remove(haplo, "A|a") == "mf"),
#           af = freq * as.numeric(str_remove(haplo, "M|m") == "af"),
#           a  = freq * as.numeric(grepl(pattern = "a", haplo)), 
#           m  = freq * as.numeric(grepl(pattern = "m", haplo)), 
#           f  = freq * as.numeric(grepl(pattern = "f", haplo) )) %>%
#    group_by(pop, gen) %>%
#    summarise(am = sum(am),
#              mf = sum(mf),
#              af = sum(af),
#              a  = sum(a),
#              m  = sum(m),
#              f  = sum(f)) %>%
#    group_by(pop, gen)
#}
#findLD    <- function(tidied_data){
#  findFreqs(tidied_data) %>%
#    summarise(D_am = am - (a*m),
#              D_mf = mf - (m*f),
#              D_af = af - (a*f),
#              r_am = D_am / sqrt( a*(1-a) * m *(1-m) ), 
#              r_mf = D_mf / sqrt( m*(1-m) * f *(1-f) ),
#              r_af = D_af / sqrt( a*(1-a) * f *(1-f) )) %>%
#    ungroup() 
#}
