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




tidying   <- function(output, blank.template = NULL){
  if(is.null(blank.template )){
    recover()
  blank.template <- gather(data = output[,-grep("reinf",names(output))], key = genopop, value = freq, - gen)
  tmp <- do.call(rbind,strsplit(blank.template$genopop,""))[, -1]
  tmp[,1] <- ifelse( tmp[,1] == "0", "a","A")
  tmp[,2] <- ifelse( tmp[,2] == "0", "m","M")
  tmp[,3] <- ifelse( tmp[,3] == "0", "f","F")
  tmp[,4] <- ifelse( tmp[,4] == "0", "a","A")
  tmp[,5] <- ifelse( tmp[,5] == "0", "m","M")
  tmp[,6] <- ifelse( tmp[,6] == "0", "f","F")
  blank.template <- blank.template %>% 
    mutate(pop = tmp[,9],
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
  tmp.output <- tmp.output %>%
    mutate(pollen_style = str_remove_all(string = haplo ,pattern = "A|a"),
           pollen_adapt = str_remove_all(string = haplo ,pattern = "F|f"),
           style_adapt = str_remove_all(string = haplo ,pattern = "M|m"),
           local_adapt  = as.numeric(grepl(pattern = "A", x = haplo)),
           style_pref   = as.numeric(grepl(pattern = "F", x = haplo)),
           pollen_sig   = as.numeric(grepl(pattern = "M", x = haplo))) 
  return(tmp.output)
}


tidyReinforcement <- function(output){
  reinf <- output[,c(1,grep("reinf",names(output)))]
  colnames(reinf) <- c("gen","p0","p1")
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
