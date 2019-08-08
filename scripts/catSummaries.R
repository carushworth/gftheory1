library("tidyverse")
files.to.summarize <- str_subset(list.files(),"outs") %>% str_subset(".txt")
do.call(rbind, lapply(files.to.summarize, function(X){t(read.table(X))}) ) %>% 
  as_tibble()  %>% 
  mutate_at(vars(-this.order), as.numeric) %>%
  write_csv(path = "summarystats.csv")
