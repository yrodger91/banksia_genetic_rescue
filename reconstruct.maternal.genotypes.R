#### Code created with the help of Dr Paul Harrison from the Monash Bioinformatics Platform ####

library(tidyverse)

#first I created a csv file containing the genotypes for the selfed offspring from the mother of interest
#individual ID in the first column followed by two columns for each locus (allele 1 and allele 2), each row is a different individual

df <- read_csv("M88.csv", col_names=FALSE)

n <- 2035 #total number of loci
colnames(df) <- c("id", paste0( rep(1:2035, each=2), "_", rep(1:2, times=n) )) #creating column names in the form locus1_allele1 and locus1_allele2, etc

#the following code merges the genotypes across the individuals: if all homozygous, then resulting maternal genotype must be homozygous
#if one of the individuals is heterozygous at a locus, the mother must be heterozygous
result <- df %>% 
  pivot_longer(!id, names_to="name", values_to="code")%>%
  separate(name, c("locus","allele"), convert=TRUE) %>%
  filter(code != 0) %>%
  group_by(locus) %>%
  summarize(code = paste(sort(unique(code)),collapse=",") ) %>%
  arrange(locus)

result %>% pivot_wider(names_from="locus", values_from="code")

#the following code pivots the above results into the original format and separates out the 2 alleles per locus again
result2 <- result %>%
    group_by(locus) %>%
    summarize(
        allele=c(1,2),
        code = 
            if (code == "1") c(1,1) 
            else if (code == "2") c(2,2) 
            else c(1,2)) %>%
    ungroup() %>%
    mutate(name=paste0(locus,"_",allele)) %>%
    select(-locus, -allele) %>%
    pivot_wider(names_from="name", values_from="code")

write.table(result2, "M88 mother genotype 1row.txt", quote = F, row.names = F, col.names = F) # outputs 1 row structure format

#use pgdspider to convert to genepop for use in hybridlab simulations
