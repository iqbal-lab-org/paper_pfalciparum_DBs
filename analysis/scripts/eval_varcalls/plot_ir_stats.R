library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(forcats)
library(argparser, quietly=TRUE)

df = read_tsv("/home/brice/Desktop/main_PhD/analyses/plasmo_surfants/tmp_work/ir_stats.tsv")

means <- df %>% group_by(gene, tool) %>% summarise(
  mean_bases_mapped=mean(bases_mapped_cigar),
  mean_error_rate=mean(base_error_rate),
  mean_ir_scaled_eddist=mean(induced_ref_scaled_eddist),
  )
ggplot(means, aes(y=mean_bases_mapped,x=gene,fill=tool)) + geom_col(position="dodge")
ggplot(means, aes(y=mean_error_rate,x=gene,fill=tool)) + geom_col(position="dodge")
ggplot(means, aes(y=mean_ir_scaled_eddist,x=gene,fill=tool)) + geom_col(position="dodge")

df_myo_7 = df %>% filter(gene == "LSA1" | gene == "DBLMSP" | gene == "DBLMSP2"
                           | gene == "LSA3" | gene == "CSP" | gene == "MSP1" | gene == "MSP2")
means_myo_7 <- df_myo_7 %>% group_by(gene, tool) %>% summarise(
  mean_bases_mapped=mean(bases_mapped_cigar),
  mean_error_rate=mean(base_error_rate),
  mean_ir_scaled_eddist=mean(induced_ref_scaled_eddist),
  )
ggplot(means_myo_7, aes(y=mean_bases_mapped,x=gene,fill=tool)) + geom_col(position="dodge")
ggplot(means_myo_7, aes(y=mean_error_rate,x=gene,fill=tool)) + geom_col(position="dodge")
ggplot(means_myo_7, aes(y=mean_ir_scaled_eddist,x=gene,fill=tool)) + geom_col(position="dodge")

ggplot(df, aes(y=induced_ref_scaled_eddist,x=base_error_rate,colour=tool)) + geom_point()
