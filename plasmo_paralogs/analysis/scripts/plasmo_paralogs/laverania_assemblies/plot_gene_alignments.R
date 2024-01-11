library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(ggsci)
library(viridis)
# library(cowplot)
# library(ggrepel)
# library(ggExtra)
# library(argparser, quietly=TRUE)

setwd("/home/adminbrice/Desktop/research/PhD/main_PhD/analyses/plasmo_paralogs")
df <- read_tsv("analysis/outputs/laverania_assemblies/alignments/alignment_stats.tsv")
# unique(df$assembly_species)
# colnames(df)

# Ensure one assembly per species, and remove dimorphic genes
df_filtered <- df %>% filter(! assembly_ID %in% c("PGABG01", "PREICH001")) %>%
  filter(! gene_ID %in% c("EBA175", "MSP1"))

df_filtered$assembly_species <- sub("Plasmodium_","P. ",df_filtered$assembly_species)
# Order species by divergence level to P. falciparum, as estimated in https://doi.org/10.1038/s41564-018-0162-2
df_filtered$assembly_species <- factor(df_filtered$assembly_species, levels=c("P. praefalciparum", "P. reichenowi", "P. billcollinsi", "P. blacklocki", "P. gaboni", "P. adleri"))
# Order genes by 5'-3' occurrence on P. falciparum 3D7 genome
df_filtered$gene_ID <- factor(df_filtered$gene_ID, levels=c("GLURP", "MSP3", "MSP6","DBLMSP", "MSP11", "DBLMSP2", "LSA1","AMA1"))
# df$gene_ID <- factor(df$gene_ID, levels=c("EBA175","MSP1", "GLURP", "MSP3", "MSP6","DBLMSP", "MSP11", "DBLMSP2", "LSA1","AMA1"))

ggplot(df_filtered, aes(x=assembly_species, y=hit_length_as_frac_of_3d7_length, colour=hit_percent_identity)) + 
  geom_point(size=2) + 
  facet_wrap(vars(gene_ID),ncol=3) + 
  scale_colour_viridis_c() + 
  theme(axis.text.x = element_text(angle=20,size=11), text=element_text(size=14)) + 
  ylab("Hit length\n(as fraction of 3D7 sequence length)") +
  xlab("Plasmodium species\n(Ordered by increasing level of divergence to P. falciparum)") + 
  labs(title="Gene presence/absence in Laverania\n(Genes ordered by occurrence on P. falciparum 3D7 genome)", colour="Hit edit distance to 3D7 sequence\n(as fraction of hit length)")

ggsave("tmp_work/gene_presence_absence_laverania.pdf", width=16,height=9)
