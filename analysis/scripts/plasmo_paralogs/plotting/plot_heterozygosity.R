rm(list=ls())

library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(dbplyr)
library(RSQLite)
library(DBI)
library(argparser, quietly=TRUE)

args <- arg_parser("Plot homozygosity stats")
args <- add_argument(args, "sqlite_db", help="")
args <- add_argument(args, "plotting_utils", help="Rscript containing plotting utilities")
args <- add_argument(args, "base_dir_for_output_plots", help="output directory for plots")
argv <- parse_args(args)

BASE_PLOT_DIR <- argv$base_dir_for_output_plots
dir.create(BASE_PLOT_DIR, recursive=TRUE)
source(argv$plotting_utils)

con <- dbConnect(RSQLite::SQLite(), argv$sqlite_db)
#con <- dbConnect(RSQLite::SQLite(), "tmp_work/seq_stats.sqlite")

df <- tibble()
gene_names <- c("AMA1","MSP3","MSP6","DBs","MSP2","EBA175")
tool_names <- c("malariaGEN","gramtools")
seqtypes <- c("dna","protein")
seqtypes <- c("protein")

# There are two gramtools DBs seq (protein MSAs) files to be aware of:
# md5sum: 1ed85d032361f267bc447a7ebada8e2a is the one WITH single-snp-resolved
# md5sum: d04b38ef5ec0acd51eae1022491ac0e1 is the one WITHOUT single-snp-resolved
table_names <- dbListTables(con)
for (gene_name in gene_names){
  for (used_seqtype in seqtypes){
    for (tool_name in tool_names){
    table_name <- table_names[grepl(paste(tool_name,gene_name,used_seqtype,"homozygosity",sep="_"), table_names)]
    if (length(table_name) > 1){
      table_name <- grep("1ed85d032361f267bc447a7ebada8e2a",table_name,value=TRUE)
    }
    if (length(table_name) > 0){
      res <- dbSendQuery(con, paste("SELECT * FROM ",table_name))
      t <- dbFetch(res)
      t <- mutate(t,seqtype=rep(used_seqtype))
      df <- rbind(df, t)
      dbClearResult(res)
    }
    }
  }
}

DBs <- c("DBLMSP","DBLMSP2","DBLMSP_DBLMSP2")

## domain definitions from Pfam
# These were computed by taking the Pfam coordinates of the domain in 3D7 gene sequence 
# and transposing them to MSA space
DBL_domain_location_protein <- c(217,337)
DBL_domain_location_dna <- c(626,989)
SPAM_domain_location <- c(595,873)

## Operational definition of DBL domain: shared sequences either side of Pfam DBL
df_diverged <- filter(df, homozygosity == 0 & gene_ID == DBs[3] & tool_name == tool_names[2])
domains <- tibble()
for (used_seqtype in seqtypes){
  df_diverged_seqtype <- filter(df_diverged, seqtype == used_seqtype) %>%
    mutate(pos_diff = pos - lag(pos)) %>% filter(pos_diff < 3)
  used_domain_location <- DBL_domain_location_protein
  if (used_seqtype == "dna") used_domain_location <- DBL_domain_location_dna
  # The +1 are to express in 1-based
  Operational_DBL_domain_location <- c(
    max(df_diverged_seqtype$pos[df_diverged_seqtype$pos < 
                                  used_domain_location[1]]) + 1,
    min(df_diverged_seqtype$pos[df_diverged_seqtype$pos > 
                                  used_domain_location[2]]) + 1
  )
  t <- tibble(
    intercept=c(used_domain_location,Operational_DBL_domain_location),
    domain=c(rep("DBL(Pfam)",2),rep("DBL(Operational)",2)),
    seqtype=rep(used_seqtype)
  )
  domains <- rbind(domains, t)
}


positions <- ggplot(filter(df,seqtype=="protein"), aes(x=pos,y=1-homozygosity)) + 
  geom_point() + 
  facet_wrap(vars(gene_ID),ncol=1) +
  ylab("Heterozygosity (prob that two random amino acids differ)") +
  xlab("Position in protein") + 
  theme(text=element_text(size=15))
plot_save("positional_aa_heterozygosity.pdf",width=16,height=13)

df_means <- df %>% group_by(gene_ID,seqtype) %>% summarise(mean_heterozygosity=mean(1-homozygosity)) %>%
  arrange(mean_heterozygosity)

plot_positions <- function(df, domains, use_domains = FALSE){
  p <- ggplot(df, aes(x=pos,y=1-homozygosity)) + 
    geom_point()
  if (use_domains){
    p <- p+ geom_vline(data=filter(domains,seqtype=="protein"),aes(xintercept = intercept - 1,colour=domain), linetype=2)
  }
  p <- p + facet_grid(cols=vars(tool_name),rows=vars(gene_ID)) +
  ylab("Heterozygosity (prob that two random amino acids differ)") +
  xlab("Position in protein") + 
  theme(text=element_text(size=17))
  p
  return(p)
}

df_DBs <-filter(df,seqtype=="protein" & gene_ID %in% DBs) %>% mutate(
  gene_ID = factor(gene_ID, levels=DBs),
  tool_name = factor(tool_name, levels=tool_names))
plot_positions(df_DBs,domains,TRUE)
plot_save("DBs_aa_heterozygosity.pdf",width=16,height=13)

## Plot other genes 
non_DB_genes <- grep(DBs,unique(df$gene_ID), value=TRUE,invert=TRUE)
for (GID in non_DB_genes){
  df_gene <-filter(df,seqtype=="protein" & gene_ID == GID) %>% mutate(
    tool_name = factor(tool_name, levels=tool_names))
  plot_positions(df_gene, domains)
  plot_save(paste(GID,"_aa_heterozygosity.pdf",sep=""),width=15,height=10)
}
