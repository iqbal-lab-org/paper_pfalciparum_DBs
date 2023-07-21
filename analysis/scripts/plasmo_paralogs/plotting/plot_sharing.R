rm(list=ls())

library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(dbplyr)
library(RSQLite)
library(DBI)
library(ggsci)
library(ggrepel)
library(scales)
library(cowplot)
library(argparser, quietly=TRUE)

args <- arg_parser("Plot sharing-related stats computed by seq_stats workflow")
args <- add_argument(args, "sqlite_db", help="")
args <- add_argument(args, "plotting_utils", help="Rscript containing plotting utilities")
args <- add_argument(args, "base_dir_for_output_plots", help="output directory for plots")
argv <- parse_args(args)

BASE_PLOT_DIR <- argv$base_dir_for_output_plots
dir.create(BASE_PLOT_DIR, recursive=TRUE)
source(argv$plotting_utils)

## Set up SQLite connection + load common tables ##
#con <- dbConnect(RSQLite::SQLite(), "tmp_work/seq_stats.sqlite")
con <- dbConnect(RSQLite::SQLite(), argv$sqlite_db)
table_names <- dbListTables(con)

# Domain coordinates
query <- dbSendQuery(con, "SELECT * FROM domain_coordinates")
df_domains <- dbFetch(query)
dbClearResult(query)
t <- df_domains %>% filter(Feature == "DBL_operational") %>% select(protein_MSA_coords)
target_pos_range <- as.integer(unlist(strsplit(t[1,1],":")))
target_pos_range[1] = target_pos_range[1] - 10
target_pos_range[2] = target_pos_range[2] + 10

# Pf6 metadata
res <- dbSendQuery(con, "SELECT * FROM pf6_metadata")
df_pf6 <- dbFetch(res)
dbClearResult(res)

RETAINED_GEO_SHARING <- "all_countries_min1_samples"


#################################
####Definition of sharing #######
#################################

query_table <- grep("DBs_protein_sharing_def_by_geo",table_names,value = TRUE)
query <- dbSendQuery(con, paste("SELECT * FROM",query_table))
df_geo_sharing <- dbFetch(query) 
dbClearResult(query)
df_geo_sharing <- df_geo_sharing %>% mutate(sharing_ID=sub("@","_",sharing_ID))

query_table <- grep("DBs_protein_num_used_samples",table_names,value = TRUE)
query <- dbSendQuery(con, paste("SELECT * FROM",query_table))
df_used_sample_nums <- dbFetch(query)
dbClearResult(query)

# How many different shared peptides per geo grouping?
num_shared_peptides <- df_geo_sharing %>% group_by(geo_definition) %>% summarise(
  num=sum(grepl("_",sharing_ID)))

# Define a shared peptide as is seen on both genes in one of the countries
# Then, plot, whether each peptide is shared in each country.
df_geo_sharing_wide <- pivot_wider(df_geo_sharing,
                                   id_cols=c("pos","kmer_seq"),
                                   names_from=geo_definition,
                                   values_from=sharing_ID, 
                                   values_fill = list(sharing_ID="missing")) %>%
  filter(pos >= target_pos_range[1] + 10 & pos <= target_pos_range[2] - 10)

df_geo_sharing_wide_filtered <- df_geo_sharing_wide %>% filter(
  grepl("_",get(RETAINED_GEO_SHARING))
) %>% arrange(pos) %>% select(! starts_with("all_") & ! starts_with("shared_"))
df_heatmap <- pivot_longer(df_geo_sharing_wide_filtered, !c(pos,kmer_seq), names_to="country",values_to="sharing_ID") %>%
 mutate(id = paste(pos,kmer_seq,sep=":")) %>%
  arrange(id)
ggplot(df_heatmap,aes(x=country,y=id,fill=sharing_ID)) + geom_tile() + scale_fill_lancet() +
  theme(text=element_text(size=15),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(angle=15)
        ) +
  labs(y="Shared peptides (ordered by position)",x="Country",fill="Sharing status")
plot_save("sharing_definition/peptide_sharing_heatmap.pdf",width=15,height=11)

# Are shared peptides seen across many countries?
sum_over_countries <- group_by(df_heatmap, id) %>% 
  summarise(num_countries=sum(grepl("_",sharing_ID)))  %>% 
  group_by(num_countries) %>% summarise(num_peptides=n()) %>%
  mutate(frac_peptides = round(num_peptides / sum(num_peptides),2))
p_shared_num<-ggplot(sum_over_countries,aes(x=num_countries,y=num_peptides,label=frac_peptides)) + 
  geom_col() + geom_text(nudge_y=4,size=5.5) +
  labs(x="Number of high-sampled countries where\npeptide is seen across both genes",y="Number of peptides")+
  theme(text=element_text(size=18))
plot_save("sharing_definition/peptide_sharing_freq_histo.pdf",plot=p_shared_num,width=15,height=11)

# Does number of shared peptides increase with sampling?
total_num_peptides <- sum(grepl("_",df_geo_sharing_wide$all_countries_min50_samples))
df_num_peptides <- filter(df_geo_sharing,!grepl("all_|shared_",geo_definition)) %>%
  mutate(is_shared=grepl("_",sharing_ID)) %>%
  filter(is_shared) %>%
  group_by(geo_definition) %>% 
  summarise(num_peptides=sum(is_shared)) %>%
  mutate(country=geo_definition) %>%
  merge(df_used_sample_nums) %>%
  mutate(frac_peptides = num_peptides / total_num_peptides)
ggplot(df_num_peptides,aes(x=num_samples,y=num_peptides,label=country)) + geom_point(size=2) + 
  geom_text_repel(size=5) + labs(x="Number of analysed samples",y="Number of shared peptides found") + 
  theme(text=element_text(size=17)) + 
  geom_text(label=paste("Total num shared peptides:",total_num_peptides),x=550,y=250,size=5)
plot_save("sharing_definition/shared_peptides_by_num_samples.pdf",width=15,height=11)

#################################
####END: Definition of sharing #######
#################################

####################################
####Num shared/private kmers #######
####################################
DBs <- c("DBLMSP","DBLMSP2","Both")
t <- df_geo_sharing %>% 
  filter(geo_definition == RETAINED_GEO_SHARING & seqtype == "protein") %>%
  group_by(pos,sharing_ID) %>% 
  filter(pos>=(target_pos_range[1]) & pos <=(target_pos_range[2]))
df_counts <- df_geo_sharing %>% 
  filter(geo_definition == RETAINED_GEO_SHARING & seqtype == "protein") %>%
  group_by(pos,sharing_ID) %>% summarise(num_kmers=length(unique(kmer_seq))) %>%
  filter(pos>=(target_pos_range[1]) & pos <=(target_pos_range[2]))
df_counts$sharing_ID <- sub(".*_.*","Both",df_counts$sharing_ID)
df_counts$sharing_ID <- factor(df_counts$sharing_ID,levels=DBs)

ggplot(df_counts, aes(x=pos,y=num_kmers)) + 
  geom_line() +
  facet_wrap(vars(sharing_ID),ncol=1) + 
  scale_y_continuous(breaks=c(1,2,4,6,8)) +
  ylab("Number of distinct peptides") +
  xlab("Position in protein") + 
  theme(text=element_text(size=17))
plot_save("counts_and_freqs/DBs_peptide_counts.pdf",width=15,height=11)

####################################
####END: Num shared/private kmers #######
####################################

##########################################
####Per-country sharing frequencies#######
##########################################
query_table <- grep("DBs_protein_sharing_assignment_by_sample",table_names,value = TRUE)
query <- dbSendQuery(con, paste("SELECT * FROM",query_table))
df_shared_kmers <- dbFetch(query)
dbClearResult(query)

df_sharing_levels <- df_shared_kmers %>% group_by(pos, country, gene) %>%
  summarise(frac_shared = sum(is_shared)/n(),num_seqs=n())

high_sampled_countries<-grep("all_|shared_",unique(df_geo_sharing$geo_definition),invert=TRUE,value=TRUE)

p_shared_freqs<-ggplot(filter(df_sharing_levels, country %in% high_sampled_countries &
                pos >= target_pos_range[1]-10 &
                pos <= target_pos_range[2] + 10)
              ,aes(x=pos,y=frac_shared,colour=gene)) + 
  geom_line() + 
  facet_wrap(vars(country)) +
  theme(text=element_text(size=15)) + 
  labs(x="Position in protein alignment",y="Fraction of population with shared peptides",colour="Gene")
plot_save("counts_and_freqs/per_gene_shared_peptide_freqs.pdf", plot=p_shared_freqs,width=14, height=10)
plot_grid(p_shared_num,p_shared_freqs,ncol=1,label_size=20)

##########################################
####END: Per-country sharing frequencies#######
##########################################

#################################
####Diversity and divergence#######
#################################
## For only the shared region
#for (margin in c(10, 0)){
for (margin in c(0)){
  start <- target_pos_range[1] + margin
  end <- target_pos_range[2] - margin
  table_grep <- paste(start,end,"sample_identities",sep="_")
  file_name <- paste("DBs_diversity_divergence_",start,"_",end,".pdf",sep="")
  query_table <- grep(table_grep,table_names,value = TRUE)
  query <- dbSendQuery(con, paste("SELECT * FROM",query_table))
  df_sample_identities <- dbFetch(query)
  df_sample_identities <- filter(df_sample_identities,!grepl("random",gene_ID))
  dbClearResult(query)
  
  ggplot(df_sample_identities, aes(x=percent_identity)) + geom_histogram(bins=50) + 
    facet_wrap(vars(gene_ID),ncol=1) + labs(y="Sample count",x="Fraction of codons shared in the population") +
    theme(text=element_text(size=16))
  plot_save(file_name,width=14,height=12)
}

#################################
###### Positional matching#######
#################################
## This was used to test for whether sequences were more or less shared than
## expected at random, at each position in alignment. Now done at the whole-gene level,
## as above.
#query_table <- grep("DBs_protein_positional_matching",table_names,value = TRUE)
#query <- dbSendQuery(con, paste("SELECT * FROM",query_table))
#df_pos_matching <- dbFetch(query)
#dbClearResult(query)
#
#ggplot(filter(df_pos_matching,pos >= target_pos_range[1] & pos <= target_pos_range[2]),
#       aes(x=pos,y=FPMK, colour = model)) + 
#  geom_line() + 
#  scale_y_log10() +
#  facet_wrap(vars(kmer_size),ncol=1)
#ggsave("positional_linkage.pdf", width=14, height=9)


#################################
###### END: Positional matching#######
#################################

##### Sharedness ####
#df_shared <- read_tsv("/home/brice/Desktop/main_PhD/analyses/plasmo_paralogs/tmp_work/gene_seqs/sharedness/169_408__DBs_thresh0.9_cdhit_k10.tsv")
#df_shared <- read_tsv("/home/brice/Desktop/main_PhD/analyses/plasmo_paralogs/tmp_work/gene_seqs/sharedness/DBs_protein_131_397_thresh96_regfull_sharedness.tsv")
#df_shared <- read_tsv("/home/brice/Desktop/main_PhD/analyses/plasmo_paralogs/tmp_work/gene_seqs/sharedness/DBs_protein_131_397_thresh999_regfull_sharedness.tsv")
#
#plot_sharedness <- function(df_shared, k, num_colours){
#  df_shared_long <- pivot_wider(
#    df_shared %>% filter(kmer_size == k), 
#    id_cols="sample_id",
#    names_from="pos",
#    values_from="colour_mapping")
#num_cols <- dim(df_shared_long)[2]
#data <- as.matrix(df_shared_long[2:num_cols])
#rownames(data) <- df_shared_long$sample_id
#heatmap(data,Colv=NA,scale="none",col=hcl.colors(num_colours,"viridis"))
#}
#
#pdf("/home/brice/Desktop/cur_results/DB_sharedness/positional_sharing_k5_3colours.pdf",width=14,height=9)
#plot_sharedness(df_shared, 5, 3)
#dev.off()
#plot_sharedness(df_shared, 5, 10)
#hmap <- plot_sharedness(df_shared, 5, 10)
#pdf("/home/brice/Desktop/cur_results/DB_sharedness/positional_sharing_k5_100colours.pdf",width=16,height=10)
#plot_sharedness(df_shared, 5, 100)
#dev.off()
#pdf("/home/brice/Desktop/cur_results/DB_sharedness/positional_sharing_k10_3colours.pdf",width=16,height=10)
#plot_sharedness(df_shared, 10, 3)
#dev.off()
#plot_sharedness(df_shared, 10, 10)
#plot_sharedness(df_shared, 10, 20)
#pdf("/home/brice/Desktop/cur_results/DB_sharedness/positional_sharing_k10_50colours.pdf",width=16,height=10)
#plot_sharedness(df_shared, 10, 50)
#dev.off()
#plot_sharedness(df_shared, 10, 1000)
#plot_sharedness(df_shared, 15, 20)
#plot_sharedness(df_shared, 15, 40)
#plot_sharedness(df_shared, 15, 100)
#
#
## Plot specific samples, with genes from same sample sorted together
#df_shared <- read_tsv("/home/brice/Desktop/main_PhD/analyses/plasmo_paralogs/tmp_work/gene_seqs/sharedness/DBs_protein_131_397_all_reg131-397_sharedness.tsv")
#plot_sharedness_for_samples <- function(df_shared, k, num_colours,sample_names, sample_countries = FALSE){
#  df_shared_long <- pivot_wider(
#    df_shared %>% filter(kmer_size == k), 
#    id_cols="sample_id",
#    names_from="pos",
#    values_from="colour_mapping") %>%
#    mutate(sample_id_only=sub(".*_","",sample_id),.before=1) %>% arrange(sample_id_only,sample_id) %>%
#    filter(sample_id_only %in% sample_names)
#  num_cols <- dim(df_shared_long)[2]
#  data <- as.matrix(df_shared_long[3:num_cols])
#  rownames(data) <- df_shared_long$sample_id
#  if (!isFALSE(sample_countries)){
#    matched_countries <- sample_countries %>% filter(Sample %in% sample_names)
#    merged_df <- merge(df_shared_long, matched_countries,by.x="sample_id_only",by.y="Sample")
#    row_names <- paste(merged_df$sample_id,merged_df$Country,sep=":")
#    rownames(data) <- row_names
#  }
#  heatmap(data,Colv=NA,Rowv=NA,scale="none",col=hcl.colors(num_colours,"viridis"))
#}
#
### Plot pacb_ilmn samples
##df_shared <- read_tsv("/home/brice/Desktop/main_PhD/analyses/plasmo_paralogs/tmp_work/gene_seqs/sharedness/DBs_protein_otto_seqs_reg131-397_sharedness.tsv")
##otto_samples <- unique(grep("3D7|otto",df_shared$sample_id, value=TRUE))
##otto_samples<-unique(sub(".*_","",otto_samples))
##plot_sharedness_for_samples(df_shared,10,20,otto_samples)
#
#
#high_shared<- df_frac_linkage %>% filter(kmer_size == 10 & model == "empirical_linkage" & FPMK > 0.5)
#high_shared_samples <- sub(".*:","",high_shared$sample_ID)
#pdf(file="/home/brice/Desktop/cur_results/DB_sharedness/high_shared_samples.pdf",width=16,height=10)
#plot_sharedness_for_samples(df_shared, 10, 20, high_shared_samples, sample_countries)
#dev.off()
#low_shared<- df_frac_linkage %>% filter(kmer_size == 10 & model == "empirical_linkage" & FPMK < 0.01)
#low_shared_samples <- sub(".*:","",low_shared$sample_ID)
#pdf(file="/home/brice/Desktop/cur_results/DB_sharedness/low_shared_samples.pdf",width=16,height=10)
#plot_sharedness_for_samples(df_shared, 10, 20, low_shared_samples, sample_countries)
#dev.off()
#mid_shared<- df_frac_linkage %>% filter(kmer_size == 10 & model == "empirical_linkage" & FPMK > 0.05 & FPMK < 0.4)
#mid_shared_samples <- sub(".*:","",mid_shared$sample_ID)[1:40]
#pdf(file="/home/brice/Desktop/cur_results/DB_sharedness/mid_shared_samples.pdf",width=16,height=10)
#plot_sharedness_for_samples(df_shared, 10, 20, mid_shared_samples, sample_countries)
#dev.off()
#
#
#
#