library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(forcats)
library(cowplot)
library(ggsci)
library(argparser, quietly=TRUE)

this_dir <- dirname(sys.frame(1)$ofile)
setwd("/home/brice/Desktop/research/PhD/main_PhD/analyses/plasmo_surfants")
source("analysis/scripts/common_utils/ir_stat_functions.R")

df_ir = read_tsv("tmp_work/ir_stats.tsv")

df_ir$tool <- sub("gram_jointgeno.*gapfiller.*","gram_joint_geno_gapfiller",df_ir$tool)
df_ir$tool <- sub("gram_jointgeno.*gram_adju.*","gram_joint_geno_gram_adju",df_ir$tool)
df_ir$tool <- sub("gram_adju_.*","gram_adju",df_ir$tool)
df_ir <- df_ir %>% filter(!grepl("fws95",tool))
means_ir <- get_means_ir(df_ir)

## Filter to gene subset
filtered_out_genes <- c("MSP4","MSP5","MSP7","P48_45","Pfs25","RH5","SURF4_2","TRAP","CelTOS","HRP2")
filtered_in_tools <- c("baseline","pf6","gram_joint_geno_gapfiller")
means_ir_filtered <- means_ir %>% filter(! gene %in% filtered_out_genes & tool %in% filtered_in_tools)
means_ir_filtered <- means_ir %>% filter(! gene %in% filtered_out_genes)


plot_means(means_ir_filtered, "base_error_rate")
plot_means(means_ir_filtered, "fraction_disagreeing_pileup_min5x")
plot_means(means_ir_filtered, "fraction_positions_0x_or_less")
plot_means(means_ir_filtered, "fraction_reads_paired_one_unmapped")
plot_means(means_ir_filtered, "fraction_reads_properly_paired_aligner")
plot_means(means_ir_filtered, "fraction_reads_above_mean_ins_size_plus_two_std")
plot_means(means_ir_filtered, "max_insert_size")

ggsave("/home/brice/Desktop/cur_results/eval_varcalls/fraction_reads_above_mean_ins_size_plus_two_std.pdf", width=14, height=9)

plot_gene_stat_tool_histo(df_ir,"fraction_positions_0x_or_less","DBLMSP_DBL","gapfiller")
plot_gene_stat_histos(df_ir,"fraction_positions_0x_or_less","DBLMSP_DBL")
plot_gene_stat_histos(df_ir,"fraction_positions_4x_or_less","DBLMSP_DBL")
plot_gene_stat_histos(df_ir,"fraction_positions_0x_or_less","DBLMSP2_DBL")
plot_gene_stat_histos(df_ir,"fraction_disagreeing_pileup_min5x","DBLMSP_DBL")
plot_gene_stat_histos(df_ir,"fraction_disagreeing_pileup_min5x","DBLMSP2_DBL")
plot_gene_stat_histos(df_ir,"fraction_reads_above_mean_ins_size_plus_two_std","DBLMSP_DBL")
plot_gene_stat_histos(df_ir,"fraction_reads_above_mean_ins_size_plus_two_std","DBLMSP2_DBL")
plot_gene_stat_histos(df_ir,"fraction_reads_above_mean_ins_size_plus_two_std","EBA175")
plot_gene_stat_histos(df_ir,"fraction_reads_above_mean_ins_size_plus_two_std","AMA1")
plot_gene_stat_histos(df_ir,"fraction_reads_properly_paired_aligner","DBLMSP_DBL")

plot_gene_stat_histos(df_ir,"fraction_positions_0x_or_less","DBLMSP_DBL")

ggsave(file.path("/home/brice/Desktop/main_PhD/analyses/plasmo_surfants/tmp_work/plots",
                 "mean_err_rate_per_tool_gene_subset.pdf"), width=13, height=9)

## Filtering samples based on above stats
filter_by <- function(df_ir, gene_name, tool_name, stat_names, values){
  initial_size<-length(unique(df_ir$sample))
 gene_df <- df_ir %>% filter(gene == gene_name & grepl(tool_name,tool))
 for (i in seq(length(stat_names))){
   if (grepl("positions_[0-9]x|disagreeing_pileup|reads_above_mean_ins",stat_names[i])){
     gene_df <- gene_df %>% filter(get(stat_names[i]) < values[i])
   } else if (grepl("properly_paired",stat_names[i])) {
     gene_df <- gene_df %>% filter(get(stat_names[i]) > values[i])
   } else {
     print("Unrecognised stat")
     return()
   }
   new_size <-length(unique(gene_df$sample))
   print(paste("Leaves ",new_size,"samples (",round(new_size/initial_size,3),"%)"))
 }
 return(gene_df)
}

filter_by(df_ir, "DBLMSP", "joint_geno_gapfiller",
          c("fraction_positions_4x_or_less","fraction_disagreeing_pileup_min5x","fraction_reads_above_mean_ins_size_plus_two_std"),
          c(0.00001,0.00001,0.15))
filter_by(df_ir, "DBLMSP", "pf6",
          c("fraction_positions_4x_or_less","fraction_disagreeing_pileup_min5x","fraction_reads_above_mean_ins_size_plus_two_std"),
          c(0.00001,0.00001,0.15))
DB1_pass <- filter_by(df_ir, "DBLMSP", "joint_geno_gapfiller",
          c("fraction_positions_4x_or_less","fraction_disagreeing_pileup_min5x","fraction_reads_above_mean_ins_size_plus_two_std"),
          c(0.00001,0.00001,0.15))
DB2_pass <- filter_by(df_ir, "DBLMSP2", "joint_geno_gapfiller",
          c("fraction_positions_4x_or_less","fraction_disagreeing_pileup_min5x","fraction_reads_above_mean_ins_size_plus_two_std"),
          c(0.00001,0.00001,0.15))

## Get samples that are improved in gramtools pipeline
pivoted <- pivot_wider(df_ir, 
                       id_cols=c("sample","gene"),
                       names_from="tool",
                       values_from="fraction_positions_0x_or_less")
pivoted <- mutate(pivoted,gram_diff=gram_joint_geno_gapfiller - pf6) %>% arrange(sample)

## Get genes that are 'perfectly resolved'
means_ir %>% group_by(tool) %>% summarise(sum(fraction_positions_0x_or_less < 0.01 & fraction_disagreeing_pileup_min5x < 0.01))

## Plot the truth assembly based validation
pacb_truth_df = read_tsv("/home/brice/Desktop/main_PhD/analyses/plasmo_surfants/tmp_work/pacb_ilmn_stats.tsv")
pacb_truth_df$tool <- sub("gram_jointgeno.*gapfiller.*","gram_joint_geno_gapfiller",pacb_truth_df$tool)
pacb_truth_df$tool <- sub("gram_jointgeno.*gram_adju.*","gram_joint_geno_gram_adju",pacb_truth_df$tool)
pacb_truth_df$tool <- sub("gram_adju_.*","gram_adju",pacb_truth_df$tool)
pacb_truth_df <- pacb_truth_df %>% filter(! grepl("fws95",tool))
pacb_truth_means <- pacb_truth_df %>% drop_na() %>% group_by(gene, tool) %>% summarise(mean_NM = mean(NM))
filtered_out_genes <- c("MSP4","MSP5","MSP7","P48_45","Pfs25","RH5","SURF4_2","TRAP","CelTOS","HRP2")
pacb_truth_means_filtered <- pacb_truth_means %>% filter(! gene %in% filtered_out_genes)
#kept_genes <- c("DBLMSP","DBLMSP2","EBA175","MSP1","HRP2","LSA1")
#kept_tools <- c("baseline","pf7","gramtools_joint_genotyping")
#pacb_truth_means_filtered <- pacb_truth_means %>% filter( gene %in% kept_genes & tool %in% kept_tools)
#pacb_truth_means_filtered$tool <- factor(pacb_truth_means_filtered$tool, levels=c("baseline","pf7","gramtools"))
                                                           

ggplot(pacb_truth_means_filtered, aes(y=mean_NM,x=gene,fill=tool)) + geom_col(position="dodge") +
  theme(text = element_text(size=14), axis.text.x = element_text(angle=0)) + 
  ylab("Mean distance to truth assembly") + xlab("Gene") + scale_fill_lancet()
ggsave("/home/brice/Desktop/cur_results/R_mean_NM.pdf", width=14, height=9)
ggsave("/home/brice/Desktop/cur_results/R_mean_NM_subset_of_genes_2.pdf", width=11, height=8)


mean_truth_per_tool <-pacb_truth_means_filtered %>% group_by(tool) %>% summarise(mean=mean(mean_NM))
ggplot(pacb_truth_means_filtered, aes(x=mean_NM)) + geom_histogram(bins=40) + 
  facet_grid(rows=vars(tool)) + geom_vline(data=mean_truth_per_tool,aes(xintercept=mean),linetype=3)  +
  theme(text = element_text(size=15)) + ylab("Gene count") + xlab("Mean distance to truth assembly")
ggsave(file.path("/home/brice/Desktop/main_PhD/analyses/plasmo_surfants/tmp_work/plots",
                 "mean_truth_dist_per_tool_histo.pdf"), width=13, height=9)

## How well does mapping based stat correlate with distance to truth assembly?
merged_df <- means_ir %>% inner_join(pacb_truth_means,by=c("gene","tool"))

plot_fit<-function(merged_df, stat1, stat2){
  lm_fit<-lm(get(stat1)~get(stat2),merged_df)
  predicted_df <- data.frame(prediction = predict(lm_fit, merged_df), merged_df[stat2])
  r_squared <-round(summary(lm_fit)$r.squared,3)

  ggplot(merged_df,aes_string(x=stat2, y=stat1)) + geom_point(aes(colour=tool))+
    geom_line(data=predicted_df,aes(x=get(stat2),y=prediction),linetype="dotted") +
    annotate("text",x=0.021,y=0.1,label=paste("R-squared: ",toString(r_squared),sep="")) +
    xlab(paste("Mean",stat2)) + ylab(paste("Mean",stat1)) + 
    theme(text = element_text(size=15)) 
}
plot_fit(merged_df,"mean_NM","base_error_rate")
plot_fit(merged_df,"mean_NM","fraction_disagreeing_pileup_min5x")
ggsave(file.path("/home/brice/Desktop/main_PhD/analyses/plasmo_surfants/tmp_work/plots",
                 "regression_NM_map_err.pdf"), width=13, height=9)

lm_fit<-lm(mean_NM~base_error_rate,merged_df)
summary(lm_fit)


######### Coverage-based analysis #########
library(ggExtra)
df_ir_covs <- df_ir %>% filter(tool == "gram_joint_geno_gapfiller") %>% 
  filter(! is.nan(mean_read_coverage)) %>%
  select(contains("read_coverage") | c("sample","gene","tool")) %>% 
  mutate(lower_bound = pmax(mean_read_coverage - 2*std_read_coverage,0), upper_bound = mean_read_coverage + 2*std_read_coverage)

interval_overlap <- function(ref_l, ref_h, target_l, target_h){
  # Below two conditions model containment, in which case overlap is complete
  if (target_l >= ref_l && target_h <= ref_h) return(1)
  if (ref_l >= target_l && ref_h <= target_h) return(1)
  # Below computes overlap of target to ref, divided by target length.
  len <- target_h - target_l
  if (target_l <= ref_l){
    result <- max((target_h - ref_l)/len, 0)
  }
  else {
    result <- max((ref_h - target_l)/len, 0)
  }
  return(min(result, 1))
}

cov_comparator <- function(sub_df, ignored_arg){
 ref_gene <- filter(sub_df, gene == "AMA1") 
 ref_lower <- ref_gene$lower_bound
 ref_mean <- ref_gene$mean_read_coverage
 ref_upper <- ref_gene$upper_bound
 result <- sub_df %>% rowwise() %>% 
   mutate(overlap = interval_overlap(ref_lower, ref_upper, lower_bound, upper_bound), mean_ratio=mean_read_coverage / ref_mean)
 return(result)
}

df_ir_covs_with_comps <- df_ir_covs %>% group_by(sample) %>% group_modify(cov_comparator)
DBs <- c("DBLMSP", "DBLMSP2")
genes <- c("DBLMSP", "DBLMSP2","CelTOS","Pfs25")
DBs_DBLs <- c("DBLMSP_DBL", "DBLMSP2_DBL")
xlab_text <- "Ratio of mean coverage in gene to mean coverage in AMA1"
ylab_text <- "Coverage overlap with AMA1"

plot_putative_CNVs <- function(df,dot_size=FALSE){
  if (dot_size){
    p<-ggplot(df %>% filter(gene %in% DBs),aes(x=mean_ratio,y=overlap,colour=gene,size=mean_read_coverage)) +
    scale_size(name="Mean read \ncoverage in gene")
  }
  else p<-ggplot(df %>% filter(gene %in% DBs),aes(x=mean_ratio,y=overlap,colour=gene))
  p<- p + geom_point() + 
    scale_x_continuous(breaks=c(0,1,2,5,10)) +
    xlab(xlab_text) + ylab(ylab_text) +
    scale_colour_discrete(name="Gene")+ theme(text = element_text(size=15))
  p<-ggMarginal(p,type="histogram") 
  p
  return(p)
}

df_ir_covs_with_comps_filterpass <- filter(df_ir_covs_with_comps, 
                                           (gene == "DBLMSP" & sample %in% unique(DB1_pass$sample)) | 
                                            (gene == "DBLMSP2" & sample %in% unique(DB2_pass$sample)) )
p<-plot_putative_CNVs(df_ir_covs_with_comps,dot_size=TRUE)
ggsave("~/Desktop/cur_results/paper_figs/fold_coverages.pdf",plot=p,height=10, width=15)
plot_putative_CNVs(df_ir_covs_with_comps_filterpass)

p2<-ggplot(df_ir_covs_with_comps %>% filter(gene %in% DBs),aes(x=mean_ratio,y=mean_read_coverage,colour=gene)) + 
  geom_point() + 
  scale_x_continuous(breaks=c(0,1,2,5,10)) +
  scale_y_continuous(breaks=c(0,10,20,50,100,150,300)) +
  xlab(xlab_text) + ylab("Mean coverage")
ggMarginal(p2,type="histogram")

p3<-ggplot(df_ir_covs_with_comps %>% filter(gene %in% DBs),aes(x=mean_read_coverage,y=overlap,colour=gene)) + 
  geom_point() + 
  xlab("Mean read coverage") + ylab(ylab_text)
ggMarginal(p3,type="histogram")

high_ratios <- df_ir_covs_with_comps_f %>% filter(gene %in% DBs & mean_ratio > 2)
t <- df_ir_covs_with_comps %>% filter(sample %in% unique(high_ratios$sample)) %>% filter(gene %in% c(DBs,"AMA1"))
low_overlaps <- df_ir_covs_with_comps %>% filter(gene %in% DBs & overlap < 0.1)
t2 <- df_ir_covs_with_comps %>% filter(sample %in% unique(low_overlaps$sample)) %>% filter(gene %in% c(DBs,"AMA1"))

write_tsv(df_ir_covs_with_comps,"tmp_work/ir_stats_fold_coverages.tsv")
