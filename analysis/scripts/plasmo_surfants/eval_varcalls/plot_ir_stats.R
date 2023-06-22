library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(forcats)
library(cowplot)
library(ggsci)
library(ggrepel)
library(ggExtra)
library(argparser, quietly=TRUE)

#this_dir <- dirname(sys.frame(1)$ofile)
#source("analysis/scripts/common_utils/ir_stat_functions.R")

args <- arg_parser("Plot stats computed by eval_varcalls workflow")
args <- add_argument(args, "ir_stats_all_tsv", help="stats tsv - all analysis-set samples")
args <- add_argument(args, "ir_stats_500_tsv", help="stats tsv - 500 validation samples")
args <- add_argument(args, "pacb_ilmn_tsv", help="stats tsv - truth-assembly-based evaluation")
args <- add_argument(args, "ir_stat_functions", help="Rscript containing plotting utilities")
args <- add_argument(args, "base_dir_for_output_plots", help="output directory for plots")
args <- add_argument(args, "base_dir_for_output_tsvs", help="output directory for tsvs: fold-coverages, 'single-SNP'-solvable samples")
argv <- parse_args(args)

ir_stats_all_tsv <- argv$ir_stats_all_tsv
#ir_stats_all_tsv <- "tmp_work/ir_stats_all.tsv"
ir_stats_500_tsv <- argv$ir_stats_500_tsv
ir_stats_pacb_ilmn <- argv$pacb_ilmn_tsv
BASE_PLOT_DIR <- argv$base_dir_for_output_plots
dir.create(BASE_PLOT_DIR, recursive=TRUE)
BASE_STATS_DIR <- argv$base_dir_for_output_tsvs
dir.create(BASE_STATS_DIR, recursive=TRUE)
source(argv$ir_stat_functions)

plot_save <- function(fname, width, height,plot=""){
  ofname <- paste(BASE_PLOT_DIR, fname, sep="/")
  directory <- dirname(ofname)
  dir.create(directory,recursive=TRUE,showWarnings = FALSE)
  if (plot == ""){
    ggsave(ofname,width=width, height=height)  
  }
  else {
    ggsave(ofname,width=width, height=height, plot=plot)  
  }
  
}


# Some of those that 'work'
filtered_in_genes <- c("AMA1","EBA175","DBLMSP","DBLMSP2","MSP2")
DBs <- c("DBLMSP","DBLMSP2")
# Some of those that don't 'work'
#filtered_in_genes <- c("HRP2","HRP3","LSA1","CLAG3_1","CLAG3_2","CSP")

filtered_in_tools <- c("baseline","cortex","octopus","gram_adju","gapfiller","gram_joint_geno","pf6","pf7","malariaGEN")
tool_colours <- c("baseline","cortex","octopus","gram_adju","gapfiller","gram_joint_geno","malariaGEN")
LANCET_COLOURS <- pal_lancet("lanonc")(length(tool_colours))
names(LANCET_COLOURS) <- tool_colours

###############################################
## Data: Ir stats on 500 validation samples ###
###############################################
df_ir_500 = read_tsv(ir_stats_500_tsv)
df_ir_500$tool <- sub("gram_jointgeno.*gapfiller.*","gram_joint_geno",df_ir_500$tool)
df_ir_500$tool <- sub("gram_adju_.*","gram_adju",df_ir_500$tool)
df_ir_500$tool <- sub("pf6","malariaGEN",df_ir_500$tool)

means_ir_500 <- get_means_ir(df_ir_500)
stat_list <- list(
  stats = c("fraction_positions_4x_or_less","fraction_disagreeing_pileup_min5x",
            "fraction_reads_above_mean_ins_size_plus_two_std", "fraction_reads_below_mean_ins_size_minus_two_std"),
  desc = c("Mean fraction of positions with \n <5x coverage",
           "Mean fraction of positions with \ndisagreeing pileup (at >=5x coverage)",
           "Mean fraction of reads with insert size \n above the mean + two standard deviations",
           "Mean fraction of reads with insert size \n below the mean - two standard deviations"
  )
)

plot_ir_stats <- function(means_ir,in_genes,in_tools,legend=TRUE){
  means_ir_filtered <- means_ir %>% filter(gene %in% in_genes & tool %in% in_tools)
  means_ir_filtered$tool <- factor(means_ir_filtered$tool,levels=filtered_in_tools)
  means_ir_filtered$gene <- factor(means_ir_filtered$gene,levels=in_genes)
  p1<-plot_means(means_ir_filtered, stat_list$stats[2]) + ylab(stat_list$desc[2])
  p2<-plot_means(means_ir_filtered, stat_list$stats[1])+ ylab(stat_list$desc[1])
  if (legend){
    legend <- get_legend(p2 + theme(legend.position = "bottom"))
    p<-plot_grid(p1 + guides(fill="none") + xlab(""),p2 + guides(fill="none"),legend, 
              ncol=1,rel_heights=c(1,1,.2), 
              vjust=1,hjust=-0.2)
  }
  else
    p<-plot_grid(p1 + guides(fill="none") + xlab(""),p2 + guides(fill="none"),
              ncol=1,
              vjust=1,hjust=-0.2)
  p
}
plot_ir_stats(means_ir_500, DBs, filtered_in_tools)
plot_save("pipeline_validation/DBs.pdf", width=13,height=14)
plot_ir_stats(means_ir_500, filtered_in_genes, filtered_in_tools)
plot_save("pipeline_validation/mean_irs_working_gene_subset.pdf", width=13,height=14)

# Plot for all genes, and each panel separately
in_genes <- grep("_DBL",unique(means_ir_500$gene),invert=TRUE,value=TRUE)
means_ir_500_filtered <- means_ir_500 %>% filter(gene %in% in_genes & tool %in% filtered_in_tools)
means_ir_500_filtered$tool <- factor(means_ir_500_filtered$tool,levels=filtered_in_tools)
means_ir_500_filtered$gene <- factor(means_ir_500_filtered$gene,levels=in_genes)
p1<-plot_means(means_ir_500_filtered, stat_list$stats[2]) + ylab(stat_list$desc[2])
plot_save("pipeline_validation/mean_irs_pileup_diffs_allgenes.pdf", width=14,height=12, plot=p1)
p2<-plot_means(means_ir_500_filtered, stat_list$stats[1])+ ylab(stat_list$desc[1])
plot_save("pipeline_validation/mean_irs_pileup_gaps_allgenes.pdf", plot=p2,
       width=14,height=12)

###############################################
## Data: Ir stats on all samples      #########
###############################################
df_ir = read_tsv(ir_stats_all_tsv)
df_ir$tool <- sub("gram_jointgeno.*gapfiller.*","gram_joint_geno",df_ir$tool)
df_ir$tool <- sub("gram_adju_.*","gram_adju",df_ir$tool)
means_ir <- get_means_ir(df_ir)


###############################################
## Filtering samples using ir stats##
###############################################
## Filtering samples based on above stats
df_ir_filtering <- mutate(df_ir, gaps=get(stat_list$stats[1]),pileup_diff=get(stat_list$stats[2]),
                          large_inserts=get(stat_list$stats[3]),small_inserts=stat_list$stats[4]) %>%
  select(sample,tool,gene,gaps,pileup_diff,large_inserts)

single_pileup_errors <- df_ir_filtering %>% filter(
  gene %in% DBs &
    tool == "gram_joint_geno" &
    gaps<=0 &
    pileup_diff >0 &
    pileup_diff < 0.0005 &
    large_inserts <= 0.15
)
write_tsv(single_pileup_errors,paste(BASE_STATS_DIR,"single_SNP_solvable_seqs.tsv",sep="/"))
DBs_pass <- df_ir_filtering %>% filter(
  gene %in% DBs &
    tool == "gram_joint_geno" &
    gaps<=0 &
    pileup_diff <=0 &
    large_inserts <= 0.15
)
  
df_ir_resolved <- df_ir_filtering %>% group_by(tool,gene) %>% 
  summarise(
    num_perfectly_resolved = sum(gaps <= 0& pileup_diff <=0 & large_inserts <=0.15,na.rm=TRUE),
    no_gaps = sum(gaps<=0,na.rm=TRUE),
    plus_no_pileup_diffs = sum(gaps <= 0 & pileup_diff <= 0,na.rm=TRUE),
    plus_no_large_inserts = num_perfectly_resolved,
    num_seqs = sum(!is.na(gaps <= 0 & pileup_diff <=0 & large_inserts <= 0.15))
    ) %>%
        mutate(fraction_perfectly_resolved = num_perfectly_resolved / num_seqs) %>%
  filter(tool %in% filtered_in_tools)
df_ir_resolved$tool <- sub("pf6","malariaGEN",df_ir_resolved$tool)
df_ir_resolved$tool <- factor(df_ir_resolved$tool)

filter_stats <- c("no_gaps","plus_no_pileup_diffs","plus_no_large_inserts")
colours <- pal_npg("nrc")(9)
df_ir_resolved_long <- pivot_longer(df_ir_resolved,
                                    filter_stats, 
                                    names_to="filtering",values_to="remaining_seqs") %>%
  mutate(frac_remaining_seqs=remaining_seqs/num_seqs,
         filtering=factor(filtering,levels=filter_stats))
df_ir_labels <- select(df_ir_resolved_long,c("tool","filtering","remaining_seqs","gene"))
ggplot(df_ir_resolved_long %>% filter(gene %in% DBs),aes(y=frac_remaining_seqs,x=gene,fill=filtering,label=remaining_seqs)) + 
  geom_col(position="dodge") + 
    theme(axis.text.x = element_text(angle=20), text=element_text(size=22),
          legend.position="bottom") + 
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5,size=4.5) + 
  facet_wrap(vars(tool)) + 
  scale_fill_manual(
    values=c(colours[4],colours[7],colours[9]),
    name="Filtering",labels=c("No pileup-based gaps (<5x)","+ No pileup-based differences (min 5x)","+ <15% large read pair inserts")) +
  xlab("Gene") + ylab("Fraction of analysis sequences \nremaining")
plot_save("sample_filtering/filtered_sequences_DBs.pdf", width=15,height=12)

# Find genes with large improvements in gram joint geno
df_ir_wide <- pivot_wider(df_ir_filtering %>% filter(gene == "DBLMSP2" & tool %in% filtered_in_tools),
                          id_cols=c("sample"),names_from="tool",
                  values_from=c("gaps","pileup_diff"))
df_ir_wide_candidates <- df_ir_wide %>% filter(gaps_gram_joint_geno == 0 & pileup_diff_gram_joint_geno == 0
                                               & pileup_diff_pf6 < pileup_diff_baseline &
                                                 pileup_diff_baseline != 0 &
                                                 gaps_pf6 > 0.1)

df_ir_wide <- pivot_wider(df_ir_filtering %>% filter(gene %in% DBs & tool %in% filtered_in_tools),
                          id_cols=c("sample","gene"),names_from="tool",
                  values_from=c("gaps","pileup_diff"))
df_ir_wide_means <- df_ir_wide %>%  group_by(gene) %>% 
  summarise(across(where(is.numeric), ~mean(.x,na.rm=TRUE))) %>%
  mutate(gaps_pf6_rel = 1-(gaps_pf6 / gaps_baseline), pileup_diff_pf6_rel = 1 - (pileup_diff_pf6 / pileup_diff_baseline))

###############################################
## Data: Evaluation on long-read assemblies  ##
###############################################
pacb_truth_df = read_tsv(ir_stats_pacb_ilmn)
pacb_truth_df$tool <- sub("gram_jointgeno.*gapfiller.*","gram_joint_geno",pacb_truth_df$tool)
pacb_truth_df$tool <- sub("gram_adju_.*","gram_adju",pacb_truth_df$tool)
pacb_truth_df$tool <- sub("pf7","malariaGEN",pacb_truth_df$tool)
num_samples <- length(unique(pacb_truth_df$sample))
pacb_truth_means <- pacb_truth_df %>% drop_na() %>% group_by(gene, tool) %>% summarise(
  mean_NM = mean(NM), 
  frac_poorly_resolved=sum(NM>0.05) / num_samples,
  num_poorly_resolved=sum(NM>0.05),
  frac_well_resolved=sum(NM<0.01) / num_samples,
  num_well_resolved=sum(NM<0.01) 
  )


plot_pacb_means <- function(pacb_truth_means,in_genes,in_tools,metric="mean_NM",metric_desc="Mean distance to truth assembly"){
  pacb_truth_means_filtered <- pacb_truth_means %>% filter(gene %in% in_genes & tool %in% in_tools)
  pacb_truth_means_filtered$tool <- factor(pacb_truth_means_filtered$tool,levels=in_tools)
  pacb_truth_means_filtered$gene <- factor(pacb_truth_means_filtered$gene,levels=in_genes)
  p<-ggplot(pacb_truth_means_filtered, aes(y=get(metric),x=gene,fill=tool)) + geom_col(position="dodge") +
    theme(text = element_text(size=18), axis.text.x = element_text(angle=20)) + 
    ylab(metric_desc) + xlab("Gene") + scale_fill_manual(breaks=in_tools,values=LANCET_COLOURS)
  p
}

plot_pacb_means(pacb_truth_means, DBs, filtered_in_tools)
plot_save("pipeline_end_comp/mean_dist_truth_assemblies_DBs.pdf",
       width=10,height=8)
plot_pacb_means(pacb_truth_means, filtered_in_genes,filtered_in_tools)
plot_save("pipeline_end_comp/mean_dist_truth_assemblies_working_gene_subset.pdf",
       width=12,height=11)
#plot_pacb_means(pacb_truth_means, unique(pacb_truth_means$gene),filtered_in_tools)

###############################################
## Plot both evaluation types together       ##
###############################################
plot_both_evals_two_cols <- function(df_ir, pacb_truth_df, in_genes, in_tools){
  p1<-plot_ir_stats(df_ir, in_genes, in_tools,legend=FALSE)
  p2<-plot_pacb_means(pacb_truth_df, in_genes, in_tools) + theme(legend.position = "bottom")
  res <- plot_grid(p2,p1,labels=c("a)","b)"),rel_heights=c(.7,1),rel_widths=c(.8,1),vjust=1,label_size=18)
  return(res)
}
## IMPORTANT: for this plot, load df_ir for all samples, not just the 500
tools_pipelines <- c("baseline","gram_joint_geno","malariaGEN")
plot_both_evals_two_cols(means_ir %>% mutate(tool=sub("pf6","malariaGEN",tool)), 
                         pacb_truth_means, 
                         DBs, tools_pipelines)
plot_save("pipeline_end_comp/eval_varcalls_DBs_pipelines.pdf",
       width=13.5,height=11)

plot_both_evals_two_cols(means_ir_500 %>% filter(tool != "malariaGEN"), 
                         pacb_truth_means %>% filter(tool!="malariaGEN"), 
                         DBs, filtered_in_tools)
plot_save("pipeline_validation/eval_varcalls_DBs_pipelines.pdf",
       width=13.5,height=11)


plot_both_evals_one_col <- function(df_ir, pacb_truth_df, in_genes, in_tools){
  means_ir_filtered <- df_ir %>% filter(gene %in% in_genes & tool %in% in_tools)
  means_ir_filtered$tool <- factor(means_ir_filtered$tool,levels=in_tools)
  means_ir_filtered$gene <- factor(means_ir_filtered$gene,levels=in_genes)
  p1<-plot_pacb_means(pacb_truth_df, in_genes, in_tools) + guides(fill="none") + xlab("") +
    theme(axis.text.x = element_blank()) 
  p2<-plot_means(means_ir_filtered, stat_list$stats[2]) + 
    ylab(stat_list$desc[2]) +
    xlab("") + guides(fill="none") + 
    theme(axis.text.x = element_blank()) 
  p3<-plot_means(means_ir_filtered, stat_list$stats[1])+ 
    ylab(stat_list$desc[1]) + 
    theme(legend.position="bottom",axis.text.x=element_text(angle=10))
  legend <- get_legend(p3)
  res <- plot_grid(p1,p2,p3 + guides(fill="none"),legend,
                   labels=c("a)","b)",""),vjust=1,hjust=-0.2,
                   rel_heights=c(1.2,1.2,1.2,0.3),label_size=18,ncol=1)
  res
  return(res)
}
plot_both_evals_one_col(means_ir, pacb_truth_means, filtered_in_genes, tools_pipelines)
plot_save("pipeline_end_comp/eval_varcalls_working_gene_subset_gram_joint_onecol.pdf",
       width=13,height=15)


###############################################
## Plot distributions instead of means       ##
###############################################
plot_list <- list()
for (gene in DBs){
  for (i in seq(1,2)){
    p <- plot_gene_stat_histos(
      mutate(filter(df_ir, tool %in% filtered_in_tools),tool=sub("pf6","malariaGEN",tool)) %>%
      mutate(tool=factor(tool,levels=filtered_in_tools))
      ,stat_list$stats[i],gene) + 
      xlab(sub("Mean fraction","Fraction",stat_list$desc[i])) +
      #ggtitle(paste("Gene: ",gene)) +
      ylab("Sample count")
  plot_list[[paste(gene,i)]] <- p
  #ggsave(paste("Chapter_plasmodium_antigens/eddist_evals/histograms/eval_varcalls_gene_stat_histos_",gene,"_",stat_list$stats[i],".pdf"),
  #     width=13,height=11)
  }
}
plot_list$`DBLMSP 1` <-plot_list$`DBLMSP 1` + xlab("")
plot_list$`DBLMSP 2` <-plot_list$`DBLMSP 2` + xlab("") + ylab("")
plot_list$`DBLMSP2 2` <-plot_list$`DBLMSP2 2` + ylab("")
plot_grid(plotlist=plot_list, ncol=2, labels=c("a)DBLMSP","","b)DBLMSP2",""))
plot_save("pipeline_end_comp/histograms/eval_varcalls_gene_stat_histos_DBs.pdf",
     width=14,height=13)

gene <- "AMA1"
plot_gene_stat_histos(
  df_ir %>% filter(tool %in% filtered_in_tools) %>%
  mutate(tool=sub("pf6","malariaGEN",tool)) %>%
  mutate(tool=factor(tool,levels=filtered_in_tools))
  ,stat_list$stats[3],gene) + 
  xlab(sub("Mean fraction","Fraction",stat_list$desc[3])) +
  ggtitle(paste("Gene: ",gene)) +
  ylab("Sample count")
  plot_save(paste("pipeline_end_comp/histograms/eval_varcalls_gene_stat_histos_",gene,"_",stat_list$stats[3],".pdf"),
       width=13,height=11)

## Compare tools for capacity to resolve genes (e.g. cortex vs octopus)
callers_NM_df <- pivot_wider(pacb_truth_df,id_cols=c("sample","gene"),names_from=tool,values_from=NM)

plot_pairwise <- function(callers_NM_df, tool1, tool2){
  freqs <- seq(0,max(drop_na(mutate(callers_NM_df,sample=NULL,gene=NULL))),by=0.01)
  ggplot(callers_NM_df,aes(x=get(tool1),y=get(tool2))) + geom_point(aes(colour=gene)) + 
    geom_line(data=tibble(x=freqs,y=freqs),aes(x=x,y=y)) + xlab(tool1) + ylab(tool2)
}
plot_pairwise(callers_NM_df, "cortex","octopus")
plot_pairwise(filter(callers_NM_df, gene %in% filtered_in_genes), "cortex","octopus")
genes <- c("DBLMSP","DBLMSP2","SERA5")
#ggplot(filter(callers_NM_df,gene %in% genes),aes(x=cortex,y=octopus)) + geom_point(aes(colour=gene)) + 
#  geom_line(data=tibble(x=freqs,y=freqs),aes(x=x,y=y))

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
plot_save("metric_correlations/regression_NM_map_err.pdf", width=13, height=9)
plot_fit(merged_df,"mean_NM","fraction_disagreeing_pileup_min5x")

lm_fit<-lm(mean_NM~base_error_rate,merged_df)
summary(lm_fit)


###############################################
########### Coverage-based analysis ###########
###############################################
library(ggExtra)
df_ir_covs <- df_ir %>% filter(tool == "gram_joint_geno") %>% 
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

p<-plot_putative_CNVs(df_ir_covs_with_comps,dot_size=TRUE)
plot_save("coverage_analysis/fold_coverages_all_analysis_set_samples.pdf",plot=p,height=10, width=15)

DB1_pass <- filter(DBs_pass, gene == "DBLMSP")
DB2_pass <- filter(DBs_pass, gene == "DBLMSP2")
df_ir_covs_with_comps_filterpass <- filter(df_ir_covs_with_comps, 
                                           (gene == "DBLMSP" & sample %in% unique(DB1_pass$sample)) | 
                                            (gene == "DBLMSP2" & sample %in% unique(DB2_pass$sample)) )
p<-plot_putative_CNVs(df_ir_covs_with_comps_filterpass)
plot_save("coverage_analysis/fold_coverages_filterpass_analysis_set_samples.pdf",plot=p,height=10, width=15)

high_ratios <- df_ir_covs_with_comps_filterpass %>% filter(gene %in% DBs & mean_ratio > 2 & overlap < 0.5)
write_tsv(high_ratios,paste(BASE_STATS_DIR,"putative_duplications_in_filterpass.tsv",sep="/"))
#t <- df_ir_covs_with_comps %>% filter(sample %in% unique(high_ratios$sample)) %>% filter(gene %in% c(DBs,"AMA1"))
#t2 <- df_ir_covs_with_comps %>% filter(sample %in% unique(low_overlaps$sample)) %>% filter(gene %in% c(DBs,"AMA1"))

#p2<-ggplot(df_ir_covs_with_comps %>% filter(gene %in% DBs),aes(x=mean_ratio,y=mean_read_coverage,colour=gene)) + 
#  geom_point() + 
#  scale_x_continuous(breaks=c(0,1,2,5,10)) +
#  scale_y_continuous(breaks=c(0,10,20,50,100,150,300)) +
#  xlab(xlab_text) + ylab("Mean coverage")
#ggMarginal(p2,type="histogram")
#
#p3<-ggplot(df_ir_covs_with_comps %>% filter(gene %in% DBs),aes(x=mean_read_coverage,y=overlap,colour=gene)) + 
#  geom_point() + 
#  xlab("Mean read coverage") + ylab(ylab_text)
#ggMarginal(p3,type="histogram")


write_tsv(df_ir_covs_with_comps,paste(BASE_STATS_DIR,"ir_stats_fold_coverages.tsv",sep="/"))


###############################################
## Well vs poorly-resolved genes##
###############################################
used_tools <-c("baseline","gapfiller","malariaGEN")
all_genes <- unique(pacb_truth_df$gene)
plot_pacb_means(pacb_truth_means, all_genes, used_tools)
plot_save("gene_resolvedness/mean_dist_truth_assemblies_all_genes.pdf",
       width=16,height=12)
plot_pacb_means(pacb_truth_means, all_genes, used_tools,metric="num_well_resolved","Number low-distance samples (NM <0.01)")
plot_save("gene_resolvedness/num_low_dist_truth_assemblies_all_genes.pdf",
       width=16,height=12)
plot_pacb_means(pacb_truth_means, all_genes, used_tools,metric="num_poorly_resolved","Number high-distance samples (NM >0.05)")
plot_save("gene_resolvedness/num_high_dist_truth_assemblies_all_genes.pdf",
       width=16,height=12)
well_resolved <- pacb_truth_means %>% filter(num_well_resolved >= 7 & num_poorly_resolved <= 2)
t<-well_resolved %>% group_by(tool) %>% summarise(paste(gene,collapse=","))
t2 <- well_resolved %>% filter(tool == "gapfiller")
gapfiller_resolved <- t2$gene
poorly_resolved <- setdiff(all_genes,gapfiller_resolved)


###############################################
## Tandem repeat analysis ##
###############################################
if (FALSE){
repeat_df = read_tsv("TODO")

longest_repeats <- repeat_df %>% filter(tool == "truth") %>% group_by(gene) %>%
  summarise(max_repeat_len=max(total_size))
longest_repeats <- full_join(longest_repeats,filter(pacb_truth_means,tool=="gapfiller")) 
longest_repeats$max_repeat_len <- replace_na(longest_repeats$max_repeat_len,0)
lr_plt <- ggplot(longest_repeats,aes(x=mean_NM,y=max_repeat_len,label=gene)) +
  geom_point(size=3)+
  geom_text_repel(size=5.5)+
  xlab("Mean distance to truth assemblies") +
  ylab("Maximal tandem repeat\n size in gene population (bp)") +
  scale_y_log10() +
  theme(text = element_text(size=18))
# Match up repeat sequences for given gene + sample, and compare
# the repeat size inference between my pipeline and the truth assembly.
# Doesn't work for repeat sequences that differ slightly, but are 'really' the 
# same. Clustering would be required to get to that.
repeat_df_long <- pivot_wider(repeat_df, names_from=tool,values_from=total_size,
                              id_cols=c(gene,sample,repeat_consensus),
                              values_fn=function(elems){return(elems[1])})
repeat_diffs <-repeat_df_long %>% filter(truth != pipeline_inference) %>% 
  mutate(diff=abs(truth - pipeline_inference))
repeat_diffs <- full_join(repeat_diffs,filter(pacb_truth_df,tool=="gapfiller"))
gene_sizes = read_tsv("/home/brice/Desktop/main_PhD/analyses/plasmo_surfants/tmp_work/repeats/gene_sizes.tsv")
repeat_diffs <- full_join(repeat_diffs,gene_sizes) %>% mutate(diff_NM = diff / size)

rd_plt <- ggplot(repeat_diffs,aes(x=NM,y=truth,label=gene,size=diff)) + 
  geom_point() + 
  geom_text_repel(size=5) + 
  xlab("Distance to sample truth assembly") +
  ylab("Tandem repeat size\n in sample gene (bp)") +
  labs(size="Repeat size error (bp)") +
  theme(text = element_text(size=18))

plot_grid(lr_plt, rd_plt, ncol=1,rel_heights=c(0.55,0.45),labels=c("a)","b)"),label_size=18)
plot_save("repeat_analysis/tandem_repeats.pdf",height=16,width=13)
}