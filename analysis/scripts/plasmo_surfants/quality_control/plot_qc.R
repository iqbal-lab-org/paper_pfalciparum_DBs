library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(forcats)
library(argparser, quietly=TRUE)

df_rs = read_tsv("/home/brice/Desktop/main_PhD/analyses/plasmo_surfants/tmp_work/quality_control/read_stats.tsv")

ggplot(df_rs %>% filter(metric_category == "mapping_depth" & metric != "stdv"),aes(x=value,fill=metric)) + geom_histogram()
ggplot(df_rs %>% filter(metric_category == "base_quality" & metric != "stdv"),aes(x=value,fill=metric)) + geom_histogram()
ggplot(df_rs %>% filter(metric_category == "read_length" & metric != "stdv"),aes(x=value,fill=metric)) + geom_histogram()


df_vs = read_tsv("/home/brice/Desktop/main_PhD/analyses/plasmo_surfants/tmp_work/quality_control/vcf_stats.tsv")
df_vs <- df_vs %>% mutate(midpoint = start + ((end - start) / 2),
                          avg=value / (end - start + 1))
dblmsp <- filter(df_vs, gene == "DBLMSP" & tool == "gram_jointgeno_ebf7bcd5__pf6_analysis_set_fws95__7__13" &
                   metric %in% c("DP","COV"))
ggplot(dblmsp, aes(x=midpoint,y=avg,group=midpoint)) + facet_wrap(vars(metric)) + geom_boxplot(outlier.alpha = 0.1)

dblmsp2 <- filter(df_vs, gene == "DBLMSP2" & tool == "gram_jointgeno_ebf7bcd5__pf6_analysis_set_fws95__7__13" &
                   metric %in% c("DP","COV"))
ggplot(dblmsp2, aes(x=midpoint,y=avg,group=midpoint)) + facet_wrap(vars(metric)) + geom_boxplot()

eba175 <- filter(df_vs, gene == "EBA175" & tool == "gram_jointgeno_ebf7bcd5__pf6_analysis_set_fws95__7__13" &
                   metric %in% c("DP","COV"))
ggplot(eba175, aes(x=midpoint,y=avg,group=midpoint)) + facet_wrap(vars(metric)) + geom_boxplot()

