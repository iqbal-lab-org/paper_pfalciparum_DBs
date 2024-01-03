library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(forcats)
library(argparser, quietly=TRUE)
library(ggsci)
library(scales)

df_rs = read_tsv("/home/brice/Desktop/main_PhD/analyses/plasmo_surfants/tmp_work/quality_control/read_stats.tsv")
used_metrics <- c("base_quality","read_length","insert_size_below_10000","mapping_depth")
df_rs <- filter(df_rs,metric_category %in% used_metrics)
df_rs$metric_category <- factor(df_rs$metric_category,levels=used_metrics,
                                  labels=c(
                                    "Read base quality (Phred-scale)",
                                    "Read length",
                                    "Sequenced fragment lengths (<10000)",
                                    "Read fold-coverage"
                                  )
                                )
ggplot(df_rs %>% filter(metric == "mean"),
       aes(x=value)) + geom_histogram(colour="#1B1919FF",fill="#00468BFF") + facet_wrap(vars(metric_category),scales="free",ncol=2) +
  xlab("Metric value") + ylab("Number of samples") + theme(text = element_text(size=15))
ggsave("/home/brice/Desktop/main_PhD/analyses/plasmo_surfants/tmp_work/quality_control/read_qc.pdf",width=12,height=10)

df_rs_wide <- pivot_wider(filter(df_rs,metric=="mean"),names_from=metric_category,values_from=value)
ggplot(df_rs_wide,aes(x=read_length,y=insert_size_below_10000)) + geom_point()

show_col(pal_lancet("lanonc")(9))
