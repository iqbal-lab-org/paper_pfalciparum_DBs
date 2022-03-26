library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(forcats)
library(argparser, quietly=TRUE)

get_means_ir <- function(df_ir){
  result <- df_ir %>%  group_by(gene, tool) %>% 
    summarise(across(where(is.numeric), ~mean(.x,na.rm=TRUE)))
  return(result)
}

plot_means <- function (means_ir,stat_name){
  ggplot(means_ir, aes_string(y=stat_name,x="gene",fill="tool")) + geom_col(position="dodge") +
  theme(text = element_text(size=14), axis.text.x = element_text(angle=30)) + ylab(paste("Mean",stat_name))
}

plot_gene_stat_tool_histo <- function(df_ir, stat_name, gene_name,tool_name) {
  ggplot(df_ir %>% filter(gene == gene_name & tool == tool_name),aes_string(x=stat_name)) +
  geom_histogram(bins=100)
}
plot_gene_stat_histos <- function(df_ir, stat_name, gene_name) {
  means<-df_ir %>% filter(gene == gene_name) %>% group_by(tool) %>% summarise(means=mean(get(stat_name)))
  p<-ggplot(df_ir %>% filter(gene == gene_name),aes_string(x=stat_name)) + 
    geom_histogram(bins=100) +
    facet_wrap(vars(tool))
  p + geom_vline(data=means,aes(xintercept=means,colour="red")) + 
    geom_text(data=means,aes(x=means,label=round(means,3),y=200),nudge_x=0.03) +
    theme(legend.position="none")
}