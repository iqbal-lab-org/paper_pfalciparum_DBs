library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(forcats)
library(ggsci)

get_means_ir <- function(df_ir){
  result <- df_ir %>%  group_by(gene, tool) %>% 
    summarise(across(where(is.numeric), ~mean(.x,na.rm=TRUE)))
  return(result)
}

plot_means <- function (means_ir,stat_name){
  res <- ggplot(means_ir, aes_string(y=stat_name,x="gene",fill="tool")) + geom_col(position="dodge") +
    scale_fill_manual(values=LANCET_COLOURS) +
  theme(text = element_text(size=18), axis.text.x = element_text(angle=30)) + 
    ylab(paste("Mean",stat_name)) + xlab("Gene")
  res
  return(res)
}

plot_gene_stat_tool_histo <- function(df_ir, stat_name, gene_name,tool_name) {
  ggplot(df_ir %>% filter(gene == gene_name & tool == tool_name),aes_string(x=stat_name)) +
  geom_histogram(bins=100)
}
plot_gene_stat_histos <- function(df_ir, stat_name, gene_name) {
  means<-df_ir %>% filter(gene == gene_name) %>% group_by(tool) %>% summarise(
    means=mean(get(stat_name),na.rm=TRUE),
    adjust=mean(get(stat_name),na.rm=TRUE) + 0.2 * max(get(stat_name)))
  p<-ggplot(df_ir %>% filter(gene == gene_name),aes_string(x=stat_name)) + 
    geom_histogram(bins=100) + theme(text=element_text(size=17), axis.text.x=element_text(angle=20))+
    facet_wrap(vars(tool))
  p <- p + geom_vline(data=means,aes(xintercept=means,colour="red")) + 
    geom_text(data=means,aes(x=adjust,label=round(means,3),y=1000)) +
    theme(legend.position="none")
  p
  return(p)
}