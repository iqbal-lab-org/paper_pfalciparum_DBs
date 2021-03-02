library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(forcats)
library(argparser, quietly=TRUE)

p <- arg_parser("Plot alignments")
p <- add_argument(p, "stats_file", help="")
p <- add_argument(p, "outdir", help="")
argv <- parse_args(p)
df = read_tsv(argv$stats_file)

means <- df %>% drop_na() %>% group_by(gene, tool) %>% summarise(mean_NM = mean(NM))
write_tsv(means,file.path(argv$outdir, "mean_NM.tsv"))

ggplot(means, aes(y=mean_NM,x=gene,fill=tool)) + geom_col(position="dodge")
ggsave(file.path(argv$outdir, "R_mean_NM.pdf"), width=14, height=9)

p <- ggplot(df, aes(x=NM)) + geom_histogram(binwidth=0.005)
p + facet_grid(rows=vars(gene),cols=vars(tool))
ggsave(file.path(argv$outdir, "R_NM_facetgrid.pdf"), width=15, height=12)

NAs<- df %>% group_by(tool, gene) %>% summarise(num_no_mapping=sum(is.na(NM)))
ggplot(NAs, aes(y=num_no_mapping,x=gene,fill=tool)) + geom_col(position="dodge")
ggsave(file.path(argv$outdir, "R_num_no_mapping.pdf"), width=14, height=9)
