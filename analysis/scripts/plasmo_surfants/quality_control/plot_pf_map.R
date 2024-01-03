library(maps)
library(ggplot2)
library(ggsci)
library(scales)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)

my_pal = pal_lancet()(7)
show_col(my_pal)

# Map drawing code taken from:
# https://sarahleejane.github.io/learning/r/2014/09/20/plotting-beautiful-clear-maps-with-r.html
world_map <- map_data("world")
p <- ggplot() + coord_fixed() +
  xlab("") + ylab("")
#Add map to base plot
base_world_messy <- p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                     colour=my_pal[3], fill=my_pal[3])
#Strip the map down so it looks super clean (and beautiful!)
cleanup <- 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.line = element_line(colour = "white"), legend.position="none",
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank())

base_world <- base_world_messy + cleanup


setwd("/home/adminbrice/Desktop/research/PhD/main_PhD/analyses/plasmo_surfants/tmp_work")
pf_analysis_set <- read_tsv("pf6_analysis_set_fws95.tsv")
pf_analysis_set_counted <- group_by(pf_analysis_set, Lat,Long) %>% summarise(count=n())

p2 <- base_world_messy + geom_point(data=pf_analysis_set_counted, 
             aes(x=Long, y=Lat, size=count), colour="Deep Pink", 
             fill="Pink",pch=21, alpha=I(0.7)) + scale_size(range=c(3,9)) + labs(size="Number of samples")


ggsave("sample_num_counts_pf6_fws95.pdf",width=12,height=8)
