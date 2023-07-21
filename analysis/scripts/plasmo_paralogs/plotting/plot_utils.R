library(ggplot2)

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
