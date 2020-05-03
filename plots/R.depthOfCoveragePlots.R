require(cowplot) # for doing multi-graph plots

fnPlotlog2Depth <- function(data) {
  stepSize <- 1000000
  
  plot <- ggplot() + 
    geom_point(data=data, map=aes(x=Pos, y=log2depth), size=.5) +
    theme_bw() +
    scale_x_discrete(breaks=(round(seq(min(data$Pos), max(data$Pos), by = stepSize),1)), labels=round((seq(min(data$Pos), max(data$Pos), by = stepSize) / stepSize),0)) +
    scale_y_continuous( limits = c(-1.5,1.5), expand = c(0,0) ) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  return(plot)
}


stocks <- c("chrX-01-M","chrX-02-M","chrX-03-M","chrX-04-M","chrX-05-M","chrX-06-M","chrX-07-M","chrX-08-M","chrX-09-M","chrX-10-M","chrX-11-M","chrX-12-M","chrX-13-M","chrX-14-M","chrX-15-M","chrX-16-M","chrX-17-M","chrX-18-M","chrX-19-M","chrX-20-M","chrX-21-M","chrX-22-M","chrX-23-M","chrX-24-M","chrX-25-M","chrX-26-M","chrX-27-M","chrX-28-M","chrX-29-M","chrX-30-M","chrX-31-M","chrX-32-M","chrX-33-M","chrX-34-M","chrX-35-M","chrX-36-M");
stocks <- c("mcm5-01","mcm5-02","mcm5-03","mcm5-04","mcm5-05","mcm5-06","mcm5-07","mcm5-08","mcm5-09","mcm5-10","mcm5-11","mcm5-12","mcm5-13","mcm5-14","mcm5-15","mcm5-16","mcm5-17","mcm5-18","mcm5-19","mcm5-20","mcm5-21","mcm5-22","mcm5-23","mcm5-24","mcm5-25","mcm5-26","mcm5-27","mcm5-28")


chromosomes <- c("chrX","chr2L","chr2R","chr3L","chr3R","chr4")

for (chr in chromosomes) {
  pdfFile <- paste("/Users/danny/projects/mcm5/log2depth/plot_depthOfCoverage.",chr,".pdf",sep="")
  print(pdfFile)
  pdf(pdfFile)
  
  plotList <- list()
  nameList <- list()
  count <- 0
  
  for (stock in stocks) {
    count = count + 1
    
    fullFile <- paste("/Users/danny/projects/mcm5/log2depth/",stock,"/",stock,".realigned.log2depth.tsv",sep="")
    #print(fullFile)
    data <- read.csv(fullFile,sep = "\t", header = T)
    
    chrPlot <- fnPlotlog2Depth(subset(data,Chr==chr))
    plotList[[count]] = chrPlot
    nameList[[count]] = stock
    
    if (count==4) {
      n1 <- paste(chr,nameList[[1]],sep=" ")
      n2 <- paste(chr,nameList[[2]],sep=" ")
      n3 <- paste(chr,nameList[[3]],sep=" ")
      n4 <- paste(chr,nameList[[4]],sep=" ")
      
      print(plot_grid(plotList[[1]], plotList[[2]], plotList[[3]], plotList[[4]], rel_heights = c(1,1,1,1), labels = c(n1,n2,n3,n4), nrow=4, ncol=1))
      count = 0
      #plotList <- list()
      #nameList <- list()
    }
  }
  
  n1 <- paste(chr,nameList[[1]],sep=" ")
  n2 <- paste(chr,nameList[[2]],sep=" ")
  n3 <- paste(chr,nameList[[3]],sep=" ")
  n4 <- paste(chr,nameList[[4]],sep=" ")

  print(plot_grid(plotList[[1]], plotList[[2]], plotList[[3]], plotList[[4]], rel_heights = c(1,1,1,1), labels = c(n1,n2,n3,n4), nrow=4, ncol=1))
  
  dev.off()
}





