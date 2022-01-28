library(cowplot)
library(ggplot2)
library(readr)

#datadf <- "/mnt/diversity/divyaratan_popli/100arc/inbreeding/bronze_age_newSim/mergedwin_remHighdiv_fil0/pw_diff.csv"
#datatf <- "/mnt/diversity/divyaratan_popli/100arc/inbreeding/bronze_age_newSim/mergedwin_remHighdiv_fil0/pw_total.csv"
#p1f <- "/mnt/diversity/divyaratan_popli/100arc/inbreeding/bronze_age_newSim/hmm_parameters_fil0/p1_file"
#chrmf <- "/mnt/diversity/divyaratan_popli/100arc/inbreeding/bronze_age_newSim/mergedwin_fil0/merged_chrm.csv"
egplot1 <- function(p1f, datadf, datatf, chrmf, outimg){

  ind1=1907+1 #AITI62B_OTTM156
  ind2=6783+1 #UNTA5867_UNTA5868Sk1
  ind3=1029+1 #AITI40_AITI72
  ind4=244+1 #AITI2_AITI55
  p1 <- read_file(p1f)
  p1 <- as.double(gsub("[\r\n]", "", p1))
  p12=p1/2
  p34=p1*3/4

  datad=read.csv(file = datadf, header = FALSE)
  datat=read.csv(file = datatf, header = FALSE)
  chrm=read.csv(file = chrmf, header = FALSE)

  d1=datad[,ind1]/datat[,ind1]
  d2=datad[,ind2]/datat[,ind2]
  d3=datad[,ind3]/datat[,ind3]
  d4=datad[,ind4]/datat[,ind4]

  wins=seq(1,300,1)

  chrm1=chrm%%2
  df1 <- data.frame(d1,wins,chrm1)
  df2 <- data.frame(d2,wins,chrm1)
  df3 <- data.frame(d3,wins,chrm1)
  df4 <- data.frame(d4,wins,chrm1)

  plot1 <- ggplot(data=df1, aes(x=wins,y=d1)) +
    geom_rect(aes(
      xmin = wins-1,
      xmax = wins,
      ymin = -Inf,
      ymax = +Inf,
      fill = factor(V1),
      #alpha=0.8
    )) +
    scale_fill_manual(values = c("grey90", "white")) +
    theme_void() +
    geom_line() +
    geom_segment(aes(x=1,xend=300,y=p1,yend=p1),linetype='dotted') + geom_segment(aes(x=1,xend=300,y=p12,yend=p12),linetype='dotted') +
    geom_segment(aes(x=1,xend=300,y=p34,yend=p34),linetype='dotted') +theme_bw() +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text = element_text(size = 8))+
    labs(y=expression(D[w]/N[w])) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) + coord_cartesian(ylim = c(0, 0.65))

  plot1 <-plot1+ annotate(geom="text", x=-2, y=0.265, label="p0",
                          color="black", size=3) +
    annotate(geom="text", x=-2, y=0.20, label="p1",
             color="black", size=3) +
    annotate(geom="text", x=-2, y=0.14, label="p2",
             color="black", size=3)

  plot2 <- ggplot(data=df2, aes(x=wins,y=d2)) +
    geom_rect(aes(
      xmin = wins-1,
      xmax = wins,
      ymin = -Inf,
      ymax = +Inf,
      fill = factor(V1),
      #alpha=0.8
    )) +
    scale_fill_manual(values = c("grey90", "white")) +
    theme_void() +
    geom_line() +
    geom_segment(aes(x=1,xend=300,y=p1,yend=p1),linetype='dotted') + geom_segment(aes(x=1,xend=300,y=p12,yend=p12),linetype='dotted') +
    geom_segment(aes(x=1,xend=300,y=p34,yend=p34),linetype='dotted') +theme_bw() +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text = element_text(size = 8))+
    labs(y=expression(D[w]/N[w])) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) + coord_cartesian(ylim = c(0, 0.65))

  plot2 <-plot2+ annotate(geom="text", x=-2, y=0.265, label="p0",
                          color="black", size=3) +
    annotate(geom="text", x=-2, y=0.20, label="p1",
             color="black", size=3) +
    annotate(geom="text", x=-2, y=0.14, label="p2",
             color="black", size=3)


  plot3 <- ggplot(data=df3, aes(x=wins,y=d3)) +
    geom_rect(aes(
      xmin = wins-1,
      xmax = wins,
      ymin = -Inf,
      ymax = +Inf,
      fill = factor(V1),
      #alpha=0.8
    )) +
    scale_fill_manual(values = c("grey90", "white")) +
    theme_void() +
    geom_line() +
    geom_segment(aes(x=1,xend=300,y=p1,yend=p1),linetype='dotted') + geom_segment(aes(x=1,xend=300,y=p12,yend=p12),linetype='dotted') +
    geom_segment(aes(x=1,xend=300,y=p34,yend=p34),linetype='dotted') +theme_bw() +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text = element_text(size = 8))+
    labs(y=expression(D[w]/N[w])) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) + coord_cartesian(ylim = c(0, 0.65))

  plot3 <-plot3+ annotate(geom="text", x=-2, y=0.265, label="p0",
                          color="black", size=3) +
    annotate(geom="text", x=-2, y=0.20, label="p1",
             color="black", size=3) +
    annotate(geom="text", x=-2, y=0.14, label="p2",
             color="black", size=3)


  plot4 <- ggplot(data=df4, aes(x=wins,y=d4)) +
    geom_rect(aes(
      xmin = wins-1,
      xmax = wins,
      ymin = -Inf,
      ymax = +Inf,
      fill = factor(V1),
      #alpha=0.8
    )) +
    scale_fill_manual(values = c("grey90", "white")) +
    theme_void() +
    geom_line() +
    geom_segment(aes(x=1,xend=300,y=p1,yend=p1),linetype='dotted') + geom_segment(aes(x=1,xend=300,y=p12,yend=p12),linetype='dotted') +
    geom_segment(aes(x=1,xend=300,y=p34,yend=p34),linetype='dotted') +theme_bw() +
    theme(legend.position = "none", axis.ticks.x=element_blank(), strip.text = element_text(size = 8))+
    labs(y=expression(D[w]/N[w]), x = "Windows along the genome") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) + coord_cartesian(ylim = c(0, 0.65))

  plot4 <-plot4+ annotate(geom="text", x=-2, y=0.265, label="p0",
                           color="black", size=3) +
    annotate(geom="text", x=-2, y=0.20, label="p1",
             color="black", size=3) +
    annotate(geom="text", x=-2, y=0.14, label="p2",
             color="black", size=3)


  panel<-plot_grid(plot1,plot2,plot3,plot4, labels= c("A","B","C","D"), ncol=1, rel_widths = c(1,1,1,1),rel_heights = c(1,1,1,1.3), align="v")
  ggsave(outimg,
         width = 8, height = 5, dpi = 150, units = "in", device='png')
}


egplot1(p1f=snakemake@input[["p1f"]], datadf=snakemake@input[["datadf"]], datatf=snakemake@input[["datatf"]], chrmf=snakemake@input[["chrmf"]], outimg=snakemake@output[["outimg"]])
