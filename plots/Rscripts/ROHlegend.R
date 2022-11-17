library(ggplot2)
library(cowplot)



make_panel1 <- function(data,acrval,leg){
  plot1<-ggplot(data=data[data$acr==acrval,], aes(x=Wins, y=Prop)) +
    geom_rect(aes(
      xmin = Wins-1,
      xmax = Wins,
      ymin = -Inf,
      ymax = +Inf,
      fill = factor(chrm),
      #alpha=0.8
    )) +
    scale_fill_manual(values = c("grey90", "white")) +
    theme_void()
  plot1<-plot1+geom_line(aes(colour=Coverage),alpha=0.7, size=0.8)
  plot1<-plot1 + theme_bw()+ theme(legend.position = "none", axis.title.x = element_blank())+labs(y= "Heterozygosity", x = "Windows along genome") +
    scale_color_manual(values=c("0.1"="#999933", "0.2"="#117733","4"="#661100")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  return(plot1)
}

make_panel2 <- function(data,acrval,leg){
  plot2<-ggplot(data=data[data$acr==acrval,]) +
    geom_rect(aes(
      xmin = Wins-1,
      xmax = Wins,
      ymin = -Inf,
      ymax = +Inf,
      fill = factor(chrm),
      #alpha=0.8
    )) +
    scale_fill_manual(values = c("grey90", "white")) +
    theme_void()
  plot2<-plot2+geom_line(aes(x=Wins, y=True_ROH, linetype="Truth"))
  #plot2
  plot2<-plot2+geom_line(aes(x=Wins, y=Model_pred, linetype="Inferred", colour=Coverage), alpha=0.7, size=0.8)

  plot2<-plot2+scale_linetype_manual(values=c("Truth"="dashed", "Inferred"="solid"))
  plot2 <- plot2 + scale_color_manual(values=c("0.1"="#999933", "0.2"="#117733","4"="#661100"))
  if (leg==1){
    plot2<-plot2 + theme_bw() + theme(legend.position = "bottom")+labs(y= "1 - P(ROH)", x = "Windows along genome")
    plot2<- plot2+ guides(linetype = guide_legend(nrow = 1, byrow = TRUE))+ guides(colour = guide_legend(nrow = 1, byrow = TRUE))
  } else if (leg==0){
    plot2<-plot2 + theme_bw() + theme(legend.position = "none")+labs(y= "1 - P(ROH)", x = "Windows along genome")
  }
  plot2<-plot2 +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  #panel<-plot_grid(plot1,plot2, labels= c("A","B"), ncol=1, rel_widths = c(1,1))
  #panel
  return(plot2)

}


ROH_plot <- function(infile,outfile){


  data=read.csv(file = infile, header = TRUE)
  data$Coverage <- as.factor(data$Coverage)

  chrmlist=list()
  for (val in 1:264)
  {
    if (val%%2==1){
      lchr=rep(0,10)
    }
    else if (val%%2==0){
      lchr=rep(1,10)
    }

    chrmlist <- append(chrmlist, lchr)
  }


  data$chrm =as.integer(chrmlist)

  panel11 <- make_panel1(data=data, acrval='0_0_0',leg=0)
  panel12 <- make_panel1(data=data, acrval='0_0_1',leg=0)
  panel21 <- make_panel2(data=data, acrval='0_0_0',leg=0)
  panel22 <- make_panel2(data=data, acrval='0_0_1',leg=0)

  panel31 <- make_panel1(data=data, acrval='2_0_1',leg=0)
  panel32 <- make_panel1(data=data, acrval='0_1_1',leg=0)
  panel41 <- make_panel2(data=data, acrval='2_0_1',leg=1)
  #panel42 <- make_panel2(data=data, acrval='0_1_1',leg=0)

  #panel<-plot_grid(panel11,panel12,panel21,panel22,panel31,panel32,panel41,panel42, labels= c("A","B"," "," ","C","D"," "," "), ncol=2, rel_widths = c(1,1,1,1,1,1,1,1), rel_heights =  c(1,1,1,1,1,1,1,1), align="v",axis="lr", label_size=10, label_x=0, label_y=1, hjust = 0, vjust = 3)
  #panel


  ggsave(outfile,
         width = 8, height = 5, dpi = 150, units = "in", device='png')
}

ROH_plot(infile=snakemake@input[["infile_roh"]], outfile=snakemake@output[["outfile_roh"]])

#infile='/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/plot_ROH.csv'
#outfile='/mnt/diversity/divyaratan_popli/100arc/inbreeding/model_git/plots/test.png'
