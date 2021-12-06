library(ggplot2)
library(cowplot)

make_panel1 <- function(data,acrval,leg){
  plot1<-ggplot(data=data[data$acr==acrval,], aes(x=Wins, y=Prop))
  plot1<-plot1+geom_line(aes(colour=Coverage),alpha=0.4, size=0.8)
  plot1<-plot1 + theme_bw()+ theme(legend.position = "none", axis.title.x = element_blank())+labs(y= "Heterozygosity", x = "Windows along genome")

  return(plot1)
}
make_panel2 <- function(data,acrval,leg){
  plot2<-ggplot(data=data[data$acr==acrval,])
  plot2<-plot2+geom_line(aes(x=Wins, y=True_ROH, linetype="Truth"))
  #plot2
  plot2<-plot2+geom_line(aes(x=Wins, y=Model_pred, linetype="Inferred", colour=Coverage), alpha=0.4, size=0.8)

  plot2<-plot2+scale_linetype_manual(values=c("Truth"="dashed", "Inferred"="solid"))
  plot2 <- plot2 + scale_color_manual(values=c("0.2"="Green", "0.5"="Blue","4"="Red"))
  if (leg==1){
    plot2<-plot2 + theme_bw() + theme(legend.position = "bottom")+labs(y= "1 - P(ROH)", x = "Windows along genome")
    plot2<- plot2+ guides(linetype = guide_legend(nrow = 1, byrow = TRUE))+ guides(colour = guide_legend(nrow = 1, byrow = TRUE))
  } else if (leg==0){
    plot2<-plot2 + theme_bw() + theme(legend.position = "none")+labs(y= "1 - P(ROH)", x = "Windows along genome")
  }
  #panel<-plot_grid(plot1,plot2, labels= c("A","B"), ncol=1, rel_widths = c(1,1))
  #panel
  return(plot2)

}


ROH_plot <- function(infile,outfile){


  data=read.csv(file = infile, header = TRUE)
  data$Coverage <- as.factor(data$Coverage)


  panel11 <- make_panel1(data=data, acrval='0_0_0',leg=0)
  panel12 <- make_panel1(data=data, acrval='0_0_1',leg=0)
  panel21 <- make_panel2(data=data, acrval='0_0_0',leg=0)
  panel22 <- make_panel2(data=data, acrval='0_0_1',leg=0)

  panel31 <- make_panel1(data=data, acrval='2_0_1',leg=0)
  panel32 <- make_panel1(data=data, acrval='0_1_1',leg=0)
  panel41 <- make_panel2(data=data, acrval='2_0_1',leg=1)
  panel42 <- make_panel2(data=data, acrval='0_1_1',leg=0)

  panel<-plot_grid(panel11,panel12,panel21,panel22,panel31,panel32,panel41,panel42, labels= c("A","B"," "," ","C","D"," "," "), ncol=2, rel_widths = c(1,1,1,1,1,1,1,1), rel_heights =  c(1,1,1,1,1,1,1,1), align="v",axis="lr")
  panel


ggsave(outfile,
       width = 8, height = 5, dpi = 150, units = "in", device='png')
}

ROH_plot(infile=snakemake@input[["infile_roh"]], outfile=snakemake@output[["outfile_roh"]])
