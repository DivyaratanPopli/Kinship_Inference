library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
#fname='/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/plot_ibd.csv'
#p1f='/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/contam0/inbred0/run1/coverage4/asc0/inputMode_hapProbs_fil0/hmm_parameters/p1_file.gz'

ibd_plot <- function(inf,p1f,outf){

  #p1 <- read_file(p1f)
  #p1 <- as.double(gsub("[\r\n]", "", p1))
  p1=0.021
  p12=p1/2
  p34=p1*3/4

  data=read.csv(file = inf, header = TRUE)
  colnames(data) <- c("X", "Wins","Data","True","Unrelated","Degree3","Degree2","Parent-Child","Siblings","Identical")

  require(tidyverse)

  d2 = rbind(data, c(0, NA, NA, rep(0.5,7)), c(0,NA, NA, rep(1,7))) %>% as_tibble

  d3 <-d2 %>% pivot_longer(Data:Identical)
  d3$name_f = factor(d3$name, levels=c('Data','True','Identical','Parent-Child','Siblings','Degree2','Degree3','Unrelated'))

  chrmlist=list()
  for (val in 1:22)
  {
    if (val%%2==1){
      lchr=rep(0,10*8)
    }
    else if (val%%2==0){
      lchr=rep(1,10*8)
    }
    chrmlist <- append(chrmlist, lchr)

  }
  nalist=rep(NA,16)
  chrmlist <-append(chrmlist,nalist)

  d3$chrm =as.integer(chrmlist)


  d_ibd<-subset(d3, name!="Data")
  d_prop<-subset(d3, name=="Data")

  dat_text <- data.frame(
    label = as.factor(c("", "LL=-1478.24", "LL=-1259.44", "LL=-1213.61", "LL=-1236.09", "LL=-1261.19", "LL=-1281.77")),
    name_f = as.factor(c("True", "Identical", "Parent-Child", "Siblings", "Degree2", "Degree3", "Unrelated")),
    x     = c(5,5,5,5,5,5,5),
    y     = c(0.85,0.85,0.85,0.85,0.85,0.85,0.85)
  )

  datcol <-data.frame(
    from=c(1,20,40),
    to=c(10,30,50)
  )

  plot1i<-ggplot(data=d_ibd,aes(x=Wins, y=value)) +
    geom_rect(aes(
      xmin = Wins-1,
      xmax = Wins,
      ymin = -Inf,
      ymax = +Inf,
      fill = factor(chrm),
      #alpha=0.8
    )) +
    scale_fill_manual(values = c("grey90", "white")) +
    theme_void() +
    geom_line() + facet_grid(name_f~., scales='free') +
    scale_y_continuous(breaks=c(0.5,0.75,1.0), labels = c("2","1","0")) +theme_bw() + labs(y="IBD state", x = "Windows along the genome") +
    theme(strip.text = element_text(size = 10), legend.position = "none") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
    #scale_alpha_manual(values = c(0.1, 0.1))
    #geom_rect(aes(xmin = from, xmax = to, ymin = -Inf, ymax = Inf), alpha = 0.4,fill = factor(chrm))

  plot1 <- plot1i + geom_text(
    data    = dat_text,
    mapping = aes(x = x, y = y, label = label)
  )

  plot2i <-ggplot(data=d_prop,aes(x=Wins, y=value)) +
    geom_rect(aes(
      xmin = Wins-1,
      xmax = Wins,
      ymin = -Inf,
      ymax = +Inf,
      fill = factor(chrm),
      #alpha=0.8
    )) +
    scale_fill_manual(values = c("grey90", "white")) +
    theme_void() +
    geom_line() + facet_grid(name_f~., scales='free') +
    geom_segment(aes(x=0,xend=220,y=p1,yend=p1),linetype='dotted') + geom_segment(aes(x=0,xend=220,y=p12,yend=p12),linetype='dotted') +
    geom_segment(aes(x=0,xend=220,y=p34,yend=p34),linetype='dotted') +theme_bw()+
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text = element_text(size = 10))+
    labs(y=expression(D[w]/N[w])) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))

  plot2 <-plot2i+ annotate(geom="text", x=0, y=0.024, label="p0",
                           color="black", size=3) +
                  annotate(geom="text", x=0, y=0.018, label="p1",
                           color="black", size=3) +
                  annotate(geom="text", x=0, y=0.013, label="p2",
                           color="black", size=3)

  panel<-plot_grid(plot2,plot1, ncol=1, rel_widths = c(1,1),rel_heights = c(1,7), align="v")


  ggsave(outf,
         width = 8, height = 8, dpi = 150, units = "in", device='pdf')
}

ibd_plot(inf=snakemake@input[["infile_ibd"]], p1f=snakemake@input[["p1_file"]], outf=snakemake@output[["output_ibd"]])
