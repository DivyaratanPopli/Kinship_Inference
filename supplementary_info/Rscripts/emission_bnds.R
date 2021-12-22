library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

#B1file <- "/mnt/diversity/divyaratan_popli/100arc/inbreeding/model_git/supplementary_info/contam0_inbred1_run57_coverage0.2_asc0_inputMode_hapProbs_fil0_B1file_0_15_relid"
#B2file <- "/mnt/diversity/divyaratan_popli/100arc/inbreeding/model_git/supplementary_info/contam0_inbred1_run57_coverage0.2_asc0_inputMode_hapProbs_fil0_B2file_0_15_relid"
#xnew1 <- "/mnt/diversity/divyaratan_popli/100arc/inbreeding/model_git/supplementary_info/contam0_inbred1_run57_coverage0.2_asc0_inputMode_hapProbs_fil0_xnew1file_0_15_relid"
#xnew2 <- "/mnt/diversity/divyaratan_popli/100arc/inbreeding/model_git/supplementary_info/contam0_inbred1_run57_coverage0.2_asc0_inputMode_hapProbs_fil0_xnew2file_0_15_relid"
#dfile <- "/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz_highDiversity/contam0/inbred1/run57/coverage0.2/asc0/inputMode_hapProbs_fil0/remHighdiv_diff.csv.gz"
#tfile <- "/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz_highDiversity/contam0/inbred1/run57/coverage0.2/asc0/inputMode_hapProbs_fil0/remHighdiv_total.csv.gz"
#p1f <- "/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz_highDiversity/contam0/inbred1/run57/coverage0.2/asc0/inputMode_hapProbs_fil0/hmm_parameters/p1_file.gz"

plot_b01 <- function(B1file,B2file,xnew1,xnew2,dfile,tfile,p1f,outf){

  ind=15
  df=read.csv(file = dfile, header = FALSE,sep = ",")[,ind]
  tf=read.csv(file = tfile, header = FALSE,sep = ",")[,ind]
  B1=read.csv(file = B1file, header = FALSE,sep = ",")
  B2=read.csv(file = B2file, header = FALSE,sep = ",")
  x1=read.csv(file = xnew1, header = FALSE,sep = ",")
  x2=read.csv(file = xnew2, header = FALSE,sep = ",")

  B1[B1$V1==-9,'V1']=NA
  B2[B2$V1==-9,'V1']=NA

  p1 <- read_file(p1f)
  p1 <- as.double(gsub("[\r\n]", "", p1))
  p12=p1/2
  p34=p1*3/4

  prop=df/tf
  df_prop <- do.call(rbind, Map(data.frame, value=prop, Wins=c(1:length(prop)), Bname='Data'))

  df_B1 <- do.call(rbind, Map(data.frame, value=B1$V1, Wins=c(1:nrow(B1)), Bname='B1'))
  df_B2 <- do.call(rbind, Map(data.frame, value=B2$V1, Wins=c(1:nrow(B2)), Bname='B2'))

  B_all <- rbind(df_B1,df_B2)

  Blabel <- c('B1'='Emission w/o constraint', 'B2'='Emission with constraint')


  plot1<-ggplot(data=df_prop,aes(x=Wins, y=value)) + geom_line() + facet_grid(Bname~., scales='free') +
    geom_segment(aes(x=0,xend=220,y=p1,yend=p1),linetype='dotted') + geom_segment(aes(x=0,xend=220,y=p12,yend=p12),linetype='dotted') +
    geom_segment(aes(x=0,xend=220,y=p34,yend=p34),linetype='dotted') +theme_bw()+
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text = element_text(size = 6))+
    labs(y=expression(D[w]/N[w])) + theme(strip.text = element_text(size = 7))


  plot2<-ggplot(data=B_all,aes(x=Wins, y=value)) + geom_line() + facet_grid(Bname~., scales='free', labeller = as_labeller(Blabel)) +
    theme_bw() + labs(y=expression("Emission probabilities for Z"[w]*"=2"), x = "Windows along the genome (w)") +
    theme(strip.text = element_text(size = 7))

  panel <- plot_grid(plot1,plot2, ncol=1, rel_widths = c(1,1),rel_heights = c(1,2), align="v",axis='b')
  ggsave(outf,
         width = 8, height = 5, dpi = 150, units = "in", device='png')

}

plot_b01(dfile=snakemake@input[["dfile"]], tfile=snakemake@input[["tfile"]], p1f=snakemake@input[["p1f"]],
            xnew1=snakemake@input[["xnew1"]], xnew2=snakemake@input[["xnew2"]],
            B1file=snakemake@input[["B1file"]], B2file=snakemake@input[["B2file"]],
            outf=snakemake@output[["outf"]])
