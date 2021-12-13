library(oompaBase)
library(TailRank)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(readr)

#dfile <-'/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz_highDiversity/contam0/inbred0/run1/coverage4/asc0/inputMode_hapProbs_fil0/merged_wind.csv.gz'
#tfile <-'/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz_highDiversity/contam0/inbred0/run1/coverage4/asc0/inputMode_hapProbs_fil0/merged_wint.csv.gz'
#p1file<-'/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz_highDiversity/contam0/inbred0/run1/coverage4/asc0/inputMode_hapProbs_fil0/hmm_parameters/p1_file.gz'
#deltafile0<-'/mnt/diversity/divyaratan_popli/100arc/inbreeding/model_git/supplementary_info/datafiles/contam0_inbred0_run1_coverage4_asc0_inputMode_hapProbs_fil0_delta_0_1_un.csv.gz'
#deltafile1<-'/mnt/diversity/divyaratan_popli/100arc/inbreeding/model_git/supplementary_info/datafiles/contam0_inbred0_run1_coverage4_asc0_inputMode_hapProbs_fil0_delta_0_8_pc.csv.gz'
#deltafile2<-'/mnt/diversity/divyaratan_popli/100arc/inbreeding/model_git/supplementary_info/datafiles/contam0_inbred0_run1_coverage4_asc0_inputMode_hapProbs_fil0_delta_0_15_id.csv.gz'
#mergePanels(dfile, tfile, p1file, deltafile0, deltafile1, deltafile2, '/mnt/diversity/divyaratan_popli/100arc/inbreeding/model_git/supplementary_info/test.png')

makepanelVar <- function(dfile, tfile, p1file, ind, pratio, deltaind,deltafile){

  diff=read.csv(file = dfile, header = FALSE)[,ind]
  total=read.csv(file = tfile, header = FALSE)[,ind]

  deltalist=read.csv(file = deltafile, header = FALSE)

  prop=diff/total
  propdf=as.data.frame(prop)

  p <- read_file(p1file)
  p <- as.double(gsub("[\r\n]", "", p))
  p1=p*pratio

  alpha=deltalist[deltaind,]
  beta=deltalist[deltaind+1,]

  diff_binom=rbinom(length(prop),total,p1)
  prop_binom=diff_binom/total
  propdf$binom=prop_binom

  diff_betabinom=rbb(length(prop),total,alpha,beta)
  prop_betabinom=diff_betabinom/total
  propdf$betabinom=prop_betabinom

  data <- propdf %>%
    pivot_longer(c(prop,binom,betabinom), names_to = "dist", values_to = "proportion")#Plot
  data$dist <- factor(data$dist, levels= c('prop','binom','betabinom'))

  return(data)

}


mergePanels <- function(dfile, tfile, p1file, deltafile0, deltafile1, deltafile2, outf){



  Var0=makepanelVar(dfile=dfile, tfile=tfile, p1file=p1file, ind=1, pratio=1, deltaind=3,deltafile=deltafile0)
  Var1=makepanelVar(dfile=dfile, tfile=tfile, p1file=p1file, ind=8, pratio=3/4, deltaind=1,deltafile=deltafile1)
  Var2=makepanelVar(dfile=dfile, tfile=tfile, p1file=p1file, ind=15, pratio=1/2, deltaind=5,deltafile=deltafile2)

  colors <- c("prop" = "grey", "betabinom" = "blue", "binom" = "red")
  fills <- c("prop" = "grey", "betabinom" = "blue", "binom" = "red")
  alphas <- c("prop" = 0.8, "betabinom" = 0.1, "binom" = 0.1)
  colnames <- c("prop" = "Heterozygosity", "betabinom" = "Betabinomial", "binom" = "Binomial")

  panel0binom <- ggplot(data=Var0, aes(x=proportion, color=dist, fill=dist, alpha=dist)) +
    geom_histogram(binwidth = 0.002,position="identity") +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = fills) +
    scale_alpha_manual(values = alphas) +
    theme_bw() + theme(axis.title.x = element_blank(), legend.position = "none") +
    labs(y= "Count")

  panel1binom <- ggplot(data=Var1, aes(x=proportion, color=dist, fill=dist, alpha=dist)) +
    geom_histogram(binwidth = 0.002,position="identity") +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = fills) +
    scale_alpha_manual(values = alphas) +
    theme_bw() + theme(axis.title.y = element_blank(), legend.position = "none") +
    labs(x= "Proportion of differences")

  panel2binom <- ggplot(data=Var2, aes(x=proportion, color=dist, fill=dist, alpha=dist)) +
    geom_histogram(binwidth = 0.002,position="identity") +
    scale_color_manual(values = colors, labels=colnames, name='Distribution') +
    scale_fill_manual(values = fills, labels=colnames, name='Distribution') +
    scale_alpha_manual(values = alphas, labels=colnames, name='Distribution') +
    theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())


  plot_grid(panel0binom,panel1binom,panel2binom, labels= c("A","B","C"),
            ncol=3, rel_widths = c(1,1,1), rel_heights =  c(1,1,1), align = "h", axis="b")
  ggsave(outf,
         width = 8, height = 5, dpi = 150, units = "in", device='png')
}

mergePanels(dfile=snakemake@input[["dfile"]], tfile=snakemake@input[["tfile"]], p1file=snakemake@input[["p1file"]],
              deltafile0=snakemake@input[["del0"]], deltafile1=snakemake@input[["del1"]],
              deltafile2=snakemake@input[["del2"]], outf=snakemake@output[["binomplot"]])
