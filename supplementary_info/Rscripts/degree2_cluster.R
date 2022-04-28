library(ggplot2)
library(tidyr)
library(tidyverse)
library(plyr)

hsib=c('8_9','8_10')
avu='9_12'
gr=c('1_12','2_12','3_12','4_12','10_13','11_13','12_14','12_16','5_14')


#fnames=c("/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/hmm_numba_corr_IBD/contam0/inbred0/run1/coverage4/filtered0/asc0/IBDstates_allLikelihoods_inphapProbs.csv.gz", "/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/hmm_numba_corr_IBD/contam0/inbred0/run2/coverage4/filtered0/asc0/IBDstates_allLikelihoods_inphapProbs.csv.gz")
#outplot="/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/test.png"
deg2_plot <- function(fnames, outplot){

  list_all <<- NULL
  for (r in 1:60){


    deg2=read.csv(file = fnames[r], header = TRUE,sep = ",")

    list_hsib=deg2[deg2$pair %in% hsib,]
    list_avu=deg2[deg2$pair %in% avu,]
    list_gr=deg2[deg2$pair %in% gr,]

    list_run <- rbind(list_avu,list_gr,list_hsib)
    list_run$rel='gr'
    list_run[list_run$pair %in% hsib,'rel']= 'hsib'
    list_run[list_run$pair %in% avu,'rel']= 'avu'

    drops <- c("pair","rel","IBD_len","IBD_num")
    list_run <-list_run[, (names(list_run) %in% drops)]

    list_all <- rbind(list_all, list_run)


    }


  colors <- c("avu" = "cadetblue3", "hsib" = "blueviolet", "gr" = "chocolate1")
  #alphas <-c("0"=0.01,"1"=1)

  ggplot(data=list_all, aes(x=IBD_len, y=IBD_num, color=rel)) +
    geom_point(alpha=0.5) + theme_bw() + labs(x="Total length of IBD fragments (Mb)", y= "Total number of IBD fragments") +
    scale_color_manual("Relatedness",values = colors, labels = c("avu" = "Avuncular", "hsib" = "Half-siblings", "gr" = "Grandparent-Grandchild"))

  ggsave(outplot,
         width = 8, height = 5, dpi = 150, units = "in", device='png')

}



deg2_plot(fnames=snakemake@input[["fnames"]],outplot=snakemake@output[["outplot"]])
