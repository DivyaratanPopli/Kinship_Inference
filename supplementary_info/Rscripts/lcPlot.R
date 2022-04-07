library(oompaBase)
library(TailRank)
library(ggplot2)
library(cowplot)
library(tidyr)
library(tidyverse)
library(plyr)

#lcres <- "/mnt/diversity/divyaratan_popli/100arc/inbreeding/laurits_writeup/Ch9_rem/neavar/noCh12_9mar_newmodel/lcMLkin/output.relate"
#likf <-"/mnt/diversity/divyaratan_popli/100arc/inbreeding/laurits_writeup/Ch9_rem/neavar/newSim_ibd/hmm_numba/filtered0/merged_relatable_allLikelihoods.csv"
#kinf <- "/mnt/diversity/divyaratan_popli/100arc/inbreeding/laurits_writeup/Ch9_rem/neavar/newSim_ibd/hmm_numba/filtered0/IBDstates_allLikelihoods.csv"


lcplot <- function(lcres,likf,kinf,outplot){

  lc<- read.csv(file = lcres, header = TRUE,sep = "\t")
  lik <-read.csv(file = likf, header = TRUE, sep = "\t")

  lc %>%
    separate(Ind1, c("name1", "bam1"), "\\.") %>% separate(Ind2, c("name2", "bam2"), "\\.") ->lc
  lc <-unite(lc, pair, c(name1, name2), remove=FALSE)

  drops <- c("name1","name2","bam1","bam2","nbSNP","pi_HAT")
  lc <-lc[ , !(names(lc) %in% drops)]

  lc <-merge(lik[c('pair','loglik_ratio')], lc, by = c("pair"))

  lc$r_hat <- lc$k1_hat/2 + lc$k2_hat
  lc$Relatedness <- 'Unrelated'
  lc[lc$pair=='Chagyrskaya07_Chagyrskaya17','Relatedness'] <- 'Parent-Child'
  lc[lc$pair=='Chagyrskaya13_Chagyrskaya19','Relatedness'] <- 'Identical'
  lc[lc$pair=='Chagyrskaya1141_Chagyrskaya13','Relatedness'] <- 'Identical'
  lc[lc$pair=='Chagyrskaya1141_Chagyrskaya19','Relatedness'] <- 'Identical'
  lc[lc$pair=='Chagyrskaya01_Chagyrskaya60','Relatedness'] <- 'Degree2'
  lc[lc$pair=='Chagyrskaya06_Chagyrskaya14','Relatedness'] <- 'Identical'

  lc=lc[lc$pair!='Chagyrskaya06_Chagyrskaya14',]

  colors1 <- c("Unrelated" = "grey80", "Degree2" = "darkorchid1", "Parent-Child" = "orangered3", "Identical" = "darkgoldenrod1")

  panel1<-ggplot(data=lc[which(lc$loglik_ratio >=1),], mapping=aes(x=k0_hat, y=r_hat, color=Relatedness)) +
  geom_point() + theme_bw() + labs(x=expression(k["0"]), y= "Coefficient of Relatedness") +
  scale_color_manual(values = colors1) + ylim(0, 1) + xlim(0, 1) +
  theme(legend.position = "none")



  ###kin

  kin<- read.csv(file = kinf, header = TRUE,sep = ",")
  newdf <-merge(lik[c('pair','loglik_ratio')], kin, by = c("pair"))
  #newdf$good <-"NA"
  #newdf[newdf$loglik_ratio>=1,"good"]=1

  newdf$Relatedness="Unrelated"
  newdf[newdf$pair=='Chagyrskaya07_Chagyrskaya17','Relatedness'] <- 'Parent-Child'
  newdf[newdf$pair=='Chagyrskaya13_Chagyrskaya19','Relatedness'] <- 'Identical'
  newdf[newdf$pair=='Chagyrskaya1141_Chagyrskaya13','Relatedness'] <- 'Identical'
  newdf[newdf$pair=='Chagyrskaya1141_Chagyrskaya19','Relatedness'] <- 'Identical'
  newdf[newdf$pair=='Chagyrskaya01_Chagyrskaya60','Relatedness'] <- 'Degree2'
  newdf[newdf$pair=='Chagyrskaya06_Chagyrskaya14','Relatedness'] <- 'Identical'

  newdf=newdf[newdf$pair!='Chagyrskaya06_Chagyrskaya14',]


  newdf$r_prop <- newdf$k1_prop/2 + newdf$k2_prop

  #newdf[newdf$loglik_ratio<1,]=NA

  colors <- c("Unrelated" = "grey80", "Degree2" = "darkorchid1", "Parent-Child" = "orangered3", "Identical" = "darkgoldenrod1")
  #alphas <-c("0"=0.01,"1"=1)
  panel2<-ggplot(data=newdf[which(newdf$loglik_ratio >=1),], mapping=aes(x=k0_prop, y=r_prop, color=Relatedness)) +
    geom_point() + theme_bw() + labs(x=expression(k["0"]), y= "Coefficient of Relatedness") +
    scale_color_manual(values = colors) +  ylim(0, 1) + xlim(0, 1) +
    #scale_alpha_manual(values = alphas) +
    theme(axis.title.y = element_blank())

  plot_grid(panel1,panel2, labels= c("A","B"),
            ncol=2, rel_widths = c(1,1.3), align = "hv", axis="bt")

  ggsave(outplot,
         width = 8, height = 5, dpi = 150, units = "in", device='png')
}

lcplot(lcres = snakemake@input[["lcres"]],likf = snakemake@input[["likf"]],kinf = snakemake@input[["kinf"]],outplot = snakemake@output[["outplot"]])




#########combine 2 plots

# lc %>%
#   separate(Ind1, c("name1", "bam1"), "\\.") %>% separate(Ind2, c("name2", "bam2"), "\\.") ->lc
# lc$pair <- paste(lc$name1,lc$name2)
#
# lc <-unite(lc, pair, c(name1, name2), remove=FALSE)
# drops <- c("name1","name2","bam1","bam2","nbSNP","pi_HAT")
# lc <-lc[ , !(names(lc) %in% drops)]
# lc$method <- "lcMLkin"
#
# dropk <- c("rel","k0_abs","k1_abs","k2_abs","IBD_len","IBD_num","X")
# kin <- kin[ , !(names(kin) %in% dropk)]
# kin$method <- "KIn"
# colnames(kin)<- c("pair","k0_hat","k1_hat","k2_hat","Relatedness","r_hat","method")
#
# alldata <- rbind(kin,lc)
#
# colors <- c("Unrelated" = "grey80", "Degree2" = "darkorchid1", "Parent-Child" = "orangered3", "Identical" = "darkgoldenrod1")
#
# ggplot(data=alldata, mapping=aes(x=k0_hat, y=r_hat, color=Relatedness, shape=method)) +
#   geom_point() + theme_bw() + labs(x=expression(k["0"]), y= "Coefficient of Relatedness") +
#   scale_color_manual(values = colors)
