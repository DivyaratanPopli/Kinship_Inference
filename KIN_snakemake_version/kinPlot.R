library(ggplot2)
library(tidyr)
library(tidyverse)
library(plyr)


#overlap_table <-"/mnt/diversity/divyaratan_popli/100arc/inbreeding/laurits_writeup/Ch9_rem/neavar/newSim_probCorr/overlap_fil0"
#relatable <- "/mnt/diversity/divyaratan_popli/100arc/inbreeding/laurits_writeup/Ch9_rem/neavar/newSim_probCorr/hmm_numba/filtered0/merged_relatable_allLikelihoods.csv"
#out_rel <- "aaaaaa.txt"
#out_over <-"aaaaaaaaaa.txt"
#tarf="/mnt/diversity/alba_bm/Belgian/Nuclear_Belgian/KIN_Indiv_contam_newSim/targets.txt

bubbleplot_kin <- function(overlap_table, relatable, out_rel, out_over, tarf) {

  allnames <- read.csv(tarf, header = FALSE, sep = ",")
  over=read.csv(file=overlap_table, header=TRUE, row.names=1, sep=',')
  rel=read.csv(file=relatable, row.names = 1, header=TRUE, sep='\t')
  info=merge(over, rel, by=c("pair"),all.x=TRUE)

  info$pair <- gsub("[()]", "", info$pair)
  matr=separate(info, col = pair, into = c("lib1","lib2"), sep = "_")
  matr1 <- matr %>%
    select(lib2,lib1,overlap,relatedness,second_guess,loglik_ratio,withinDeg_second_guess,withinDeg_ll)
  colnames(matr1) <- colnames(matr)
  matr1 <- rbind(matr, matr1)
  matr1$siz <- cut(matr1$loglik_ratio,
                   breaks = c(0 , 2 , 10 , Inf),
                   labels = c("<=2" , "<=10" , ">10" ))



  matr1$relatedness<- revalue(matr1$relatedness, c("id"="Identical","sib"="Siblings","pc"="Parent-Child","sec"="Second Degree","un"="Unrelated", "deg3"="Third Degree"))
  matr1$relatedness<-factor(matr1$relatedness, levels=c("Identical", "Siblings", "Parent-Child","Second Degree","Third Degree","Unrelated"))
  matr1$overlap <- cut(matr1$overlap,breaks = c(-Inf,0,10,20,50,100,500,1000,5000,10000,Inf),right = FALSE,
                       label=c("0","0-10","10-20","21-50","50-100","101-500","501-1000","1000-5000","5000-10000",">10000"))
  matr1$loglik_ratio <- round(matr1$loglik_ratio,1)
  #data[, c(1, 2)] <- as.data.frame(t(apply(data[, c(1, 2)], 1, sort)))data$V1 <- factor(data$V1, levels=c(“Fonds”,“Spy8",“Spy1”,“Spy94a”,“GoyetQ57-3”,“GoyetQ57-2",“GoyetQ57-1”,“GoyetQ305-4",“GoyetQ56-1”,“GoyetQ56-1-lowCov”,“Goyet374a-1”,“Goyet305-7",“GoyetQ55-4”,“Goyet1424-3D”,“GoyetC5-1”,“GoyetQ119-2",“GoyetQ376-25”))
  matr1$lib2 <- factor(matr1$lib2, levels=allnames$V1)#data$SP_V1 <- ifelse(data$V1 == “A29253”, “red”, “black”)
  matr1$lib1 <- factor(matr1$lib1, levels=allnames$V1)#data$SP_V1 <- ifelse(data$V1 == “A29253”, “red”, “black”)

  colors <- c("seagreen4","steelblue2","mediumblue","steelblue3","lightpink3","thistle3","gray45","lightsalmon3","darkred","firebrick","tan","black","gray70","lightcyan3","darkorchid","pink3","pink4")#data[, c(1, 2)] <- as.data.frame(t(apply(data[, c(1, 2)], 1, sort)))#Plot all the infered relatedness

  ##########trying heatmap.2
  #i1=unique(matr1$lib1)
  #i2=unique(matr1$lib2)

  #xx <- data.frame(row.names = i1)
  #colnames(xx) <- i2

  #colnames(xx) <- i1
  #row.names(xx) <- i2
  #filter(matr1, lib1 == i1[j] & lib2 == i2[k])$relatedness
  #write.table(x=matr1,file = "/mnt/diversity/divyaratan_popli/100arc/inbreeding/Chagyrskaya_remapped24sept/noCh12/plot_inbred/input_heatmap.csv", sep = ',')
  #abc=read.csv(file="/mnt/diversity/divyaratan_popli/100arc/inbreeding/Chagyrskaya_remapped24sept/noCh12/plot_inbred/pyinput_heatmap.csv", row.names=1,header=TRUE, sep=',',)
  #m1 <- as.matrix(abc)
  #heatmap(m1)

  gfig_rel=ggplot(matr1, aes(lib1, lib2), col= "white") +
    geom_tile(aes(fill=relatedness), col="white")  +
    scale_fill_manual(values = c("Identical"="khaki1","Parent-Child"="orangered3", "Siblings"="orange2","Unrelated"="gray80","Second Degree"="darkorchid1","Third Degree"="aquamarine2")) +

    labs(title = "relatedness", x = "lib_ID", y = "lib_ID")+
    theme(
      panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white",colour = "white", size = 0.5, linetype = "solid"))+
    theme(axis.text.x = element_text(angle = 90,hjust = 1, colour = colors, size = 14)) +
    theme(axis.text.y = element_text(angle = 0,hjust = 1, colour = colors, size = 14)) +
    geom_text(aes(lib1, lib2, label = matr1$loglik_ratio), color = "black", size = 5)#Plot overlap +
  theme(legend.text = element_text(size=14),
        legend.title=element_text(size=15))

  ggsave(out_rel,height=12,width=14, plot = gfig_rel)

  gfig_over=ggplot(matr1, aes(lib1, lib2), col= "white") +
    geom_tile(aes(fill=overlap), col="white")  +
    scale_fill_manual(values = c("gray90","darkslateblue","dodgerblue3","mediumseagreen","palegreen","yellowgreen","skyblue4","khaki1","firebrick")) +
    labs(title = "Overlap filter", x = "lib_ID", y = "lib_ID")+
    theme(
      panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white",colour = "white", size = 0.5, linetype = "solid"))+
    theme(axis.text.x = element_text(angle = 90,hjust = 1, colour = colors, size = 14)) +
    theme(axis.text.y = element_text(angle = 0,hjust = 1, colour = colors, size = 14)) +
    theme(legend.text = element_text(size=14),
          legend.title=element_text(size=15))


  ggsave(out_over,height=12,width=14, plot = gfig_over)
}


bubbleplot_kin(overlap_table=snakemake@input[["over"]], relatable = snakemake@input[["rel"]], out_rel=snakemake@output[["outfile_rel"]], out_over=snakemake@output[["outfile_over"]], tarf=snakemake@input[["tarf"]])
