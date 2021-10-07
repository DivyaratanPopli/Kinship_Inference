library(ggplot2)
library(tidyr)
library(tidyverse)
library(plyr)

bubbleplot_read <- function(overlap_table, relatable, out_rel, out_over) {


  print(overlap_table)
  print(relatable)
  print(out_rel)
  print(out_over)

  over=read.csv(file=overlap_table, header=TRUE, row.names=1, sep=',')
  row.names(over) <-seq(1, nrow(over), by=1)

  rel=read.csv(file=relatable, header=TRUE, sep='\t')

  colnames(rel) <- c("pair", "relatedness","zu","zl")
  rel$pair <- gsub('^.|.$', '', rel$pair)
  rel$pair <- gsub('__', '_', rel$pair)
  info=merge(over, rel, by=c("pair"),all.x=TRUE)

  info$pair <- gsub("[()]", "", info$pair)
  matr=separate(info, col = pair, into = c("lib1","lib2"), sep = "_")

  matr$dist_thresh <- matr$zu
  for (i in 1:nrow(matr)) {
    if (is.na(matr$dist_thresh[i])) {
      matr$dist_thresh[i] <- abs(matr$zl[i])
    }
  }


  matr1 <- matr %>%
    select(lib2,lib1,overlap,relatedness,zu,zl,dist_thresh)
  colnames(matr1) <- colnames(matr)
  matr1 <- rbind(matr, matr1)

  matr1$relatedness<- revalue(matr1$relatedness, c("IdenticalTwins/SameIndividual"="Identical","First Degree"="First Degree","Second Degree"="Second Degree","Unrelated"="Unrelated"))
  matr1$relatedness<-factor(matr1$relatedness, levels=c("Identical", "First Degree", "Second Degree","Unrelated"))
  matr1$overlap <- cut(matr1$overlap,breaks = c(-Inf,0,10,20,50,100,500,1000,5000,10000,Inf),right = FALSE,
                      label=c("0","0-10","10-20","21-50","50-100","101-500","501-1000","1000-5000","5000-10000",">10000"))
  matr1$dist_thresh <- round(matr1$dist_thresh,1)
  #data[, c(1, 2)] <- as.data.frame(t(apply(data[, c(1, 2)], 1, sort)))data$V1 <- factor(data$V1, levels=c(“Fonds”,“Spy8",“Spy1”,“Spy94a”,“GoyetQ57-3”,“GoyetQ57-2",“GoyetQ57-1”,“GoyetQ305-4",“GoyetQ56-1”,“GoyetQ56-1-lowCov”,“Goyet374a-1”,“Goyet305-7",“GoyetQ55-4”,“Goyet1424-3D”,“GoyetC5-1”,“GoyetQ119-2",“GoyetQ376-25”))
  matr1$lib2 <- factor(matr1$lib2, levels=c("Fonds","Spy1","Spy8","Spy-LowCov-deam","GoyetQ305-7","GoyetQ374a-1","Goyet-LowCov","GoyetQ56-1","GoyetQ305-4","GoyetQ57-1","GoyetQ57-2","GoyetQ57-3"))#data$SP_V1 <- ifelse(data$V1 == “A29253”, “red”, “black”)
  matr1$lib1 <- factor(matr1$lib1, levels=c("Fonds","Spy1","Spy8","Spy-LowCov-deam","GoyetQ305-7","GoyetQ374a-1","Goyet-LowCov","GoyetQ56-1","GoyetQ305-4","GoyetQ57-1","GoyetQ57-2","GoyetQ57-3"))#data$SP_V1 <- ifelse(data$V1 == “A29253”, “red”, “black”)

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
    scale_fill_manual(values = c("Identical"="khaki1","First Degree"="orangered3","Unrelated"="gray80","Second Degree"="darkorchid1")) +
    labs(title = "relatedness", x = "lib_ID", y = "lib_ID")+
    theme(
      panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white",colour = "white", size = 0.5, linetype = "solid"))+
    theme(axis.text.x = element_text(angle = 90,hjust = 1, colour = colors, size = 14)) +
    theme(axis.text.y = element_text(angle = 0,hjust = 1, colour = colors, size = 14)) +
    geom_text(aes(lib1, lib2, label = matr1$dist_thresh), color = "black", size = 5) + #Plot overlap
    theme(legend.text = element_text(size=14),
          legend.title=element_text(size=15))

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


  ggsave(out_rel,height=12,width=14, plot = gfig_rel)
  ggsave(out_over,height=12,width=14, plot = gfig_over)
}

bubbleplot_read(overlap_table=snakemake@input[["over"]], relatable = snakemake@input[["rel"]], out_rel=snakemake@output[["outfile_rel"]], out_over=snakemake@output[["outfile_over"]])

#datafolder="/mnt/diversity/divyaratan_popli/100arc/inbreeding/Chagyrskaya_remapped24sept/curated/"
#overlap_table=paste(datafolder,"overlap_table_fil1.csv",sep="")
#relatable0=paste(datafolder,"hmm_inbred/filtered0/relatable_allLikelihoods.csv", sep = "")
#relatable=paste(datafolder,"hmm_inbred/filtered1/read/READ_results", sep = "")
