
library(ggplot2)
#install.packages("magrittr")
#install.packages("dplyr")
library(magrittr)
library(dplyr)
library(tidyr)

#inf='/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/hmm_numba_incorrect_roh_i1_i1_insteadOf_i1_i2/roc/contam0_inbred0_model_performance_allroc_asc0.csv.gz'

cutoff_plot <- function(infile,outfile){


  compinp=read.csv(file = infile, header = TRUE)
  compinp$Relatedness <- factor(compinp$Relatedness, levels = c("id", "fir", "sec", 'un_all','pc','sib','deg3','un_all_withdeg3'))

  data <- compinp %>%
    pivot_longer(c(True_positive,False_positive), names_to = "T_F", values_to = "classification_proportion")#Plot
  data$T_F <- factor(data$T_F, levels= c('True_positive','False_positive'))

  data<-data[!(data$method=="read_inppshap" | data$method=="allLikelihoods_inppshap"),]
  #data = subset(data, select = -c(method) )

  data$Coverage <- as.factor(data$Coverage)

  #changing some column names:

  Relatedness.labs <- c("Unrelated", "Unrelated w/o \n 3rd Degree", "3rd Degree", "2nd Degree", "1st Degree", "Siblings", "Parent-Child", "Identical")
  names(Relatedness.labs) <- c("un_all", "un_all_withdeg3", "deg3", "sec", "fir", "sib", "pc", "id")

  T_F.labs <- c("False Positive", "True Positive")
  names(T_F.labs) <- c("False_positive", "True_positive")

  plot<-ggplot(data=data, aes(x=cutoff, y=classification_proportion, color=Coverage))
  plot<-plot+geom_line()
  plot<- plot + facet_grid(T_F ~ Relatedness, labeller = labeller(T_F = T_F.labs, Relatedness = Relatedness.labs))
  #plot + geom_hline(yintercept=0.05, linetype="dashed")
  plot <- plot + geom_hline(data = data.frame(T_F='False_positive'), aes(yintercept = 0.05), linetype = "dotted")
  plot <- plot+ labs(y="Classification Proportions", x = "Log likelihood Ratio Cutoff")
  plot + theme_bw() + theme(
    strip.background = element_rect(
      color="transparent", fill="transparent", size=1.5, linetype="solid"
    ), legend.position="bottom", text = element_text(size=10)
  )
  ggsave(outfile, width = 8, height = 5, dpi = 150, units = "in", device='png')

}


cutoff_plot(infile=snakemake@input[["infile_cut"]], outfile=snakemake@output[["outfile_cut"]])
#cutoff_plot(infile=inf,outfile="/mnt/diversity/divyaratan_popli/test.png")
