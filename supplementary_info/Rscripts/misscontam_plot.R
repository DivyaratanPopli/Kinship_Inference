
library(ggplot2)
#install.packages("magrittr")
#install.packages("dplyr")
library(magrittr)
library(dplyr)
library(tidyr)

#inf='/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/hmm_numba_incorrect_roh_i1_i1_insteadOf_i1_i2/roc/contam0_inbred0_model_performance_allroc_asc0.csv.gz'
#infile='/mnt/diversity/divyaratan_popli/review_sim/hmm_numba/miss_contam/output/missC10/contam1_inbred0_model_performance_allLikelihoods_inphapProbs_asc0.csv'

contamination_plot <- function(infile,outfile){


  compinp=read.csv(file = infile, header = TRUE)
  compinp$Relatedness <- factor(compinp$Relatedness, levels = c("fir", "sec", 'un_all','pc','sib','deg3'))

  data <- compinp %>%
    pivot_longer(c(True_positive,False_positive), names_to = "T_F", values_to = "classification_proportion")#Plot
  data$T_F <- factor(data$T_F, levels= c('True_positive','False_positive'))

  #data<-data[!(data$method=="read_inppshap" | data$method=="allLikelihoods_inppshap"),]
  #data = subset(data, select = -c(method) )

  data$Coverage <- as.factor(data$Coverage)

  #changing some column names:

  Relatedness.labs <- c("Unrelated", "3rd Degree", "2nd Degree", "1st Degree", "Siblings", "Parent-Child")
  names(Relatedness.labs) <- c("un_all", "deg3", "sec", "fir", "sib", "pc")

  T_F.labs <- c("False Positive", "True Positive")
  names(T_F.labs) <- c("False_positive", "True_positive")

  colors1 <- c("0.1"="#999933", "0.2"="#117733","4"="#661100", "0.05"="#0072B2", "0.5"="#E69500")

  plot<-ggplot(data=data, aes(x=contamination, y=classification_proportion, color=Coverage))
  plot<-plot+geom_line()
  plot<- plot + facet_grid(T_F ~ Relatedness, labeller = labeller(T_F = T_F.labs, Relatedness = Relatedness.labs), scales="free") +
  scale_color_manual(values = colors1)
  #plot + geom_hline(yintercept=0.05, linetype="dashed")
  plot <- plot + geom_hline(data = data.frame(T_F='False_positive'), aes(yintercept = 0.05), linetype = "dotted")
  plot <- plot + geom_vline(data = data, aes(xintercept = 2.5), linetype = "dashed")

  plot <- plot+ labs(y="Classification Proportions", x = "miss-specified contamination level(%)")
  plot + theme_bw() + theme(
    strip.background = element_rect(
      color="transparent", fill="transparent", size=1.5, linetype="solid"
    ), legend.position="bottom", text = element_text(size=12)
  )
  ggsave(outfile, width = 8, height = 5, dpi = 150, units = "in", device='png')

}


contamination_plot(infile=snakemake@input[["infile"]], outfile=snakemake@output[["outfile"]])
