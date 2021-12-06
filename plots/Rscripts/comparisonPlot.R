library(ggplot2)
#install.packages("magrittr")
#install.packages("dplyr")
library(magrittr)
library(dplyr)
library(tidyr)

#inf='/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/hmm_numba_incorrect_roh_i1_i1_insteadOf_i1_i2/roc/comparison_plot_data.csv.gz'

comparison_plot <- function(inf,outf){


  compinp=read.csv(file = inf, header = TRUE)
  compinp$Relatedness <- factor(compinp$Relatedness, levels = c("id", "fir", "sec", 'un_all','pc','sib','deg3','un_all_withdeg3'))

  data <- compinp %>%
    pivot_longer(c(True_positive,False_positive), names_to = "T_F", values_to = "classification_proportion")#Plot
  data$T_F <- factor(data$T_F, levels= c('True_positive','False_positive'))

  colnames(data)[4] <- "Method"
  data$Method=as.character(data$Method)
  data[data$Method=='allLikelihoods_inphapProbs','Method'] = 'KIn'
  data[data$Method=='read_inppshap','Method'] = 'READ'




  colnames(data)[9] <- "Simulation"
  data$Simulation=as.character(data$Simulation)
  data[data$Simulation=='0_0_0','Simulation'] = 'Control'
  data[data$Simulation=='0_0_1','Simulation'] = 'With ROH'
  data[data$Simulation=='0_2_0','Simulation'] = 'With Ascertainment'
  data[data$Simulation=='1_0_0','Simulation'] = 'With Contamination'

  data$Coverage <- as.factor(data$Coverage)
  data$Simulation <-as.factor(data$Simulation)

  #changing some column names:

  Relatedness.labs <- c("Unrelated", "Unrelated w/o \n 3rd Degree", "3rd Degree", "2nd Degree", "1st Degree", "Siblings", "Parent-Child", "Identical")
  names(Relatedness.labs) <- c("un_all", "un_all_withdeg3", "deg3", "sec", "fir", "sib", "pc", "id")

  T_F.labs <- c("False Positive", "True Positive")
  names(T_F.labs) <- c("False_positive", "True_positive")


  plot<-ggplot(data=data, aes(x=Coverage, y=classification_proportion, group=interaction(Simulation,Method), color=Simulation, lty=Method))
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
  ggsave(outf,
         width = 8, height = 5, dpi = 150, units = "in", device='png')

}

comparison_plot(inf=snakemake@input[["infile_comp"]], outf=snakemake@output[["outfile_comp"]])
#comparison_plot(inf,"/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/hmm_numba_incorrect_roh_i1_i1_insteadOf_i1_i2/roc/comparison_plot_test.png")
