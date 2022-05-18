library(ggplot2)


#tabf <- "/mnt/diversity/divyaratan_popli/100arc/inbreeding/laurits_writeup/Ch9_rem/neavar/curated/hmm_inbred/filtered0/read/output_table.csv"

readout <- function(tabf,out){

  tab <- read.csv(tabf, header = TRUE,sep = ',')
  tab$ind12 <- paste(tab$Individual1, tab$Individual2,sep='_')

  tabs <- tab[order(tab$NonNormalizedP0),]

  normalization_val = median(tab$NonNormalizedP0, na.rm = TRUE)

  row.names(tabs) <- NULL


  ggplot(
    data = tabs,
    aes(x = reorder(ind12, NonNormalizedP0), y = NonNormalizedP0, ymin = NonNormalizedP0 - (2*NonNormalizedStandardError), ymax = NonNormalizedP0 + (2*NonNormalizedStandardError) ))+
    geom_point() + geom_errorbar(width = 0.2) +
    geom_hline(aes(yintercept = normalization_val,linetype = "median"), color = "blue", size = 0.5) +
    geom_hline(aes(yintercept = normalization_val*0.90625,linetype = "2nd Degree"), color = "blue", size = 0.5) +
    geom_hline(aes(yintercept = normalization_val*0.8125,linetype = "1st Degree"), color = "blue", size = 0.5) +
    geom_hline(aes(yintercept = normalization_val*0.625,linetype = "Identical"), color = "blue", size = 0.5) +
    scale_linetype_manual(name="Upper Threshold:", values = c("median"=1,"2nd Degree"=2,"1st Degree"=4,"Identical"=3)) +
    labs(x = "Pair of individuals", y = "Pairwise differences") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.9))

  ggsave(filename=out, width = 20, height = 10)
}

readout(tabf=snakemake@input[["table"]], out=snakemake@output[["png"]])

#guide = guide_legend(override.aes = list(color = c("blue", "blue","blue","blue")))
#graph.1 <- ggplot(data = data.frame(x = 1:1000, y = cumulative.sum.of.means), aes(x, y)) + geom_line(size = 1) +
#  geom_hline(aes(yintercept = mean(means.random.numbers),linetype = "NameOfLine1"), color = "blue", size = 1) +
#  geom_hline(aes(yintercept = 1/0.2,linetype = "NameOfLine2"), color = "red", size = 1) +
#  scale_linetype_manual(name="Title of legend", values = c(1,2), guide = guide_legend(override.aes = list(color = c("blue", "red"))))


#code for legend:https://www.reddit.com/r/rstats/comments/8q5a2n/adding_a_custom_ggplot_legends_while_using_geom/
