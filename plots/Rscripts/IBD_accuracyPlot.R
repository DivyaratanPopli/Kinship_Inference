library(ggplot2)

#infile='/mnt/diversity/divyaratan_popli/100arc/inbreeding/model_git/plots/plotdata/plot_IBDaccuracy.csv'
IBD_accuracy_plot <- function(infile,outfile){

  data=read.csv(file = infile, header = TRUE)


  data[data$rel=='deg4','rel']='Unrelated'
  data[data$rel=='deg5','rel']='Unrelated'
  data[data$rel=='un','rel']='Unrelated'

  data[data$rel=='avu','rel']='2nd Degree'
  data[data$rel=='gr','rel']='2nd Degree'
  data[data$rel=='hsib','rel']='2nd Degree'

  data[data$rel=='id','rel']='Identical'
  data[data$rel=='pc','rel']='Parent-Child'
  data[data$rel=='sib','rel']='Siblings'

  data[data$rel=='deg3','rel']='3rd Degree'

  allrel=c("Unrelated", "3rd Degree", "2nd Degree", "Siblings", "Parent-Child", "Identical")
  allcov=c(4,0.5,0.2,0.1,0.05,0.03)

  accuracy=c()
  stdev=c()
  rel=c()
  cov=c()
  for(cv in allcov) {

    datacv=data[data$cov==cv,]
    for(i in allrel) {
      temp=datacv[datacv$rel==i,]
      accuracy=c(accuracy,mean(temp$accuracy))
      stdev=c(stdev,sd(temp$accuracy))
      rel=c(rel,i)
      cov=c(cov,cv)
    }
  }
  df=data.frame(accuracy,stdev,rel,cov)

  df$cov=as.factor(df$cov)

  colors1 <- c("Unrelated" = "grey30", "3rd Degree" = "blue", "2nd Degree" = "darkorchid1", "Siblings"= "Green", "Parent-Child" = "orangered3", "Identical" = "darkgoldenrod1")

  f <- ggplot(
    df,
    aes(x = cov, y = accuracy, ymin = accuracy-stdev, ymax = accuracy+stdev, colour=rel)
  )
  f <- f + geom_errorbar(width = 0.2) +
    geom_point(size = 1.5)
  f <- f + labs(y="Accuracy", x="Coverage") +
    scale_color_manual(values = colors1, name="Relatedness") + theme_bw() +geom_path(aes(x=cov, y=accuracy, group=rel)) +
    theme(legend.position = "bottom")

  ggsave(outfile,
         width = 10, height = 4, dpi = 150, units = "in", device='png')
  }

IBD_accuracy_plot(infile=snakemake@input[["inf"]], outfile=snakemake@output[["outf"]])
