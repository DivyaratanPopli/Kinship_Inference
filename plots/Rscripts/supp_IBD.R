library(ggplot2)
#infile1='/mnt/diversity/divyaratan_popli/100arc/inbreeding/model_git/plots/plotdata/plot_IBDaccuracy.csv'
#infile2='/mnt/diversity/divyaratan_popli/100arc/inbreeding/model_git/plots/plotdata/plot_IBDaccuracy_cont.csv'
#infile3='/mnt/diversity/divyaratan_popli/100arc/inbreeding/model_git/plots/plotdata/plot_IBDaccuracy_asc.csv'
#infile4='/mnt/diversity/divyaratan_popli/100arc/inbreeding/model_git/plots/plotdata/plot_IBDaccuracy_roh.csv'

IBD_accuracy_plot <- function(infile1,infile2,infile3,infile4,outfile){
  bigd <-data.frame()
  for(j in c(1,2,3,4)) {
    if (j==1){
      infile=infile1
    }else if (j==2){
      infile=infile2
    }else if (j==3){
      infile=infile3
    }else if (j==4){
      infile=infile4
    }
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
    df$Scenario=j
    bigd=rbind(bigd,df)
  }
  
  #bigd_un=bigd[bigd$rel=='Siblings',]
  #bigd_un_control=bigd_un[bigd_un$Scenario==1,]
  #bigd_un_c=bigd_un[bigd_un$Scenario==2,]
  #bigd_un_a=bigd_un[bigd_un$Scenario==3,]
  #bigd_un_r=bigd_un[bigd_un$Scenario==4,]
  
  bigd$Scenario=as.factor(bigd$Scenario)
  colors1 <- c("Unrelated" = "grey30", "3rd Degree" = "blue", "2nd Degree" = "darkorchid1", "Siblings"= "Green", "Parent-Child" = "orangered3", "Identical" = "darkgoldenrod1")
  line1 <- c("1" = "solid", "2" = "dashed", "3" = "dotted","4" = "12345678")
  
  f <- ggplot(
    bigd,
    aes(x = cov, y = accuracy, ymin = accuracy-stdev, ymax = accuracy+stdev, group = interaction(rel,Scenario) ,colour=rel, linetype = Scenario)
  )
  f <- f + geom_errorbar(width = 0.2) +
    geom_point(size = 1.5)
  f <- f + labs(y="Accuracy", x="Coverage") +
    scale_color_manual(values = colors1, name="Relatedness") + theme_bw() +geom_path(aes(x=cov, y=accuracy, group=interaction(rel,Scenario),colour=rel, linetype = Scenario)) +
    scale_linetype_manual(values = line1, name="Scenario", labels=c("Control","C","A","R")) +
    theme(legend.position = "bottom", text = element_text(size=15))
  f
  ggsave(outfile,
         width = 10, height = 4, dpi = 150, units = "in", device='png')
}

IBD_accuracy_plot(infile1=snakemake@input[["inf"]],infile2=snakemake@input[["inf_c"]], infile3=snakemake@input[["inf_a"]], infile4=snakemake@input[["inf_r"]], outfile=snakemake@output[["outf"]])
