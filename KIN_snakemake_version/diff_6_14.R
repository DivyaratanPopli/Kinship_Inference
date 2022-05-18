#allid="/mnt/diversity/divyaratan_popli/100arc/inbreeding/laurits_writeup/Ch9_rem/neavar/noCh12/ch8libs_fil1/hmm_parameters_fil1/p_all"
#x="/mnt/diversity/divyaratan_popli/100arc/inbreeding/laurits_writeup/Ch9_rem/neavar/noCh12/mergedwin_fil1/pw_Chagyrskaya06_Chagyrskaya14.csv"
#lines="/mnt/diversity/divyaratan_popli/100arc/inbreeding/laurits_writeup/Ch9_rem/neavar/noCh12/hmm_inbred/filtered1/read/output_table.csv"
#fout="/mnt/diversity/divyaratan_popli/100arc/inbreeding/laurits_writeup/Ch9_rem/neavar/noCh12/hmm_inbred/test1.png"



plot6_14 <- function(p1Ch8, f6_14, linesf, fout) {
  
  x <- read.csv(f6_14)
  allid <- read.csv(p1Ch8)

  xprop=sum(x$dis)/sum(x$count)
  
  lines=read.csv(linesf,sep=',')
 
  ch1319prop=lines$NonNormalizedP0[which(lines$Individual1 == 'Chagyrskaya13' &  lines$Individual2 =='Chagyrskaya19')]
  ch131141prop=lines$NonNormalizedP0[which(lines$Individual1 == 'Chagyrskaya1141' &  lines$Individual2 =='Chagyrskaya13')]
  ch191141prop=lines$NonNormalizedP0[which(lines$Individual1 == 'Chagyrskaya1141' &  lines$Individual2 =='Chagyrskaya19')]


  dprop <- vector()

  for (i in 1:22) {
    #print(i)
    d <- x[x$chrom !=i,]
    dprop[i] <-sum(d$dis)/sum(d$count)
  }

  png(file=fout,
      width=750, height=550)


  c1 <- rgb(133,186,255,max = 255, alpha = 100, names = "lt.blue")
  c2 <- rgb(255,162,203, max = 255, alpha = 30, names = "lt.pink")

  hist(allid$prop,breaks = 50, col=c1, ylim = c(0,80), freq = TRUE, xlab = "Proportion of differences", main = "Histogram of differences in identical individuals" )
  hist(dprop,col=c2, ylim=c(0,80), freq=TRUE, add=T)
  box()
  #c11 <- rgb(250,12,233, max = 255, alpha = 100, names = "lt.pink")theme(axis.text.x = element_text(angle = 90, vjust = 0.9)) +
  
  segments(x0=xprop,y0=0,x1=xprop,y1=80,col="red")
  segments(x0=ch1319prop,y0=0,x1=ch1319prop,y1=80,col="green")
  segments(x0=ch131141prop,y0=0,x1=ch131141prop,y1=80,col="orange")
  segments(x0=ch191141prop,y0=0,x1=ch191141prop,y1=80,col="purple")
  
  segments(x0=quantile(allid$prop,0.975),y0=0,x1=quantile(allid$prop,0.975),y1=80,col="black",lty = 2)
  segments(x0=quantile(allid$prop,0.025),y0=0,x1=quantile(allid$prop,0.025),y1=80,col="black",lty = 2)
  

  legend(title="Vertical lines",0.05, 60, legend=c("Ch13_Ch19", "Ch13_Ch1141", "Ch1141_Ch19","Avg diff for Ch06_Ch14", "2.5% & 97.5% quantiles"),
         col=c("green", "orange", "purple","red","black"), lty=c(1,1,1,1,2), cex=0.8)

  legend(title="Histograms","topright", c("Identical_differences", "Ch06_Ch14"), col=c(c1,c2), lwd=10)

  dev.off()
  
  print("Ch08")
  print(allid)
  print("6_14")
  print(dprop)
  print("ch1319")
  print(ch1319prop)
  print("ch131141")
  print(ch131141prop)
  print("ch191141")
  print(ch191141prop)
  print("614")
  print(xprop)
  
}


  plot6_14(p1Ch8=snakemake@input[["phalf"]], f6_14 = snakemake@input[["f6_14"]], linesf= snakemake@input[["linesf"]], fout=snakemake@output[["fout"]])

  #plot6_14(p1Ch8=allid, f6_14 =x, linesf=lines, fout=fout)

  
  