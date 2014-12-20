library(ggplot2)
library(reshape2)

args=commandArgs(trailingOnly = TRUE)

infile= args[1]
outfile= args[2]

data=read.table(infile,sep="\t",header=TRUE)

data$Gene=data$TargetAmplicon

data$Gene=gsub("\\.\\S+","",data$Gene,perl=T)
data$Gene=gsub("_\\S+","",data$Gene,perl=T)

melted=melt(data,id.vars=c("TargetAmplicon","Gene"))
pointPlot=ggplot(melted,aes(TargetAmplicon,log10(value),color=variable)) + 
geom_point() +
geom_smooth(aes(group=variable),se=FALSE) +
theme(axis.text.x = element_text(angle=90,size=5)) +
ylab("Coverage (based on full read overlap, log10 scaled") + xlab("Gene ")
ggsave(filename=outfile,width=30)


