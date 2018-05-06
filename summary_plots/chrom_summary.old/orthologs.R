orthologsarms <- read.table("orthologs_arms.dat",header=T)
orthologsctr <- read.table("orthologs_center.dat",header=T)
pdf("orthologs.pdf")
boxplot(orthologsarms$TOTAL,orthologsctr$TOTAL,main="orthologs Density BoxPlot", outline=FALSE, names=c("Arms","Center"))
