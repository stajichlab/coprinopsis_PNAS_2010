genesarms <- read.table("genes_arms.dat",header=T)
genesctr <- read.table("genes_center.dat",header=T)
pdf("genes.pdf")
boxplot(genesarms$TOTAL,genesctr$TOTAL,main="genes Density BoxPlot", outline=FALSE, names=c("Arms","Center"))
