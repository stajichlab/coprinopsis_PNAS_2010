paralogsarms <- read.table("paralogs_arms.dat",header=T)
paralogsctr <- read.table("paralogs_center.dat",header=T)
pdf("paralogs.pdf")
boxplot(paralogsarms$TOTAL,paralogsctr$TOTAL,main="paralogs Density BoxPlot", outline=FALSE, names=c("Arms","Center"))
ks.test(paralogsarms$TOTAL,paralogsctr$TOTAL)
