repeatsarms <- read.table("repeats_arms.dat",header=T)
repeatsctr <- read.table("repeats_center.dat",header=T)
pdf("repeats.pdf")
boxplot(repeatsarms$TOTAL,repeatsctr$TOTAL,main="repeats Density BoxPlot", outline=FALSE, names=c("Arms","Center"))
