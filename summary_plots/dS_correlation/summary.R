cold = read.table("dS_COLD.tab",header=F);
hot = read.table("dS_HOT.tab",header=F);
neut = read.table("dS_NEUTRAL.tab",header=F);
unscored = read.table("dS_UNSCORED.tab",header=F);
cold = subset(cold$V1,cold$V1 < 3 & cold$V1 > 0)
hot = subset(hot$V1,hot$V1 < 3 & hot$V1 > 0)
neut = subset(neut$V1,neut$V1 < 3 & neut$V1 > 0)
unscored = subset(unscored$V1,unscored$V1 < 3 & unscored$V1 > 0)

pdf("dS_correlation_boxplot.pdf")

boxplot(cold,hot,neut,unscored,main="Paralog dS comparison among recombination", outline=FALSE, names=c("COLD","HOT","NEUTRAL","UNSCORED"))
pdf("dS_correlation_histogram.pdf")
hist(cold,100);
hist(hot,100);
hist(neut,100);
hist(unscored,100);
summary(cold)
summary(hot)
summary(neut)
summary(unscored)
ks.test(cold,hot)
ks.test(cold,hot)
ks.test(cold,neut)
ks.test(unscored,neut)
ks.test(unscored,hot)
