cold = read.table("dS_COLD.tab",header=F);
notcold= read.table("dS_NOTCOLD.tab",header=F);
cold = subset(cold$V1,cold$V1 < 3 & cold$V1 > 0)
notcold= subset(notcold$V1,notcold$V1 < 3 & notcold$V1 > 0)

pdf("dS_correlation_boxplot_notcold.pdf")

boxplot(cold,notcold,main="Paralog dS comparison among recombination", outline=FALSE, names=c("COLD","NOT-COLD"))
pdf("dS_correlation_histogram_notcold.pdf")
hist(cold,100);
hist(notcold,100);
summary(cold)
summary(notcold)
ks.test(cold,notcold)
t.test(cold,notcold)
