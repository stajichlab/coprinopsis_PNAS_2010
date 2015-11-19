A2Bmut <- read.table("A2Bmutr_onlydups.exp",header=T,sep="\t")
AmutB2 <- read.table("AmutB2r_onlydups.exp",header=T,sep="\t")
AmutBmut <- read.table("AmutBmutr_onlydups.exp",header=T,sep="\t")

# A2Bmut

HH <- subset(A2Bmut,A2Bmut$geneconfig == "HeadtoHead")
HT <- subset(A2Bmut,A2Bmut$geneconfig == "HeadtoTail")
TT <- subset(A2Bmut,A2Bmut$geneconfig == "TailtoTail")
NADJ <- subset(A2Bmut,A2Bmut$geneconfig == "NotAdjacent")
length(HH$model1_start);
length(HT$model1_start);
length(TT$model1_start)
length(NADJ$model1_start);
pdf("adjacent_duplicated.pdf")

HH$model1_cat_exp == HH$model2_cat_exp

boxplot(HH$model1_cat_exp,HH$model2_exp_cat,
        HT$model1_cat_exp,HT$model2_exp_cat,
        TT$model1_cat_exp,TT$model2_exp_cat,
        NADJ$model1_cat_exp,NADJ$model2_cat_exp,
        names=c("HH 1","HH 2","HT 1","HT 2","TT 1","TT 2","NA 1","NA 2"),
        main="A2-Bmut model1/model2 ratio");


# AmutB2
HH <- subset(AmutB2,AmutB2$geneconfig == "HeadtoHead")
HT <- subset(AmutB2,AmutB2$geneconfig == "HeadtoTail")
TT <- subset(AmutB2,AmutB2$geneconfig == "TailtoTail")
NADJ <- subset(AmutB2,AmutB2$geneconfig == "NotAdjacent")

length(HH$model1_start);
length(HT$model1_start);
length(TT$model1_start)
length(NADJ$model1_start);

plot(HH$model1_expression,HH$model2_expression,main="Head to Head Amut-B2")
plot(HT$model1_expression,HT$model2_expression,main="Head to Tail Amut-B2")
plot(TT$model1_expression,TT$model2_expression,main="Tail to Tail Amut-B2")
plot(NADJ$model1_expression,NADJ$model2_expression,main="Not Adjacent Amut-B2")
cor(HH$model1_expression,HH$model2_expression)
cor(HT$model1_expression,HT$model2_expression)
cor(TT$model1_expression,TT$model2_expression)
boxplot(HH$model1_expression,HH$model2_expression,
        HT$model1_expression,HT$model2_expression,
        TT$model1_expression,TT$model2_expression,
        NADJ$model1_expression,NADJ$model2_expression,
        names=c("HH 1","HH 2","HT 1","HT 2","TT 1","TT 2","NA 1","NA 2"),
        main="Amut-B2 model1/model2 ratio");


# AmutBmut
HH <- subset(AmutBmut,AmutBmut$geneconfig == "HeadtoHead")
HT <- subset(AmutBmut,AmutBmut$geneconfig == "HeadtoTail")
TT <- subset(AmutBmut,AmutBmut$geneconfig == "TailtoTail")
NADJ <- subset(AmutBmut,AmutBmut$geneconfig == "NotAdjacent")

length(HH$model1_start);
length(HT$model1_start);
length(TT$model1_start)
length(NADJ$model1_start);

plot(HH$model1_expression,HH$model2_expression,main="Head to Head Amut-Bmut")
plot(HT$model1_expression,HT$model2_expression,main="Head to Tail Amut-Bmut")
plot(TT$model1_expression,TT$model2_expression,main="Tail to Tail Amut-Bmut")
plot(NADJ$model1_expression,NADJ$model2_expression,main="Not Adjacent Amut-Bmut")
cor(HH$model1_expression,HH$model2_expression)
cor(HT$model1_expression,HT$model2_expression)
cor(TT$model1_expression,TT$model2_expression)
boxplot(HH$model1_expression,HH$model2_expression,
        HT$model1_expression,HT$model2_expression,
        TT$model1_expression,TT$model2_expression,
        NADJ$model1_expression,NADJ$model2_expression,
        names=c("HH 1","HH 2","HT 1","HT 2","TT 1","TT 2","NA 1","NA 2"),
        main="Amut-Bmut model1/model2 ratio");
