#################################################
# Statistical Analysis with BioMedical Data (I) #
#################################################

library(MASS)
attach(Pima.tr)
shapiro.test(bmi)
qqnorm(bmi)
qqline(bmi)

t.test(bmi, mu=30)
bmi.ttest <- t.test(bmi, mu=30)
names(bmi.ttest)
bmi.ttest$p.value
t.test(bmi, mu=30, alternative="greater")
t.test(bmi, mu=30, alternative="less")

var.test(bmi ~ type)

t.test(bmi ~ type)

data(anorexia)

FT <- subset(anorexia, Treat=='FT')
head(FT)
shapiro.test(FT$Prewt - FT$Postwt)
t.test(FT$Prewt, FT$Postwt, paired=TRUE)

CBT <- subset(anorexia, Treat=='CBT')
shapiro.test(CBT$Prewt - CBT$Postwt)
wilcox.test(CBT$Prewt, CBT$Postwt, paired=TRUE)

placebo <- c(7, 5, 6, 4, 12)
new_drug <- c(3, 6, 4, 2, 1)
wilcox.test(placebo, new_drug, exact=FALSE)

attach(iris)
cor(Sepal.Length, Petal.Width)
cor.test(Sepal.Length, Petal.Width)
cor(iris[, 1:4])
pairs(iris[, 1:4])

install.packages("Hmisc")
library(Hmisc)
rcorr(as.matrix(iris[, 1:4]))

iris.na.test <- iris[, 1:4]
iris.na.test[1,1] <- NA
iris.na.test[3,2] <- NA
iris.na.test[4,3] <- NA
head(iris.na.test)
cor(iris.na.test)

head(iris.na.test)
cor(iris.na.test, use="complete.obs")

head(iris.na.test)
cor(iris.na.test, use="pairwise.complete.obs")

doctorA <- c(4,1,3,2,6,5,8,7)
doctorB <- c(5,3,1,2,6,4,7,8)
cor.test(doctorA, doctorB, method="spearman")

newTest <- c(50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)
standardTest <- c(61, 61, 59, 71, 80, 76, 90, 106, 98, 100, 114)
dat <- data.frame(newTest, standardTest)
plot(standardTest ~ newTest, data=dat, xlim=c(0, 110), ylim=c(0, 120))
dat.lm <- lm(standardTest ~ newTest, data=dat)
summary(dat.lm)
abline(dat.lm, col="red")

coef(dat.lm)
predict(dat.lm, newdata=data.frame(newTest=80))

par(mfrow=c(1,2))
plot(resid(dat.lm))
abline(1, 0, col="red")
qqnorm(resid(dat.lm))
qqline(resid(dat.lm), col="red")
shapiro.test(resid(dat.lm))
install.packages("lmtest")
library(lmtest)
dwtest(standardTest ~ newTest, data=dat)
