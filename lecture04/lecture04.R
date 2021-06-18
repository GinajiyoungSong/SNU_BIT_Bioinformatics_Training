#################################################
# Statistical Analysis with BioMedical Data (II)#
#################################################

# Chi-Square Test (적합성검정)
chisq.test(c(24,16), p=c(0.7,0.3))

# (독립성검정) 두 변수가 독립인지 검정
countTable <- matrix( c(10845, 180, 10933, 104), nrow=2, byrow=TRUE)
rownames(countTable) <- c("Placebo", "Aspirin")
colnames(countTable) <- c("No Heart Attack", "Heart Attack")
# 아스피린 처방과 심장마비가 상호 독립적인지 아닌지 검정
chisq.test(countTable)


# 신생아 저체중과 산모데이터 로드
library(MASS)
head(birthwt) # 데이터 확인

birthwt$smoke <- factor(birthwt$smoke, labels = c("Non Smoker", "Smoker"))
birthwt$low <- factor( birthwt$low, label=c("No", "Yes"))
smoke_low_tb <- table( birthwt$smoke, birthwt$low)

# 2x2 분할표를 input 하여 독립성 검정
chisq.test(smoke_low_tb)



# fisher test
TeaTasting <- matrix( c(3, 1, 1, 3), nrow=2 )
fisher.test( TeaTasting )


# Cocran-Armitage Tend Test
prop.trend.test( c(13,7,21), c(42,14,28))

chisq.test(matrix(c(13, 29, 10, 4, 10, 18), ncol=3))


# McNemar's Test
mcnemar.test(matrix(c(5,15,5,7), ncol=2))



# CrossTable()함수
install.packages("gmodels")
library(gmodels)
CrossTable(countTable, fisher = TRUE, chisq = TRUE, mcnemar = TRUE, expected = TRUE, sresid = TRUE, format="SPSS")





# ANOVA Test
# 3개 이상 다수의 집단을 비교하고자 할 떄
attach(anorexia)
Change <- Postwt - Prewt
boxplot( Change ~ Treat, col=rainbow(3))
# 3개의 집단, 대조군, 가족치료(FT), CBT 치료 결과 변화 시각화

aov.out <- aov( Change ~ Treat )
summary(aov.out)
# ANOVA 결과 pr<0.001 이므로 세가지 치료방법에 따라 몸무게 평균변화가 있음 


# 사후검정
# 3 집단간에 어떤 집단의 모평균이 차이가 나는지 알기 위해서 다중비교


# Tukey 사후검정
TukeyHSD( aov.out )
# 2집단간 pairwise 하게 비교하는 결과가 나옴, FT-Cont p = 0.0045127 유의한것을 알수 있음
plot(TukeyHSD(aov.out))

# LSD 검정 함수
pairwise.t.test( Change, Treat) #본페로니 보정 결과 유의한 값 도출 p=0.0048



getwd()
setwd("/Users/SJY/Documents/유전체 분석과정/lecture04")
teaching_time <- read.table("teaching_time.txt", header = TRUE, sep=" ")


#반복이 없는 이원분산분석 ( two-way ANOVA )
aov.out <- aov(days ~ ageGroup + method, data=teaching_time)
summary(aov.out)


# 반복이 있는 이원분산분석
summary(ToothGrowth)
ToothGrowth$dose <- factor(ToothGrowth$dose)
summary(ToothGrowth)

aov.out <- aov( len ~ supp*dose, data=ToothGrowth )
summary(aov.out)

# 교호작용
interaction.plot(ToothGrowth$dose, ToothGrowth$supp, ToothGrowth$len)
# interaction.plot( pred1, pred2, resp ) pred1, pred2 : 범주형 예측 변수, resp:반응변수 이빨의 길이

