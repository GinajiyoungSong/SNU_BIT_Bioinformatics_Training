#################################################
# Statistical Analysis with BioMedical Data (I) #
#################################################

library(MASS)

'''이 실습에서는 R을 이용하여 정규성 검정, 분산비 검정, 두 군 간의 평균 비교,
상관관계 분석 및 기본적인 회귀분석 절차에 관해 기술하였다. 두 군 간의 평균비
교에서는 모수적/비모수적 검정을 사용하며 상관관계 에서는 Pearson,
Spearman, Kendall 상관 계수 구하기 및 통계적 검정을 수행한다. 여러 차원의
변수가 있는 경우의 all pair-wise 한 상관계수 구하는 법에 대해 설명한다. 회귀
분석은 단순선형회귀분석을 다루고 있으며 주어진 데이터를 이용하여 선형 회귀식
을 적합하고 선형식의 기울기 및 절편 그리고 잔차 분석을 수행한다'''

attach(Pima.tr) # 데이터 로드

# Normality Test (정규성 검정)
# 체지방 지수를 나타내는 BMI(body mass index)변수가 정규분포를 따르는지 확인
shapiro.test(bmi) 
# p < 0.05 이면 귀무가설을 기각 = 정규분포를 따르지 않는다
# p > 0.05 이면 귀무가설 = 정규분포를 따른다. (정규성 확인)

# 시각화
qqnorm(bmi)  # 정규 분포의 Q-Q plot을 그리는 것으로 정규분포에 얼마나 근접한지
qqline(bmi)  # 1Q와 3Q를 지나는 선을 그려준다



# One sample t-test
# Pima 여성들의 평균 bmi 값을 mu 이라고 했을 때, mu=30 -> 귀무가설
t.test(bmi, mu=30)

''' p-vlaue가 0.05 이하이므로 귀무가설을 기각한다. 즉 mu=30 이라는 가설은 틀렸다/
당뇨의 유무에 따라 bmi의 차이가 있다고 결론 내릴수 있다.'''


bmi.ttest <- t.test(bmi, mu=30)
names(bmi.ttest) # 확인 할수 있는 parameter name
bmi.ttest$p.value
bmi.ttest$statistic



# 양측 검정(two-tailed test) opt 추가
# 양측검정을 수행할 경우에는 alternative= 옵션을 사용하면 된다.

t.test(bmi, mu=30, alternative="greater") 
# 대립 가설이 mu > 30 경우
# 평균이 30이라는 귀무가설을 기각하고 평균이 30보다 크다는 대립가설을 채택.

t.test(bmi, mu=30, alternative="less") #  mu < 30인 경우






# two sample t-test & F test
'''서로 독립인 두 집단의 평균을 비교( 두 그룹간의 평균이 같다/같지 않다)하기 위
해서는 two-sample t-test를 사용할 수 있다. 예를 들어 실험군과 대조군을 설정
하여 특정 약물이나 치료의 효과를 비교하는 것이 대표적인 예일 것이다. 이 때
두 집단의 변수가 정규분포를 하고 있다면 두 집단의 분산이 같은지에 따라 각기
다른 검정법을 사용해야 된다. '''

# F-test 를 이용하여 등분산 검정을 먼저 수행한다
# 두 집단의 등분산 검정을 위해서 : var.test()
var.test(bmi ~ type)


# 등분산 검정 결과 p-value < 0.05 귀무가설 기각. 등분산이 아니다



# 두 집단의 분산이 같은 경우, pooled variance를 이용 -> two-sample test
# 등분산 확인했을 때 분산이 같지 않음 확인 -> Welch t-test 수행
t.test(bmi ~ type)

t.test( value ~ group ), t.test( sample1, sample2 ) : Welch t-test
t.test( value ~ group, var.equal=TRUE) or
t.test( sample1, sample2, var.equal=TRUE): Pooled variance t-test


# p-vlaue가 0.05 이하이므로 당뇨의 유무에 따라 bmi의 차이가 있다



# MASS 패키지 anorexia 데이터를 이용
# Paired t-test 

data(anorexia)

'''거식증 환자에 대해 가족치료(Family Treatment, FT)를 수행한 전 후의 몸무게
차이변화가 유의한지 검정

귀무가설 Mu1 - Mu2 = 0 *평균이 같음
anorexia 데이터에는 73명의 거식증 환자의 인지행동치료(CBT, Cognitive Behavior
treatement), 가족치료(FT, Family Treatment) 및 컨트롤군(Cont, control) 총 3
군의 데이터로 이루어져 있으며 각 치료법을 사용하기 전 후의 환자의 몸무게를
측정한 결과이다.

'''
FT <- subset(anorexia, Treat=='FT')
head(FT)

#몸무게 전 - 몸무게 후
shapiro.test(FT$Prewt - FT$Postwt)
# p-value> 0.05 : 유의수준(significance level)을 0.05로 할 경우에 귀무가설을 기각시키지 못하므로 이 분포는 정규분포를 따른다고 할 수 있다.

# paired opt 추가
t.test(FT$Prewt, FT$Postwt, paired=TRUE)
# p-value <0.05 이므로 귀무가설 기각 = 치료전후 체중차이있다. 가족치료가 효과가 있었다







# 비모수적 검정방법

CBT <- subset(anorexia, Treat=='CBT')
shapiro.test(CBT$Prewt - CBT$Postwt)

wilcox.test(CBT$Prewt, CBT$Postwt, paired=TRUE)
# 경고메세지 = 데이터의 순위 중 동률이 존재하기 때문에 출력되는 결과

wilcox.test(CBT$Prewt, CBT$Postwt, paired=TRUE, exact=FALSE)
# error 제거 exact=FALSE
# p-value>0.05 = 귀무가설 = 유의하지않음 = CBT치료에 의한 몸무게 변화의 차이가 없다.




# Wilcoxon rank sum test( 또는 Mann-Whitney U Test )
placebo <- c(7, 5, 6, 4, 12)
new_drug <- c(3, 6, 4, 2, 1)
wilcox.test(placebo, new_drug, exact=FALSE)
# p-value 0.05를 넘어 기각 = 신약 효능에 별 차이가 없다




attach(iris)
cor(Sepal.Length, Petal.Width)
# p-value, CI 결과return 해주는 상관계수 함수 사용
cor.test(Sepal.Length, Petal.Width)



# col 열 변수 모두 상관계수 분석 = heatmap 시각화 용도로 많이 사용
cor(iris[, 1:4])
pairs(iris[, 1:4]) # 데이터 시각화 산점도 그래프



install.packages("Hmisc")
library(Hmisc)
# 상관계수와 P-value 값을 return 하는 함수를 사용
# 패키지 설치 후 사용해야함 - default corr = pearson 모수 검정
rcorr(as.matrix(iris[, 1:4]))



# 결측치가 있는 경우 상관계수 도출 방법
iris.na.test <- iris[, 1:4] # 결측데이터 만들기위해 복사
iris.na.test[1,1] <- NA
iris.na.test[3,2] <- NA
iris.na.test[4,3] <- NA
head(iris.na.test)
cor(iris.na.test) # 결측치가 있는 경우 cor 함수가 작동하지 않음
#데이터 전처리의 중요성 - 결측치 처리를 해야 검정이 가능하다

#결측치 있는 행을 버리고 상관계수 뽑아냄 = List-wise deletion
cor(iris.na.test, use="complete.obs")


# 결측치가 있는 pair-wise 값만 제거함 : pair-wise deletion
cor(iris.na.test, use="pairwise.complete.obs")




# Spearman corr 계산 => 비모수적 검정

doctorA <- c(4,1,3,2,6,5,8,7)
doctorB <- c(5,3,1,2,6,4,7,8)
cor.test(doctorA, doctorB, method="spearman")
# p-value <0.05 이므로 귀무가설을 기각한다 (두 의사의 질환 점수 평가 r=0 이다)
# rho 값이 1에 가가운 것으로 보아 좋은 일치도를 나타낸다고 볼 수 있다







# 선형 회귀 lm() 함수
# 반응변수 ~ 설명변수

newTest <- c(50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)
standardTest <- c(61, 61, 59, 71, 80, 76, 90, 106, 98, 100, 114)
dat <- data.frame(newTest, standardTest)

plot(standardTest ~ newTest, data=dat, xlim=c(0, 110), ylim=c(0, 120))

dat.lm <- lm(standardTest ~ newTest, data=dat)
# 회귀식 lm 함수 -> 1차식의 계수와 상수


summary(dat.lm)
# summary R square 계산결과를 보임. 상세 정보, R2설명력, F 통계량 등등등

abline(dat.lm, col="red")
#회귀선 그려주는 함수

coef(dat.lm)
#계수. 상수 호출 함수


predict(dat.lm, newdata=data.frame(newTest=80))
# 80 이 회귀식에 들어가서 계산된 값을 도출함

# 예측구간 시각화
pp <- predict(dat.lm, interval="prediction", newdata=dat)
matlines(dat$newTest, pp, lty=c(1,2,3), col="blue")

# 신뢰구간 시각화
pc <-predict( dat.lm, interval = "confidence", newdata = dat)
matlines( dat$newTest, pc, lty=c(1,2,3), col="red")


# 잔차 시각화 => 잔차가 등분산성을 만족하는 지 확인해야 함
par(mfrow=c(1,2))
plot(resid(dat.lm))
# 잔차 값이 0을 중심으로 random 하게 퍼져 정규분포를 따름
# Patten 은 언제나 random 분산임

abline(1, 0, col="red")
qqnorm(resid(dat.lm))
qqline(resid(dat.lm), col="red")

shapiro.test(resid(dat.lm))
# 정규성 확인 =  p-value >0.05 


# 오차항의 등분산성을 만족하지 않음. 자기 상관관계(autocorrelation)가 존재
'''오차가 서로 독립이 아닌 경우 주로 독립변수가 시간을 나타내거나 
관측값이 관측순서에 영향을 받는 경우가 되므로, 오차항이 서로 독립이 아니면서
상관관계를 가지는 것을 ‘자기상관(autocorrelation)’이 존재한다고 한다.
'''


# 자기 상관관계를 갖는 지 확인 = Durbin-Waston 검정통계량
install.packages("lmtest")
library(lmtest)
# 귀무가설 = 자기상관관계는 0 이다.
# D의 값이 0과 근접하면 높은 양의 자기상관 관계를 가지며,
# D의 값이 4에 근접하면 높은 음의 자기 상관관계를 , 
# D의 값이 2와 근접하면 자기 상관관계가 존재하지 않는다
dwtest(standardTest ~ newTest, data=dat)

