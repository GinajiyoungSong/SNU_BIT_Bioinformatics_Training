# 295페이지
install.packages("survminer")
install.packages("JM")
install.packages("ISwR")
library(survminer)
library(survival)
library(splines)
library(KMsurv)
library(lattice)
library(JM)
library(ISwR)

# survival analysis 분석 강의
''' 
시간에 따른 '사망'이나 '재발' 등의 변화를 관찰하는 분석
사건 - 생존분석에서 사망이나 재발과 같이 연구자가 관심을 가지고 있는 변화

Censored(중도절단)된 자료가 있는 것이 특징
Kaplan-Meier분석 - 특정 집단의 '생존율'을 추정
Log-rank test - 두 집단의 생존율이 같은지 아닌지를 '검정'
Cox 비례위험모형 - 생존율에 영향을 미치는' 위험분자'를 분석

가장 기본 : uncensored / censored data 를 구분해야함'''

# 시간값을 종속변수로 회귀분석 하면 중도절단된 데이터에 대한 고려가 안됨
# 생존여부를 종속변수로 분류분석 하면 생존기간이 다른 자료들이 동일하게 취급됨(시간 고려 없음)
# 즉 \시간 \절단데이터 고려가 필요


# Actuarial, Follow-up life table Approach
time.interval <- 0:5
Total.nSample <- 100
nDeath <- c(20,10,10,10,5)
nCensored <- c(0,18,6,4,17) 

lifetab(time.interval, Total.nSample, nCensored, nDeath)

''' 
lifetab() 함수가 데이터프레임 형태로 만들어줌
nsubs : 생존자수, nlost : 중도절단수, nrisk: 유효 인원수, 
nevent : 구간 이벤트 발생수, surv : 구간에서의 생존확률'''

fit.lifetab <-lifetab(time.interval, Total.nSample, nCensored, nDeath)
plot(time.interval[1:5], fit.lifetab[,5], type="l", xlab="Duration", ylab="Survival function life Table")





#### Kaplan-Meier Approach 비모수 생존함수 추정방법
'''
이는 매 시간마다 event를 측정해서 누적 확률을 계산한다. 
-> 기간 유효인원수를 구할 필요 없이 매 시점의 전체 인원으로 사망률을 계산한다.
사망률 = 사망자수 / 해당 시점 전체 인원수

중도절단을 그냥 이탈수로 전체수에서 제거해줌, 중도절단수는 사망률에 고려하지 않음'''

Status <- c(rep(1,55), rep(0,45)) # 1 = 사망 55명, 0 = 생존 45명
Time <- c(rep(1,20),rep(2,10),rep(3,10),rep(4,10),rep(5,5),
          rep(2,18), rep(3,6), rep(4,4),rep(5,17)) 
          # 각 달의 사망자, 중도절단수 

Sample.data1<-data.frame(Status,Time)
Surv(Sample.data1$Time,Sample.data1$Status)#censored data에는 +로 구분됨


fit1<-survfit(Surv(Time,Status)~1,data=Sample.data1)
summary(fit1) 

''' 
lifetab 과 같이 데이터프레임 형태
time : 시점, n.risk : 해당 시점의 전체 인원수
n.event : 이벤트 발생 인원수, survival : 생존 확률'''

plot(fit1,xlab="Time", ylab="Survival funtion") #standard error도 점섬으로 표시됨








#### 로그순위 검정 d log-rank test
''' 두 그룹으로 관찰된 자료가 있을떄 두 그룹의 생존함수 검정 - 검정통계량 계산

library(survival)
library(ISwR) 데이터 제공 '''


head(melanom)
?melanom # 데이터 기본정보 확인 - status 3: dead from other cause: censored중도절단

'''
Surv(time, event) 
- time : the follow up time
- event : the status indicator - 0 = alive / - 1 = dead / - Others.'''

head(Surv(melanom$days, melanom$status==1))

'''
10+ : the patient did not die from melanoma within 10 days
35+ : alive / 185 : the patient died from the disease
'''
surv.bysex <- survfit( Surv(days, status==1) ~ sex, data=melanom)
summary(surv.bysex) # sex = 1/2 성별로 나눠서 테이블


# Kaplan-Meier curve
plot(surv.bysex) # 생존률 그래프
plot(surv.bysex, conf.int=TRUE, col=c("black","gray"))

# whether two or more survival curves are identical
survdiff(Surv(days,Status==1)~sex, data=melanom)

'''
주어진 표에서 중도절단 데이터를 놓치지 않는게 중요하다
각 EVENT 발생 외에 다른 시점에서 데이터가 줄어들었다면 censored data
'''

Chemo.Group <- c(rep(1,10), rep(2,10))# Group vector 생성aChemo.Status <- c(rep(1,6), rep(0,4), rep(0,7), rep(1,3))# Status vector 생성.Month <- c(8,12,14,21,26,27,8,14,28,33,
                21,21,33,rep(41,4),28,33,41)
Chemo.data<- data.frame(Time.Month, Chemo.Status, Chemo.Group)
Surv(Chemo.diff<-survdiff(Surv(Time.Month,Chemo.Status)~Chemo.Group, data=Chemo.data)
Chemo.diff
'''
whether two or more survival curves are identical
p = 0.009 귀무가설 기각=> H1 대립가설 채택
H0 : 두 군의 생존 곡선은 같다
H1 : 두 군의 생존 곡선은 같지 않다'''


#서바이벌 데이터 제대로 만들었는지 확인
Surv(Chemo.data$ Time.Month, Chemo.data$ Chemo.Status)

.plot <- survfit(Surv(Time.Month, Chemo.Status)~Chemo.Group, data=Chemo.data)
plot(Chemo.plot, xlab="Time(Month)", ylab="%survival", col=c("red","blue"), main = "Survival Chemotherapy (n=20)")
legend("bottomleft", legend=c("Group1","Group2"), lty=c(1,1), col=c("red","blue")).S
### �








��중 로그순위 검정
surv (Weighted log-rank test)d'''
앞에서 로그 순위검정은 두 생존함수를 비교할 때, 모든 시간대에 동일한 가중치(=1)를 부여하였다. 
즉, 모든 시간대의 영향력이 동일하게 반영됨을 의미. 
실제 자료분석 시 두 생존함수의 차이가 초반에 두드러지는 경우, 이러한 차이를 
부각하기 위해서는 초반에 큰 가중치를 부여하는 것이 높은 검정력을 가져올 수 있다.

Survdiff 함수 내에서는 parameter rho가 p의 역할을 한다. 
default 값이 p=0으로 잡혀 있으므로 로그 순위 검정을 시행하게 된다.
 p = 1일 땐 Peto-Peto Gehan-Wilcoxon test를 시행하게 된다.
 p = -1일 땐 뒤로 갈수록 가중치를 주는 역할을 하게 된다

'''
iff(Surv(Time.Month, Chemo.Status)~Chemo.Group, data=Chemo.data, rho=-1) # -1 로그 > 시간이 지날수록 가중치 계속 올라감

Chemo.diff <- survdiff(Surv(Time.Month, Chemo.Status)~Chemo.Group, data=Chemo.data)
Chemo.diff # rho 값을 주기전보다 Chisq 값이 높음. p_value 가 작은것도 볼 수 있음



# JM패키지에서 pbc2.id 데이터 로드

str('''
해당 데이터는 Mayo Clinic에서 rare autoimmune liver disease인 원발성 담
즙성 간경병증(primary biliary cirrhosis) 환자 312명을 random하게 뽑아서 추
적 연구한 데이터이다.
'''pdata("pbc2.id",package="JM")

bc2.id)
pbc2.fit <- survreg(Surv(years,statu status2) ~ (drug + sex) * (age + I(age^2)), data = pbc2.id.pbbc2.fit

'''
모델에서는 age의 linear effect와 quadratic한 변수의 각각의 해석과
그 교호작용등의 해석이 어렵다는 단점이 있다. 이것들을 시각화 하는 것이 Effect plot이다.
'''
Combi.dataset <- with(pbc2.id, expand.grid( age = seq(30, 65, length.out 
                                                      = 25), sex = levels(sex), drug = levels(drug)))

pc2 <- p<- predict(pbc2.fit, Combi.dataset, se.fit=TRUE,type="lp"l
Combi.dataset$pred <- prs.pbc2[[1]]
Combi.dataset$se <- prs.pbc2[[2]]
Combi.dataset$lo <- Combi.dataset$pred - 1.96 * Combi.dataset$se
Combi.dataset$up <- Combi.dataset$pred + 1.96 * Combi.dataset$se
head(Combi.dataset)

xyplot(pred + lo + up ~ age | drug*sex, data = Combi.dataset, type = "l", 
       lty = c(1,2,2), col = "black", lwd = 2, xlab = "Age", 
       ylab = "Log Survival Time") #lattice 패키지의 xyplot





### Cox proportional hazards model
'''k사건 발생에 관계되는 인자가 하나일 때는 test도 충분하지만 실제로
는 둘 이상인 경우가 많아
생존�  영향을 미치는 다른 위험 인자들도 보정해야 한다.

시간�
치료방법에 따른 생존의 차이를 보고자 할 때 치료방법 뿐 아니라 나이나 성별, 
환자가 가지고 있는 질환이나 다른 위험 인자 들도 직간접적으로 생존에 영향을 미치
므로 이러한 변수(potential confounders)들을 보정해야 한다.

생존분석은 odds ratio와 유사한 hazard ratio를 다루므로 로지스틱 회귀분석의 알고리즘을
차용할 수 있으며 이러한 분석 방법을 Cox regression이라고 한다

시간에 관계없이 hazard ratio가 일정하다는 가정'''

coxph(Surv(years,status=="dead")~sex+age,data=pbc2.id)
'''앞선 모수회귀모형은 추정된 회귀계수가 공변량과 '생존시간'과의 관계라면 
비례위험모형에서는 공변량과 '위험함수'에 대한것이므로 
추정값이 '서로 반대부호'를 가지게 된다. 

- actual hazards ratio : exp(coef)

결과값을 해석하면 위험함수 즉, 사망률에 대한 것이므로 female의 사망 위험률이
0.6189 expection 감소 (coef - 마이너스) age 1살 증가시 1.04, 즉 4% 증가'''
fit.coxph<-coxph(Surv(years,status=="dead")~sex+age+drug,data=pbc2.id)
ggforest(fit.coxph,data=pbc2.id) #hazard ratio plot




### Stratified analysis
#  Compute the log-rank test for a gender effect stratified by ulceration(궤양)
coxph(Surv(days, status==1) ~ sex ,data=melanom)
str.plot<-survfit( Surv(days, status==1) ~ sex + strata(ulc),data=melanom )
plot(str.plot,col=c("red","blue"),lty=c(1,1,2,2))
legend("bottomleft",c("Female with present of ulceration",
                      "Female with absent of ulceration","Male with present of ulceration",
                      "Male with absent of ulceration"), lty=c(1,1,2,2))

'''남자는 병이 심각하게 진행되었을 때야 방문하는 경향이 있을 수 있다. 
따라서 병의 진행도를 통한 보정을 해주면 성별간의 차이가 줄어들 수 있을
것이다.'''