# 250 페이지 - Evaluation and Validation

# 수업에 필요한 패키지 install
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DAAG")
library(MASS)
library(DAAG)


# permutation test - r'''A permutation test is a type of statistical significance test in which the distribution of the 
test statistic under the null hypothesis is obtained by calculating all possible values of the 
test statistic under all possible rearrangements of the observed data points'''
 
# 실습 - 저체중아를 출산한 산모의 체중이 그렇지 않은 산모의 체중과 차이가 나는지를 permutation test를 통해 알아본다.
= birthwt[(birthwt[,1] == "0"), 3] # 3??� 3번째 col = 출산한 산모의 몸무게 = birthwt[(birthwt[,1] == "1"), 3]

to #  normal vs low wtbs# 두 그룹의 산모 체중에 대한 t 검정 값
 = t.test(nor, low)
tobs <- tobs$statistic

'''
????오픈 소스코드.biostat.wisc.edu/~kbroman/teaching/stat371/permfunc.R
'''
getwd()
source("C:/Users/woguSJYuments/perm유전체 분석과정/func.R.txt")

tpe# 두 그룹간의 permutation testrm = perm.test
pva(nor, low, n.perm=1000)lhist(tperm)

abline(v=abs(tobs), lty=2, col=2)

ue = mean(abs(tperm >= abs(tobs)))


#pvalue ??��Permutation test in R code ength(nor)
nlow = length(low)
n = nnor + nlow

tot = c(# Data re-compiled
tot = c(nor, low)
z = rep(0:1,c(nnor,nlow))

# do 1000 permutations of the data
n.perm = 1000
tobs = 1:n.perm
for(i in 1:n.perm) {
  z = sample(z)
  xn = tot[z==1]
  yn = tot[z==0]
  tobs[i] = t.test(xn,yn,var.equal=TRUE)$statistic
}
tobs



#### Cross validation
# 출생아 체중 단위 바꾸기 ( g -> kg ) 
bwt_kg = birthwt$bwt / 1000

# 산모 체중 단위 바꾸기 (lb -> kg)
lwt_kg = birthwt$lwt * 0.45 

str(birthwt)
# data type 바꾸기 (int -> factor)
smk = as.factor(birthwt$smoke)

mydata = data.frame(bwt_kg, lwt_kg, smk)
head(mydata) # 새로 생성된 data 확인하기

# 훈련용 데이터 만들기
train = mydata[-seq(1,180,by=20),]
attach(train)

# 시험용 입/출력 데이터 만들기
test_in = mydata[seq(1,180,by=20),c(2,3)]
test_out = mydata[seq(1,180,by=20),1]

# 출생아 체중 예측 모델 만들기 
model = lm(bwt_kg ~ lwt_kg + smk, data=train)
summary(model)


cross = CVlm(data = train, m = 10, form.lm = model)

# 시험용 입력 데이터로 신생아 체중 예측하기
p_out = predict(model, test_in)
p_out

# 시험용 데이터 실제 결과
test_out

plot(p_out, test_out) 
# 예측 결과와 실제 결과 비교
abline(0, 1, col='red', lwd = 0.5)

# bonferroni 보정
pbon = p.adjust(pvalue, method="bonferroni")

# FDR 보정
pfdr = p.adjust(pvalue, method="fdr")
