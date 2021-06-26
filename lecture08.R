# 141
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("golubEsets")
BiocManager::install("genefilter")
BiocManager::install("e1071")


# http://bioconductor.org/packages/release/data/experiment/html/golubEsets.html
# exprSets for golub leukemia data
library(golubEsets)
library(Biobase)
library(genefilter)

# 핵심데이터 microarray exp 수치, ALL/AML 진단데이터
data("Golub_Train")
data("Golub_Test")
data("Golub_Merge")
ls()

View(Golub_Merge)
head(Golub_Test@assayData[["exprs"]])

# Train - Test Split 패키지에서 이미 split 시켜서 제공

# 전처리

# 1. 자료 수치 변환 : 발현 수치 log화 
# 2. Feature selection : anova test를 통해 유의성이 높은 gene들만 뽑는다
# 3. Normalization : expression 수치 정규화


# 1. exp 수치 변경 함수 생성
GolubTrans <- function(eSet){
  x <- exprs(eSet)
  x[x<100] <- 100
  x[x>16000] <- 16000
  x <- log2(x)}

gTrn <- GolubTrans(Golub_Train)
gTest <- GolubTrans(Golub_Test)

gMerge <- GolubTrans(Golub_Merge)

# 2. 데이터 filtering
mmfilt <- function(r=5, d=500, na.rm=TRUE){
  function(x){
    minval <- min(2^x, na.rm=na.rm)
    maxval <- max(2^x, na.rm=na.rm)
    (maxval/minval > r) && (maxval-minval > d)
  }
}
mmfun <- mmfilt() # 함수 instance

ffun <- filterfun(mmfun)
sub <- genefilter(gTrn, ffun)
sum(sub)
'''Filter 되고 남은 유전자의 개수는 3,054개 이다. 
sub 변수에 TRUE만 되는 유전자만 따로 뽑아 자료행렬을 만든다.'''

gTrnS <- gTrn[sub,]
gTrnS <- gTrn[sub,]
gTestS <- gTest[sub,]
gMergeS <- gMerge[sub,]


head(sub) # 결과보면 False 값들 확인됨

# output 변수 y 정리
Ytr <- Golub_Train$ALL.AML
Ytest <- Golub_Test$ALL.AML
Ymerge <- Golub_Merge$ALL.AML


summary(Ytr)


# 2. feature selection - Anova test
af <- Anova(c(rep(1,27), rep(2,11)), 0.00001) # gene filtering 유의확률 le-5
anovasub <- sub[sub==TRUE] # True인 값들만 추출해냄

for (i in 1:sum(anovasub)) {
  anovasub[i] <- af(gTrnS[i,])
}

sum(anovasub) # 68개 유의한 유전자 추출완료


gTrA <- gTrnS[anovasub,]
gTestA <- gTestS[anovasub,]
gMeA <- gMergeS[anovasub,]

# Median Absolute Deviation (MAD) 값이 0인 유전자를 제거하는 함수를 만들고
# apply함수를 통해 전체 유전자를 필터링하여 새로운 변수에 저장한다.
whBad1 <- (apply(gTrA, 1, mad) == 0)
whBad2 <- (apply(gTestA, 1, mad) == 0)
whBad <- whBad1 | whBad2
sum(whBad)

gTrA <- gTrA[!whBad, ]
gTeA <- gTestA[!whBad, ]
gMeA <- gMeA[!whBad, ]

dim(gTrA)
dim(gTeA)


# 3. normalization

star <- function(x)( x- median(x))/mad(x)
TrExprs <- t(apply(gTrA, 1, star))
TeExprs <- t(apply(gTeA, 1, star))
MeExprs <- t(apply(gMeA, 1, star))

# Modeling
library(class)
knn1 <- knn(t(TrExprs), t(TeExprs), pData(Golub_Train)$ALL.AML, k=1)
table(knn1, Golub_Test$ALL.AML)

knn3 <- knn(t(TrExprs), t(TeExprs), pData(Golub_Train)$ALL.AML, k=3)
table(knn3, Golub_Test$ALL.AML)

knn5 <- knn(t(TrExprs), t(TeExprs), pData(Golub_Train)$ALL.AML, k=5)
table(knn5, Golub_Test$ALL.AML)




knn3.cvpreds <-  knn.cv(t(MeExprs), Ymerge, k=3) # CV  시행
table(knn3.cvpreds, Ymerge)

knn5.cvpreds <-  knn.cv(t(MeExprs), Ymerge, k=5) 
table(knn5.cvpreds, Ymerge)






# Model - SVM 
install.packages("e1071")
library(e1071)
model <- svm(t(gTrA), Golub_Train$ALL.AML, type = "C-classification", kernel="linear")

trpred <- predict(model, t(gTrA))
sum(trpred != Golub_Train$ALL.AML)
table(trpred, Golub_Train$ALL.AML)

'''C-SVM 방법으로 training dataset을 10 fold cross validation 하여 나은 classifier를 찾는다. 
summary 함수를 이용하여 다음의 모델링된 결과를 볼 수 있다.'''
cv_model <- svm(t(gTrA), Golub_Train$ALL.AML, type = "C-classification", 
                kernel="linear", cross=10)
summary(cv_model)

tepred <- predict(cv_model, t(gTeA))
sum(tepred != Golub_Test$ALL.AML)
table(tepred, Golub_Test$ALL.AML)



# LDA - Linear Discriminant Analaysis 
'''관찰 개체가 가지고 있는 여러 특징들을 모아 그 값을 보고 어떤 집단에 속하는 것이지 
통계적으로 판단하는 기법으로 연속형인 독립변인을 통해 범주형 자료를 판별 혹은 예측하기 위한 통계기법이다.
LDA는 집단 구성원 개개인의 적합 여부를 판별하고 분류해 주기 때문에 
일정한 예측변인에 따라 개인의 특성을 파악 할 수 있는 장점이 있다.'''

library(MASS)
gTr.lda <- lda(t(TrExprs), Ytr)
plot(gTr.lda)
preds.lda <- predict(gTr.lda, t(TeExprs))
table(preds.lda$class, Ytest)
