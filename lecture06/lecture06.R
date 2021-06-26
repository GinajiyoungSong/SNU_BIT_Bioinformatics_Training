### 1. BioManager 교재 118페이지

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("affy", force = TRUE) # for analysis affy chip

# https://www.bioconductor.org/packages/release/bioc/html/affy.html


# Microarray 데이터 분석에 필요한 패키지 설치
library(affy) # 분석에 필요한 전처리 패키지


# 파일 위치 setting
getwd()
setwd("C:/Users/SJY/Documents/유전체 분석과정/lecture06") 


# Microarray 실험파일은 CEL 파일로 제공됨

'''실습에서 사용할 데이터는 TCGA(The Cancer Genome Atlas)의 데이터 중 Glioblastoma Multiforme 샘플을 
Affymetrix HT_HG-U133A 칩을 이용하여 만들어진 데이터로써 전체 샘플은 6개로 구성되어 있다.'''

d = ReadAffy()			# .CEL 파일 읽기

dl = log2(exprs(d))	# Expression table 불러오기
image(dl)	



### 2. Image plot
par(mfrow=c(2,3))
image(d)




####### Affy 패키지의 preprocessing - Quantile normalization 
q = as.matrix(read.table("quantile.txt", header=TRUE, row.names=1, sep='\t'))
## sample 데이터 불러오기

ncol(q)
dim(q)


# sort된 유전자발현량을 담을 매트릭스 만들기
qs = matrix( ncol=ncol(q), nrow=nrow(q))
dim(qs)

qr = matrix( ncol=ncol(q), nrow=nrow(q))  
# sample 내에 유전자 발현량 순위를 담을 벡터 만들기


# for문으로 sorting 매트릭스 및 rank로 매핑된 매트릭스 만들기
for(i in 1:ncol(q)){
  qs[,i] = sort(q[,i])
  qr[,i] = rank(q[,i])
}


# mean 값 만들기
qm = apply(qs, 1, mean) # row기준 = 1
qm

# quantile normalization 결과 값을 담을 벡터 만들기
qn = matrix( ncol=ncol(q), nrow=nrow(q) ) 

# r 인스턴스 = ranking 값 을 mean에 매칭
for(i in 1:length(qr)){
  r = qr[i]
  qn[i] = qm[r]
}

# sample 데이터 column명, row명 할당
rownames(qn) = rownames(q);colnames(qn) = colnames(q)
qn

### image
par(mfrow=c(1,2))
boxplot(q, main="Before normalization")
boxplot(qn, main="After normalization")


 
# quantile normalization 이외에 RMA, MAS5 등 다양한 preprocessing 함수들을 제공

####### Affy 패키지의 preprocessing - RMA
BiocManager::install("DBI",force = TRUE)
BiocManager::install("RSQLite", force = TRUE)
BiocManager::install("hthgu133acdf", force = TRUE)


# RMA 로 background correction
d_rma = rma(d)
dr = exprs(d_rma) # RMA 결과 테이블 뽑기



####### Affy 패키지의 preprocessing - mas5
# MAS5 background correction & expression table
d_mas5 = mas5(d)
dm = exprs(d_mas5)
head(dm) # sample file G1 ~ G6 발현량 확인


####### Affy 패키지의 preprocessing - Expresso 실행
d_es = expresso(d, bgcorrect.method = "none",
                normalize.method = "quantiles",
                pmcorrect.method = "pmonly",
                summary.method="medianpolish")

des = exprs(d_es)  #Expresso expression table
head(des)



### plot
par(mfrow=c(1,2))
plotDensity(dl)# density pl
ot for before normalization
plotDensity(des)# density plot for after normalization



###### box plot ############
boxplot(dl)  # boxplot for before normalization
boxplot(des) # boxplot for after normalization

boxplot(dl ,col=rainbow(15), main="before norm") # boxplot for before normalization
boxplot(des ,col=rainbow(15), main="after norm") # boxplot for after normalization



###### scatter plot ############
par(mfrow=c(1,1))
pairs(dl)		# scatter plot before normalization
pairs(des) 	# scatter plot after normalization

pairs(dl, lower.panel=panel.smooth,upper.panel = panel.smooth)
pairs(des, lower.panel=panel.smooth,upper.panel = panel.smooth)

###### MA plot ############
mva.pairs(dl)  # MA plot  before normalization
mva.pairs(des) # MA plot  after normalization
