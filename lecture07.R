## 1. package installation ###
# 133
# https://www.bioconductor.org/packages/release/bioc/html/affy.html

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("affy","hgu133a.db","samr"))
BiocManager::install("affy", force=TRUE)


## 2. package loading ###
library(hgu133a.db) #  probe를 gene symbol로 변환하기 위한 “hgu133a.db” package,
library(affy) # Affymetrix사의 chip을 가공하기 위한 “affy” package
library(samr) # SAM을 수행하기 위한 “samr” package


'''실습에서 사용할 데이터는 TCGA (The Cancer Genome Atlas)의 데이터 중 Glioblastoma
Multiforme 샘플을 Affymetrix HT_HG-U133A 칩을 이용해 만든 데이터로써, 전체 샘플은 6개로
Glioblastoma Multiforme 3개와 Normal 3개로 구성되어 있다.'''


## 3. set working directory ###

getwd()
setwd("C:/Users/SJY/Documents/유전체 분석과정/Lecture07/Data")

## 4. file loading ###
# after-normalization data
no = as.matrix(read.table("af_no.txt", header = TRUE, sep='\t'))
a = c(1, 1, 1, 2, 2, 2)
g = factor(a)
t.test(no[1,] ~g) # sample t-testing of 1 row


tp = vector()# p-value 담을 백터
tf = vector()# 추정값 담을 백터


#  for 문을 이용하여 전체 유전자에 대한 T-test를 수행
'''for 문이 한번 수행될 때 마다 한 개의 유전자에 대해서 T-test가 수행되며 
이 결과 값으로 p-value는 tp vector에, 추정 값은 tf vector에 저장'''

for (i in 1:nrow(no)) {
  tmp = t.test(no[i,] ~g, paired=FALSE)# t-test 수행
  tf[i] = (log2(tmp$estimate[1]/tmp$estimate[2]))# fold 차이 값
  tp[i] = tmp$p.value
}


id = which(tp < 0.05) # p-value < 0.05 값을 가진 probe의 위치
head(tp[id])   # 조건에 맞는 결과만 추출
length(tp[id]) # p-value < 0.05 값을 가진 probe 갯수 확인인





### Multiple testing correction (FDR)
'''Microarray 데이터에서 차별 발현 유전자 추출 시, 한 번에 많은 수의 유전자 집합에 대한 검사를 
시행하기 때문에 다중 검사 보정 (multiple testing correction)이 반드시 필요하다. 
따라서 이번 실습에서는 앞서 수행한 T-test의 결과로부터 FDR을 수행해보자. 
FDR을 먼저 수행한 후 adjusted p-value값이 0.05 미만인 유의한 유전자를 찾고 수행결과를 확인해보겠다.
'''
pfdr = p.adjust(tp, method="fdr") 
f.id = which(pfdr < 0.05) # adjusted p-value < 0.05 값을 가진 probe의 위치
head(tf[f.id])
length(tf[f.id]) # adjusted p-value < 0.05 값을 가진 probe의 갯수 확인


# t-test와 FDR 을 통해 찾은 유전자 = 차별 유전자를 도식화함
# volcano plot

plot(tf, -log10(tp))
points(tf[f.id], -log10(tp)[f.id], col="red")


###### DEG 결과 파일로 저장
f.no = cbind(no, tp, tf, pfdr)
  # expression matrix, T-test p_value, 
  # T-test estimation value, adjusted p-value 테이블 생성

dep = f.no[f.id,] # adjusted p-value < 0.05 조건의 probe 추출
head(dep)
dim(dep) # DEG 갯수 확인

write.table(dep, "ttest_fdr.txt", sep = '\t', row.names = T)




####### Probe를 Gene symbol로 변환하기

# hgu133a.db는 Affymetrix HT_HG-U133A chip 전용 annotation package

pid = row.names(no)[f.id] # adjusted p-value < 0.05 값을 가진 probe의 명칭만 추출

# annotation package로부터 gene symbol-probe 이름 관계 추출
gn = unlist(mget(pid[!is.na(pid)], hgu133aSYMBOL))

deg = cbind(dep, gn) # gene symbol 붙이기
head(deg)

degs = deg[order(deg[,9]),] # adjusted p-value가 낮은 순으로 배열
write.table(degs, "deg.txt", sep = '\t', row.names = T) # 결과 파일 작성





##### Significance Analysis of Microarray (SAM)

'''SAM은 데이터 배열에 기초한 DNA microarray 실험에서 위양성을 평가하기 위해 만들어졌다. 
SAM에서 q-value는 마치 T-test에서 p-value와 같은 의미를 갖게 된다. 
다시 말해 q-value가 작을수록 유전자 발현의 차이가 유의하다는 것을 의미하게 된다. 
먼저 "samr" packae를 로딩한 후 원하는 그룹에 따라 각 그룹에 라벨을 붙여준다.'''

# sam input matrix 를 만들고 sam test 수행
sm = list(x=no, y=a, logged2=TRUE)

# SAM test
st = samr(sm, resp.type = "Two class unpaired", nperms = 100)

# delta table 작성
dt = samr.compute.delta.table(st)
head(dt)


# FDT < 0.05 로 만들기 위해 delta 1.30 설정 ( 파라미터 조절)

d = 1.30 # 역치값
samr.plot(st, d) 

st = samr.compute.siggenes.table(st, d, sm, dt) # result table 작성
names(st) # 테이블 결과 element 확인

head(st$genes.up)
write.table(st$genes.up, "sem.txt", sep='\t', row.names = T)
