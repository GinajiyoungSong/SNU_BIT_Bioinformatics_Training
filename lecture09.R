# 153 page
install.packages("ggplot2")
install.packages("reshape2")
library(ggplot2)
library(reshape2)

setwd("C:/Users/SJY/Documents/유전체 분석과정/lecture09")
getwd()
dir()


Samples <- read.table("Samples_detail.txt", header=T, stringsAsFactors = F, sep="\t")
dim(Samples)

head(Samples)

Tissue <- table(Samples$Tissue1)
head(Samples$Tissue1); Tissue # Table로 만들어서 어느 조직에서 나왔는지 빈도 정리

par(mar=c(10,3,3,3)) # margin 설정
barplot(Tissue, col=rainbow(19), las=2) # las 는 x축에 colname 방향


# 데이터프레임 생성 + 빈도 시각화
DF_tissue <- data.frame(Tissue)
colnames(DF_tissue) <- c("Tissue", "Freq")
head(DF_tissue)


p <- ggplot(DF_tissue, aes(x=Tissue, y=Freq, fill=Tissue)) + geom_bar(stat="identity")
p <- p + theme(axis.text.x = element_text(angle=-30, hjust=0))
print(p)



###
TCGA <- table(Samples$TCGA)
DF_TCGA <- data.frame(TCGA)
colnames(DF_TCGA) <- c("TCGA", "Freq")
head(DF_TCGA)

p <- ggplot(DF_TCGA, aes(x=TCGA, y=Freq, fill=TCGA)) + geom_bar(stat="identity")
p <- p + theme(axis.text.x = element_text(angle=-30, hjust=0))
p <- p + geom_text(data=DF_TCGA, aes(TCGA, Freq, label=Freq), position=position_dodge(width=0.9))
print(p)



# Cosmic Cancer Cell Line Project & Drug
Drugs <- read.table("DrugInfo.tsv", header = T, stringsAsFactors = F, sep="\t")
head(Drugs)

# 약물반응성, IC50, AUC 값 확인 = 낮을 수록 cell이 약물에 민감하다 판단
GDSC <- read.table("DoseResponse.txt", header=T, sep="\t")
dim(GDSC)
head(GDSC,5)

'''R에서는 reshape2라는 패키지에서 Data frame구조에서 matrix로 casting해주는 함수인 
acast라는 함수가 있기 때문에 쉽게 구조를 변경할수 있다'''
# 값(value)이 rowname - colname 으로 들어가게 만들어주는 함수

#        데이터, 행에 넣는 값 ~ 열 ,값으로 넣고자 하는 열,  같은 행열에 값이 다르면 mean으로
AUC <- acast(GDSC,COSMIC_ID~DRUG_ID, value.var = "AUC", fun.aggregate = mean)
IC50 <- acast(GDSC, COSMIC_ID~DRUG_ID, value.var = "LN_IC50", fun.aggregate = mean)
AUC[1:5, 1:5]
IC50[1:5, 1:5]



# ID 는 넘버이므로 실제 name 으로 변환시켜줘야 함
idx <- match(rownames(IC50), Samples$COSMIC_ID)
rownames(IC50) <- Samples[idx, "CellLines"]
idx <- match(colnames(IC50), Drugs$DRUG_ID)
colnames(IC50) <-Drugs[idx, "DRUG_NAME"]
IC50[1:5, 1:5]



# Significant Mutation Comparison By Drug Sensitivity
Mutation <- read.table("MutatedMatrix.txt", sep="\t", check.names = FALSE)
Expression <- read.table("GeneExpression.tsv", sep="\t", check.names = FALSE)

# missing value 없애기 위해 intersect 교집합 데이터만 추출
Names <- intersect(colnames(Mutation), colnames(Expression))
Names <- intersect(Names,rownames(IC50))
head(Names) # 3개의 column이 다 정보를 가지는 Cell line name

idx_Mut <- match(Names, colnames(Mutation))
idx_Exp <- match(Names, colnames(Expression))
idx_IC50 <- match(Names, rownames(IC50))
idx_Sams <- match(Names, Samples$CellLines)


# 원하는 Mutation idx 를 가지는 데이터만 추출
Mutation <- Mutation[,idx_Mut]
head(Mutation)
Expression <- Expression[,idx_Exp]
IC50 <- IC50[idx_IC50,]
CancerType <- factor(Samples$TCGA[idx_Sams])
TissueType <- factor(Samples$Tissue1[idx_Sams])
Mutation <- as.matrix(Mutation)
Expression <- as.matrix(Expression)
head(Mutation)





# Analysis
res <- lm(IC50[,"Afatinib"]~ Mutation["EGFR",])
summary(res)

res <- lm(IC50[,"Afatinib"]~ Mutation["EGFR",] +  TissueType) # 보정변수 추가
summary(res)

res <- lm(IC50[,"Gefitinib"]~ Mutation["EGFR",])
summary(res)

res <- lm(IC50[,"Gefitinib"]~ Mutation["EGFR",] +  TissueType) # 보정변수 추가
summary(res)

# data.frame 으로 만들어서 plot 화
DF <- data.frame(Sensitivity=IC50[,"Afatinib"], EGFR=factor(Mutation["EGFR",]), TissueType)
p = ggplot(DF, aes(x=EGFR, y=Sensitivity, fill=EGFR)) + geom_boxplot()+facet_wrap(~TissueType,ncol=5)
print(p)

DF <- data.frame(Sensitivity=IC50[,"Gefitinib"], EGFR=factor(Mutation["EGFR",]), TissueType)
p = ggplot(DF, aes(x=EGFR, y=Sensitivity, fill=EGFR)) + geom_boxplot()+facet_wrap(~TissueType,ncol=5)
print(p)

table(Mutation["EGFR",])
# plot 을 보면 Lung cancer에서 굉장히 큰효과를 보이는 것을 알 수 있다.

idx_Lung <- grepl("lung", TissueType)
res <- lm(Sensitivity ~ EGFR, data=DF[idx_Lung,])
summary(res)



# example
EGFR <- factor(Mutation["EGFR",]) # 유전자 돌연변이 0,1 (유/무)
result <- c()
for (i in 1:100) {
  Exp = Expression[i,] # 발현량
  res <- t.test(Exp ~ EGFR)
  result <- rbind(result, c(rownames(Expression)[i], res$p.value))
}
q<- p.adjust(result[,2], method='fdr')
result <- cbind(result, q)
result.sig <- result[which(result[,3] <= 0.05),]
