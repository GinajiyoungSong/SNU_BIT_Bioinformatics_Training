# install 설치 
''' Case study
I. Association of BRCA1 and BRCA2 mutations with survival in ovarian cancer, 2011, JAMA
II. Integrated genomic analyses of ovarian carcinoma, 2012, Nature - 333 페이지'''


##########   CASE I  ##########
##### 1. CGDS-R package install

install.packages("cgdsr")
library(cgdsr)

''' 데이터 http://www.snubi.org/R_semi
nar/case_study.txt
nature 논문데이터 - TCGA 데이터들은 GDC 데이터포톨(http://portal.gdc.cancer.gov/)에서 다운로드가능 
실제 TCGA data access 와 preprocessing에 사용되는 CGDS-R package를 이용하여 
TCGA데이터를 손쉽게 이용할 수 있는 방법을 익힌다. (CGDS: Cancer Genomic Data Server)
'''

##### 2. Get list of cancer studies at server
mycgds = CGDS("http://www.cbioportal.org/")
test(mycgds)

cancerstudy <-getCancerStudies(mycgds)
head(cancerstudy)
cancerstudy$name
# Ovarian cancer 저널 데이터 로드
mycancerstudy <- cancerstudy[210,1]



##### 3. Extract samples and features

getCaseLists(mycgds,mycancerstudy)[,1] # 환자리스트 로드

mycaselist <- getCaseLists(mycgds, mycancerstudy)[2,1]
mycaselist

# mutation, expression, methylation profile 각각 불러오기
mymutationprofile <- getGeneticProfiles(mycgds,mycancerstudy)[6,1] 
mymethylationprofile <- getGeneticProfiles(mycgds,mycancerstudy)[5,1] 

'''
Ovarian cancer의 예후에 BRCA1, BRCA2 gene의 mutation과 methylation 여부가 
영향을 미치는지를 확인하기 위해 overall survival과 비교하고자 한다.
따라서 전체 환자군을 BRCA1 mutation group, BRCA2 muation group, 
BRCA1 methylation group, Wild type으로 나누어 분석을 수행한다.
'''


##### 4. Get mutation profiles for BRCA1 and BRCA2
brca_mutation = getMutationData(mycgds, mycaselist, mymutationprofile, c('BRCA1','BRCA2') ) 
head(brca_mutation) #  BRCA1과 BRCA2 에 mutation 이 있는 환자들을 먼저 추출

table(brca_mutation$gene_symbol)  #BRCA1,2 mutation을 가진 환자 수 확인

# mutation 그룹 분류
brca1_mutated_cases <- brca_mutation[which(brca_mutation$gene_symbol=='BRCA1'),3]
brca2_mutated_cases <- brca_mutation[which(brca_mutation$gene_symbol=='BRCA2'),3]
head(brca1_mutated_cases)
head(brca2_mutated_cases)

# Extract samples with BRCA1, BRCA2 methylation

'''BRCA1과 BRCA2 gene에 methylation이 된 경우 
mutation과 마찬가지로 해당 gene의 expression이 낮아지게 되고, 
이것이 생존률에 영향을 줄 수 있다. 
따라서 두 gene에 methylation 된 환자들을 각각 추출한다.
'''
brca_methylation = getProfileData(mycgds, c('BRCA1','BRCA2'),mymethylationprofile, mycaselist)
head(brca_methylation)

# #BRCA1 gene 의 expression과 anti correlation이 0.8 이상인 methylation이 있는 경우만 추출
brca1_methylation_cases = rownames(brca_methylation[which(brca_methylation$BRCA1>0.8),])
brca2_methylation_cases = rownames(brca_methylation[which(brca_methylation$BRCA2>0.8),])


##### 5. Clinical data integration
'''BRCA1, 2 gene에 각각 mutation이 있거나 methylation이 있는 환자들을 뽑았다. 
다음 단계로 생존분석을 하기 위해 임상 변수 중 생존 기간 변수를 추출한다.
'''
myclinicaldata = getClinicalData(mycgds,mycaselist)
head(myclinicaldata)

myclinicaldata$OS_STATUS[myclinicaldata$OS_STATUS == ""] <- NA

head(myclinicaldata$OS_STATUS)

myclinicaldata$OS_MONTHS


# BRCA gene 에 mutation이나 methylation이 있는 환자와 그렇지 않은 환자를 구분한다.
total_sample <- rownames(myclinicaldata) 
type <- rep('Wild', length(total_sample))  #모든환자를 wild type간주. base type table을 만든다.
names(type) <- total_sample
head(type)

'''환자 이름 일치시키기: mutation profile 의 환자명과 clinical, methylation
profile 의 환자명 표기가 다르기 때문에 이를 일치시킨다.'''
brca1_mutated_cases = gsub("-",".",brca1_mutated_cases) 
brca2_mutated_cases = gsub("-",".",brca2_mutated_cases)
head(brca1_mutated_cases)

# 각 mutation profile 에 따라 환자를 분류한다
type[brca1_methylation_cases] <- "BRCA1_methylation"
type[brca2_methylation_cases] <- "BRCA2_methylation"
type[brca1_mutated_cases] <- "BRCA1_mutation"
type[brca2_mutated_cases] <- "BRCA2_mutation"

type <- type[names(type) %in% rownames(myclinicaldata)]
type <- factor(type, levels= c("Wild","BRCA1_mutation","BRCA2_mutation","BRCA2_methylation") )




##### 6. Survival Analysis
'''BRCA mutation과 methylation profile 에 따라 분류된 4군에 대해 생존분석을 수행한다'''
install.packages('survival')
library(survival)


out <- survfit(Surv(OS_MONTHS, OS_STATUS=="1:DECEASED")~type, data=myclinicaldata)
# type별 분류 count

survdiff(Surv(OS_MONTHS, OS_STATUS=="1:DECEASED") ~ type, data=myclinicaldata)
coxph(Surv(OS_MONTHS, OS_STATUS=="1:DECEASED") ~ type, data=myclinicaldata)

#시각화
color <- c("black","Skyblue","Blue", "Red") # 분류 color로 나타냄
plot(out, col=color , main="Association of BRCA1/2 Mutations with Survival", xlab="Time,days", ylab="Proportion", lty=1:4, lwd=2)
legend("topright", levels(type), col=color, lty=1:4, lwd=3) 


# survminer ggplot 

install.packages("survminer")
library(survminer)
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/survminer")

ggsurvplot_res <- ggsurvplot(out, data=myclinicaldata, pval=T, risk.table=T, 
                             conf.int=F, break.time.by = 30, 
                             legend.title = "Patient types", 
                             risk.table.fontsize = 2.5,
                             surv.plot.height = 0.45, legend="right")
ggsurvplot_res




# microarray를 수행한 후 이를 z-normalization한 데이터
### 데이터 로드
getwd()
setwd('C:/Users/SJY/Documents/유전체 분석과정/lecture15/')
expression <- read.csv(header=TRUE,'expression.csv')
dim(expression)
length(unique(rownames(expression)))


# Check variablility of gene expression across patients, using the MAD
''' 환자들간에 유전자 발현양의 변화가 비교적 큰 유전자를 고르기 위해
Median Absolute Deviation (MAD) 방법을 이용한다.'''

myMad <- function(x){mad(x,na.rm=T)}    
result <- apply(expression, 1, myMad)

# Get the highest variability genes
'''MAD에서 상위에 랭크되는 유전자를 추출한다. 
논문에서는 1500개의 유전자를 추출했으나, 
본 실습에서는 시간상 100개로 줄여서 분석을 수행한다.'''
result2 <- sort(result, decreasing = T)[1:100]
exp_result <- expression[which(rownames(expression)%in%names(result2)),]



# NMF clustering  - Non-negative matrix factorization clustering package 불러오기
# 음수 미포함 행렬 분해 -> 정확하지 않지만 대략적인 해를 구하게 된다
install.packages("NMF")
library(NMF)

# Make the values to non-negative
# NMF의 경우 negative value를 허용하지 않는다. 모든 값을 non-negative 로 전환한다.

nonNegativeFx <- function(x){a <-min(x)
                          if(a<0){
                              a <- a*(-1)
                              b <- x + a}
                          else {
                              b <- x  }
                return(b)}

res2 <- apply(exp_result, 1, nonNegativeFx)
res2 <- as.data.frame(res2)

# NMF 분석을 수행하여 유의한 clustering 개수를 찾는다.
nmf_res <- nmf(res2, 2:6, nrun=10, seed=123456)
plot(nmf_res)

consensusmap(nmf_res, annCol = colnames(res), labCol = NA, labRow = NA)

# Draw heatmap and NMF analysis
'''NMF package에서 제공하는 heatmap 기능으로 
유전자 발현량에 따라 환자군이 clustering되는 패턴이 있는지 확인한다.'''
aheatmap(res2, Rowv="correlation", hclustfun = "complete")

