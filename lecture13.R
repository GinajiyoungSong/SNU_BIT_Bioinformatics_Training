# 오후 수업 269
# Pathway analysis and Network Visualization in Cancers

# 수업에 필요한 패키지 install
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

install.packages("fansi") # dependencies
install.packages("tidytree")# dependencies

BiocManager::install(c("TCGAbiolinks","ReactomePA","mygene","org.Hs.eg.db","clusterProfiler","pathfindR","EDASeq","edgeR"),force=TRUE)
install.packages(c("tibble","rlang"))


library(TCGAbiolinks)
library(ReactomePA)
library(mygene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tibble)
library(rlang)




'''이번 실습 : Breast cancer 환자의 예시 데이터를 통해 
DEG 테스트를 수행하여 gene set을 뽑아내고, 
추출한 유의한 유전자 군을 pathway 분석을 통해 해석하는 방법론을 실습해보고,
웹 기반의 툴을 이용하여 network을 효과적으로 visualization하는 것 까지 수행'''

#환자 10명의 RNA-Seq expression 데이터
head(dataBRCA)
str(dataBRCA) # factor 나 기타 변수타입에 대해 확인해야 함수실행이 가능함
dim(dataBRCA)




d_nor <- TCGAanalyze_Normalization(tabDF = dataBRCA, geneInfo=geneInfo)
''' normalization 결과 4step
첫 번째 스텝에서는 분석 가능한 expression set을 뽑고, 
두 번째, 세 번째 스텝을 통해 within, between normalization을 수행'''

head(d_nor)  # 발현량 정규화
'''TCGAanalyze_Normalization 함수를 이용하여 normalization을 수행한다.
이 함수는 추가적으로 분석 가능한 유전자군을 culturing 해주는 과정도 수행'''


# quantile filtering of genes
qf <- TCGAanalyze_Filtering(tabDF = d_nor, method = "quantile", qnt.cut = 0.25)


# Patient grouping - 환자 그룹핑 후 유의미한 유전자 확인
pat_normal <- TCGAquery_SampleTypes(barcode = colnames(qf), typesample = c("NT"))
pat_tumor <- TCGAquery_SampleTypes(barcode = colnames(qf), typesample = c("TP"))
pat_normal ; pat_tumor
# 정상과 암환자군의 그룹 구분 후 데이터 확인 - 두 집단간의 발현량 체크



#### DEG analysis 
'''옵션 중 mat1과 mat2에는 각각 filter된 데이터의 normal과 tumor 환자 바코드를 넣어 
개별적인 matrix를 인풋으로 넣는다. '''
dataDEG <- TCGAanalyze_DEA(mat1 = qf[,pat_normal], mat2 = qf[,pat_tumor],
                           Cond1type = "Normal", Cond2type = "Tumor", fdr.cut=0.01, logFC.cut = 1, method="glmLRT")
head(dataDEG)
''' 유전자 발현량의 유의성 확인
각각 5명의 normal과 tumor 샘플이 분석되었고, 
15,236개의 분석결과가 추출된 것을 확인할 수 있다'''


# DEG table with expression values in normal and tumor samples
''' 분석결과를 테이블로 저장하기 위해서 tumor, normal 라벨을 붙인다
개별 유전자들의 이름과 fold-change, FDR값 등이 차례로 나열되어 있는 것을 확인'''
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEG, "Tumor", "Normal",qf[,pat_tumor], qf[,pat_normal])
head(dataDEGsFiltLevel)


# 유의한 10개 유전자 뽑아내기
# order 명령어를 통해 특정 column에 대해 내림차순 또는 오름차순으로 조건
genelist <- head(dataDEGsFiltLevel[order(dataDEGsFiltLevel$FDR),], 10)[,1]
head(genelist) #  FDR 기준값으로 가장 높은 유전자list


### Pathway Analysis

# Gene symbol conversion
'''우리가 사용하는 package의 함수에 따라 인풋의 형태가 달라질 수 있는데 
이럴 경우 gene symbol conversion을 해주어야 한다. 
따라서 mygene이라는 패키지를 이용하여 변환하는 것을 실습한다'''
genelist <- as.vector(genelist)

## conversion from hgnc to entrez 
# mygene 패키지 내 queryMany 함수를 이용하여 변환
a_gene <-queryMany(genelist, scopes = "symbol", fields = "entrezgene", species = "human")
a_gene # data.frame 형태



# 변환된 entrez gene id 추출 
sig_gene <- a_gene$entrezgene



## pathway enrichment test ###
'''가장 기본적으로 해볼 수 있는 분석은 유의하게 뽑혀진 유전자군이 
어떤 pathway에 enrich 되어 있는가에 대한 분석일 수 있다. 
패키지에서 제공하는 함수를 통해 enrichment test를 시행한다.
결과는 패킹되어있으므로,
as.data.frame과 head 명령어를 통해 대략적으로 살펴볼 수 있다.'''

enrich_result <- enrichPathway(gene = sig_gene, pvalueCutoff = 0.05, readable = T)
head(as.data.frame(enrich_result))
# 맵핑된 pathway에 대한 정보와 내가 가진 유전자들 중 어떤 유전자들이 맵핑되었는지에 대한 통계적 정보들
# pathway visualization 
barplot(enrich_result, showCategory = 5)
dotplot(enrich_result, showCategory = 5)



###  Gene set - cluster analysis
# cluster visualization
require(clusterProfiler)


data(gcSample)
str(gcSample) # 유전자 cluster set
# 8개의 클러스터에 200여개에서 900여개까지의 유전자들이 모여 있는 dataset
# gene id는 entrez gene id로 맵핑
res <- compareCluster(gcSample, fun="enrichPathway")
dotplot(res)


data(geneList, package = "DOSE")


'''GSEA분석은 gsePathway라는 함수를 통해 수행한다. 
여러 파라미터 값들을 조정하여 결과를 얻을 수 있다.'''
gs_enrich_result <- gsePathway(geneList, minGSSize = 120, eps = 0, 
                               pvalueCutoff = 0.2, pAdjustMethod = "BH", verbose = FALSE)
gs_enrich_result <- as.data.frame(gs_enrich_result)
head(gs_enrich_result)







##### Network visualization
# Pathway cluster visualization - pathfindR 패키지를 사용
library(pathfindR) # java 버전이 맞지 않으면 실행이 안됨
data("RA_output")
head(RA_output)

RA_clustered <- cluster_enriched_terms(RA_output)
RA_clustered

# 데이터 프레임을 2줄만 정렬하여 출력 ###
knitr::kable(head(RA_clustered, 2))


# http://genemania.org/  - Web-based network analysis
'''네트워크 기반의 분석 및 시각화와 관련해서는 많은 툴들이 있는데, 
그 중에서도 많은 정보를 제공해주는 웹 기반의 genemania 툴을 소개한다.

검색 창에 www.genemania.org를 입력하여 접속
사용방법은 왼쪽 상단에 있는 검색창에 유전자리스트들을 넣으면 가능하다
예제 292페이지 참조'''