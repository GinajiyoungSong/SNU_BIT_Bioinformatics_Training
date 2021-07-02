# 6월 26일 토요일 강의
install.packages('reshape2')
install.packages('dplyr')

# wide & long format 수정
''' long format을 수정해야 하는 이유는 model function 을 활용할 때 long format 읽음

ANOVA test 실습
분산분석: 3개 이상 다수의 집단을 비교할 때 사용하는 가설 검정 방법
집단간 분산/집단내 분신 기반의 F분포 이용'''

library(reshape2)
dat.wide <- read.table(header=TRUE, text='
                      age sex condA condB condC
                      20s M 7 9 10
                      20s F 8 9 10
                      30s F 8 9 10
                      30s F 7 5 8')
melt(dat.wide, id.var=c("age","sex")) # melt()함수 사용

# 2개이상factor가 있는 변수를 하나로 합칠 때 variable.name
dat.long <- melt(dat.wide, id.var=c("age","sex"), variable.name = "condition", value.name = "score")

result <- aov(score ~ sex + condition, dat.long) # longformat으로 anova test
summary(result)
'''score에 어떤 columns이 영향을 미치는 지 확인
sex 와 condition factor의 영향을 확인하기 위해 ANOVA test에서 p_value확인
0.05 미만의 p_value가 없으므로 score 에 두 factor 모두 영향을 끼치지 않음'''

df.wide <- read.table(header = TRUE, text='
                      ageGroup methodA methodB methodC
                      <20 7 9 10
                      20-29 8 9 10
                      30-39 9 9 12
                      40-49 10 9 12
                      >50 11 12 14')
df.long <- melt(df.wide, id.var='ageGroup', variable.name = 'method', value.name = 'score')

result <- aov(score ~ ageGroup + method, df.long)
summary(result) # p_value 확인




# 반대로 long format -> wide format으로 바꿔주는 함수: dcast()
dcast(dat.long, age + sex ~ condition, value.var="score")

# funtion 옵션 사용시 -> 

#age(성별통합 -명시하지 않은 col이 통합됨). condition은 function 적용.  
dcast(dat.long, age ~ condition, mean, value.var = "score")
#sex(나이통합)
dcast(dat.long, sex ~ condition, mean, value.var = "score")
#condition
dcast(dat.long, age ~ sex, mean, value.var = "score")



# data.frame - col selection, filterting, groupby etc
library(dplyr)
search()
install.packages('MASS')
library(MASS)

head(birthwt)
''' Risk factors associated with low infant birth weight,
the birthwt dataframe has 189 rows and 10 columns.
the data were collected at Baystate Medical Center'''
birthwt_df <- as.data.frame(birthwt)

#dplyr 패키지에서 selection 옵션션
head(dplyr::select(birthwt_df,low,age,smoke,bwt))

test1 <- dplyr::select(birthwt_df, -c(ptl,ht,ui,ftv))
head(test1)

test2 <- dplyr::select(birthwt_df, low:smoke)
head(test2)


#filtering
test3 <- filter(birthwt_df, age==21, smoke==1)
head(test3)

#mutation - 새로운 계산식으로 새로운 column만듦
test4 <- dplyr::mutate(birthwt_df, bwt_kg=bwt/1000)
head(test4)

head(dplyr::arrange(birthwt_df,age, bwt))
test5 <- dplyr::arrange(birthwt_df, desc(age), desc(bwt))
head(test5)

# Group_by, summarise() / Chain operator %>%
summarise(birthwt_df, bwt_mean=mean(bwt, na.rm=TRUE), bwt_median=median(bwt))

group_smoke <- group_by(birthwt_df, smoke)
head(group_smoke)

summarise(group_smoke, count=n(), bwt_mean=mean(bwt), bwt_median=median(bwt))
birthwt_df %>% filter( age>=21 & age<=30) %>% group_by(race,smoke) %>% summarise(count=n(), bwt_mean=mean(bwt, na.rm=TRUE))


# preprocess data
bwt <- with(birthwt, {
    race <- factor(race, labels=c('white','black','other'))
    ptd <- factor(ptl>0)
    ftv <- factor(ftv)
    levels(ftv)[-(1:2)] <- "2+"
    data.frame(low=factor(low), age, lwt, race, smoke=(smoke>0),
               ptd, ht=(ht>0), ui=(ui>0), ftv)})
head(bwt)

# glm() 모델 사용
bw.glm <- glm(low~., family=binomial, data=bwt)
round(summary(bw.glm)$coef,2)


# 실습 문제 1.
dat.wide <- read.table(header = TRUE, text='
                  id population1 population2 population3
                  1 30 29 33
                  2 27 27 32
                  3 28 24 29
                  4 26 26 37
                  5 29 25 30
                  6 27 26 38
                  7 29 27 35')

dat.long <- melt(dat.wide, id.var='id', variable.name = 'population', value.name = 'bmi')
result <- aov(bmi ~ population, dat.long)
summary(result)
# p_value 확인 -population이 유의미함 - 어떻게 다른 군인지 사후분석 필요

# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/TukeyHSD
tukey.test <- TukeyHSD(result, ordered=FALSE, conf.level =0.95)
plot(tukey.test)
