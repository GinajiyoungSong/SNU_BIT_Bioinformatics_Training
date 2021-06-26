# EDA 실습 R Graphic
y <- 10:50
plot(y)

y <- 1:5
# main = 그래프 제목, sub= 바닥에 추가 문구
plot(y, main="This is\nmain", sub="This is sub", xlab="This is xlab", ylab="This is ylab",
     type="o", lwd=2, col="blue", pch=24, bg="yellow", cex=1.5)

# xlab, ylab 축의 label
# lwd = line width
# lty = line type
# col = line color
# pch = plotting character or symbol 
# bg = symbor color
# cex = symbor expansion(size)


# type 변화
par( mfrow=c(3,3))
plot(y, type = 'p', main="type='p'")
plot(y, type = 'l', main="type='l'")
plot(y, type = 'b', main="type='b'")
plot(y, type = 'c', main="type='c'")
plot(y, type = 'o', main="type='o'")
plot(y, type = 'h', main="type='h'") # 수직막대
plot(y, type = 's', main="type='s'") # 계단 그래프
plot(y, type = 'S', main="type='S'") # 계단 그래프 반대
plot(y, type = 'n', main="type='n'") # no plotting


# Line type
par( mfrow=c(2,3))
for( i in 1:6 ) plot(y, type='o', lty=i, main=paste("line type=",i,sep=''))

colors() # 컬러 이름 다 호출
# 컬러 이름의 숫자 명시도 가능 col=552
which( colors() == "red") 



# pch ( plot character or symbol )
par( mfrow=c(3,3))
for( i in 1:9 ) plot(y, type='o', pch=i, main=paste('symbol type=',i,sep=''))



# lines, axis, box, title, legend, matplot 함수
par( mfrow=c(1, 1) )
cars <- data.frame(
  standard=c(1, 3, 6, 4, 9),
  truck=c(2, 5, 4, 5, 12),
  suv=c(4, 4, 6, 6, 16)
)
cars 
  

rownames(cars) <- c("Mon", "Tue", "Wed", "Thu", "Fri")
plot( cars$standard, type="o")
plot( cars$truck, type="p")
plot( cars$suv, type="l")

max_y <- max( cars )
# 최대값으로 y range 범위 제한
# axes=False 축 제거
plot( cars$standard, type="o", col="blue", ylim=c(0, max_y), axes=TRUE, ann=FALSE )
plot( cars$standard, type="o", col="blue", ylim=c(0, max_y), axes=FALSE, ann=TRUE )


# at = break 그래프 간격
axis(1, at=1:5, lab=rownames(cars) ) # 축변호 1 = 가로축
axis(2, at=seq(0, max_y, by=4), las=1 ) # 축번호 2 = 세로축
# box() 그래프 박스 그려줌

# plot 하위엔 lines() 함수로 선추가
lines( cars$truck, type="o", pch=22, lty=2, col="red" )
lines( cars$suv, type="o", pch=23, lty=3, col="green" )

# 그래프 main 이름 삽입
title(main="Car Rental", xlab="Weekday", ylab="The Number of Cars")
# 축 label 덮어씌워줌, reset 해서 다시 그려야 함
legend( "topleft", colnames(cars), col=c("blue", "red", "green"), pch=21:23, lty=1:3 )
legend("center", colnames(cars))
# "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" 그리고 "center"






library(MASS)
head(birthwt)
factor(birthwt$race)
birthwt$race <- factor( birthwt$race, levels=c(1,2,3), labels=c("white","black","other"))
factor(birthwt$race)


head(birthwt$bwt)
attach(birthwt)
head(bwt)


par(mfrow=c(2,3))
hist(age)
hist(age,freq=FALSE) # y축의 값을 frequencey 가 아닌 probability density로 계산
# x축의 break point를 사용자가 직접 지정
hist(age, breaks = c(10,15,17,19,21,24,25,30,35,45))
# bar에 라벨을 입력ㅊ(python의annotation)
hist(age, labels=c("1~14", "15~19", "20~24", "25~29", "30~34", "35~39", "40~44") )
hist(age, col="red")
# histogram을 채울 사선의 밀도, 1 inch 당 line의 수.
hist(age, col="red", density=10)


head( birthwt[, c("bwt", "race")] )


# density plot

par(mfrow=c(1,1))
plot(density(bwt), main = 'bwt')
plot(density(bwt[race=='black']), col="black" )
lines(density(bwt[race=="white"]), col="red")
lines(density(bwt[race=="other"]), col="blue")


# box plot

par(mfrow=c(1,2))
boxplot(bwt)
boxplot(bwt ~ race)
boxplot(bwt~smoke, col=rainbow(2), xlab="Smoke", ylab="Birth Weight")
boxplot(bwt~race, col=rainbow(3), xlab="Race", ylab="Birth Weight")


detach(birthwt)


# iris 데이터 로드

par(mfrow=c(1,1))
str(iris)
attach(iris)
iris.mean <- aggregate(iris[,1:4], iris["Species"], mean)
iris.mean
iris.sd <- aggregate(iris[,1:4], iris["Species"], sd)
iris.sd

iris.sd.upper <- iris.mean[,2:5] + iris.sd[,2:5]
iris.sd.lower <- iris.mean[,2:5] - iris.sd[,2:5]

b <- barplot( as.matrix(iris.mean[,2:5]), beside=TRUE, col=c("red","blue","purple"), ylim=c(0,8))

# 잔차표시
arrows(b, as.matrix(iris.sd.upper), b, as.matrix(iris.sd.lower), angle=90, length=0.05, code=3)
# legend 표시
legend("topright", legend=iris.mean$Species, fill=c("Red","blue","purple"))
