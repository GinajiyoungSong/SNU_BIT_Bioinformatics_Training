# ggplot2 EDA 
library(ggplot2)
head(iris)
pairs(iris[,1:4]) # pair plot 으로 전체 데이터의 분포를 간단하게 확인
attach(iris)

plot(Sepal.Length, Petal.Length, pch=c(21,21,21), bg=c("black", "blue", "red")[unclass(Species)])
legend('topleft', c('septosa','versicolor','virginica'), pch=c(21,21,21), pt.bg = c("black", "blue", "red"))
factor(Species)


# ggplot2 이용
myGraph <- ggplot(iris, aes(Sepal.Length, Petal.Length))+ geom_point()
# value 를 표시하기 위해, object type 을 골라서 geom_type() layer를 쌓아줌
myGraph <- ggplot(iris, aes(Sepal.Length, Petal.Length))+ geom_point(aes(color=Species))
print(myGraph)

# point size 는 가변적인 내용이므로 aes() 함수 안에 있음
graph <- ggplot(iris, aes(Sepal.Length, Petal.Length))+ geom_point(aes(color=Species, shape=Species, size=Sepal.Length/Sepal.Width))
print(graph)


# geom_smooth()는 데이터의 회귀라인을 추가
Graph <- myGraph + geom_point(aes(color=Species), size=3)+ geom_smooth(aes(color=Species)) 
print(Graph) # Loess lines

lm_graph <- myGraph + geom_point(aes(color=Species), size=3)+ geom_smooth(aes(color=Species), method = "lm") 
print(lm_graph) # line regression lines


### facet_grid()
graph <- ggplot(iris, aes(Sepal.Length, Petal.Length))
graph + geom_point(aes(color=Species), size=3)+ geom_smooth(aes(color=Species), method = "lm") + facet_grid(.~Species)

graph + geom_point(aes(color=Species), size=3)+ geom_smooth(aes(color=Species), method = "lm") + facet_grid(Species~.)


### histogram
graph <- ggplot(iris, aes(Sepal.Length))
graph + geom_histogram(alpha=0.7, color="red", fill="blue")

# density 선 추가
graph + geom_histogram(aes(y=..density..),alpha=0.7, color="red", fill="blue") + geom_density()

# facet_grid() 사용해서 Species 별로 각각 히스토그램 생성
graph + geom_histogram(aes(y=..density.., fill=Species),alpha=0.7) + geom_density(aes(color=Species))+ facet_grid(Species~.)

graph <- ggplot(iris, aes(Species,Sepal.Length, fill=Species))
graph + geom_boxplot()
graph + geom_boxplot() + geom_point()
graph + geom_boxplot() + geom_jitter() # 점들이 서로 겹치지 않게 random하게 뿌림



# 데이터 형변환 -> melt() 함수 사용
library(reshape2)
iris.melt <- melt(iris, id.var="Species")
head(iris.melt)

graph <- ggplot(iris.melt, aes(variable, value, fill=Species))
graph + geom_boxplot()

graph <- ggplot(iris.melt, aes(Species, value, fill=variable))
graph + geom_boxplot()


# bar_plot()

graph <- ggplot(iris.melt, aes(variable, value, color=Species, fill=Species))
graph <- graph + stat_summary(fun.y=mean, geom="bar", position = "dodge")
# Error bar를 그리기 위해서는 평균값과 편차값( lower, upper 모두)이 필요하다.
graph + stat_summary(fun.data = mean_cl_normal, position = position_dodge(width=0.9), geom="errorbar", width=0.3, color="black") +
  ggtitle(("IRIS data")) + xlab("Sepal/Petal") + ylab("Length/Width")

