setwd("C:/Users/wangn/Downloads/zuk_et_al-master/data")
library(data.table)
library(qqman)
library(foreach)
library(doSNOW)
cl<-makeCluster(6) #change the 2 to your number of CPU cores

registerDoSNOW(cl)
load('data.RData')
xpca<-prcomp(X,center = TRUE,scale. = TRUE)

plot(xpca$x[,1],xpca$x[,2])

cbind(xpca$x[,1],pc[,3])
z<-scale(xpca$x[,1:10])
summary(lm(y~z))


z<-X[,which(snp2set[,1]>0)]
zpca<-prcomp(z,center = TRUE,scale. = TRUE)
summary(lm(y~zpca$x[,3]))


z<-X[,which(snp2set[,4]>0)]
dim(z)
prob<-bdmcmc(z,y,1,100,100,1)
plot(prob)
index<-which(prob>0.5)
index
summary(lm(y~z[,index]))


require(fields)
# Make a 10x10 matrix

image.plot(Ka)
for (x in 1:965)
  for (y in 1:965)
    text((x-1)/9, (y-1)/9, sprintf("%0.2f", m[x,y]))