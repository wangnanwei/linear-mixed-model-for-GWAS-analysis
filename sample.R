setwd("C:/Users/wangn/Downloads/zuk_et_al-master/data")
library(data.table)
library(qqman)
library(foreach)
library(doSNOW)
cl<-makeCluster(6) #change the 2 to your number of CPU cores

registerDoSNOW(cl)

#reading
data = fread("ypgly.nomissing.mat",header = T)
pheno = fread("y.txt",header = F)
pc = fread("pc.txt")
grm = fread("grm.standardized.txt") # additive linear kernel
snp2set = fread("snp2set.data", sep=",") # n by m matrix n being number of snps and m being number of pathways
pathways = fread("snp2set.pathway_names", sep="\n") # col name for snp2set
snps = fread("snp2set.snp_names", sep=",") #row names for snp2set -- equal to col names of data matrix

# scaling and preping the X
j = ncol(data)
i = nrow(data)
x.nomiss = as.matrix(data[2:i, 3:j])
rownames(x.nomiss) <- data$IID[2:i]
class(x.nomiss) <- "numeric"
scale = apply(x.nomiss, 2, sd, na.rm = TRUE)
center = colMeans(x.nomiss)
X = scale(x.nomiss, center=center, scale=scale)

# data cast / transform
y = as.numeric(pheno$V1)
Ka = as.matrix(grm)
#X = as.matrix()

# sample indices
n = 10^6   # just a few pairs to test
sample1 = sample(1:ncol(X), n, replace = T)
sample2 = sample(1:ncol(X), n, replace = T)

I = diag(nrow(Ka))
Vadditive = 0.689431*Ka + 0.599742*I
Vinv.additive = solve(Vadditive)
residual.additive = Vinv.additive %*% y

v <- foreach(i=1:n) %dopar%{
  column = scale(X[,sample1[i]]*X[,sample2[i]])[,1]
  (column %*% residual.additive)**2/(column %*% Vinv.additive %*% column)
}
v <- unlist(v)
pval.grm.add = pchisq(v, df=1, lower.tail = F)
qq(pval.grm.add, main="Pairwise LMM with additive GRM (Linear Kernel)")
# compare with uniform random
qq(runif(200))




v <- foreach(i=1:ncol(X)) %dopar%{
  column = scale(X[,i])[,1]
  (column %*% residual.additive)**2/(column %*% Vinv.additive %*% column)
}
v <- unlist(v)
pval.grm.add = pchisq(v, df=1, lower.tail = F)
qq(pval.grm.add, main="Single SNP LMM with additive GRM (Linear Kernel)")
# compare with uniform random
qq(runif(200))
snp_name=snp_name=names(data)[3:ncol(data)]
temp=order(pval.grm.add)[1:5]
snp_name[temp]
pval.grm.add[temp]

beta=rep(0,82597)
for(j in 1:82597){
  beta[j]<-(0.689431/82597)*X[,j]%*%residual.additive
  #beta[j]<-(X[,j]%*%Vinv.additive %*% y)/
}

v<-foreach(i=1:338) %dopar%{
  index<-which(snp2set[,i]>0)

  p1<-scale(X[,index]%*% beta[index])[,1]
  #p1<-scale(X[,index]%*% beta[index])[,1]
  (t(p1) %*% residual.additive)**2/(t(p1) %*% Vinv.additive %*% p1) 
}
v<-unlist(v)

pval = pchisq(v, df=1, lower.tail = F)

qq(pval, main="WPM LMM with additive GRM (Linear Kernel)")




v<-foreach(i=1:338) %dopar%{
  index<-which(snp2set[,i]>0)
  index2<-combn(index,2)
  z<-as.data.frame(X[,index])
  z=model.matrix(~.^2-.-1,data=z)
  
  zpca<-prcomp(z,center = TRUE,scale. = TRUE)
  p1=scale(zpca$x[,1])
  #p1<-scale(X[,index]%*% beta[index])[,1]
  #p1<-scale(X[,index]%*% beta[index])[,1]
  (t(p1) %*% residual.additive)**2/(t(p1) %*% Vinv.additive %*% p1) 
}
v<-unlist(v)

pval = pchisq(v, df=1, lower.tail = F)

qq(pval, main="WPM LMM with first principle component")








cb<-combn(1:338,2)


v<-foreach(i = 1:ncol(cb)) %dopar% {
  snp2set<-as.matrix(snp2set)
  index1<-which(snp2set[,cb[1,i]]>0)
  index2<-which(snp2set[,cb[2,i]]>0)
  z1<-X[,index1]
  zpca1<-prcomp(z1,center = TRUE,scale. = TRUE)
  p1=scale(zpca1$x[,1])
  z2<-X[,index2]
  zpca2<-prcomp(z2,center = TRUE,scale. = TRUE)
  p2=scale(zpca2$x[,1])
  
  pp<-scale(p1*p2)[,1]
  #p1<-scale(X[,index]%*% beta[index])[,1]
  (pp %*% residual.additive)**2/(pp %*% Vinv.additive %*% pp)
}

v<-unlist(v)
pval = pchisq(v, df=1, lower.tail = F)
qq(pval, main=" BPM LMM with additive GRM (Linear Kernel)")












v<-foreach(i = 1:338) %dopar% {
  snp2set<-as.matrix(snp2set)
  index1<-which(snp2set[,i]>0)
  p1<-scale(X[,index1]%*% beta[index1])[,1]
  
  pp<-scale(p1)[,1]
  #p1<-scale(X[,index]%*% beta[index])[,1]
  (pp %*% residual.additive)**2/(pp %*% Vinv.additive %*% pp)
}

v<-unlist(v)
pval = pchisq(v, df=1, lower.tail = F)
qq(pval, main=" BMP LMM with additive GRM (Linear Kernel)")




h=0.689431
e<-0.599742
h2 = h / (h+e)
lambda<-(1/h2)-1
library(glmnet)

model = glmnet(X,y, lambda=lambda, alpha = 0)
beta = model$beta[,1]


da
MODEL<-lm.ridge(y~.,as.data.frame(X),lambda = lambda)
model<-linearRidge(y~., as.data.frame(X),lambda=lambda)



detv<-det(Vadditive)
