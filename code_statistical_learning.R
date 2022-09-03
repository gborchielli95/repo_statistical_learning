library(datasets)
library(leaps)

dataset=longley

fit_BSS <- regsubsets(Employed~.,data = dataset)

fit_BSS$

sfit <- summary(fit_BSS)
id=which.min(sfit$bic)
sfit$outmat

coef(fit_BSS,i=id)['GNP']


L=seq(1,16,2)
Learn <- dataset[L,]
Inference <- dataset[-L,]


fit_Lear <- regsubsets(Employed~.,data = Learn)
id_max_lear=which.min(summary(fit_Lear)$bic)
S=names(coef(fit_Lear,id_max_lear)[coef(fit_Lear,id_max_lear)!=0])



fit_Infern <- lm(Employed~Unemployed+Armed.Forces+Year,data = Inference)
sum_Infern <- summary(fit_Infern)
pvalueInfer <- sum_Infern$coefficients[-1,4]
p.adjust(pvalueInfer,'bonferroni')['Unemployed']


sLear <- summary(fit_Lear)



fit_Learn_bictop <-lm(Employed~Unemployed+Armed.Forces+Year,data = Learn)
predictEmI <- predict(fit_Learn_bictop,Inference[,-1])
sort(abs(Inference$Employed-predictEmI))[6]

m=8
a=1/3

ceiling()
is.infinite(9/3)
is.integer((1-a)*(m+1))

ceiling(((1-a)*(m+1)))
?ceiling
sort(abs(residuals(fit_Learn_bictop)),decreasing = F)[6]


library(MASS)
dataset <- MASS::mcycle
lamdas=10^seq(from=-4,to=2,by=0.1)
#B-splines#


tpower <- function(x, t, deg){
  (x - t) ^ deg * (x > t)
}
bbase <- function(x, xl, xr, ndx, deg){
  dx <- (xr - xl) / ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  Delta <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ 
                                              deg)
  B <- (-1) ^ (deg + 1) * P %*% t(Delta)
  B
}
x <- dataset$times
xl=min(dataset$times)
xr=max(dataset$times)
ndx=50
bdeg=3
B <- bbase(x, xl, xr, ndx, bdeg)


LOO <- vector(length = length(lamdas))
y <- dataset$accel
O <- 2
D <- diag(ncol(B))
for (k in 1:O) D <- diff(D)
dim(B)
dof <- vector(length = length(lamdas))
for (j in 1:length(lamdas)) {
  lamda=lamdas[j]
  y.hat <- vector(length = length(y))
  S <- B%*%solve(crossprod(B)+lamda*crossprod(D))%*%t(B)
  dof[j] <- sum(diag(S)) 
  for (i in 1:length(y)){
    B.i <- B[-i,]
    beta_hat <- solve(t(B.i) %*% B.i + lamda * t(D) %*% D, t(B.i) %*% y[-i])
    y.hat[i] <- (B[i,]%*%beta_hat)
  }
  
  LOO[j] <- 1/nrow(dataset)*sum((y-y.hat)^2)
  
}
lamdas[which.min(LOO)]
LOO[which.min(LOO)]
dof[which.min(LOO)]
beta_hat
