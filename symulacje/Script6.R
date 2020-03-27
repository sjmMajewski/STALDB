X<-rnorm(40000,0,0.1);
X<-matrix(X,nrow=200);
beta<-rep(0,200);
beta[1:20]<-4;
Y<-X%*%beta+rnorm(200);
library('bigstep');
d<-prepare_data(Y,X);


wyn<-fast_forward(d,crit=aic);
wyn<-stepwise(d,crit=bic);
wyn<-stepwise(d,crit=mbic);
wyn<-stepwise(d,crit=mbic2);

beta[1:20]<-10;


Y<-X%*%beta+rcauchy(200);

d<-prepare_data(Y,X);

wyn<-fast_forward(d,crit=aic);
wyn<-stepwise(d,crit=bic);
wyn<-stepwise(d,crit=mbic);
wyn<-stepwise(d,crit=mbic2);


rY<-rank(Y);

d<-prepare_data(rY,X);
wyn<-fast_forward(d,crit=aic);
wyn<-stepwise(d,crit=bic);
wyn<-stepwise(d,crit=mbic);
wyn<-stepwise(d,crit=mbic2);

library(MASS)
wyn0<-lm(Y~X[,1:40]);
s0<-summary(wyn0);
t0<-s0$coefficients;
plot(t0[2:41,3]~beta[1:40])
wyn1<-rlm(X[,1:40],Y,psi=psi.huber)
s1<-summary(wyn1);
t1<-s1$coefficients;
plot(t1[,3]~beta[1:40])

wyn2<-rlm(X[,1:40],Y,psi=psi.bisquare)
s2<-summary(wyn2);
t2<-s2$coefficients;
plot(t2[,3]~beta[1:40])






