

p<-1000;
k<-50;
rep<-1000
se<-rep(0,rep*5);
se<-matrix(se,nrow=rep);
mu<-rep(0,p);
m1<-sqrt(2*log(p));
TP<-rep(0,rep*2);
FD<-rep(0,rep*2);
FDP<-rep(0,rep*2);
FDP<-matrix(FDP, nrow=rep);
TP<-matrix(TP,nrow=rep);
FD<-matrix(FD,nrow=rep);


mu[1:k]<-m1;
%mu[1:k]<-rnorm(k,0,2*m1);
for (i in 1:rep)
{X<-mu+rnorm(p);
 se[i,1]<-sum((X-mu)^2);
 JS1<-(1-(p-2)/sum(X^2))*X;
 se[i,2]<-sum((JS1-mu)^2);
 d<-(p-3)/(p-1)/var(X);
 JS2<-(1-d)*X+d*mean(X);
 se[i,3]<-sum((JS2-mu)^2);
 X1<-X;
 ADB<-(abs(X)<qnorm(1-0.1/2/p));
 X1[ADB]<-0;
 TP[i,1]<-k-sum(ADB[1:k]);
 FD[i,1]<-(p-k)-sum(ADB[(k+1):p]);
 FDP[i,1]<-FD[i,1]/max(1,(FD[i,1]+TP[i,1]));
 se[i,4]<-sum((X1-mu)^2);
 XS<-sort(abs(X), decreasing=TRUE);
 ss1<-seq(1:p);
 thrBH<-qnorm(1-ss1*0.1/2/p);
 ind2<-which((XS-thrBH)>0);
 if (length(ind2)>0)
 {indmax<-max(ind2);
 thr<-XS[indmax];
 X2<-X;
 ADBH<-(abs(X)<thr)
 X2[ADBH]<-0;
 TP[i,2]<-k-sum(ADBH[1:k]);
 FD[i,2]<-(p-k)-sum(ADBH[(k+1):p]);
 FDP[i,2]<-FD[i,2]/(FD[i,2]+TP[i,2])}
 else {X2=rep(0,p)}
 se[i,5]<-sum((X2-mu)^2);
}
mse<-apply(se,2,'mean');
power<-apply(TP,2,'mean')/k;
FDR<-apply(FDP,2,'mean')
boxplot(se, names=c('MLE','JS1','JS2','Bon','BH'))






