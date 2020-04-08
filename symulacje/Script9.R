rep1=10;
n<-500;
p<-c(50,100,250,450);
k<-50;
l<-rep(0,rep1);
msels<-matrix(rep(0,rep1*4),nrow=rep1);
surels<-matrix(rep(0,rep1*4),nrow=rep1);
mseridge<-matrix(rep(0,rep1*4),nrow=rep1);
sureridge<-matrix(rep(0,rep1*4),nrow=rep1);
msebls<-matrix(rep(0,rep1*4),nrow=rep1);
msebridge<-matrix(rep(0,rep1*4),nrow=rep1);

library(mvtnorm);





up<-100/(3.5^2*50-50);
gamma=0.15;

for (j in 1:4)
   {
     beta<-rep(0,p[j]);
     beta[1:k]<-3.5; 
     
     sigma<-matrix(0.6,p[j],p[j])
     diag(sigma)<-1

         for (i in 1:rep1)
{ 

         X<-(rmvnorm(n,numeric(p[j]),sigma))/sqrt(n);
         Y=X%*%beta+rnorm(n);
         hbeta<-coefficients(lm(Y~X-1));
         msels[i,j]<-sum((X%*%(hbeta-beta))^2);
         msebls[i,j]<-sum((hbeta-beta)^2);

         surels[i,j]<-sum((X%*%hbeta-Y)^2)+2*p[j]-n;

         
         sridge<-solve(t(X)%*%X+gamma*diag(p[j]));
         trM<-sum(diag(X%*%sridge%*%t(X)));
         hridge<-sridge%*%t(X)%*%Y;
         mseridge[i,j]<-sum((X%*%(hridge-beta))^2);
         msebridge[i,j]<-sum((hridge-beta)^2);
         sureridge[i,j]<-sum((X%*%hridge-Y)^2)+2*trM-n;
      }
    }
par(mfrow=c(1,2))
boxplot(msels,ylim=c(0,600))
boxplot(surels,ylim=c(0,600))

par(mfrow=c(1,2))
boxplot(msebls,ylim=c(0,15000))
boxplot(msebridge,ylim=c(0,15000))


par(mfrow=c(1,2))
boxplot(msebls,ylim=c(0,1500))
boxplot(msebridge,ylim=c(0,1500))

diff<-matrix(rep(0,4*rep1),nrow=rep1);
diff<-msebls-msebridge;

par(mfrow=c(1,2))
boxplot(diff);
boxplot(diff[,1:3])


par(mfrow=c(1,2))
boxplot(msels,ylim=c(0,600))
boxplot(mseridge,ylim=c(0,600))

diff<-matrix(rep(0,4*rep1),nrow=rep1);
diff<-msels-mseridge;
boxplot(diff);


diffls<-matrix(rep(0,4*rep1),nrow=rep1);
diffls<-surels-msels;
diffridge<-matrix(rep(0,4*rep1),nrow=rep1);
diffridge<-sureridge-mseridge;

par(mfrow=c(1,2))
boxplot(diffls,ylim=c(-70,70))
boxplot(diffridge,ylim=c(-70,70))

         
         
         
