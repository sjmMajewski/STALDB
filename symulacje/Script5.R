rep1=100;
n<-500;
p<-500
k<-10;
X<-rnorm(n*p,0,1/sqrt(n));
X<-matrix(X,nrow=n);
beta<-rep(0,p);
beta[1:k]<-3.5;
l<-rep(0,rep1)

for (i in 1:rep1)
   {
 
     Y=X%*%beta+rnorm(n);
     dat<-prepare_data(Y,X); 
     obj1<-fast_forward(dat, crit='bic', maxf=200);
     res<-as.numeric(obj1$model);
     l[i]<-sum(res>k)
  }
hist(l, main='False Discoveries by BIC', xlab='FD');
mean(l)     
 



  
     
     
   
