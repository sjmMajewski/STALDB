

p<-seq(from=100,by=500, to=5000);
l1<-length(p);
k<-1;

rep<-1000
DB<-rep(0,rep*l1);
DB<-matrix(DB, nrow=rep);
DC<-rep(0,rep*l1);
DC<-matrix(DC, nrow=rep);

  for(j in 1:l1)
    { mu<-rep(0, p[j]);
      #k<-floor((p[j])^(3/4));
      #mu[1:k]<-0.7;
      mu[1:k]<-0.75*log(p[j]);
       for (i in 1:rep)
      {
          X<-mu+rnorm(p[j]);
         DB[i,j]<-(max(abs(X))>qnorm(1-0.1/2/p[j]));
         DC[i,j]<-(sum(X^2)>qchisq(1-0.1,p[j]))
     }
}
     
PowerB<-apply(DB,2,'mean');
PowerC<-apply(DC,2,'mean');
plot(PowerB, type='l', ylim=c(0,1), col='red');
points(PowerC, type='l', col='blue')


