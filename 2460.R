
#read data
setwd("C://Users//zhang//Desktop")
aa=read.table("AA.F1_F2.dat")
aa=as.matrix(aa)
ao=read.table("AO.F1_F2.dat")
ao=as.matrix(ao)
test1=read.table("test_1.dat")
test1=as.matrix(test1[,2:3])

test2=read.table("test_2.dat")
test2=as.matrix(test2[,-1])
#visualize and analyze data
dim(aa)
dim(ao)
dim(test1)
dim(test2)
par(mfrow=c(1,2))
plot(aa[,2],aa[,3],xlab = "f1",ylab = "f2",main = "Formants of AA")
plot(ao[,2],ao[,3],xlab = "f1",ylab = "f2",main = "Formants of AO")
#build the model
aaf1=as.numeric(aa[,2])
aaf2=as.numeric(aa[,3])

aof1=as.numeric(ao[,2])
aof2=as.numeric(ao[,3])
#define a function, then doing cross validation would be easier below
successrate=function(aaf1,aaf2,aof1,aof2,test1,test2){
  #comput the det of sigma, which is the covariance matrix
mu=c(mean(aaf1),mean(aaf2))
a=cbind(aaf1,aaf2)
a=as.matrix(a)
det1=det(cov(a))


b=as.matrix(cbind(aof1,aof2))
det2=det(cov(b))


#test 
  #compute prior
p1=length(aaf1)/(length(aaf1)+length(aof1))
p2=1-p1

#test aa

AAcorrectrate=0
AArightsample=NULL

for (i in 1: nrow(test1)){c=c((test1[i,1]-mean(aaf1)),(test1[i,2]-mean(aaf2)))
d=c((test1[i,1]-mean(aof1)),(test1[i,2]-mean(aof2)))
  if (p1/det1^0.5*exp(-1/2*t(c)%*%solve(cov(a))%*%c) >
      p2/det2^0.5*exp(-1/2*t(d)%*%solve(cov(b))%*%d ))
    {AAcorrectrate=AAcorrectrate+1
     AArightsample=c(AArightsample,i)}
}
AAcorrectrate=AAcorrectrate/nrow(test1)

#test ao
AOcorrectrate=0
AOrightsample=NULL
for (i in 1: nrow(test2)){c=c((test2[i,1]-mean(aaf1)),(test2[i,2]-mean(aaf2)))
d=c((test2[i,1]-mean(aof1)),(test2[i,2]-mean(aof2)))
if (p1/det1^0.5*exp(-1/2*t(c)%*%solve(cov(a))%*%c) <
    p2/det2^0.5*exp(-1/2*t(d)%*%solve(cov(b))%*%d )){AOcorrectrate=AOcorrectrate+1
    AOrightsample=c(AOrightsample,i)}
}
AOcorrectrate=AOcorrectrate/nrow(test2)

return(c(AAcorrectrate,AOcorrectrate,AArightsample,AOrightsample))
}
#view the output of the classifier
successrate(aaf1,aaf2,aof1,aof2,test1,test2)
#correct rate for testing AA is 87.5% and only the sixth prediction is wrong
par(mfrow=c(1,1))
plot(test1[-6,1],test1[-6,2],xlab = "f1",ylab = "f2",xlim=c(50,300),ylim=c(650,1000),main = "test data")
abline(lm(test1[-6,2]~test1[-6,1]))
par(new=TRUE)
plot(test1[6,1],test1[6,2],pch=4,xlab = "f1",ylab = "f2",xlim=c(50,300),ylim=c(650,1000),main = "test data")
#correct rate for testing AA is 100%
#The classifier is effient 

##cross validation
totalaaf1=c(aaf1,test1[,1])
totalaaf2=c(aaf2,test1[,2])
totalaof1=c(aof1,test2[,1])
totalaof2=c(aof2,test2[,2])
length(totalaaf1)
length(totalaof1)
#partition into two sets
list1=sample(1:87,43)
list2=sample(1:108,54)
newaaf1=totalaaf1[list1 ]
newaaf2=totalaaf2[list1]
newaof2=totalaof2[list2]
newaof1=totalaof1[list2]
newtestaa=cbind(totalaaf1[-list1],totalaaf2[-list1])
newtestao=cbind(totalaof1[-list2],totalaof2[-list2 ])

successrate(newaaf1,newaaf2,newaof1,newaof2,newtestaa,newtestao)

# swap the test set and trainning set
swapaaf1=totalaaf1[-list1 ]  
swapaaf2=totalaaf2[-list1]
swapaof2=totalaof2[-list2]
swapaof1=totalaof1[-list2]
swaptestaa=cbind(totalaaf1[list1],totalaaf2[list1])
swaptestao=rbind(totalaof1[list2],totalaof2[list2 ])

successrate(swapaaf1,swapaaf2,swapaof1,swapaof2,swaptestaa,swaptestao)
#THE FIRST TWO ENTRY RETURNED ARE THE CORRECT RATE FOR AA AND AO RESPECTIVELY

#compute boundary
#a=cbind(aaf1,aaf2)
 #a=as.matrix(a)
#det1=det(cov(a))
 #b=as.matrix(cbind(aof1,aof2))
 #solve(cov(a))-solve(cov(b))
 #mu1=apply(a,2,mean)
 #mu2=apply(b,2,mean)
 #t(mu2)%*%solve(cov(b))%*%mu2-t(mu1)%*%solve(cov(a))%*%mu1
 #log(det2/det1)
l#log(p2/p1)