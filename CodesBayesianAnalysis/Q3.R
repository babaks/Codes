getPostLambda=function(lambda,eta,sigma,yObs,n,alpha1,beta1){
prior=log(dbeta(lambda,alpha1,beta1));
likelihood=sum(log(lambda*dnorm(yObs,0,sigma)+(1-lambda)*dnorm(yObs,0,sqrt(sigma^2+eta^2))))
}

getPostEta=function(lambda,eta,sigma,yObs,n,alpha2,beta2){
prior=log(dbeta(eta,alpha1,beta1));
likelihood=sum(log(lambda*dnorm(yObs,0,sigma)+(1-lambda)*dnorm(yObs,0,sqrt(sigma^2+eta^2))))
}

getPostSigma=function(lambda,eta,sigma,yObs,n,alpha3,beta3){
prior=log(dbeta(sigma,alpha3,beta3));
likelihood=sum(log(lambda*dnorm(yObs,0,sigma)+(1-lambda)*dnorm(yObs,0,sqrt(sigma^2+eta^2))))
}



mySliceSampler = function(nIterations=100){
# this is the observed dataset
#setwd("C:/Documents and Settings/sjavanma.UCI-ICS/Desktop/HW2/")
yObs= scan("tests.txt",what=numeric());
n=length(yObs);

#these are the verctors for the sample points
lambdas=numeric();
etas=numeric();
sigmas=numeric()

# these are parameters for lambda
alpha1=0.7; beta1=0.1;
lambdas[1]=0.7
etas[1]=2;
sigmas[1]=2;

w=0.1;m=20;
for (i in 2:nIterations){
#starting points should be in range

z=getPostLambda(lambdas[i-1],etas[i-1],sigmas[i-1],yObs,n,alpha1,beta1)-rexp(1);
currentParam=lambdas[i-1]
#Stepping Out
u=runif(1);
L=currentParam-w*u;
R=L+w;
v=runif(1);
J=floor(m*v);
K=(m-1)-J;

browser()
while(J>0 && L>0 && z<getPostLambda(L,etas[i-1],sigmas[i-1],yObs,n,alpha1,beta1)){
 L= L-w;
 L=max(0,L)
 J= J-1;
}#while

while(K>0 && R<1 && z<getPostLambda(R,etas[i-1],sigmas[i-1],yObs,n,alpha1,beta1)){
 R = R+w
 R = min(1,R);
 K = K-1;
}#while


# Shrinkage to obtain a sample
u=runif(1);
newParam=L+u*(R-L);
while(z>getPostLambda(newParam,etas[i-1],sigmas[i-1],yObs,n,alpha1,beta1)){
if(newParam<currentParam){
 L=newParam;}
else{R=newParam;}
 u = runif(1);
 newParam=L+u*(R-L)

}#while

newParam
lambdas[i]= newParam;
print(c(i, newParam))
browser()
}#for nIterations

 return (list(lambdas=lambdas,etas=etas,sigmas=sigmas));
}#function


# Running the program
res = mySliceSampler(100);