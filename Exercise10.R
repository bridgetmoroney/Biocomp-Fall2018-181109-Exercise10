setwd("/Users/bridgetmoroney/Desktop/Biocomp-Fall2018-181109-Exercise10/")
data=read.table(file="data.txt",header = TRUE,sep = ",",stringsAsFactors = FALSE)

#using custom liklihood function for quadratic function
nlllike=function(p,x,y){
  a=p[1]
  b=p[2]
  c=p[3]
  sigma=exp(p[4])
  
  expected=a+b*x+c*x^2
  
  nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll)
}

#estimating parameters by minimizing negative log likelihood
initialGuess=c(1,1,1,1)
fit=optim(par=initialGuess,fn=nlllike,x=data$x,y=data$y)
print(fit)

#for simpler model
nlllike2=function(p,x,y){
  d=p[1]
  e=p[2]
  sigma2=exp(p[3])
  
  expected2=e+d*x
  
  nll2=-sum(dnorm(x=y,mean=expected2,sd=sigma2,log=TRUE))
  return(nll2)
}

#estimating parameters for simpler model
initialGuess2=c(1,1,1)
fit2=optim(par=initialGuess2,fn=nlllike2,x=data$x,y=data$y)
print(fit2)

#finding teststat
teststat=2*(fit2$value-fit$value)
#finding degrees of freedom
df=length(fit$par)-length(fit2$par)
#comparing to chi-squared value
1-pchisq(q=teststat,df=df)
#this value is 1 meaning that pchisq is so close to zero that the computer
#rounds to one. This means that the simpler model has a much better fit than the
#complex quadratic model

#part 2
library(deSolve)
library(ggplot2)

#defining function to write competition models
ddSim=function(t,y,p){
  N1=y[1]
  N2=y[2]
  R1=p[1]
  R2=p[2]
  alpha11=p[3]
  alpha12=p[4]
  alpha22=p[5]
  alpha21=p[6]
  
  model1=R1*(1-N1*alpha11-N2*alpha12)*N1
  model2=R2*(1-N2*alpha22-N1*alpha21)*N2
  
  return(list(c(model1,model2)))
}
#R should be less than one and alphas less than 0.1
#make sure alpha12 < alpha11 and alpha21 < alpha22
#model simulation one
params=c(0.8,0.5,0.09,0.03,0.08,0.04)
NO=c(0.5,1)
times=1:100

#simulating model using ode and following parameters
modelSim=ode(y=NO,times=times,func=ddSim,parms=params)
#converting to dataframe
modelOutput=data.frame(time=modelSim[,1],N1=modelSim[,2],N2=modelSim[,3])
#plotting simulation output
ggplot(modelOutput,aes(x=times,y=N1,color="blue"))+
  geom_line()+theme_classic()+geom_line(aes(x=times,y=N2,color="red"))+
  ylab("Species")+xlab("Time")+theme(legend.position = "none")
#plot shows that both species have reached coexistence and are causing
#each other to go extinct

#model simulation two, following parameters in article
params2=c(0.5,0.7,0.08,0.02,0.05,0.025)
NO2=c(0.4,0.9)
times2=1:100

modelSim2=ode(y=NO2,times=times2,func=ddSim,parms=params2)
modelOutput2=data.frame(time=modelSim2[,1],N1=modelSim2[,2],N2=modelSim2[,3])
ggplot(modelOutput2,aes(x=time,y=N1,color="blue"))+geom_line()+theme_classic()+
  geom_line(aes(x=times,y=N2,color="red"))+ylab("Species")+xlab("Time")+
  theme(legend.position = "none")
#shows same concept as above

#model simulation three, going against parameters stated in article
params3=c(0.2,0.9,0.01,0.02,0.02,0.03)
NO3=c(1,2)
times3=1:200

modelSim3=ode(y=NO3,times=times3,func=ddSim,parms=params3)
modelOutput3=data.frame(time=modelSim3[,1],N1=modelSim3[,2],N2=modelSim3[,3])              
ggplot(modelOutput3,aes(x=time,y=N1,color="blue"))+geom_line()+theme_classic()+
  geom_line(aes(x=times3,y=N2,color="red"))+ylab("Species")+xlab("Time")+
  theme(legend.position = "none")
#this model shows that one species is driven to extinction by the other
#when parameters for coexistence are not followed