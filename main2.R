source('GCVsmoothing.R')
#install.packages('fda')
#install.packages('fields')
#install.packages('clusteval')
library(fda)
library(fields)
#library(clusteval)
#######################################################################
########Generate data for simulation study (part 1)####################
#######################################################################
ch = 1000
T = 1000
J = floor((T+1)/2)-1

kk=100 # number of replicates
span = seq(from = 3, to = 15, by = 1)
result1 = matrix(0,25,kk)
result2 = matrix(0,25,kk)
result3 = matrix(0,25,kk)
result4 = matrix(0,25,kk)
result5 = matrix(0,25,kk)
result6 = matrix(0,25,kk)
sim_1 = array(0,c(40,499,25))
for(k in 1:kk){
  
  set.seed(k)
  a = rnorm(5,0,0.01)#random difference
  a1 = a[1]
  a2 = a[2]
  a3 = a[3]
  a4 = a[4]
  a5 = a[5]
  data1 = matrix(nrow=200,ncol=1000);   #cluster 1 
  for(i in 1:40){
    data1[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.8+a1,0.1)),n=1000)}
  for(i in 41:80){
    data1[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.8+a2,0.1)),n=1000)}
  for(i in 81:120){
    data1[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.8+a3,0.1)),n=1000)}
  for(i in 121:160){
    data1[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.8+a4,0.1)),n=1000)}
  for(i in 161:200){
    data1[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.8+a5,0.1)),n=1000)}
  
  
  data2 = matrix(nrow=200,ncol=1000);  #cluster 2
  for(i in 1:40){
    data2[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.9+a1,-0.9)),n=1000)}
  for(i in 41:80){
    data2[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.9+a2,-0.9)),n=1000)}
  for(i in 81:120){
    data2[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.9+a3,-0.9)),n=1000)}
  for(i in 121:160){
    data2[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.9+a4,-0.9)),n=1000)}
  for(i in 161:200){
    data2[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.9+a5,-0.9)),n=1000)}
  
  data3 = matrix(nrow=200,ncol=1000);  #cluster 3
  for(i in 1:40){
    data3[i,]=arima.sim(list(order=c(2,0,0),ar=c(-0.1+a1,-0.9)),n=1000)}
  for(i in 41:80){
    data3[i,]=arima.sim(list(order=c(2,0,0),ar=c(-0.1+a2,-0.9)),n=1000)}
  for(i in 81:120){
    data3[i,]=arima.sim(list(order=c(2,0,0),ar=c(-0.1+a3,-0.9)),n=1000)}
  for(i in 121:160){
    data3[i,]=arima.sim(list(order=c(2,0,0),ar=c(-0.1+a4,-0.9)),n=1000)}
  for(i in 161:200){
    data3[i,]=arima.sim(list(order=c(2,0,0),ar=c(-0.1+a5,-0.9)),n=1000)}
  
  data4=matrix(nrow=200,ncol=1000); #cluster 4
  for(i in 1:40){
    data4[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.9+a1,-0.9)),n=1000)}
  for(i in 41:80){
    data4[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.9+a2,-0.9)),n=1000)}
  for(i in 81:120){
    data4[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.9+a3,-0.9)),n=1000)}
  for(i in 121:160){
    data4[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.9+a4,-0.9)),n=1000)}
  for(i in 161:200){
    data4[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.9+a5,-0.9)),n=1000)}
  
  data5 = matrix(nrow=200,ncol=1000);  #cluster 5
  for(i in 1:40){
    data5[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.8+a1,0.1)),n=1000)}
  for(i in 41:80){
    data5[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.8+a2,0.1)),n=1000)}
  for(i in 81:120){
    data5[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.8+a3,0.1)),n=1000)}
  for(i in 121:160){
    data5[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.8+a4,0.1)),n=1000)}
  for(i in 161:200){
    data5[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.8+a5,0.1)),n=1000)}
  #data=rbind(data1,data2,data3,data4,data5)  # 1000 time series
  
  data11 = data1
  data22 = (4/5)*data1+(1/10)*data2
  data33 = (3/5)*data1+(1/10)*data3
  data44 = (2/5)*data1+(1/10)*data4
  data55 = (1/5)*data1+(1/10)*data5
  data=rbind(data11,data22,data33,data44,data55)
  dataout=addgamma(data,0.33)
  
  
  ### GCV smoothing ###
  span = seq(from = 3, to = 15, by = 1)
  simm = LogSmoothGVC(inputdata =data,span,T,J,ch)
  simmout = LogSmoothGVC(inputdata =dataout,span,T,J,ch)
  
  sim = array(0,c(40,499,25))
  simout = sim
  
  for(i in 1:25){
    sim[,,i] = simm[((i-1)*40+1):(i*40),]
    simout[,,i] = simmout[((i-1)*40+1):(i*40),]
  }
  if(k == 1){sim_1 = sim}
  resultfm = robustcluster(sim,5,'fm')
  resultcr = robustcluster(sim,5,'cr')
  resultm = robustcluster(sim,5,'mean')
  resultfmout = robustcluster(simout,5,'fm')
  resultcrout = robustcluster(simout,5,'cr')
  resultmout = robustcluster(simout,5,'mean')
  

  result1[,k] = resultfm$res
  result2[,k] = resultcr$res
  result3[,k] = resultm$res
  result4[,k] = resultfmout$res
  result5[,k] = resultcrout$res
  result6[,k] = resultmout$res
  print(k)
}
####Rand index
tru = c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5) #true cluster
rifm = 0
for(i in 1:100){
  rifm = rifm + cluster_similarity(result1[,i],tru)/100
}
ricr = 0
for(i in 1:100){
  ricr = ricr + cluster_similarity(result2[,i],tru)/100
}
rim=0
for(i in 1:100){
  rim = rim + cluster_similarity(result3[,i],tru)/100
}
rifmout = 0
for(i in 1:100){
  rifmout = rifmout+cluster_similarity(result4[,i],tru)/100
}
ricrout = 0
for(i in 1:100){
  ricrout = ricrout+cluster_similarity(result5[,i],tru)/100
}
rimout = 0
for(i in 1:100){
  rimout = rimout+cluster_similarity(result6[,i],tru)/100
}
########add constant#######################
kkk = 100
res1 = matrix(0,25,kkk)
res2 = res1
res3 = res1
set.seed(100)
for(i in 1:kkk){
  simout2 = array(0,c(40,499,25))
  for(j in 1:25){  
    simout2[,,j] = addcon(sim_1[,,j],0.15)
  }  
  resultfmout2 = robustcluster(simout2,5,'fm')
  resultcrout2 = robustcluster(simout2,5,'cr')
  resultmout2 = robustcluster(simout2,'mean')
  res1[,i] = resultfmout2$t
  res2[,i] = resultcrout2$t
  res3[,i] = resultmout2$t
}
#####Rand index
rifmout1 = 0
for(i in 1:100){
  rifmout1 = rifmout1 + cluster_similarity(res1[,i],tru)/100
}
ricrout1 = 0
for(i in 1:100){
  ricrout1 = ricrout1 + cluster_similarity(res2[,i],tru)/100
}
rimout1 = 0
for(i in 1:100){
  rimout1 = rimout1 + cluster_similarity(res3[,i],tru)/100
}


##########reporduce Figure 8 and 11##
example = sim_1
exampleout = example
for(i in 1:25){  
  exampleout[,,i] = addcon(example[,,i],0.15)
}  

median = matrix(0,25,499)
mean = median
medianout = median
meanout = median
for(i in 1:25){
  median[i,] = fmed(example[,,i])
  mean[i,] = colMeans(example[,,i])
  medianout[i,] = fmed(exampleout[,,i])
  meanout[i,] = colMeans(exampleout[,,i])
}

xl = seq(from=0.1,to=49.9,by=0.1)
set.panel(2,2)
plot(xl,example[1,,1],type='l',ylim=c(-3,10),xlab='Freq in Hertz',ylab='Log-periodogram',main='(a) Data without contaminations')
for(i in 2:40){
  lines(xl,example[i,,1])
}
fbplot(t(example[,,1]),ylim=c(-3,10),xlab='Freq in Hertz',ylab='Log-periodogram',main='(b) Fbplot of the data')
plot(xl,exampleout[1,,1],type='l',ylim=c(-3,10),xlab='Freq in Hertz',ylab='Log-periodogram',main='(c) Data with contaminations')
for(i in 2:40){
  lines(xl,exampleout[i,,1])
}
fbplot(t(exampleout[,,1]),ylim=c(-3,10),xlab='Freq in Hertz',ylab='Log-periodogram',main='(d) Fbplot of the data')


set.panel(2,2)
plot(xl,median[1,],type='l',ylim=c(-3.5,6),xlab='Freq in Hertz',ylab='Log-periodogram',main='(a) Medians without contaminations')
for(i in 2:25){
  lines(xl,median[i,])
}
plot(xl,medianout[1,],type='l',ylim=c(-3.5,6),xlab='Freq in Hertz',ylab='Log-periodogram',main='(b) Medians with contaminations')
for(i in 2:25){
  lines(xl,medianout[i,])
}
plot(xl,mean[1,],type='l',ylim=c(-3.5,6),xlab='Freq in Hertz',ylab='Log-periodogram',main='(c) Means without contaminations')
for(i in 2:25){
  lines(xl,mean[i,])
}

plot(xl,meanout[1,],type='l',ylim=c(-3.5,6),xlab='Freq in Hertz',ylab='Log-periodogram',main='(d) Means with contaminations')
for(i in 2:25){
  lines(xl,meanout[i,])
}
set.panel(1,1)

