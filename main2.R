source('GCVsmoothing.R')
library(fda)
library(fields)
library(clusteval)
library(TSclust)
library(mclust)
#######################################################################
########Generate data for simulation study (part 1)####################
#######################################################################
ch = 1000
T = 1000;
J = floor((T+1)/2)-1;

kk=100 # number of replicates
span = seq(from = 3, to = 15, by = 1)
resultt1 = matrix(0, 25, kk)
resultt2 = matrix(0, 25, kk)
resultt3 = matrix(0, 25, kk)
resultt_llrdls3 = matrix(0, 25, kk)
resultt4 = matrix(0,25*4, kk)
resultt5 = matrix(0,25*4, kk)
resultt6 = matrix(0,25*4, kk)
resultt7 = matrix(0,25*4, kk)
res1 = matrix(0, 100, kk)
res2 = res1
res3 = res1
res_llrdls3 = res1
sim_1 = array(0, c(40, 499, 25))
for(k in 1:100){
  set.seed(k)
  a = rnorm(5,0,0.01)#random difference
  a1 = a[1]
  a2 = a[2]
  a3 = a[3]
  a4 = a[4]
  a5 = a[5]
  data1 = matrix(nrow = 200, ncol = 1000);   #cluster 1 
  for(i in 1:40){
    data1[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(0.8+a1, 0.1)), n=1000)}
  for(i in 41:80){
    data1[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(0.8+a2, 0.1)), n=1000)}
  for(i in 81:120){
    data1[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(0.8+a3, 0.1)), n=1000)}
  for(i in 121:160){
    data1[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(0.8+a4, 0.1)), n=1000)}
  for(i in 161:200){
    data1[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(0.8+a5, 0.1)), n=1000)}
  
  
  data2=matrix(nrow=200,ncol=1000);  #cluster 2
  for(i in 1:40){
    data2[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(0.9+a1, -0.9)), n=1000)}
  for(i in 41:80){
    data2[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(0.9+a2, -0.9)), n=1000)}
  for(i in 81:120){
    data2[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(0.9+a3, -0.9)), n=1000)}
  for(i in 121:160){
    data2[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(0.9+a4, -0.9)), n=1000)}
  for(i in 161:200){
    data2[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(0.9+a5, -0.9)), n=1000)}
  
  data3=matrix(nrow=200,ncol=1000);  #cluster 3
  for(i in 1:40){
    data3[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(-0.1+a1, -0.9)), n=1000)}
  for(i in 41:80){
    data3[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(-0.1+a2, -0.9)), n=1000)}
  for(i in 81:120){
    data3[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(-0.1+a3, -0.9)), n=1000)}
  for(i in 121:160){
    data3[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(-0.1+a4, -0.9)), n=1000)}
  for(i in 161:200){
    data3[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(-0.1+a5, -0.9)), n=1000)}
  
  data4=matrix(nrow=200,ncol=1000); #cluster 4
  for(i in 1:40){
    data4[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(-0.9+a1, -0.9)), n=1000)}
  for(i in 41:80){
    data4[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(-0.9+a2, -0.9)), n=1000)}
  for(i in 81:120){
    data4[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(-0.9+a3, -0.9)), n=1000)}
  for(i in 121:160){
    data4[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(-0.9+a4, -0.9)), n=1000)}
  for(i in 161:200){
    data4[i,] = arima.sim(list(order = c(2, 0, 0), ar = c(-0.9+a5, -0.9)), n=1000)}
  
  data5=matrix(nrow=200,ncol=1000);  #cluster 5
  for(i in 1:40){
    data5[i,] = arima.sim(list(order=c(2, 0, 0), ar=c(-0.8+a1, 0.1)), n=1000)}
  for(i in 41:80){
    data5[i,] = arima.sim(list(order=c(2, 0, 0), ar=c(-0.8+a2, 0.1)), n=1000)}
  for(i in 81:120){ 
    data5[i,] = arima.sim(list(order=c(2, 0, 0), ar=c(-0.8+a3, 0.1)), n=1000)}
  for(i in 121:160){
    data5[i,] = arima.sim(list(order=c(2, 0, 0), ar=c(-0.8+a4, 0.1)), n=1000)}
  for(i in 161:200){
    data5[i,] = arima.sim(list(order=c(2, 0, 0), ar=c(-0.8+a5, 0.1)), n=1000)}
  #data=rbind(data1,data2,data3,data4,data5)  # 1000 time series
  
  data11 = data1
  data22 = (4/5) * data1 + (1/10) * data2
  data33 = (3/5) * data1 + (1/10) * data3
  data44 = (2/5) * data1 + (1/10) * data4
  data55 = (1/5) * data1 + (1/10) * data5
  data = rbind(data11, data22, data33, data44, data55)
  dataout1 = addgamma(data, 0.25)
  dataout2 = addgamma(data, 0.3)
  dataout3 = addgamma(data, 0.35)
  # dataout4=addgamma(data,0.4)
  
  ### GCV smoothing ###
  span = seq(from = 3, to = 15, by = 1)
  simm = LogSmoothGVC(inputdata=data,span,T,J,ch)
  simmout1 = LogSmoothGVC(inputdata=dataout1, span, T, J, ch)
  simmout2 = LogSmoothGVC(inputdata=dataout2, span, T, J, ch)
  simmout3 = LogSmoothGVC(inputdata=dataout3, span, T, J, ch)
  # simmout4=LogSmoothGVC(inputdata =dataout4,span,T,J,ch)
  
  sim = array(0, c(40, 499, 25))
  simout2_0.25 = sim
  simout2_0.3 = sim
  simout2_0.35 = sim
  
  simout_0.1 = sim 
  simout_0.2 = sim
  simout_0.3 = sim

  for(i in 1:25){
    sim[,,i] = simm[((i-1)*40+1):(i*40),]
    simout2_0.25[ , , i] = simmout1[((i-1)*40+1):(i*40), ]
    simout2_0.3[ , , i] = simmout2[((i-1)*40+1):(i*40), ] 
    simout2_0.35[ , , i] = simmout3[((i-1)*40+1):(i*40), ] 

   
    simout_0.1[,,i] = addcon(sim[ , , i], 0.1)
    simout_0.2[,,i] = addcon(sim[ , , i], 0.2)
    simout_0.3[,,i] = addcon(sim[ , , i], 0.3)
  }
  
  if(k == 1){sim_1 = sim}
  resultfm = fmclust(sim, 5)
  resultcr = crclust(sim, 5)
  resultm = meanclust(sim, 5)
  resultt_llrdls3[,k] = meanclust_llrdls(sim,5)$t
  
  resultt1[, k] = resultfm$t
  resultt2[, k] = resultcr$t
  resultt3[, k] = resultm$t
  resultt4[1:25, k] = fmclust(simout_0.1, 5)$t
  resultt4[26:50, k] = fmclust(simout_0.2, 5)$t
  resultt4[51:75, k] = fmclust(simout_0.3, 5)$t

  resultt5[1:25, k] = crclust(simout_0.1, 5)$t
  resultt5[26:50, k] = crclust(simout_0.2, 5)$t
  resultt5[51:75, k] = crclust(simout_0.3, 5)$t

  resultt6[1:25, k] = meanclust(simout_0.1, 5)$t
  resultt6[26:50, k] = meanclust(simout_0.2, 5)$t  
  resultt6[51:75, k] = meanclust(simout_0.3, 5)$t

  resultt7[1:25,k] = meanclust_llrdls(simout_0.1, 5)$t
  resultt7[26:50,k] = meanclust_llrdls(simout_0.2, 5)$t  
  resultt7[51:75,k] = meanclust_llrdls(simout_0.3, 5)$t

  res1[1:25, k] = fmclust(simout2_0.25, 5)$t
  res1[26:50, k] = fmclust(simout2_0.3, 5)$t
  res1[51:75, k] = fmclust(simout2_0.35, 5)$t

  res2[1:25, k] = crclust(simout2_0.25, 5)$t
  res2[26:50, k] = crclust(simout2_0.3, 5)$t
  res2[51:75, k] = crclust(simout2_0.35, 5)$t
 
  res3[1:25, k] = meanclust(simout2_0.25, 5)$t
  res3[26:50, k] = meanclust(simout2_0.3, 5)$t
  res3[51:75, k] = meanclust(simout2_0.35, 5)$t

  res_llrdls3[1:25, k] = meanclust_llrdls(simout2_0.25, 5)$t
  res_llrdls3[26:50, k] = meanclust_llrdls(simout2_0.3, 5)$t
  res_llrdls3[51:75, k] = meanclust_llrdls(simout2_0.35, 5)$t
 
  print(k)
}


index = 1:100
tru = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5)
rifm = 0
for(i in index){
  rifm = rifm + adjustedRandIndex(resultt1[1:25, i], tru) / 100
}
ricr = 0
for(i in index){
  ricr = ricr + adjustedRandIndex(resultt2[1:25, i], tru) / 100
}
rim = 0
for(i in index){
  rim = rim + adjustedRandIndex(resultt3[1:25, i], tru) / 100
}
rim_llrdls = 0
for(i in index){
  rim_llrdls = rim_llrdls + adjustedRandIndex(resultt_llrdls3[1:25, i], tru) / 100
}

###type 2
rifmout2_0.25 = 0
rifmout2_0.3 = 0
rifmout2_0.35 = 0
for(i in index){
  rifmout2_0.25 = rifmout2_0.25 + adjustedRandIndex(res1[1:25, i], tru) / 100
  rifmout2_0.3 = rifmout2_0.3 + adjustedRandIndex(res1[26:50, i], tru) / 100
  rifmout2_0.35 = rifmout2_0.35 + adjustedRandIndex(res1[51:75, i], tru) / 100
}
ricrout2_0.25=0
ricrout2_0.3=0
ricrout2_0.35=0
for(i in index){
  ricrout2_0.25 = ricrout2_0.25 + adjustedRandIndex(res2[1:25, i], tru) / 100
  ricrout2_0.3 = ricrout2_0.3 + adjustedRandIndex(res2[26:50, i], tru) / 100
  ricrout2_0.35 = ricrout2_0.35 + adjustedRandIndex(res2[51:75, i], tru) / 100
}
rimout2_0.25 = 0
rimout2_0.3 = 0
rimout2_0.35 = 0

rimout2_llr_0.25 = 0
rimout2_llr_0.3 = 0
rimout2_llr_0.35 = 0

for(i in index){
  rimout2_0.25 = rimout2_0.25 + adjustedRandIndex(res3[1:25, i],tru) / 100
  rimout2_0.3 = rimout2_0.3 + adjustedRandIndex(res3[26:50, i], tru) / 100
  rimout2_0.35 = rimout2_0.35 + adjustedRandIndex(res3[51:75, i], tru) / 100
  
  rimout2_llr_0.25 =rimout2_llr_0.25 + adjustedRandIndex(res_llrdls3[1:25, i], tru) / 100
  rimout2_llr_0.3 = rimout2_llr_0.3 + adjustedRandIndex(res_llrdls3[26:50, i], tru) / 100
  rimout2_llr_0.35 = rimout2_llr_0.35 + adjustedRandIndex(res_llrdls3[51:75, i], tru) / 100
}
rifmout1_0.1 = 0
rifmout1_0.2 = 0
rifmout1_0.3 = 0
for(i in index){
  rifmout1_0.1 = rifmout1_0.1 + adjustedRandIndex(resultt4[1:25, i], tru) / 100
  rifmout1_0.2 = rifmout1_0.2 + adjustedRandIndex(resultt4[26:50, i], tru) / 100
  rifmout1_0.3 = rifmout1_0.3 + adjustedRandIndex(resultt4[51:75, i], tru) / 100
  }
ricrout1_0.1 = 0
ricrout1_0.2 = 0
ricrout1_0.3 = 0
for(i in index){
  ricrout1_0.1 = ricrout1_0.1 + adjustedRandIndex(resultt5[1:25, i], tru) / 100
  ricrout1_0.2 = ricrout1_0.2 + adjustedRandIndex(resultt5[26:50, i], tru) / 100
  ricrout1_0.3 = ricrout1_0.3 + adjustedRandIndex(resultt5[51:75, i], tru) / 100
  }
rimout1_0.1 = 0
rimout1_0.2 = 0
rimout1_0.3 = 0
rimout1_llr_0.1 = 0
rimout1_llr_0.2 = 0
rimout1_llr_0.3 = 0
for(i in index){
  rimout1_0.1 = rimout1_0.1 + adjustedRandIndex(resultt6[1:25, i], tru) / 100
  rimout1_0.2 = rimout1_0.2 + adjustedRandIndex(resultt6[26:50, i], tru) / 100
  rimout1_0.3 = rimout1_0.3 + adjustedRandIndex(resultt6[51:75, i], tru) / 100

  rimout1_llr_0.1 = rimout1_llr_0.1 + adjustedRandIndex(resultt7[1:25, i], tru) / 100
  rimout1_llr_0.2 = rimout1_llr_0.2 + adjustedRandIndex(resultt7[26:50, i], tru) / 100
  rimout1_llr_0.3 = rimout1_llr_0.3 + adjustedRandIndex(resultt7[51:75, i], tru) / 100
}
####################results############################
rifmout1_0.1
rifmout1_0.2
rifmout1_0.3

ricrout1_0.1
ricrout1_0.2
ricrout1_0.3

rimout1_0.1
rimout1_0.2
rimout1_0.3

rimout1_llr_0.1
rimout1_llr_0.2
rimout1_llr_0.3

rifmout2_0.25
rifmout2_0.3
rifmout2_0.35

ricrout2_0.25
ricrout2_0.3
ricrout2_0.35

rimout2_0.25
rimout2_0.3
rimout2_0.35

rimout2_llr_0.25
rimout2_llr_0.3
rimout2_llr_0.35

#########time#######
# ti = Sys.time()
# t = fmclust(simout2_0.25,5)$t
# time = Sys.time() - ti


resfm = matrix(0,25,100)
rescr = matrix(0,25,100)
resfm1_10 = matrix(0,25,100)
resfm1_20 = matrix(0,25,100)
resfm1_30 = matrix(0,25,100)

resfm2_25 = matrix(0,25,100)
resfm2_30 = matrix(0,25,100)
resfm2_35 = matrix(0,25,100)

rescr1_10 = matrix(0,25,100)
rescr1_20 = matrix(0,25,100)
rescr1_30 = matrix(0,25,100)

rescr2_25 = matrix(0,25,100)
rescr2_30 = matrix(0,25,100)
rescr2_35 = matrix(0,25,100)





ch=1000
T=1000;
J = floor((T+1)/2)-1;
span = seq(from = 3, to = 15, by = 1)
set.seed(1)
for(k in 1:100){
  
  
  data=matrix(0,1000,1000)
  for(i in 1:100){
    data[i,]=arima.sim(list(order=c(2,0,0),ar=c(0.8,0.1)),n=1000,sd=exp(0))
    data[i+100,]=arima.sim(list(order=c(2,0,0),ar=c(0.8,0.1)),n=1000,sd=exp(1))
    data[i+200,]=arima.sim(list(order=c(2,0,0),ar=c(0.8,0.1)),n=1000,sd=exp(2.2/2))
    data[i+300,]=arima.sim(list(order=c(2,0,0),ar=c(0.8,0.1)),n=1000,sd=exp(4.2/2))
    data[i+400,]=arima.sim(list(order=c(2,0,0),ar=c(0.8,0.1)),n=1000,sd=exp(4.4/2))
    data[i+500,]=arima.sim(list(order=c(2,0,0),ar=c(0.8,0.1)),n=1000,sd=exp(6.4/2))
    data[i+600,]=arima.sim(list(order=c(2,0,0),ar=c(0.8,0.1)),n=1000,sd=exp(6.6/2))
    data[i+700,]=arima.sim(list(order=c(2,0,0),ar=c(0.8,0.1)),n=1000,sd=exp(8.6/2))
    data[i+800,]=arima.sim(list(order=c(2,0,0),ar=c(0.8,0.1)),n=1000,sd=exp(8.8/2))
    data[i+900,]=arima.sim(list(order=c(2,0,0),ar=c(0.8,0.1)),n=1000,sd=exp(10.8/2))
  }
  dataout2_25 = addgamma(data,0.25)
  dataout2_30 = addgamma(data,0.3)
  dataout2_35 = addgamma(data,0.35)

  data_pg = LogSmoothGVC(inputdata =data,span,T,J,ch)
  simmout2_25 = LogSmoothGVC(inputdata =dataout2_25,span,T,J,ch)
  simmout2_30 = LogSmoothGVC(inputdata =dataout2_30,span,T,J,ch)
  simmout2_35 = LogSmoothGVC(inputdata =dataout2_35,span,T,J,ch)
  
  
  data_array = array(0,c(40,499,25))
  data2_25 = array(0,c(40,499,25))
  data2_30 = array(0,c(40,499,25))
  data2_35 = array(0,c(40,499,25))
  
  for(i in 1:5){
    data_array[1:20,,i]=data_pg[((i-1)*20+1):(i*20),]
    data_array[21:40,,i]=data_pg[100+(((i-1)*20+1):(i*20)),]
    data_array[1:20,,5+i]=data_pg[200+(((i-1)*20+1):(i*20)),]
    data_array[21:40,,5+i]=data_pg[300+(((i-1)*20+1):(i*20)),]
    data_array[1:20,,10+i]=data_pg[400+(((i-1)*20+1):(i*20)),]
    data_array[21:40,,10+i]=data_pg[500+(((i-1)*20+1):(i*20)),]
    data_array[1:20,,15+i]=data_pg[600+(((i-1)*20+1):(i*20)),]
    data_array[21:40,,15+i]=data_pg[700+(((i-1)*20+1):(i*20)),]
    data_array[1:20,,20+i]=data_pg[800+(((i-1)*20+1):(i*20)),]
    data_array[21:40,,20+i]=data_pg[900+(((i-1)*20+1):(i*20)),]
    

    data2_25[1:20,,i]=simmout2_25[((i-1)*20+1):(i*20),]
    data2_25[21:40,,i]=simmout2_25[100+(((i-1)*20+1):(i*20)),]
    data2_25[1:20,,5+i]=simmout2_25[200+(((i-1)*20+1):(i*20)),]
    data2_25[21:40,,5+i]=simmout2_25[300+(((i-1)*20+1):(i*20)),]
    data2_25[1:20,,10+i]=simmout2_25[400+(((i-1)*20+1):(i*20)),]
    data2_25[21:40,,10+i]=simmout2_25[500+(((i-1)*20+1):(i*20)),]
    data2_25[1:20,,15+i]=simmout2_25[600+(((i-1)*20+1):(i*20)),]
    data2_25[21:40,,15+i]=simmout2_25[700+(((i-1)*20+1):(i*20)),]
    data2_25[1:20,,20+i]=simmout2_25[800+(((i-1)*20+1):(i*20)),]
    data2_25[21:40,,20+i]=simmout2_25[900+(((i-1)*20+1):(i*20)),]
  
    data2_30[1:20,,i]=simmout2_30[((i-1)*20+1):(i*20),]
    data2_30[21:40,,i]=simmout2_30[100+(((i-1)*20+1):(i*20)),]
    data2_30[1:20,,5+i]=simmout2_30[200+(((i-1)*20+1):(i*20)),]
    data2_30[21:40,,5+i]=simmout2_30[300+(((i-1)*20+1):(i*20)),]
    data2_30[1:20,,10+i]=simmout2_30[400+(((i-1)*20+1):(i*20)),]
    data2_30[21:40,,10+i]=simmout2_30[500+(((i-1)*20+1):(i*20)),]
    data2_30[1:20,,15+i]=simmout2_30[600+(((i-1)*20+1):(i*20)),]
    data2_30[21:40,,15+i]=simmout2_30[700+(((i-1)*20+1):(i*20)),]
    data2_30[1:20,,20+i]=simmout2_30[800+(((i-1)*20+1):(i*20)),]
    data2_30[21:40,,20+i]=simmout2_30[900+(((i-1)*20+1):(i*20)),]
    
    data2_35[1:20,,i]=simmout2_35[((i-1)*20+1):(i*20),]
    data2_35[21:40,,i]=simmout2_35[100+(((i-1)*20+1):(i*20)),]
    data2_35[1:20,,5+i]=simmout2_35[200+(((i-1)*20+1):(i*20)),]
    data2_35[21:40,,5+i]=simmout2_35[300+(((i-1)*20+1):(i*20)),]
    data2_35[1:20,,10+i]=simmout2_35[400+(((i-1)*20+1):(i*20)),]
    data2_35[21:40,,10+i]=simmout2_35[500+(((i-1)*20+1):(i*20)),]
    data2_35[1:20,,15+i]=simmout2_35[600+(((i-1)*20+1):(i*20)),]
    data2_35[21:40,,15+i]=simmout2_35[700+(((i-1)*20+1):(i*20)),]
    data2_35[1:20,,20+i]=simmout2_35[800+(((i-1)*20+1):(i*20)),]
    data2_35[21:40,,20+i]=simmout2_35[900+(((i-1)*20+1):(i*20)),]
  }
  
  data1_10 = array(0,c(40,499,25))
  data1_20 = array(0,c(40,499,25))
  data1_30 = array(0,c(40,499,25))
  
  for(i in 1:25){
    data1_10[,,i]=addcon(data_array[,,i],0.1)
    data1_20[,,i]=addcon(data_array[,,i],0.2)
    data1_30[,,i]=addcon(data_array[,,i],0.3)
  }

  resfm[,k] = fmclust(data_array,5)$t
  rescr[,k] = crclust(data_array,5)$t
  resfm1_10[,k] = fmclust(data1_10,5)$t
  resfm1_20[,k] = fmclust(data1_20,5)$t
  resfm1_30[,k] = fmclust(data1_30,5)$t
  
  resfm2_25[,k] = fmclust(data2_25,5)$t
  resfm2_30[,k] = fmclust(data2_30,5)$t
  resfm2_35[,k] = fmclust(data2_35,5)$t
  
  
  
  rescr1_10[,k] = crclust(data1_10,5)$t
  rescr1_20[,k] = crclust(data1_20,5)$t
  rescr1_30[,k] = crclust(data1_30,5)$t
  
  rescr2_25[,k] = crclust(data2_25,5)$t
  rescr2_30[,k] = crclust(data2_30,5)$t
  rescr2_35[,k] = crclust(data2_35,5)$t
  
  print(k)
}


tru=c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5)
rifm=0
ricr=0
rifm1_10 = 0
rifm1_20 = 0
rifm1_30 = 0
rifm2_25 = 0
rifm2_30 = 0
rifm2_35 = 0
ricr1_10 = 0
ricr1_20 = 0
ricr1_30 = 0
ricr2_25 = 0
ricr2_30 = 0
ricr2_35 = 0


for(i in 1:100){
  rifm = rifm + adjustedRandIndex(resfm[1:25,i],tru)/100
  ricr = ricr + adjustedRandIndex(rescr[1:25,i],tru)/100
  
  rifm1_10 = rifm1_10 + adjustedRandIndex(resfm1_10[1:25,i],tru)/100
  rifm1_20 = rifm1_20 + adjustedRandIndex(resfm1_20[1:25,i],tru)/100
  rifm1_30 = rifm1_30 + adjustedRandIndex(resfm1_30[1:25,i],tru)/100
  rifm2_25 = rifm2_25 + adjustedRandIndex(resfm2_25[1:25,i],tru)/100
  rifm2_30 = rifm2_30 + adjustedRandIndex(resfm2_30[1:25,i],tru)/100
  rifm2_35 = rifm2_35 + adjustedRandIndex(resfm2_35[1:25,i],tru)/100
  
  ricr1_10 = ricr1_10 + adjustedRandIndex(rescr1_10[1:25,i],tru)/100
  ricr1_20 = ricr1_20 + adjustedRandIndex(rescr1_20[1:25,i],tru)/100
  ricr1_30 = ricr1_30 + adjustedRandIndex(rescr1_30[1:25,i],tru)/100
  ricr2_25 = ricr2_25 + adjustedRandIndex(rescr2_25[1:25,i],tru)/100
  ricr2_30 = ricr2_30 + adjustedRandIndex(rescr2_30[1:25,i],tru)/100
  ricr2_35 = ricr2_35 + adjustedRandIndex(rescr2_35[1:25,i],tru)/100
  
}

rifm
ricr

rifm1_10 
rifm1_20 
rifm1_30 
ricr1_10 
ricr1_20 
ricr1_30 


rifm2_25 
rifm2_30 
rifm2_35 
ricr2_25 
ricr2_30 
ricr2_35 

