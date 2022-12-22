### Return Bias-Corrected Log Periodograms
LogPdg = function(inputData){
  data1 = as.matrix(inputData)
  
  ch = dim(data1)[1];
  T = dim(data1)[2];
  J = floor((T+1)/2)-1;
  
  ##Calculate log periodogram
  logp = rep(0,ch*J);
  for(i in 1:ch) {
    logp[((i-1)*J+1):(i*J)] = log((1/T)*(abs(fft(data1[i,])[1:J]))^2);
  }
  
  return(logp+0.57721);
}

### Periodogram not log periodograms
RawPdg = function(inputData){
  data1 = as.matrix(inputData)
  ch=dim(data1)[1];
  T=dim(data1)[2];
  J = floor((T+1)/2)-1;
  Rawp = rep(0,ch*J);
  for(i in 1:ch) {
    Rawp[((i-1)*J+1):(i*J)] =(1/T)*(abs(fft(data1[i,])[1:J]))^2
  }
  return(Rawp);
}

### Convert Vector (Log Periodogram) to Matrix
### Columns are objects, with time series are stacked in rows.
VecToMatrix = function(inputData, ch, T, J){
  matrixp=matrix(0,T/2-1,ch)
  J=T/2-1
  for(i in 1:ch){
    matrixp[,i]= inputData[((i-1)*J+1):(i*J)]
  }
  return(matrixp)
}

#define neighborhood
neighborhood = function(x,x.star,y,bandwidth){
  index=c(1:length(x))
  indexn = index[(x<(x.star+bandwidth+1)) & (x>(x.star-bandwidth-1))]
  xn = x[indexn]
  yn = y[indexn]
  out=list(xn,yn)
  names(out) = c("xn", "yn")
  return(out)
}

#Calculate local estimate
localEstimates=function(xn,yn,bandwith){
  #boxcar smoother W = 1/(2p+1)
  w=rep(1/(2*bandwith + 1),length(xn))
  f.hat = as.numeric((t(w)%*%yn)/sum(w) )
  #f.hat = sum(yn)/(2*bandwith + 1)
  return(f.hat)
}

### Average Channel by using local channel around it.
AvgChannel = function(inputData,channel){
  for(i in 1:160){
    inputData[i,channel] = mean(c(inputData[i,(channel-5):(channel-1)],inputData[i,(channel+1):(channel+5)]))
  }
  return(inputData)
}

SmoothGVC = function(matrix.Rawpdg,span){
  GVCp = rep(0,length(span))
  span.min = 0
  n.trial = dim(matrix.Rawpdg)[1]
  matrix.smooth = matrix(0,dim(matrix.Rawpdg)[1],dim(matrix.Rawpdg)[2])
  y.estimate = matrix(0,length(span),dim(matrix.Rawpdg)[2])
  for(k in 1:n.trial){
    for(i in 1:(length(span))){
      bandwidth = span[i]
      first.trial.process = ProcessRawPdg(matrix.Rawpdg[k,],bandwidth)
      x.feq = seq(1,length(first.trial.process),1)
      for(j in 1:length(matrix.Rawpdg[k,])){
        window = neighborhood(x=x.feq,x.star=x.feq[j]+bandwidth,y=first.trial.process,bandwidth)
        y.estimate[i,j]=localEstimates(window$xn,window$yn,bandwidth)
      }
      GVCp[i] = CalculateGVC(matrix.Rawpdg[k,], y.estimate[i,], bandwidth)
    }
    matrix.smooth[k,] = y.estimate[which(GVCp == min(GVCp)),]
  }  
  return(matrix.smooth)
}

CalculateGVC = function(f, f.hat, bandwidth){
  M = length(f)
  sum = 0
  q = c(0.5,rep(1,M-2),0.5)
  for(i in 1:M){
    num = -log(f[i]/f.hat[i])+(f[i]-f.hat[i])/f.hat[i]
    dem = (1 - (1/(2*bandwidth + 1)))^2
    sum = sum + q[i]*(num/dem)
  }
  return(sum/M)
}

ProcessRawPdg = function(pdg,bandwidth){
  temp = rep(0,bandwidth)
  end = length(pdg)
  temp = rev(pdg[2:(bandwidth+1)])
  pdg.final = c(temp,pdg)
  temp = rev(pdg[(end-bandwidth-1):(end-1)])
  pdg.final = c(pdg.final,temp)
  return(pdg.final)
}

LogSmoothGVC = function(inputdata,span,T,J,ch){
  Channel.RawPdg = RawPdg(inputdata)
  Channel.RawPdg=t(VecToMatrix(Channel.RawPdg,ch,T,J))
  #Channel.RawPdg=AvgChannel(Channel.RawPdg,61,ch) #filter out 61 Hz, not needed
  return(log(SmoothGVC(Channel.RawPdg,span)) + 0.57721)
}

########################clustering##################

# add a constant to the log-periodogram
addcon=function(x,rate){
  dimension=dim(x)
  m = dimension[1]
  n = dimension[2] 
  ran_v = runif(m)
  noise = rep(0,m)
for(i in 1:m){
  if (ran_v[i]<rate){noise[i]=1}
}
for(i in 1:m){
  x[i,] = x[i,]+4*noise[i]
}
return(x)
}

## add an eye blink effect with gamma function
addgamma = function(x, rate){
  dimension=dim(x)
  m = dimension[1]
  n = dimension[2]
  c = (31:270)/150
  ga = gamma(c)-gamma(1.8)
  ga_1 = 15*ga[240:1]
  ga_2 = -7.5*ga
  ga_3 = t(cbind(t(ga_1),t(ga_2)))
  
  inter = seq(from=ga_3[240], to = ga_3[241], by = ((ga_3[241] - ga_3[240])/(31)))
  gam = c()
  gam[1:240] = ga_3[1:240]
  gam[241:270] = inter[2:31]
  gam[271:510] = ga_3[241:480]
  for (i in 1:m){
    wn = arima.sim(list(order=c(0,0,0)),n=510)
    eye = gam+wn
    p = runif(1)
    if(p<=rate){
      x[i,251:760] = x[i,251:760]+eye
    }
  }
  return(x)
}

##calculate the  central region

central = function(x){
  dimension = dim(x)
  m = dimension[1]
  n = dimension[2]
  depth = fMBD(x)
  dep_sorted = sort(depth)
  p = floor(n/2)
  dep50 = dep_sorted[p]
  num = matrix(0,m,p)
  flag=1
  for(i in 1:n){
    if(depth[i]>dep50){
      num[,flag] = x[,i]
      flag = flag+1
    }
  }
  area = rep(0,m)
  for (i in 1:m){
    area[i] = max(num[i,])-min(num[i,])
  }
  return(sum(area))
}


##distance matrix using central region
distcr = function(x){
  dimension = dim(x)
  p = dimension[1]
  q = dimension[2]
  r = dimension[3]
  dcr = matrix(0,r,r)
  for (i in 1:(r-1)){
    for(j in (i+1):r){
       dcr[i,j] = central(cbind(t(x[,,i]),t(x[,,j])))
       dcr[j,i] = dcr[i,j]####
    }
  }
for(i in 1:r){dcr[i,i] = 10000}  
  return(dcr)
}

##how many columns that is not empty
notzero = function(x){
  dimension = dim(x)
  m = dimension[1]
  n = dimension[2]
  for(i in 1:m){
    if(sum(x[i,]==rep(0,n))>0){
      z = i-1 ;break
    }
  }
  return(z)
}

###get the functional median curve
fmed = function(x){
  dep = fMBD(t(x))
  index = which.max(dep)
  s = x[index,]
  return(s)
}

##robust clustering algorithms
robustcluster = function(x, n, method){
  dimension = dim(x)
  n_rep = dimension[1]
  n_freq = dimension[2]
  n_clu = dimension[3]
  min_dis = rep(0, n_clu)
  res = 1:n_clu
  n_element_in_each_clu = matrix(1, 1, n_clu)
  clu_mtx = array(0,c(n_clu*n_rep, n_freq,n_clu))
  for (i in 1:n_clu){
    clu_mtx[1:n_rep, , i] = x[, , i]
  }
  if(method == 'cr'){dis = distcr(x)}
  if(method == 'fm'){
    dis=matrix(0,n_clu,n_clu)
    for(i in 1:(n_clu-1)){
      for(j in (i+1):n_clu){
        dis[i,j]=diss(fmed(x[,,i]),fmed(x[,,j]))
      }
    }
    for(i in 1:n_clu){
      dis[i,i] = 10000
    }
  }
  if(method == 'mean'){
    dis=matrix(0,n_clu,n_clu)
    for(i in 1:n_clu){
      for(j in 1:n_clu){
        ms = rbind(colMeans(x[,,i]), colMeans(x[,,j]))
        dis[i,j] = as.matrix(dist(ms))[1,2]
      }
    }
    for(i in 1:n_clu){
      dis[i,i] = 10000
    }
  }
  rcd = 1
  for(k in n_clu:(n+1)){
    index_min = which(dis==min(dis), arr.ind=TRUE)
    index_min = index_min[1,]
    if(index_min[2] <= index_min[1]){
      temp = index_min[1]
      index_min[1] = index_min[2]
      index_min[2] = temp}
    min_dis[rcd] = dis[index_min[1], index_min[2]]
    rcd = rcd+1
    clu_mtx[(n_rep*n_element_in_each_clu[index_min[1]]+1):(n_rep*(n_element_in_each_clu[index_min[1]]+n_element_in_each_clu[index_min[2]])), , index_min[1]]=clu_mtx[1:(n_rep*n_element_in_each_clu[index_min[2]]), , index_min[2]]
    clu_mtx[ , ,index_min[2]] = matrix(0,n_clu*n_rep, n_freq)
    n_element_in_each_clu[index_min[1]] = n_element_in_each_clu[index_min[1]]+n_element_in_each_clu[index_min[2]]
    n_element_in_each_clu[index_min[2]] = 0
    index = which(res==res[index_min[2]])
    res[index] = res[index_min[1]]
    for(i in 1:(n_clu-1)){
      for (j in (i+1):n_clu){
        a = notzero(clu_mtx[, , i])
        b = notzero(clu_mtx[, , j])
        if (a==0 || b==0){
          dis[i,j] = 10000
          dis[j,i] = 10000}
        else{
          if(method == 'cr'){
            dis[i,j] = central(cbind(t(clu_mtx[1:a, , i]), t(clu_mtx[1:b, , j])))
            dis[j,i] = dis[i,j]}
          if(method == 'fm'){
            dis[i,j]=diss(fmed(clu_mtx[1:a,,i]), fmed(clu_mtx[1:b,,j]))
            dis[j,i] = dis[i,j]}
          if(method == 'mean'){
            mss = rbind(colMeans(s4[1:a,,i]), colMeans(s4[1:b,,j]))
            dis[i,j] = as.matrix(dist(mss))[1,2]
            dis[j,i] = dis[i,j]
          }
        }
      }
    }
  }
  return(list(n_element_in_each_clu=n_element_in_each_clu, res=res, min_dis=min_dis))
}

#standerize the clusters
stdclu = function(x){
   len = length(x)
   s = rep(0,len)
   cl = rep(0,len)
   for(i in 1:len){
      cl[i] = (x[i] == i)
   }
   num = sum(cl)
   t = which(cl == 1)
   for(i in 1:num){
      p = which(x == t[i])
      s[p] = i
   }
return(s)
}
                
####rank sum test
ranktest = function(dataX, dataY, dataRef)
{
  n = dim(dataX)[2];
  m = dim(dataY)[2];
  r = dim(dataRef)[2];
  order = integer(n+m);
  for(i in 1:m)
  {
    sample = cbind(dataRef, dataY[,i]);
    result = fbplot(sample, plot=F);
    order[i] = sum(result$depth[1:r] <= result$depth[r+1])
  }
  for(i in 1:n)
  {
    sample = cbind(dataRef,dataX[,i]);
    result = fbplot(sample, plot=F);
    order[i+m] = sum(result$depth[1:r] <= result$depth[r+1])
  }
  rk = sort.int(order, index.return = T)$ix;
  w = sum(rk[1:m]);
  return(w)
}

#combination
combinat = function(n,p){
  
  if (n<p){combinat=0}
  else {combinat = exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))}
  
}


fBD2 = function(data){
  p = dim(data)[1]
  n = dim(data)[2]
  rmat = apply(data, 1, rank)
  down = apply(rmat, 1, min)-1
  up = n-apply(rmat, 1, max)
  (up*down+n-1) / combinat(n,2)
  }

########MBD###########################################################
fMBD = function(data){
  p = dim(data)[1]
  n = dim(data)[2]
  rmat = apply(data,1,rank)
  down = rmat-1
  up = n-rmat
  (rowSums(up*down)/p+n-1)/combinat(n,2)
}

diss=function(x1, x2){ return(sqrt(sum((x1 - x2) ^ 2))) }




