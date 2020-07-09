
library(mclust)
library(fda)
library(fields)
library(clusteval)
### Return Bias-Corrected Log Periodograms
LogPdg = function(inputData){
  data1 = as.matrix(inputData)
  
  ch=dim(data1)[1];
  T=dim(data1)[2];
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
    Rawp[((i-1)*J+1):(i*J)] = (1/T)*(abs(fft(data1[i,])[1:J]))^2
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
    if(min(matrix.Rawpdg[k,])==0){ind=which(matrix.Rawpdg[k,]==0)
    matrix.Rawpdg[k,ind]=1e-30
    }
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
## add an constant to a periodogram
addcon=function(x, rate){
  d = dim(x)
  m = d[1]
  n = d[2] 
  c = runif(m)
  cc = rep(0, m)
  for(i in 1:m){
    if (c[i] < rate){cc[i] = 1}
  }
  for(i in 1:m){
    x[i, ] = x[i, ] + 4 * cc[i]
  }
  return(x)
}

## add an eye blink effect with gamma function
addgamma = function(x, rate){
  d = dim(x)
  m = d[1]
  n = d[2]
  c = (31:270) / 150
  ga = gamma(c) - gamma(1.8)
  ga1 = 15 * ga[240:1]
  ga2 = -7.5 * ga
  ga3 = t(cbind(t(ga1), t(ga2)))
  inter = seq(from =ga3[240], to = ga3[241], by = ((ga3[241] - ga3[240])/(31)))
  gam = c()
  gam[1:240] = ga3[1:240]
  gam[241:270] = inter[2:31]
  gam[271:510] = ga3[241:480]
  for (i in 1:m){
    wn = arima.sim(list(order=c(0, 0, 0)), n = 510)
    gamm = gam + wn
    p = runif(1)
    if(p <= rate){
      x[i, 251:760] = x[i, 251:760] + gamm
    }
  }
  return(x)
}

##calculate the central region
central = function(x){
  d = dim(x)
  m = d[1]
  n = d[2]
  dep = fMBD(x)
  deps = sort(dep)
  p = floor(n / 2)
  dep50 = deps[p]
  xx = matrix(0, m, p)
  flag = 1
  for(i in 1:n){
    if(dep[i] > dep50){
      xx[, flag] = x[, i]
      flag = flag + 1
    }
  }
  xxx = rep(0, m)
  for (i in 1:m){
    xxx[i] = max(xx[i,]) - min(xx[i, ])
  }
  return(sum(xxx))
}

##distance matrix using central region
distcr=function(x){
  d = dim(x)
  p = d[1]
  q = d[2]
  r = d[3]
  dcr = matrix(0,r,r)
  for (i in 1:(r - 1)){
    for(j in (i + 1):r){
       dcr[i, j] = central(cbind(t(x[ , , i]), t(x[ , , j])))
       dcr[j, i] = dcr[i, j]####
    }
  }
for(i in 1:r){dcr[i, i] = 10000}  
  return(dcr)
}

##distance matrix using fm
distfm= function(x){
  d=dim(x)
  p=d[1]
  q=d[2]
  r=d[3]

  dis=matrix(0,r,r)
  for(i in 1:(r-1)){
    for(j in (i+1):r){
      dis[i,j] = diss(fmed(x[,,i]),fmed(x[,,j]))
      dis[j,i] = dis[i,j]
    }
  }
  return(dis)
}

##how many columns that is not empty
notzero=function(x){
  d=dim(x)
  m=d[1]
  n=d[2]
  for(i in 1:m){
    if(sum(x[i,]==rep(0,n))>0){
      z=i-1 ;break
    }
  }
  return(z)
}

###get the functional median curve
fmed=function(x){
  dep=fMBD(t(x))
  n=which.max(dep)
  s=x[n,]
  return(s)
}

##central region clustering algorithm
crclust=function(x,n){
  
  d=dim(x)
  p=d[1]
  q=d[2]
  r=d[3]
  record=matrix(0,r,r-n)
  dd=rep(0,r)
  t=1:r
  sss=matrix(1,1,r)
  s4=array(0,c(r*p,q,r))
  for (i in 1:r){
    s4[1:p,,i]=x[,,i]
  }
  dis=distcr(x)
  vb=1
  for(k in r:(n+1)){
    ss=which(dis==dis[which.min(dis)],arr.ind=T)
    t1=ss[1,]
    if(t1[2]<=t1[1]){
      tt=t1[1]
      t1[1]=t1[2]
      t1[2]=tt}
    dd[vb]=dis[t1[1],t1[2]]
    vb=vb+1
    s4[ (p*sss[t1[1]]+1):(p*(sss[t1[1]]+sss[t1[2]])),,t1[1]]=s4[1:(p*sss[t1[2]]),,  t1[2]]
    s4[,,t1[2]]=matrix(0,r*p,q)
    sss[t1[1]]=sss[t1[1]]+sss[t1[2]]
    sss[t1[2]]=0
    bb=which(t==t[t1[2]])
    t[bb]=t[t1[1]]
    for(i in 1:(r-1)){
      for (j in (i+1):r){
        a=notzero(s4[,,i])
        b=notzero(s4[,,j])
        if (a==0 || b==0){
          dis[i,j]=10000
          dis[j,i]=10000}else{
            dis[i,j]=central( cbind( t(s4[1:a,,i]),t(s4[1:b,,j])  ))
            dis[j,i]=dis[i,j]
          }
        
      }
      record[,r+1-k]=t}
  }
  return(list(tt=sss,t=t,dd=dd,record=record))
  }

##functionalmedian clustering 
fmclust=function(x,n){
  d=dim(x)
  p=d[1]
  q=d[2]
  r=d[3]
  record=matrix(0,r,r-n)
  dd=rep(0,r)
  t=1:r
  sss=matrix(1,1,r)
  s4=array(0,c(r*p,q,r))
  for (i in 1:r){
    s4[1:p,,i]=x[,,i]
  }
  dis=matrix(0,r,r)
  for(i in 1:(r-1)){
    for(j in (i+1):r){
       dis[i,j]=diss(fmed(x[,,i]),fmed(x[,,j]))
       dis[j,i]=dis[i,j]
    }
  }
  
  for(i in 1:r){
    dis[i,i]=10000
  }
  vb=1
  for(k in (r:(n+1))){
    ss=which(dis==dis[which.min(dis)],arr.ind=T)
    t1=ss[,1]
    if(t1[2]<=t1[1]){
      tt=t1[1]
      t1[1]=t1[2]
      t1[2]=tt}
    dd[vb]=dis[t1[1],t1[2]]
    vb=vb+1
    s4[ (p*sss[t1[1]]+1):(p*(sss[t1[1]]+sss[t1[2]])),,t1[1]]=s4[1:(p*sss[t1[2]]),,  t1[2]]
    s4[,,t1[2]]=matrix(0,r*p,q)
    sss[t1[1]]=sss[t1[1]]+sss[t1[2]]
    sss[t1[2]]=0
    bb=which(t==t[t1[2]])
    t[bb]=t[t1[1]]
    for(i in 1:(r-1)){
      for (j in (i+1):r){
        a=notzero(s4[,,i])
        b=notzero(s4[,,j])
        if (a==0 || b==0){
          
          dis[i,j]=10000
          dis[j,i]=10000
          }else{
            
            dis[i,j]=diss(fmed(s4[1:a,,i]),fmed(s4[1:b,,j]))
          dis[j,i]=dis[i,j]
            }
        
      }
    }
    record[,r+1-k]=t
  }
  return(list(tt=sss,t=t,dd=dd,record=record))
}

##functional mean clustering
meanclust=function(x,n){
  d=dim(x)
  p=d[1]
  q=d[2]
  r=d[3]
  dd=rep(0,r)
  record=matrix(0,r,r-n)
  t=1:r
  sss=matrix(1,1,r)
  s4=array(0,c(r*p,q,r))
  for (i in 1:r){
    s4[1:p,,i]=x[,,i]
  }
  dis=matrix(0,r,r)
  for(i in 1:r){
    for(j in 1:r){
      ctb=rbind(colMeans(x[,,i]),colMeans(x[,,j]))
      dis[i,j]=as.matrix(dist(ctb))[1,2]
    }
  }
  
 for(i in 1:r){
   dis[i,i]=10000
 }
  vb=1
  for(k in (r:(n+1))){
    ss=which(dis==dis[which.min(dis)],arr.ind=T)
    t1=ss[,1]
    if(t1[2]<=t1[1]){
      tt=t1[1]
      t1[1]=t1[2]
      t1[2]=tt}
    dd[vb]=dis[t1[1],t1[2]]
    vb=vb+1
    s4[ (p*sss[t1[1]]+1):(p*(sss[t1[1]]+sss[t1[2]])),,t1[1]]=s4[1:(p*sss[t1[2]]),,  t1[2]]
    s4[,,t1[2]]=matrix(0,r*p,q)
    sss[t1[1]]=sss[t1[1]]+sss[t1[2]]
    sss[t1[2]]=0
    bb=which(t==t[t1[2]])
    t[bb]=t[t1[1]]
    for(i in 1:(r-1)){
      for (j in (i+1):r){
        a=notzero(s4[,,i])
        b=notzero(s4[,,j])
        if (a==0 || b==0){
          dis[i,j]=10000
          dis[j,i]=10000}else{
            qw=rbind(colMeans(s4[1:a,,i]),colMeans(s4[1:b,,j]))
            dis[i,j]=as.matrix(dist(qw))[1,2]
            dis[j,i]=dis[i,j]
          }
        
      }
    }
    record[,r+1-k]=t
  }
  return(list(tt=sss,t=t,dd=dd,record=record))
}

meanclust_tvd=function(x,n){
  d=dim(x)
  p=d[1]
  q=d[2]
  r=d[3]
  dd=rep(0,r)
  record=matrix(0,r,r-n)
  t=1:r
  sss=matrix(1,1,r)
  s4=array(0,c(r*p,q,r))
  for (i in 1:r){
    s4[1:p,,i]=x[,,i]
  }
  dis=matrix(0,r,r)
  for(i in 1:r){
    for(j in 1:r){
      dis[i,j]=tvd(colMeans(x[,,i]),colMeans(x[,,j]))
    }
  }
  
  for(i in 1:r){
    dis[i,i]=10000
  }
  vb=1
  for(k in (r:(n+1))){
    ss=which(dis==dis[which.min(dis)],arr.ind=T)
    t1=ss[,1]
    if(t1[2]<=t1[1]){
      tt=t1[1]
      t1[1]=t1[2]
      t1[2]=tt}
    dd[vb]=dis[t1[1],t1[2]]
    vb=vb+1
    s4[ (p*sss[t1[1]]+1):(p*(sss[t1[1]]+sss[t1[2]])),,t1[1]]=s4[1:(p*sss[t1[2]]),,  t1[2]]
    s4[,,t1[2]]=matrix(0,r*p,q)
    sss[t1[1]]=sss[t1[1]]+sss[t1[2]]
    sss[t1[2]]=0
    bb=which(t==t[t1[2]])
    t[bb]=t[t1[1]]
    for(i in 1:(r-1)){
      for (j in (i+1):r){
        a=notzero(s4[,,i])
        b=notzero(s4[,,j])
        if (a==0 || b==0){
          dis[i,j]=10000
          dis[j,i]=10000}else{
            dis[i,j]=tvd(colMeans(s4[1:a,,i]),colMeans(s4[1:b,,j]))
            dis[j,i]=dis[i,j]
          }
        
      }
    }
    record[,r+1-k]=t
  }
  return(list(tt=sss,t=t,dd=dd,record=record))
}

meanclust_llrdls=function(x,n){
  d=dim(x)
  p=d[1]
  q=d[2]
  r=d[3]
  record=matrix(0,r,r-n)
  dd=rep(0,r)
  t=1:r
  sss=matrix(1,1,r)
  s4=array(0,c(r*p,q,r))
  for (i in 1:r){
    s4[1:p,,i]=x[,,i]
  }
  dis=matrix(0,r,r)
  for(i in 1:(r-1)){
    for(j in (i+1):r){
      dis[i,j] = diss_llr(colMeans(x[,,i]),colMeans(x[,,j]))
      dis[j,i]=dis[i,j]
    }
  }
  
  for(i in 1:r){
    dis[i,i]=10000
  }
  vb=1
  for(k in (r:(n+1))){
    ss=which(dis==dis[which.min(dis)],arr.ind=T)
    t1=ss[,1]
    if(t1[2]<=t1[1]){
      tt=t1[1]
      t1[1]=t1[2]
      t1[2]=tt}
    dd[vb]=dis[t1[1],t1[2]]
    vb=vb+1
    s4[ (p*sss[t1[1]]+1):(p*(sss[t1[1]]+sss[t1[2]])),,t1[1]]=s4[1:(p*sss[t1[2]]),,  t1[2]]
    s4[,,t1[2]]=matrix(0,r*p,q)
    sss[t1[1]]=sss[t1[1]]+sss[t1[2]]
    sss[t1[2]]=0
    bb=which(t==t[t1[2]])
    t[bb]=t[t1[1]]
    for(i in 1:(r-1)){
      for (j in (i+1):r){
        a=notzero(s4[,,i])
        b=notzero(s4[,,j])
        if (a==0 || b==0){
          dis[i,j]=10000
          dis[j,i]=10000}else{
            dis[i,j]=diss_llr(colMeans(s4[1:a,,i]),colMeans(s4[1:b,,j]))
            dis[j,i]=dis[i,j]
          }
        
      }
    }
    record[,r+1-k]=t
  }
  return(list(tt=sss,t=t,dd=dd,record=record))
}

fmclust_tvd=function(x,n){
  d=dim(x)
  p=d[1]
  q=d[2]
  r=d[3]
  record=matrix(0,r,r-n)
  dd=rep(0,r)
  t=1:r
  sss=matrix(1,1,r)
  s4=array(0,c(r*p,q,r))
  for (i in 1:r){
    s4[1:p,,i]=x[,,i]
  }
  dis=matrix(0,r,r)
  for(i in 1:(r-1)){
    for(j in (i+1):r){
      dis[i,j]=tvd(fmed(x[,,i]),fmed(x[,,j]))
      dis[j,i]=dis[i,j]
    }
  }
  
  for(i in 1:r){
    dis[i,i]=10000
  }
  vb=1
  for(k in (r:(n+1))){
    ss=which(dis==dis[which.min(dis)],arr.ind=T)
    t1=ss[,1]
    if(t1[2]<=t1[1]){
      tt=t1[1]
      t1[1]=t1[2]
      t1[2]=tt}
    dd[vb]=dis[t1[1],t1[2]]
    vb=vb+1
    s4[ (p*sss[t1[1]]+1):(p*(sss[t1[1]]+sss[t1[2]])),,t1[1]]=s4[1:(p*sss[t1[2]]),,  t1[2]]
    s4[,,t1[2]]=matrix(0,r*p,q)
    sss[t1[1]]=sss[t1[1]]+sss[t1[2]]
    sss[t1[2]]=0
    bb=which(t==t[t1[2]])
    t[bb]=t[t1[1]]
    for(i in 1:(r-1)){
      for (j in (i+1):r){
        a=notzero(s4[,,i])
        b=notzero(s4[,,j])
        if (a==0 || b==0){
          
          dis[i,j]=10000
          dis[j,i]=10000
        }else{
          
          dis[i,j]=tvd(fmed(s4[1:a,,i]),fmed(s4[1:b,,j]))
          dis[j,i]=dis[i,j]
        }
        
      }
    }
    record[,r+1-k]=t
  }
  return(list(tt=sss,t=t,dd=dd,record=record))
}
#standerize the clusters
stdclu=function(x){
   len=length(x)
   s=rep(0,len)
   cl=rep(0,len)
   for(i in 1:len){
      cl[i]=(x[i]==i)
   }
   num=sum(cl)
   t=which(cl==1)
   for(i in 1:num){
      p=which(x==t[i])
      s[p]=i
   }
return(s)
}
                
#####TVD##################
mi=function(x,y){
  xymin=rep(0,length(x))
  for(i in 1:length(x)){
    xymin[i]=min(x[i],y[i])
  }
  return(xymin)
}

tvd=function(x,y){
  x1=x-min(x)
  y1=y-min(y)
  x2=x1/sum(x1)
  y2=y1/sum(y1)
  tv=1-sum(mi(x2,y2))
  return(tv)
}




ranktest=function(dataX,dataY,dataRef)
{
  n=dim(dataX)[2];
  m=dim(dataY)[2];
  r=dim(dataRef)[2];
  order=integer(n+m);
  for(i in 1:m)
  {
    sample=cbind(dataRef,dataY[,i]);
    result=fbplot(sample,plot=F);
    order[i]=sum(result$depth[1:r]<=result$depth[r+1])
  }
  for(i in 1:n)
  {
    sample=cbind(dataRef,dataX[,i]);
    result=fbplot(sample,plot=F);
    order[i+m]=sum(result$depth[1:r]<=result$depth[r+1])
  }
  rk=sort.int(order,index.return = T)$ix;
  w=sum(rk[1:m]);
  return(w)
}

#combination

combinat=function(n,p){
  
  if (n<p){combinat=0}
  
  else {combinat=exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))}
  
}







#BD2

fBD2=function(data){
  
  p=dim(data)[1]
  
  n=dim(data)[2]
  
  rmat=apply(data,1,rank)
  
  down=apply(rmat,1,min)-1
  
  up=n-apply(rmat,1,max)
  
  (up*down+n-1)/combinat(n,2)
  
  
  
}



fMBD=function(data){
  
  p=dim(data)[1]
  
  n=dim(data)[2]
  
  rmat=apply(data,1,rank)
  
  down=rmat-1
  
  up=n-rmat
  
  (rowSums(up*down)/p+n-1)/combinat(n,2)
  
}


diss=function(x1, x2){ return(sqrt(sum((x1 - x2) ^ 2))) }



diss_llr = function(x,y){
  w = function(u,v){
    return(0.5*(log((0.5*u+0.5*v)/v  - 0.5*log(u/v) ) ))
  } 
  x=exp(x)
  y=exp(y)
  return(sum(w(x,y)+w(y,x)))
}


tvd = function(a,b){
  anorm = a - min(a)
  anorm = anorm/mean(anorm)
  bnorm = b - min(b)
  bnorm = bnorm/mean(bnorm)
  int = mean((pmin(anorm,bnorm)+pmin(anorm,bnorm))/2)
  return(1-int)
}

