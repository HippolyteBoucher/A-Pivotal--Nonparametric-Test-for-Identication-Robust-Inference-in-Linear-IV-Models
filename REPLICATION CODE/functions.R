library("np")
library("expm")
library("MASS")
library("foreach")
library("doParallel")
library("ggplot2")


#####################################################################################
##### Build covariance estimator ####################################################
#####################################################################################



b_yin<-function(y,x,z,zc,h=F){
  n<-length(y)
  yy<-cbind(y,x)
  if (h==F){
    h<-rep(1.06*n^(-1/5),dim(z)[2])
  }
  kbar<-npksum(txdat=z,exdat=t(zc),bws=h)$ksum
  ybar<-npksum(txdat=z,exdat=t(zc),tydat=yy,bws=h)$ksum/kbar
  yf<-t(apply(yy,1,function(e) e-ybar))
  ome<-as.numeric(npksum(weight=yf,tydat=yf,txdat=z,exdat=t(zc),bws=h)$ksum)/kbar
  return(ome)
}


b_covar_hom<-function(y,x,z,h=F,type="lin",omet){
  
  Y<-cbind(y,x)
  l<-dim(Y)[2]-1
  n<-dim(Y)[1]
  
  if (h==F){
    h<-1.06*n^(-1/5)
  }
  
  if (type=="lin"){
    res<-lm(Y~z)$residuals
    ome<-t(res)%*%res/n
  } else if (type=="np"){
    ome<-matrix(rowSums(vapply(1:n,function(e) b_yin(y,x,z,z[e,],h),FUN.VALUE=numeric(4))),l+1,l+1)/n
  } else if (type=="gls"){
    res<-lm(Y~z+0)$residuals
    res2<-log(res^2)
    resc<-res[,1]*res[,2]
    resfit<-lm(res2~z+0)$fitted.values
    ome<-matrix(c(mean(exp(resfit[,1])),mean(resc),mean(resc),mean(exp(resfit[,2]))),2,2)
  } else if (type=="true"){
    ome<-omet
  }
  
  return(ome)
  
}



b_covar_het<-function(y,x,z,h=F,type="lin",omet){
  
  Y<-cbind(y,x)
  z<-as.matrix(z)
  n<-dim(Y)[1]
  l1<-dim(Y)[2]
  
  if (h==F){
    h<-1.06*n^(-1/5)
  }
  
  if (type=="np"){
    ome<-vapply(1:n, function(ee) b_yin(y,x,z,z[ee,],h),FUN.VALUE=numeric(4))
  } else if (type=="lin") {
    res<-lm(Y~z+0)$residuals
    ome<-vapply(1:n, function(ee) res[ee,]%*%t(res[ee,]),FUN.VALUE=numeric(4))
  } else if (type=="gls"){
    res<-stats::lm(Y~z+0)$residuals
    res2<-res^2
    resc<-res[,1]*res[,2]
    resfit<-lm(cbind(res2,resc)~poly(z,3))$fitted.values
    ome<-rbind(as.numeric(resfit[,1]),as.numeric(resfit[,3]),as.numeric(resfit[,3]),as.numeric(resfit[,2]))
  } else if (type=="true"){
    ome<-vapply(1:n, function(e) omet*(1+z[e,1]^2)/2,FUN.VALUE=numeric(4))
  } else if (type=="gls-np"){
    res<-stats::lm(Y~z+0)$residuals
    res2<-res^2
    resc<-res[,1]*res[,2]
    resfit<-lm(cbind(res2,resc)~poly(z,2))$fitted.values
    ome<-rbind(as.numeric(resfit[,1]),as.numeric(resfit[,3]),as.numeric(resfit[,3]),as.numeric(resfit[,2]))
    ind_ninv<-which(((ome[1,]+ome[4,])/2-abs(ome[2,])<0.01)|(ome[1,]*ome[4,]-ome[3,]^2<0.01))
    if (length(ind_ninv)>0){
      ome[,ind_ninv]<-vapply(ind_ninv, function(ee) b_yin(y,x,z,z[ee,],h),FUN.VALUE=numeric(4))
    }
  }
  
  
  return(ome)
  
}



#####################################################################################
##### Build P-Values given Statistics ###############################################
#####################################################################################


b_sim<-function(s_sim,s_sim_pz,s_sim_W,T_pz,T_W,Ts){
  
  clr<-t(s_sim_pz)%*%s_sim_pz-min(Re(eigen(t(cbind(s_sim_pz,T_pz))%*%cbind(s_sim_pz,T_pz))$values))
  icm<-t(s_sim_W)%*%s_sim
  cicm<-t(s_sim_W)%*%s_sim-min(Re(eigen(t(cbind(s_sim_W,T_W))%*%cbind(s_sim,Ts))$values))
  
  return(c(clr,icm,cicm))
}



b_pval<-function(stats,xnorm,T_pz,T_W,Ts,z){
  
  l<-dim(as.matrix(T_pz))[2]
  k<-dim(z)[2]
  m<-dim(xnorm)[2]
  
  pval<-rep(0,7)
  
  s_sim<-xnorm
  s_sim_pz<-lm(xnorm~z+0)$fitted.values
  s_sim_W<-npksum(tydat=xnorm,txdat=z,bws=rep(1,k))$ksum/n
  
  dstats<-vapply(1:m, function(e) b_sim(s_sim[,e],s_sim_pz[,e],s_sim_W[,e],T_pz,T_W,Ts),FUN.VALUE=numeric(3))
  
  pval[1]<-pchisq(as.numeric(stats[1]),df=k,lower.tail=F)
  pval[2]<-pchisq(as.numeric(stats[2]),df=l,lower.tail=F)
  pval[3]<-mean(as.numeric(stats[3])<dstats[1,])
  pval[4]<-mean(as.numeric(stats[4])<dstats[2,])
  pval[5]<-pchisq(as.numeric(stats[5]),df=l,lower.tail=F)
  pval[6]<-mean(as.numeric(stats[6])<dstats[3,])
  pval[7]<-pchisq(as.numeric(stats[7]),df=l,lower.tail=F)
  
  return(as.numeric(pval))
}


b_pval_para<-function(stats,xnorm,T_pz,T_W,Ts,z){
  
  l<-dim(as.matrix(T_pz))[2]
  k<-dim(z)[2]
  m<-dim(xnorm)[2]
  
  pval<-rep(0,7)
  
  s_sim<-xnorm
  s_sim_pz<-lm(xnorm~z+0)$fitted.values
  s_sim_W<-npksum(tydat=xnorm,txdat=z,bws=rep(1,k))$ksum/n
  
  dstats<-foreach::foreach(ind=1:m,.export="b_sim",.combine=rbind) %dopar% {
    b_sim(s_sim[,ind],s_sim_pz[,eind],s_sim_W[,ind],T_pz,T_W,Ts)
  }
  
  pval[1]<-pchisq(as.numeric(stats[1]),df=k,lower.tail=F)
  pval[2]<-pchisq(as.numeric(stats[2]),df=l,lower.tail=F)
  pval[3]<-mean(as.numeric(stats[3])<dstats[1,])
  pval[4]<-mean(as.numeric(stats[4])<dstats[2,])
  pval[5]<-pchisq(as.numeric(stats[5]),df=l,lower.tail=F)
  pval[6]<-mean(as.numeric(stats[6])<dstats[3,])
  pval[7]<-pchisq(as.numeric(stats[7]),df=l,lower.tail=F)
  
  return(as.numeric(pval))
}


#####################################################################################
##### Build Stats Homoskedastic Case ################################################
#####################################################################################


b_YP<-function(y,x,z){
  Y<-cbind(y,x)
  Y_pz<-lm(Y~z+0)$fitted.values
  Y_W<-npksum(txdat=z,tydat=Y,bws=rep(1,dim(z)[2]))$ksum/dim(Y)[1]
  return(cbind(Y_pz,Y_W))
}




b_pval_hom<-function(Y,Y_pz,Y_W,z,beta0,ome,b2sls,sig,xnorm,para=F){
  
  l<-length(beta0)
  
  b0<-c(1,-beta0)
  A0<-rbind(beta0,diag(l))
  
  Ss<-Y%*%b0/as.numeric(sqrt(t(b0)%*%ome%*%b0))
  S_pz<-Y_pz%*%b0/as.numeric(sqrt(t(b0)%*%ome%*%b0))
  S_W<-Y_W%*%b0/as.numeric(sqrt(t(b0)%*%ome%*%b0))
  
  
  if (l>1){
    Ts<-Y%*%solve(ome)%*%A0%*%expm::sqrtm(solve(t(A0)%*%solve(ome)%*%A0))
    T_pz<-Y_pz%*%solve(ome)%*%A0%*%expm::sqrtm(solve(t(A0)%*%solve(ome)%*%A0))
    T_W<-Y_W%*%solve(ome)%*%A0%*%expm::sqrtm(solve(t(A0)%*%solve(ome)%*%A0))
  } else {
    Ts<-Y%*%solve(ome)%*%A0/as.numeric(sqrt(t(A0)%*%solve(ome)%*%A0))
    T_pz<-Y_pz%*%solve(ome)%*%A0/as.numeric(sqrt(t(A0)%*%solve(ome)%*%A0))
    T_W<-Y_W%*%solve(ome)%*%A0/as.numeric(sqrt(t(A0)%*%solve(ome)%*%A0))
  }
  
  
  ar<-t(S_pz)%*%S_pz
  
  S_pzt<-lm(S_pz~T_pz+0)$fitted.values
  klm<-t(S_pzt)%*%S_pzt
  if (l>1){
    clr<-ar-min(Re(eigen(t(cbind(S_pz,T_pz))%*%cbind(S_pz,T_pz))$values))
  } else {
    clr<-(ar-t(T_pz)%*%T_pz+sqrt((t(S_pz)%*%S_pz-t(T_pz)%*%T_pz)^2+4*(t(S_pz)%*%T_pz)^2))/2
  }
  icm<-t(S_W)%*%Ss
  S_Wt<-lm(Ss~T_W+0)$fitted.values
  kicm<-t(S_Wt)%*%S_Wt
  if (l>1){
    cicm<-icm-min(Re(eigen(t(cbind(Ss,Ts))%*%cbind(S_W,T_W))$values))
  } else {
    cicm<-(icm-t(T_W)%*%Ts+sqrt((t(S_W)%*%Ss-t(T_W)%*%Ts)^2+4*(t(S_W)%*%Ts)^2))/2
  }
  
  A<-rbind(rep(0,l),diag(l))
  wald<-t(b2sls-beta0)%*%t(A)%*%t(Y_pz)%*%Y_pz%*%A%*%(b2sls-beta0)/sig
  
  stats<-c(ar,klm,clr,icm,kicm,cicm,wald)
  if (para==F){
    pval<-b_pval(stats,xnorm,T_pz,T_W,Ts,z)
  } else {
    pval<-b_pval_para(stats,xnorm,T_pz,T_W,Ts,z)
  }
  return(pval)
}



#####################################################################################
##### Build Pval Heteroskedastic Case ###############################################
#####################################################################################


b_pval_het<-function(y,x,z,beta0,ome,b2sls,sig,xnorm,para=F){
  
  l<-dim(as.matrix(x))[2]
  k<-dim(z)[2]
  n<-length(y)
  b0<-c(1,-beta0)
  A0<-rbind(beta0,diag(l))
  Y<-cbind(y,x)
  
  Ss<-vapply(1:n, function(e) t(Y[e,])%*%b0/as.numeric(sqrt(t(b0)%*%matrix(ome[,e],l+1,l+1)%*%b0)),FUN.VALUE=numeric(1))
  if (l>1){
    Ts<-vapply(1:n, function(e) t(Y[e,])%*%MASS::ginv(matrix(ome[,e],l+1,l+1))%*%A0%*%expm::sqrtm(solve(t(A0)%*%matrix(ome[,e],l+1,l+1)%*%A0)),FUN.VALUE=numeric(l))
  } else {
    Ts<-vapply(1:n, function(e) t(Y[e,])%*%MASS::ginv(matrix(ome[,e],l+1,l+1))%*%A0/as.numeric(sqrt(t(A0)%*%matrix(ome[,e],l+1,l+1)%*%A0)),FUN.VALUE=numeric(1))
  }
  S_pz<-lm(Ss~z+0)$fitted.values
  S_W<-npksum(txdat=z,tydat=Ss,bws=rep(1,k))$ksum/n
  T_pz<-lm(Ts~z+0)$fitted.values
  T_W<-npksum(txdat=z,tydat=Ts,bws=rep(1,k))$ksum/n
  
  ar<-t(S_pz)%*%S_pz
  S_Tpz<-lm(Ss~T_pz+0)$fitted.values
  lm<-t(S_Tpz)%*%S_Tpz
  if (l>1){
    clr<-ar-min(Re(eigen(t(cbind(S_pz,T_pz))%*%cbind(S_pz,T_pz))$values))
  } else {
    clr<-(ar-t(T_pz)%*%T_pz+sqrt((t(S_pz)%*%S_pz-t(T_pz)%*%T_pz)^2+4*(t(S_pz)%*%T_pz)^2))/2
  }
  icm<-t(S_W)%*%Ss
  S_TW<-lm(Ss~T_W+0)$fitted.values
  kicm<-t(S_TW)%*%S_TW
  if (l>1){
    cicm<-icm-min(Re(eigen(t(cbind(S_W,T_W))%*%cbind(Ss,Ts))$values))
  } else {
    cicm<-(icm-t(T_W)%*%Ts+sqrt((t(S_W)%*%Ss-t(T_W)%*%Ts)^2+4*(t(S_W)%*%Ts)^2))/2
  }
  A<-rbind(rep(0,l),diag(l))
  Y_pz<-lm(Y~z+0)$fitted.values
  Y_pz_sig<-t(vapply(1:n, function(e) Y_pz[e,]*sig[e],FUN.VALUE=numeric(l+1)))
  wald<-t(b2sls-beta0)%*%t(A)%*%t(Y_pz)%*%Y_pz%*%A%*%solve(t(A)%*%t(Y_pz)%*%Y_pz_sig%*%A)%*%t(A)%*%t(Y_pz)%*%Y_pz%*%A%*%(b2sls-beta0)
  
  stats<-c(ar,lm,clr,icm,kicm,cicm,wald)
  if (para==F){
    pval<-b_pval(stats,xnorm,T_pz,T_W,Ts,z)
  } else {
    pval<-b_pval_para(stats,xnorm,T_pz,T_W,Ts,z)
  }
  
  return(pval)
}


#####################################################################
############### Simulation Functions ################################
#####################################################################


# 1st stage condeitional mean


fpi<-function(z,np,w1){
  z<-as.matrix(z)
  n<-dim(z)[1]
  k<-dim(z)[2]
  if (np=="lin"){
    tfpi<-z[,1]
  } else if (np=="nlin"){
    #tfpi<-(exp(z[,1])-exp(1/2))/sqrt(exp(1)*(exp(1)-1))
    tfpi<-(z[,1]+z[,2]+z[,1]*z[,2]+z[,2]^2+z[,1]^2+z[,1]^2*z[,2]^2-3)/sqrt(26)
  } else if (np=="polyp"){
    tfpi<-(z[,1]^2-1)/sqrt(3)
  } else if (np=="polysp"){
    #tfpi<-z[,1]*(1-z[,2]+z[,2]^3/3)*sqrt(3/2)
    #tfpi<-z[,1]*(1-z[,2])/sqrt(2)
    tfpi<-(z[,1]+z[,2]^2-1)/sqrt(4)
    #tfpi<-(z[,1]+z[,2])/sqrt(2)
    #tfpi<-z[,1]*z[,2]
  } else if (np=="spec"){
    tfpi<-rowMeans(z)/sqrt(k)
  }
  
  if (w1=="false"){
    tfpi<-tfpi
  } else if (w1=="semi"){
    tfpi<-tfpi/n^(1/4)
  } else if (w1=="semi2"){
    tfpi<-tfpi/n^sqrt(1/5)
  } else if (w1=="weak"){
    tfpi<-tfpi/n^(1/2)
  }
  
  return(matrix(tfpi,n,1))
}

# Simulate data and recover pvalues for 7 tests for a grid of beta0

b_sim_pval<-function(gridbeta0,xnorm,n,l,k,cuv,np="lin",w1="false",type="lin"){
  
  z<-matrix(rnorm(n*k),n,k)
  if (np=="polyp"|np=="lin"){
    z<-as.matrix(z[,1])
  }
  uv<-matrix(rnorm(n*(l+1)),n,l+1)%*%expm::sqrtm(cuv)
  
  
  u<-as.matrix(uv[,1])
  v<-as.matrix(uv[,2:(l+1)])
  
  x<-as.matrix(fpi(z,np,w1))+v
  y<-as.matrix(u)
  Y<-cbind(y,x)
  
  ome<-b_covar_hom(y,x,z,type,omet=cuv)
  
  YP<-b_YP(y,x,z)
  Y_pz<-as.matrix(YP[,1:2])
  Y_W<-as.matrix(YP[,3:4])
  
  xfit<-lm(x~z+0)$fitted.values
  b2sls<-solve(t(xfit)%*%xfit)%*%t(xfit)%*%y
  sig<-sum((y-x%*%b2sls)^2)/n
  
  sapply(gridbeta0, function(e) b_pval_hom(Y,Y_pz,Y_W,z,e,ome,b2sls,sig,xnorm,para=F))
}



b_sim_pval_het<-function(gridbeta0,xnorm,n,l,k,cuv,np="lin",w1="false",type="lin",h){
  
  z<-matrix(rnorm(n*k),n,k)
  if (np=="polyp"|np=="lin"|np=="nlin"){
    z<-as.matrix(z[,1])
  }
  uv_t<-matrix(rnorm(n*(l+1)),n,l+1)
  sr_cuv<-expm::sqrtm(cuv)
  uv<-t(vapply(1:n, function(e) sqrt((1+z[e,1]^2)/2)*t(uv_t[e,])%*%sr_cuv,FUN.VALUE=numeric(2)))
  
  u<-as.matrix(uv[,1])
  v<-as.matrix(uv[,2:(l+1)])
  
  x<-as.matrix(fpi(z,np,w1))+v
  y<-as.matrix(u)
  Y<-cbind(y,x)
  
  ome<-b_covar_het(y=y,x=x,z=z[,1],h=h,type=type,omet=cuv)
  
  xfit<-lm(x~z+0)$fitted.values
  b2sls<-solve(t(xfit)%*%xfit)%*%t(xfit)%*%y
  sig<-(y-x%*%b2sls)^2
  
  vapply(gridbeta0, function(e) b_pval_het(y=y,x=x,z=z,beta0=e,ome=ome,b2sls=b2sls,sig=sig,xnorm=xnorm,para=F),FUN.VALUE=numeric(7))
}



##### Parallel simulations

sim_pow_avg<-function(gridbeta0,m,n,l,k,cuv,np,w1,type,NS){
  
  ngrid<-length(gridbeta0)
  set.seed(1)
  xnorm<-matrix(rnorm(n*m),n,m)
  
  mat_pval<-matrix(0,7*ngrid,NS)
  mat_pval<-foreach::foreach(is=1:NS,.export=c("b_pval","b_pval_hom","b_sim_pval","m"),.packages=c("np","expm"),.combine=cbind) %dopar% {
    
    set.seed(is)
    as.numeric(b_sim_pval(gridbeta0=gridbeta0,xnorm=xnorm,n=n,l=l,k=k,cuv=cuv,np=np,w1=w1,type=type))
    
  }
  return(mat_pval)
}


het_sim_pow_avg<-function(gridbeta0,m,n,l,k,cuv,np,w1,type,h,NS){
  
  ngrid<-length(gridbeta0)
  set.seed(1)
  xnorm<-matrix(rnorm(n*m),n,m)
  
  het_mat_pval<-matrix(0,7*ngrid,NS)
  het_mat_pval<-foreach::foreach(is=1:NS,.export=c("b_pval","b_pval_het","b_covar_het",
                                                   "b_sim_pval_het","m"),
                                 .packages=c("np","expm","MASS"),.combine=cbind) %dopar% {
    
    set.seed(is)
    as.numeric(b_sim_pval_het(gridbeta0=gridbeta0,xnorm=xnorm,n=n,l=l,k=k,cuv=cuv,np=np,w1=w1,
                              type=type,h=h))
    
  }
  return(het_mat_pval)
}




###############################################################
##### Plot Functions ##########################################
###############################################################

plot_pval<-function(pval,gridbeta0,nplots,alpha,pow_curv=T,cutoff){
  
  ngrid<-length(gridbeta0)
  num_plots<-length(nplots)
  power_data<-data.frame(power=pval,bnull=rep(rep(gridbeta0,each=7),num_plots),
                         test=factor(rep(rep(1:7,ngrid),num_plots)),
                         np=rep(nplots,each=7*ngrid))
  
  colgg<-c("blueviolet","springgreen4","blue","orange","red","gold3","black")
  labgg<-c("AR test","LM test","CLR test","ICM test","KICM test","CICM test","Wald-2SLS test")
  lgg<-c(1,1,1,2,4,1,3)
  sgg<-c(0.5,0.5,0.5,1,1.5,1,1)
  
  if (pow_curv==T){
    ggall<-ggplot(power_data,aes(x=bnull,y=power,col=test,linetype=test,size=test))+geom_line()+
      scale_color_manual(name="Power built from",labels=labgg,values=colgg)+xlab("Null")+
      theme(axis.title.y=element_blank())+
      scale_linetype_manual(name="Power built from",labels=labgg,values=lgg)+
      scale_size_manual(name="Power built from",labels=labgg,values=sgg)+
      geom_hline(yintercept=alpha)+coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-0.05,1.05), expand = FALSE)+
      annotate("text",x=1,y=0.12,label="10% cutoff")+
      facet_grid(rows=vars(np))
  } else {
    ggall<-ggplot(power_data,aes(x=bnull,y=power,col=test,linetype=test,size=test))+geom_line()+
      scale_color_manual(name="Average p-value built from",labels=labgg,values=colgg)+xlab("Null")+
      theme(axis.title.y=element_blank())+
      scale_linetype_manual(name="Average p-value built from",labels=labgg,values=lgg)+
      scale_size_manual(name="Average p-value built from",labels=labgg,values=sgg)+
      geom_hline(yintercept=alpha)+coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-0.05,0.6), expand = FALSE)+
      annotate("text",x=1,y=cutoff,label=paste(cutoff,"cutoff")+
      facet_grid(rows=vars(np))
  }
  
  
  return(ggall)
  
}



