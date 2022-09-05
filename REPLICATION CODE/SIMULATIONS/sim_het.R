library("pracma")
library("np")
library("expm")
library("foreach")
library("doParallel")
library("ggplot2")

source("functions.R")



#####################################################################################
##### Simulations Setup #############################################################
#####################################################################################


NS<-5000 # number of simulations
n<-400
m<-500

l<-1 # number of endogenous variables
k<-2 # number of instruments

beta<-rep(0,l) # true beta
rho<-0.81 # endogeneity coefficient
cuv<-matrix(0,l+1,l+1) # corr between errors
diag(cuv)<-rep(1,l+1)
cuv[2:(l+1),1]<-rep(rho,l)
cuv[1,2:(l+1)]<-rep(rho,l)

alpha<-0.1 # level of the test

ngrid<-21
gridbeta0<-seq(beta-1.5,beta+1.5,length.out=ngrid)




###############################################################
##### Simulations #############################################
###############################################################


ncores<-detectCores()-1
cluster<-makeCluster(mc <- getOption("cl.cores", ncores))
registerDoParallel(cluster)
clusterExport(cl=cluster,c("b_covar_het","b_covar_hom","b_pval","b_pval_het","b_pval_hom",
                           "b_sim","b_YP","b_yin","b_sim_pval_het","fpi","het_sim_pow_avg"))


##### Strong Instruments m=500, n=400

tic()
het_mat_pval<-het_sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="lin",w1="false",type="gls-np",h=F,NS=NS)
toc()
het_mat_pow_lin_str_s500_n400<-rowMeans(het_mat_pval<alpha)
het_mat_avg_lin_str_s500_n400<-rowMeans(het_mat_pval)
het_mat_size_lin_str_s500_n400<-het_mat_pow_lin_str_s500_n400[(7*10+1):(7*11)]
het_mat_size_lin_str_s500_n400


tic()
het_mat_pval<-het_sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="nlin",w1="false",type="gls-np",h=F,NS=NS)
toc()
het_mat_pow_nlin_str_s500_n400<-rowMeans(het_mat_pval<alpha)
het_mat_avg_nlin_str_s500_n400<-rowMeans(het_mat_pval)
het_mat_size_nlin_str_s500_n400<-het_mat_pow_nlin_str_s500_n400[(7*10+1):(7*11)]
het_mat_size_nlin_str_s500_n400


tic()
het_mat_pval<-het_sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="polyp",w1="false",type="gls-np",h=F,NS=NS)
toc()
het_mat_pow_polyp_str_s500_n400<-rowMeans(het_mat_pval<alpha)
het_mat_avg_polyp_str_s500_n400<-rowMeans(het_mat_pval)
het_mat_size_polyp_str_s500_n400<-het_mat_pow_polyp_str_s500_n400[(7*10+1):(7*11)]
het_mat_size_polyp_str_s500_n400



tic()
het_mat_pval<-het_sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="polysp",w1="false",type="gls-np",h=F,NS=100)
toc()
het_mat_pow_polysp_str_s500_n400<-rowMeans(het_mat_pval<alpha)
het_mat_avg_polysp_str_s500_n400<-rowMeans(het_mat_pval)
het_mat_size_polysp_str_s500_n400<-het_mat_pow_polysp_str_s500_n400[(7*10+1):(7*11)]
het_mat_size_polysp_str_s500_n400




##### Semi-strong Instruments m=500, n=400

tic()
het_mat_pval<-het_sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="lin",w1="semi",type="gls-np",h=F,NS=NS)
toc()
het_mat_pow_lin_semi_s500_n400<-rowMeans(het_mat_pval<alpha)
het_mat_avg_lin_semi_s500_n400<-rowMeans(het_mat_pval)
het_mat_size_lin_semi_s500_n400<-het_mat_pow_lin_semi_s500_n400[71:77]
het_mat_size_lin_semi_s500_n400

tic()
het_mat_pval<-het_sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="nlin",w1="semi",type="gls-np",h=F,NS=NS)
toc()
het_mat_pow_nlin_semi_s500_n400<-rowMeans(het_mat_pval<alpha)
het_mat_avg_nlin_semi_s500_n400<-rowMeans(het_mat_pval)
het_mat_size_nlin_semi_s500_n400<-het_mat_pow_nlin_semi_s500_n400[71:77]
het_mat_size_nlin_semi_s500_n400


tic()
het_mat_pval<-het_sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="polyp",w1="semi",type="gls-np",h=F,NS=NS)
toc()
het_mat_pow_polyp_semi_s500_n400<-rowMeans(het_mat_pval<alpha)
het_mat_avg_polyp_semi_s500_n400<-rowMeans(het_mat_pval)
het_mat_size_polyp_semi_s500_n400<-het_mat_pow_polyp_semi_s500_n400[71:77]
het_mat_size_polyp_semi_s500_n400


tic()
het_mat_pval<-het_sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="polysp",w1="semi",type="gls-np",h=F,NS=NS)
toc()
het_mat_pow_polysp_semi_s500_n400<-rowMeans(het_mat_pval<alpha)
het_mat_avg_polysp_semi_s500_n400<-rowMeans(het_mat_pval)
het_mat_size_polysp_semi_s500_n400<-het_mat_pow_polysp_semi_s500_n400[71:77]
het_mat_size_polysp_semi_s500_n400




##### Weak Instruments m=500, n=400

tic()
het_mat_pval<-het_sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="lin",w1="weak",type="gls-np",h=F,NS=NS)
toc()
het_mat_pow_lin_weak_s500_n400<-rowMeans(het_mat_pval<alpha)
het_mat_avg_lin_weak_s500_n400<-rowMeans(het_mat_pval)
het_mat_size_lin_weak_s500_n400<-het_mat_pow_lin_weak_s500_n400[71:77]
het_mat_size_lin_weak_s500_n400

tic()
het_mat_pval<-het_sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="nlin",w1="weak",type="gls-np",h=F,NS=NS)
toc()
het_mat_pow_nlin_weak_s500_n400<-rowMeans(het_mat_pval<alpha)
het_mat_avg_nlin_weak_s500_n400<-rowMeans(het_mat_pval)
het_mat_size_nlin_weak_s500_n400<-het_mat_pow_nlin_weak_s500_n400[71:77]
het_mat_size_nlin_weak_s500_n400


tic()
het_mat_pval<-het_sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="polyp",w1="weak",type="gls-np",h=F,NS=NS)
toc()
het_mat_pow_polyp_weak_s500_n400<-rowMeans(het_mat_pval<alpha)
het_mat_avg_polyp_weak_s500_n400<-rowMeans(het_mat_pval)
het_mat_size_polyp_weak_s500_n400<-het_mat_pow_polyp_weak_s500_n400[71:77]
het_mat_size_polyp_weak_s500_n400

tic()
het_mat_pval<-het_sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="polysp",w1="weak",type="gls-np",h=F,NS=NS)
toc()
het_mat_pow_polysp_weak_s500_n400<-rowMeans(het_mat_pval<alpha)
het_mat_avg_polysp_weak_s500_n400<-rowMeans(het_mat_pval)
het_mat_size_polysp_weak_s500_n400<-het_mat_pow_polysp_weak_s500_n400[71:77]
het_mat_size_polysp_weak_s500_n400

rm(het_mat_pval)

stop(cluster)

save.image("sim_het.RData")




###############################################################
##### Plots ###################################################
###############################################################

load("sim_het.Rdata")



##### str s500 n400

data_plot<-c(het_mat_pow_lin_str_s500_n400,het_mat_pow_nlin_str_s500_n400,het_mat_pow_polyp_str_s500_n400)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage"),
          alpha=0.1,pow_curv=T,cutoff=0.1)

data_plot<-c(het_mat_avg_lin_str_s500_n400,het_mat_avg_nlin_str_s500_n400,het_mat_avg_polyp_str_s500_n400,
             het_mat_avg_polysp_str_s500_n400)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=F,cutoff=0.2)

##### semi s500 n400

data_plot<-c(het_mat_pow_lin_semi_s500_n400,het_mat_pow_nlin_semi_s500_n400,het_mat_pow_polyp_semi_s500_n400,
             het_mat_pow_polysp_semi_s500_n400)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=T,cutoff=0.2)

data_plot<-c(het_mat_avg_lin_semi_s500_n400,het_mat_avg_nlin_semi_s500_n400,het_mat_avg_polyp_semi_s500_n400,
             het_mat_avg_polysp_semi_s500_n400)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=F,cutoff=0.2)

##### weak s500 n400

data_plot<-c(het_mat_pow_lin_weak_s500_n400,het_mat_pow_nlin_weak_s500_n400,het_mat_pow_polyp_weak_s500_n400,
             het_mat_pow_polysp_weak_s500_n400)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=T,cutoff=0.2)

data_plot<-c(het_mat_avg_lin_weak_s500_n400,het_mat_avg_nlin_weak_s500_n400,het_mat_avg_polyp_weak_s500_n400,
             het_mat_avg_polysp_weak_s500_n400)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=F,cutoff=0.2)

