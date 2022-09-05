library("pracma")
library("np")
library("expm")
library("foreach")
library("doParallel")
library("doRNG")

source("functions.R")


#####################################################################################
##### Simulations Setup #############################################################
#####################################################################################


NS<-5000 # number of simulations

l<-1 # number of endogenous variables
k<-10 # number of instruments

beta<-rep(0,l) # true beta
rho<-0.81 # endogeneity coefficient
cuv<-matrix(0,l+1,l+1) # corr between errors
diag(cuv)<-rep(1,l+1)
cuv[2:(l+1),1]<-rep(rho,l)
cuv[1,2:(l+1)]<-rep(rho,l)

omet<-cuv

alpha<-0.1 # level of the test

ngrid<-1
gridbeta0<-0


###############################################################
##### Simulations #############################################
###############################################################



##### Strong Instruments

ncores<-detectCores()-1
cluster<-makeCluster(mc <- getOption("cl.cores", ncores))
registerDoParallel(cluster)
clusterExport(cl=cluster,c("b_covar_het","b_covar_hom","b_pval","b_pval_het","b_pval_hom",
                           "b_sim","b_YP","b_yin","b_sim_pval","fpi"))

rowMeans(sim_pow_avg(gridbeta0=0,m=200,n=100,l=l,k=5,cuv=cuv,np="spec",w1="false",type="lin",NS=500)<alpha)

##### Strong Instruments m=200, n=100

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=100,l=l,k=k,cuv=cuv,np="lin",w1="false",type="lin",NS=NS)
toc()
mat_size_lin_str_s200_n100<-rowMeans(mat_pval<alpha)
mat_size_lin_str_s200_n100

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=100,l=l,k=k,cuv=cuv,np="nlin",w1="false",type="lin",NS=NS)
toc()
mat_size_nlin_str_s200_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=100,l=l,k=k,cuv=cuv,np="polyp",w1="false",type="lin",NS=NS)
toc()
mat_size_polyp_str_s200_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=100,l=l,k=k,cuv=cuv,np="polysp",w1="false",type="lin",NS=NS)
toc()
mat_size_polysp_str_s200_n100<-rowMeans(mat_pval<alpha)
mat_size_polysp_str_s200_n100



tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=100,l=l,k=k,cuv=cuv,np="lin",w1="false",type="lin",NS=NS)
toc()
mat_size_lin_str_s500_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=100,l=l,k=k,cuv=cuv,np="nlin",w1="false",type="lin",NS=NS)
toc()
mat_size_nlin_str_s500_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=100,l=l,k=k,cuv=cuv,np="polyp",w1="false",type="lin",NS=NS)
toc()
mat_size_polyp_str_s500_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=100,l=l,k=k,cuv=cuv,np="polysp",w1="false",type="lin",NS=NS)
toc()
mat_size_polysp_str_s500_n100<-rowMeans(mat_pval<alpha)
mat_size_polysp_str_s500_n100



tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=400,l=l,k=k,cuv=cuv,np="lin",w1="false",type="lin",NS=NS)
toc()
mat_size_lin_str_s200_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=400,l=l,k=k,cuv=cuv,np="nlin",w1="false",type="lin",NS=NS)
toc()
mat_size_nlin_str_s200_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=400,l=l,k=k,cuv=cuv,np="polyp",w1="false",type="lin",NS=NS)
toc()
mat_size_polyp_str_s200_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=400,l=l,k=k,cuv=cuv,np="polysp",w1="false",type="lin",NS=NS)
toc()
mat_size_polysp_str_s200_n400<-rowMeans(mat_pval<alpha)




tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=400,l=l,k=k,cuv=cuv,np="lin",w1="false",type="lin",NS=NS)
toc()
mat_size_lin_str_s500_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=400,l=l,k=k,cuv=cuv,np="nlin",w1="false",type="lin",NS=NS)
toc()
mat_size_nlin_str_s500_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=400,l=l,k=k,cuv=cuv,np="polyp",w1="false",type="lin",NS=NS)
toc()
mat_size_polyp_str_s500_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=400,l=l,k=k,cuv=cuv,np="polysp",w1="false",type="lin",NS=NS)
toc()
mat_size_polysp_str_s500_n400<-rowMeans(mat_pval<alpha)




##### Semi-Strong Instruments

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=100,l=l,k=k,cuv=cuv,np="lin",w1="semi",type="lin",NS=NS)
toc()
mat_size_lin_semi_s200_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=100,l=l,k=k,cuv=cuv,np="nlin",w1="semi",type="lin",NS=NS)
toc()
mat_size_nlin_semi_s200_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=100,l=l,k=k,cuv=cuv,np="polyp",w1="semi",type="lin",NS=NS)
toc()
mat_size_polyp_semi_s200_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=100,l=l,k=k,cuv=cuv,np="polysp",w1="semi",type="lin",NS=NS)
toc()
mat_size_polysp_semi_s200_n100<-rowMeans(mat_pval<alpha)


tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=100,l=l,k=k,cuv=cuv,np="lin",w1="semi",type="lin",NS=NS)
toc()
mat_size_lin_semi_s500_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=100,l=l,k=k,cuv=cuv,np="nlin",w1="semi",type="lin",NS=NS)
toc()
mat_size_nlin_semi_s500_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=100,l=l,k=k,cuv=cuv,np="polyp",w1="semi",type="lin",NS=NS)
toc()
mat_size_polyp_semi_s500_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=100,l=l,k=k,cuv=cuv,np="polysp",w1="semi",type="lin",NS=NS)
toc()
mat_size_polysp_semi_s500_n100<-rowMeans(mat_pval<alpha)


tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=400,l=l,k=k,cuv=cuv,np="lin",w1="semi",type="lin",NS=NS)
toc()
mat_size_lin_semi_s200_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=400,l=l,k=k,cuv=cuv,np="nlin",w1="semi",type="lin",NS=NS)
toc()
mat_size_nlin_semi_s200_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=400,l=l,k=k,cuv=cuv,np="polyp",w1="semi",type="lin",NS=NS)
toc()
mat_size_polyp_semi_s200_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=400,l=l,k=k,cuv=cuv,np="polysp",w1="semi",type="lin",NS=NS)
toc()
mat_size_polysp_semi_s200_n400<-rowMeans(mat_pval<alpha)


tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=400,l=l,k=k,cuv=cuv,np="lin",w1="semi",type="lin",NS=NS)
toc()
mat_size_lin_semi_s500_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=400,l=l,k=k,cuv=cuv,np="nlin",w1="semi",type="lin",NS=NS)
toc()
mat_size_nlin_semi_s500_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=400,l=l,k=k,cuv=cuv,np="polyp",w1="semi",type="lin",NS=NS)
toc()
mat_size_polyp_semi_s500_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=400,l=l,k=k,cuv=cuv,np="polysp",w1="semi",type="lin",NS=NS)
toc()
mat_size_polysp_semi_s500_n400<-rowMeans(mat_pval<alpha)


##### Weak Instruments

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=100,l=l,k=k,cuv=cuv,np="lin",w1="weak",type="lin",NS=NS)
toc()
mat_size_lin_weak_s200_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=100,l=l,k=k,cuv=cuv,np="nlin",w1="weak",type="lin",NS=NS)
toc()
mat_size_nlin_weak_s200_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=100,l=l,k=k,cuv=cuv,np="polyp",w1="weak",type="lin",NS=NS)
toc()
mat_size_polyp_weak_s200_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=100,l=l,k=k,cuv=cuv,np="polysp",w1="weak",type="lin",NS=NS)
toc()
mat_size_polysp_weak_s200_n100<-rowMeans(mat_pval<alpha)


tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=100,l=l,k=k,cuv=cuv,np="lin",w1="weak",type="lin",NS=NS)
toc()
mat_size_lin_weak_s500_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=100,l=l,k=k,cuv=cuv,np="nlin",w1="weak",type="lin",NS=NS)
toc()
mat_size_nlin_weak_s500_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=100,l=l,k=k,cuv=cuv,np="polyp",w1="weak",type="lin",NS=NS)
toc()
mat_size_polyp_weak_s500_n100<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=100,l=l,k=k,cuv=cuv,np="polysp",w1="weak",type="lin",NS=NS)
toc()
mat_size_polysp_weak_s500_n100<-rowMeans(mat_pval<alpha)


tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=400,l=l,k=k,cuv=cuv,np="lin",w1="weak",type="lin",NS=NS)
toc()
mat_size_lin_weak_s200_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=400,l=l,k=k,cuv=cuv,np="nlin",w1="weak",type="lin",NS=NS)
toc()
mat_size_nlin_weak_s200_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=400,l=l,k=k,cuv=cuv,np="polyp",w1="weak",type="lin",NS=NS)
toc()
mat_size_polyp_weak_s200_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=200,n=400,l=l,k=k,cuv=cuv,np="polysp",w1="weak",type="lin",NS=NS)
toc()
mat_size_polysp_weak_s200_n400<-rowMeans(mat_pval<alpha)


tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=400,l=l,k=k,cuv=cuv,np="lin",w1="weak",type="lin",NS=NS)
toc()
mat_size_lin_weak_s500_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=400,l=l,k=k,cuv=cuv,np="nlin",w1="weak",type="lin",NS=NS)
toc()
mat_size_nlin_weak_s500_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=400,l=l,k=k,cuv=cuv,np="polyp",w1="weak",type="lin",NS=NS)
toc()
mat_size_polyp_weak_s500_n400<-rowMeans(mat_pval<alpha)

tic()
mat_pval<-sim_pow_avg(gridbeta0=0,m=500,n=400,l=l,k=k,cuv=cuv,np="polysp",w1="weak",type="lin",NS=NS)
toc()
mat_size_polysp_weak_s500_n400<-rowMeans(mat_pval<alpha)



stop(cluster)

save.image("sim_hom_size.RData")
load("sim_hom_size.Rdata")

