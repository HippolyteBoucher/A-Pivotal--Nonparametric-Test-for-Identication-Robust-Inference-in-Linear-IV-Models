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
k<-2 # number of instruments

beta<-rep(0,l) # true beta
rho<-0.81 # endogeneity coefficient
cuv<-matrix(0,l+1,l+1) # corr between errors
diag(cuv)<-rep(1,l+1)
cuv[2:(l+1),1]<-rep(rho,l)
cuv[1,2:(l+1)]<-rep(rho,l)

omet<-cuv

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
                           "b_sim","b_YP","b_yin","b_sim_pval","fpi","sim_pow_avg"))


##### Strong Instruments m=200, n=100

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=200,n=100,l=l,k=k,cuv=cuv,np="lin",w1="false",type="lin",NS=NS)
toc()
mat_pow_lin_str_s200_n100<-rowMeans(mat_pval<alpha)
mat_avg_lin_str_s200_n100<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=200,n=100,l=l,k=k,cuv=cuv,np="nlin",w1="false",type="lin",NS=NS)
toc()
mat_pow_nlin_str_s200_n100<-rowMeans(mat_pval<alpha)
mat_avg_nlin_str_s200_n100<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=200,n=100,l=l,k=k,cuv=cuv,np="polyp",w1="false",type="lin",NS=NS)
toc()
mat_pow_polyp_str_s200_n100<-rowMeans(mat_pval<alpha)
mat_avg_polyp_str_s200_n100<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=200,n=100,l=l,k=k,cuv=cuv,np="polysp",w1="false",type="lin",NS=NS)
toc()
mat_pow_polysp_str_s200_n100<-rowMeans(mat_pval<alpha)
mat_avg_polysp_str_s200_n100<-rowMeans(mat_pval)


##### Strong Instruments m=500, n=400

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="lin",w1="false",type="lin",NS=NS)
toc()
mat_pow_lin_str_s500_n400<-rowMeans(mat_pval<alpha)
mat_avg_lin_str_s500_n400<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="nlin",w1="false",type="lin",NS=NS)
toc()
mat_pow_nlin_str_s500_n400<-rowMeans(mat_pval<alpha)
mat_avg_nlin_str_s500_n400<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="polyp",w1="false",type="lin",NS=NS)
toc()
mat_pow_polyp_str_s500_n400<-rowMeans(mat_pval<alpha)
mat_avg_polyp_str_s500_n400<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="polysp",w1="false",type="lin",NS=NS)
toc()
mat_pow_polysp_str_s500_n400<-rowMeans(mat_pval<alpha)
mat_avg_polysp_str_s500_n400<-rowMeans(mat_pval)



##### Semi-strong Instruments m=200, n=100

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=200,n=100,l=l,k=k,cuv=cuv,np="lin",w1="semi",type="lin",NS=NS)
toc()
mat_pow_lin_semi_s200_n100<-rowMeans(mat_pval<alpha)
mat_avg_lin_semi_s200_n100<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=200,n=100,l=l,k=k,cuv=cuv,np="nlin",w1="semi",type="lin",NS=NS)
toc()
mat_pow_nlin_semi_s200_n100<-rowMeans(mat_pval<alpha)
mat_avg_nlin_semi_s200_n100<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=200,n=100,l=l,k=k,cuv=cuv,np="polyp",w1="semi",type="lin",NS=NS)
toc()
mat_pow_polyp_semi_s200_n100<-rowMeans(mat_pval<alpha)
mat_avg_polyp_semi_s200_n100<-rowMeans(mat_pval)


tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=200,n=100,l=l,k=k,cuv=cuv,np="polysp",w1="semi",type="lin",NS=NS)
toc()
mat_pow_polysp_semi_s200_n100<-rowMeans(mat_pval<alpha)
mat_avg_polysp_semi_s200_n100<-rowMeans(mat_pval)


##### Semi-strong Instruments m=500, n=400


tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="lin",w1="semi",type="lin",NS=NS)
toc()
mat_pow_lin_semi_s500_n400<-rowMeans(mat_pval<alpha)
mat_avg_lin_semi_s500_n400<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="nlin",w1="semi",type="lin",NS=NS)
toc()
mat_pow_nlin_semi_s500_n400<-rowMeans(mat_pval<alpha)
mat_avg_nlin_semi_s500_n400<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="polyp",w1="semi",type="lin",NS=NS)
toc()
mat_pow_polyp_semi_s500_n400<-rowMeans(mat_pval<alpha)
mat_avg_polyp_semi_s500_n400<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="polysp",w1="semi",type="lin",NS=NS)
toc()
mat_pow_polysp_semi_s500_n400<-rowMeans(mat_pval<alpha)
mat_avg_polysp_semi_s500_n400<-rowMeans(mat_pval)


##### Weak Instruments m=200, n=100

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=200,n=100,l=l,k=k,cuv=cuv,np="lin",w1="weak",type="lin",NS=NS)
toc()
mat_pow_lin_weak_s200_n100<-rowMeans(mat_pval<alpha)
mat_avg_lin_weak_s200_n100<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=200,n=100,l=l,k=k,cuv=cuv,np="nlin",w1="weak",type="lin",NS=NS)
toc()
mat_pow_nlin_weak_s200_n100<-rowMeans(mat_pval<alpha)
mat_avg_nlin_weak_s200_n100<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=200,n=100,l=l,k=k,cuv=cuv,np="polyp",w1="weak",type="lin",NS=NS)
toc()
mat_pow_polyp_weak_s200_n100<-rowMeans(mat_pval<alpha)
mat_avg_polyp_weak_s200_n100<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=200,n=100,l=l,k=k,cuv=cuv,np="polysp",w1="weak",type="lin",NS=NS)
toc()
mat_pow_polysp_weak_s200_n100<-rowMeans(mat_pval<alpha)
mat_avg_polysp_weak_s200_n100<-rowMeans(mat_pval)


##### Weak Instruments m=500, n=400

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="lin",w1="weak",type="lin",NS=NS)
toc()
mat_pow_lin_weak_s500_n400<-rowMeans(mat_pval<alpha)
mat_avg_lin_weak_s500_n400<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="nlin",w1="weak",type="lin",NS=NS)
toc()
mat_pow_nlin_weak_s500_n400<-rowMeans(mat_pval<alpha)
mat_avg_nlin_weak_s500_n400<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="polyp",w1="weak",type="lin",NS=NS)
toc()
mat_pow_polyp_weak_s500_n400<-rowMeans(mat_pval<alpha)
mat_avg_polyp_weak_s500_n400<-rowMeans(mat_pval)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=500,n=400,l=l,k=k,cuv=cuv,np="polysp",w1="weak",type="lin",NS=NS)
toc()
mat_pow_polysp_weak_s500_n400<-rowMeans(mat_pval<alpha)
mat_avg_polysp_weak_s500_n400<-rowMeans(mat_pval)


# Special 

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=200,n=100,l=l,k=4,cuv=cuv,np="spec",w1="false",type="lin",NS=NS)
toc()
mat_spec_str<-mat_pval
mat_pow_spec_str<-rowMeans(mat_spec_str<alpha)
mat_avg_spec_str<-rowMeans(mat_spec_str)


tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=200,n=100,l=l,k=4,cuv=cuv,np="spec",w1="semi",type="lin",NS=NS)
toc()
mat_spec_semi<-mat_pval
mat_pow_spec_semi<-rowMeans(mat_spec_semi<alpha)
mat_avg_spec_semi<-rowMeans(mat_spec_semi)

tic()
mat_pval<-sim_pow_avg(gridbeta0=gridbeta0,m=200,n=100,l=l,k=4,cuv=cuv,np="spec",w1="weak",type="lin",NS=NS)
toc()
mat_spec_weak<-mat_pval
mat_pow_spec_weak<-rowMeans(mat_spec_weak<alpha)
mat_avg_spec_weak<-rowMeans(mat_spec_weak)

rm(mat_pval)

stop(cluster)

save.image("sim_hom_power.RData")



###############################################################
##### Plots ###################################################
###############################################################

load("sim_hom_power.Rdata")

##### str s200 n100

data_plot<-c(mat_pow_lin_str_s200_n100,mat_pow_nlin_str_s200_n100,mat_pow_polyp_str_s200_n100,
             mat_pow_polysp_str_s200_n100)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                       "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=T,cutoff=0.2)

data_plot<-c(mat_avg_lin_str_s200_n100,mat_avg_nlin_str_s200_n100,mat_avg_polyp_str_s200_n100,
             mat_avg_polysp_str_s200_n100)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=F,cutoff=0.2)

##### semi s200 n100

data_plot<-c(mat_pow_lin_semi_s200_n100,mat_pow_nlin_semi_s200_n100,mat_pow_polyp_semi_s200_n100,
             mat_pow_polysp_semi_s200_n100)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=T,cutoff=0.2)

data_plot<-c(mat_avg_lin_semi_s200_n100,mat_avg_nlin_semi_s200_n100,mat_avg_polyp_semi_s200_n100,
             mat_avg_polysp_semi_s200_n100)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=F,cutoff=0.2)

##### weak s200 n100

data_plot<-c(mat_pow_lin_weak_s200_n100,mat_pow_nlin_weak_s200_n100,mat_pow_polyp_weak_s200_n100,
             mat_pow_polysp_weak_s200_n100)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=T,cutoff=0.2)

data_plot<-c(mat_avg_lin_weak_s200_n100,mat_avg_nlin_weak_s200_n100,mat_avg_polyp_weak_s200_n100,
             mat_avg_polysp_weak_s200_n100)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=F,cutoff=0.2)


##### str s500 n400

data_plot<-c(mat_pow_lin_str_s500_n400,mat_pow_nlin_str_s500_n400,mat_pow_polyp_str_s500_n400,
             mat_pow_polysp_str_s500_n400)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=T,cutoff=0.2)

data_plot<-c(mat_avg_lin_str_s500_n400,mat_avg_nlin_str_s500_n400,mat_avg_polyp_str_s500_n400,
             mat_avg_polysp_str_s500_n400)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=F,cutoff=0.2)

##### semi s500 n400

data_plot<-c(mat_pow_lin_semi_s500_n400,mat_pow_nlin_semi_s500_n400,mat_pow_polyp_semi_s500_n400,
             mat_pow_polysp_semi_s500_n400)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=T,cutoff=0.2)

data_plot<-c(mat_avg_lin_semi_s500_n400,mat_avg_nlin_semi_s500_n400,mat_avg_polyp_semi_s500_n400,
             mat_avg_polysp_semi_s500_n400)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=F,cutoff=0.2)

##### weak s500 n400

data_plot<-c(mat_pow_lin_weak_s500_n400,mat_pow_nlin_weak_s500_n400,mat_pow_polyp_weak_s500_n400,
             mat_pow_polysp_weak_s500_n400)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=T,cutoff=0.2)

data_plot<-c(mat_avg_lin_weak_s500_n400,mat_avg_nlin_weak_s500_n400,mat_avg_polyp_weak_s500_n400,
             mat_avg_polysp_weak_s500_n400)
plot_pval(data_plot,gridbeta0=gridbeta0,nplots=c("Linear 1st stage","Non-linear 1st stage",
                                                 "Polar polynomial 1st stage","Semi-polar polynomial 1st stage"),
          alpha=0.1,pow_curv=F,cutoff=0.2)

