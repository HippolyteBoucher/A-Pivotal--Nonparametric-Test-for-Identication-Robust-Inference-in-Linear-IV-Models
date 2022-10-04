library("rARPACK")
library("foreign")
library("expm")
library("emdbook")
library("foreach")
library("doParallel")
library("Rcpp")
library("ggplot2")
library("ivmodel")
library("nortest")
library("dplyr")
library("AER")
library("pracma")

source("functions.R")


#####################################################################################
###### Obtaining the Data ###########################################################
#####################################################################################

# The data from Angrist and Krueger (1991) is freely available at
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/ENLGZX
# only the final cleaned dataset is uploaded


#####################################################################################
###### Transforming the Data ########################################################
#####################################################################################

df<-read.dta("NEW7080.dta")
names(df)<-c("age","ageq","v3","educ","enocent","esocent",
             "v7","v8","lwklywge","married","midatl","mt",
             "neweng","v14","v15","census","v17","qob",
             "race","smsa","soatl","v22","v23","wnocent",
             "wsocent","v26","yob")

##### Cohort and other variables

df$cohort<-20.29
df$cohort[df$yob>=30 & df$yob<=39]<-30.39
df$cohort[df$yob>=40 & df$yob<=49]<-40.49
df$ageq[df$census==80]<-df$ageq[df$census==80]-1900
df$ageqsq<-df$ageq^2


##### Year of birth dummies

df$yr20<-0
df$yr20[df$yob==1920]<-1
df$yr20[df$yob==30]<-1
df$yr20[df$yob==40]<-1

df$yr21<-0
df$yr21[df$yob==1921]<-1
df$yr21[df$yob==31]<-1
df$yr21[df$yob==41]<-1

df$yr22<-0
df$yr22[df$yob==1922]<-1
df$yr22[df$yob==32]<-1
df$yr22[df$yob==42]<-1

df$yr23<-0
df$yr23[df$yob==1923]<-1
df$yr23[df$yob==33]<-1
df$yr23[df$yob==43]<-1

df$yr24<-0
df$yr24[df$yob==1924]<-1
df$yr24[df$yob==34]<-1
df$yr24[df$yob==44]<-1

df$yr25<-0
df$yr25[df$yob==1925]<-1
df$yr25[df$yob==35]<-1
df$yr25[df$yob==45]<-1

df$yr26<-0
df$yr26[df$yob==1926]<-1
df$yr26[df$yob==36]<-1
df$yr26[df$yob==46]<-1

df$yr27<-0
df$yr27[df$yob==1927]<-1
df$yr27[df$yob==37]<-1
df$yr27[df$yob==47]<-1

df$yr28<-0
df$yr28[df$yob==1928]<-1
df$yr28[df$yob==38]<-1
df$yr28[df$yob==48]<-1

df$yr29<-0
df$yr29[df$yob==1929]<-1
df$yr29[df$yob==39]<-1
df$yr29[df$yob==49]<-1


##### Quarter of birth dummies

df$qtr1<-0
df$qtr1[df$qob==1]<-1

df$qtr2<-0
df$qtr2[df$qob==2]<-1

df$qtr3<-0
df$qtr3[df$qob==3]<-1

df$qtr4<-0
df$qtr4[df$qob==4]<-1


##### YOB*QOB dummies

df$qtr120<-df$qtr1*df$yr20
df$qtr121<-df$qtr1*df$yr21
df$qtr122<-df$qtr1*df$yr22
df$qtr123<-df$qtr1*df$yr23
df$qtr124<-df$qtr1*df$yr24
df$qtr125<-df$qtr1*df$yr25
df$qtr126<-df$qtr1*df$yr26
df$qtr127<-df$qtr1*df$yr27
df$qtr128<-df$qtr1*df$yr28
df$qtr129<-df$qtr1*df$yr29

df$qtr220<-df$qtr2*df$yr20
df$qtr221<-df$qtr2*df$yr21
df$qtr222<-df$qtr2*df$yr22
df$qtr223<-df$qtr2*df$yr23
df$qtr224<-df$qtr2*df$yr24
df$qtr225<-df$qtr2*df$yr25
df$qtr226<-df$qtr2*df$yr26
df$qtr227<-df$qtr2*df$yr27
df$qtr228<-df$qtr2*df$yr28
df$qtr229<-df$qtr2*df$yr29

df$qtr320<-df$qtr3*df$yr20
df$qtr321<-df$qtr3*df$yr21
df$qtr322<-df$qtr3*df$yr22
df$qtr323<-df$qtr3*df$yr23
df$qtr324<-df$qtr3*df$yr24
df$qtr325<-df$qtr3*df$yr25
df$qtr326<-df$qtr3*df$yr26
df$qtr327<-df$qtr3*df$yr27
df$qtr328<-df$qtr3*df$yr28
df$qtr329<-df$qtr3*df$yr29


##### Cleaned data

df2<-df[df$cohort==20.29,]
rm(df)

df3<-df2
rm(df2)
df3$v3<-NULL
df3$v7<-NULL
df3$v8<-NULL
df3$v15<-NULL
df3$v17<-NULL
df3$v22<-NULL
df3$v23<-NULL
df3$v26<-NULL
df3$v3<-NULL
df3$ed<-df3$educ
df3$educ<-NULL
df3$lwage<-df3$lwklywge
df3$lwklywge<-NULL
df3$v14<-NULL
df3$yob<-NULL
df3$cohort<-NULL
df3$age<-NULL
df3$qob<-NULL
df3$census<-NULL



#####. Final Dataset

save(df3,file="final_dataset.Rdata")
load("final_dataset.Rdata")


##### OLS

modelln_1<-lm(lwage~ed+yr20+yr21+yr22+yr23+yr24+yr25+yr26
              +yr27+yr28,data=df3)
modelln_2<-lm(lwage~ed+ageq+ageqsq+yr20+yr21+yr22+yr23+yr24
              +yr25+yr26+yr27+yr28,data=df3)
modelln_3<-lm(lwage~ed+ageq+ageqsq+married+race+smsa+yr20
              +yr21+yr22+yr23+yr24+yr25+yr26+yr27+yr28,data=df3)
modelln_4<-lm(lwage~ed+ageq+ageqsq+enocent+esocent+married
              +midatl+mt+neweng+race+smsa+soatl+wnocent+wsocent
              +yr20+yr21+yr22+yr23+yr24+yr25+yr26+yr27+yr28,
              data=df3)

summary(modelln_1);summary(modelln_2);summary(modelln_3);summary(modelln_4)

rm(modelln_1);rm(modelln_2);rm(modelln_3);rm(modelln_4)


##### 2SLS

y<-df3$lwage
d<-df3$ed

x_1<-subset(df3,select=14:22)
x_2<-subset(df3,select=c(1,13:22))
x_3<-subset(df3,select=c(1,4,8:9,13:22))
x_4<-subset(df3,select=1:22)

z_1<-subset(df3,select=c(14:22,28:56))
z_2<-subset(df3,select=c(1,13:22,28:56))
z_3<-subset(df3,select=c(1,4,8:9,13:22,28:56))
z_4<-subset(df3,select=c(1:22,28:56))

z<-subset(df3,select=28:56)

rm(df3)

modeliv_1<-ivmodel(y,d,z,x_1)
modeliv_2<-ivmodel(y,d,z,x_2)
modeliv_3<-ivmodel(y,d,z,x_3)
modeliv_4<-ivmodel(y,d,z,x_4)

summary(modeliv_1)
summary(modeliv_2)
summary(modeliv_3)
summary(modeliv_4)


##### Recover IV estimators and standard-deviation

est_1<-matrix(0,4,2)
est_1[1:2,1]<-modeliv_1$kClass$point.est
est_1[1:2,2]<-modeliv_1$kClass$std.err
est_1[3,1]<-modeliv_1$LIML$point.est
est_1[3,2]<-modeliv_1$LIML$std.err
est_1[4,1]<-modeliv_1$Fuller$point.est
est_1[4,2]<-modeliv_1$Fuller$std.err

est_2<-matrix(0,4,2)
est_2[1:2,1]<-modeliv_2$kClass$point.est
est_2[1:2,2]<-modeliv_2$kClass$std.err
est_2[3,1]<-modeliv_2$LIML$point.est
est_2[3,2]<-modeliv_2$LIML$std.err
est_2[4,1]<-modeliv_2$Fuller$point.est
est_2[4,2]<-modeliv_2$Fuller$std.err

est_3<-matrix(0,4,2)
est_3[1:2,1]<-modeliv_3$kClass$point.est
est_3[1:2,2]<-modeliv_3$kClass$std.err
est_3[3,1]<-modeliv_3$LIML$point.est
est_3[3,2]<-modeliv_3$LIML$std.err
est_3[4,1]<-modeliv_3$Fuller$point.est
est_3[4,2]<-modeliv_3$Fuller$std.err

est_4<-matrix(0,4,2)
est_4[1:2,1]<-modeliv_4$kClass$point.est
est_4[1:2,2]<-modeliv_4$kClass$std.err
est_4[3,1]<-modeliv_4$LIML$point.est
est_4[3,2]<-modeliv_4$LIML$std.err
est_4[4,1]<-modeliv_4$Fuller$point.est
est_4[4,2]<-modeliv_4$Fuller$std.err

rm(modeliv_1);rm(modeliv_2);rm(modeliv_3);rm(modeliv_4)


##########################################################################
##### Inference via Test Inversion #######################################
##########################################################################

##### Preparing the data to use KICM

z<-as.matrix(z)
z_1<-as.matrix(z_1)
z_2<-as.matrix(z_2)
z_3<-as.matrix(z_3)
z_4<-as.matrix(z_4)
n<-dim(z)[1]

### Project out the controls

mxd_1<-lm(d~as.matrix(x_1))$residuals
mxy_1<-lm(y~as.matrix(x_1))$residuals
mxd_2<-lm(d~as.matrix(x_2))$residuals
mxy_2<-lm(y~as.matrix(x_2))$residuals
mxd_3<-lm(d~as.matrix(x_3))$residuals
mxy_3<-lm(y~as.matrix(x_3))$residuals
mxd_4<-lm(d~as.matrix(x_4))$residuals
mxy_4<-lm(y~as.matrix(x_4))$residuals

### Get the residuals

res_1<-lm(cbind(mxy_1,mxd_1)~z_1+0)$residuals
res_2<-lm(cbind(mxy_2,mxd_2)~z_2+0)$residuals
res_3<-lm(cbind(mxy_3,mxd_3)~z_3+0)$residuals
res_4<-lm(cbind(mxy_4,mxd_4)~z_4+0)$residuals

### Compute F-statistics

F_1<-sum(mxd_1^2-res_1[,2]^2)/sum(res_1[,2]^2)*n/28
F_2<-sum(mxd_2^2-res_2[,2]^2)/sum(res_2[,2]^2)*n/28
F_3<-sum(mxd_3^2-res_3[,2]^2)/sum(res_3[,2]^2)*n/28
F_4<-sum(mxd_4^2-res_4[,2]^2)/sum(res_4[,2]^2)*n/28
F_1;F_2;F_3;F_4

rm(res_1);rm(res_2);rm(res_3);rm(res_4)
rm(y);rm(d);rm(x_1);rm(x_2);rm(x_3);rm(x_4);rm(z)

### Compute covariances

hom_1<-b_covar_hom(mxy_1,mxd_1,z_1,h=F,type="lin")
hom_2<-b_covar_hom(mxy_2,mxd_2,z_2,h=F,type="lin")
hom_3<-b_covar_hom(mxy_3,mxd_3,z_3,h=F,type="lin")
hom_4<-b_covar_hom(mxy_4,mxd_4,z_4,h=F,type="lin")

gridbeta0<-seq(0,0.75,by=0.01)
ngrid<-length(gridbeta0)

### Build the vector S and matrix T

YP_1<-b_YP(mxy_1,mxd_1,z_1)
YP_2<-b_YP(mxy_2,mxd_2,z_2)
YP_3<-b_YP(mxy_3,mxd_3,z_3)
YP_4<-b_YP(mxy_4,mxd_4,z_4)

Y_pz_1<-YP_1[,1:2]
Y_W_1<-YP_1[,3:4]
Y_1<-cbind(mxy_1,mxd_1)

Y_pz_2<-YP_2[,1:2]
Y_W_2<-YP_2[,3:4]
Y_2<-cbind(mxy_2,mxd_2)

Y_pz_3<-YP_3[,1:2]
Y_W_3<-YP_3[,3:4]
Y_3<-cbind(mxy_3,mxd_3)

Y_pz_4<-YP_4[,1:2]
Y_W_4<-YP_4[,3:4]
Y_4<-cbind(mxy_4,mxd_4)

rm(YP_1);rm(YP_2);rm(YP_3);rm(YP_4);rm(z_1);rm(z_2);rm(z_3);rm(z_4)


##### P-value for specification (1) for all methods over grid

pval_1<-sapply(gridbeta0, function(e) b_pval_hom(Y_1,Y_pz_1,Y_W_1,z_1,
                                                 e,hom_1,est_1[2,1],
                                                 est_1[2,2]^2,xnorm=0,
                                                 para=F,sim=F))
fpval_1<-matrix(0,8,ngrid)
fpval_1[1,]<-pval_1[1,]
fpval_1[2,]<-pval_1[2,]
fpval_1[3,]<-sapply(gridbeta0, function(e) CLR(modeliv_1,beta0=e)$p.value)
fpval_1[4,]<-pval_1[3,]
fpval_1[5,]<-sapply(gridbeta0, function(e) pchisq((est_1[1,1]-e)^2/est_1[1,2]^2,
                                                  df=1,lower.tail=F))
fpval_1[6,]<-pval_1[4,]
fpval_1[7,]<-sapply(gridbeta0, function(e) pchisq((est_1[3,1]-e)^2/est_1[3,2]^2,
                                                  df=1,lower.tail=F))
fpval_1[8,]<-sapply(gridbeta0, function(e) pchisq((est_1[4,1]-e)^2/est_1[4,2]^2,
                                                  df=1,lower.tail=F))

plot_pval_ak(fpval_1,gridbeta0,alpha=0.1,cutoff=0.12)


##### P-values for specification (2) for all methods over grid

pval_2<-sapply(gridbeta0, function(e) b_pval_hom(Y_2,Y_pz_2,Y_W_2,z_2,
                                                 e,hom_2,est_2[2,1],est_2[2,2]^2,xnorm=0,para=F,sim=F))
fpval_2<-matrix(0,8,ngrid)
fpval_2[1,]<-pval_2[1,]
fpval_2[2,]<-pval_2[2,]
fpval_2[3,]<-sapply(gridbeta0, function(e) CLR(modeliv_2,beta0=e)$p.value)
fpval_2[4,]<-pval_2[3,]
fpval_2[5,]<-sapply(gridbeta0, function(e) pchisq((est_2[1,1]-e)^2/est_2[1,2]^2,
                                                  df=1,lower.tail=F))
fpval_2[6,]<-pval_2[4,]
fpval_2[7,]<-sapply(gridbeta0, function(e) pchisq((est_2[3,1]-e)^2/est_2[3,2]^2,
                                                  df=1,lower.tail=F))
fpval_2[8,]<-sapply(gridbeta0, function(e) pchisq((est_2[4,1]-e)^2/est_2[4,2]^2,
                                                  df=1,lower.tail=F))

plot_pval_ak(fpval_2,gridbeta0,alpha=0.1,cutoff=0.12)


##### P-values for specification (3) for all methods over grid

pval_3<-sapply(gridbeta0, function(e) b_pval_hom(Y_3,Y_pz_3,Y_W_3,z_3,
                                                 e,hom_3,est_3[2,1],
                                                 est_3[2,2]^2,xnorm=0,
                                                 para=F,sim=F))
fpval_3<-matrix(0,8,ngrid)
fpval_3[1,]<-pval_3[1,]
fpval_3[2,]<-pval_3[2,]
fpval_3[3,]<-sapply(gridbeta0, function(e) CLR(modeliv_3,beta0=e)$p.value)
fpval_3[4,]<-pval_3[3,]
fpval_3[5,]<-sapply(gridbeta0, function(e) pchisq((est_3[1,1]-e)^2/est_3[1,2]^2,
                                                  df=1,lower.tail=F))
fpval_3[6,]<-pval_3[4,]
fpval_3[7,]<-sapply(gridbeta0, function(e) pchisq((est_3[3,1]-e)^2/est_3[3,2]^2,
                                                  df=1,lower.tail=F))
fpval_3[8,]<-sapply(gridbeta0, function(e) pchisq((est_3[4,1]-e)^2/est_3[4,2]^2,
                                                  df=1,lower.tail=F))

plot_pval_ak(fpval_3,gridbeta0,alpha=0.1,cutoff=0.12)


##### P-values for specification (4) for all methods over grid

pval_4<-sapply(gridbeta0, function(e) b_pval_hom(Y_4,Y_pz_4,Y_W_4,z_4,
                                                 e,hom_4,est_4[2,1],
                                                 est_4[2,2]^2,xnorm=0,
                                                 para=F,sim=F))
fpval_4<-matrix(0,8,ngrid)
fpval_4[1,]<-pval_4[1,]
fpval_4[2,]<-pval_4[2,]
fpval_4[3,]<-sapply(gridbeta0, function(e) CLR(modeliv_4,beta0=e)$p.value)
fpval_4[4,]<-pval_4[3,]
fpval_4[5,]<-sapply(gridbeta0, function(e) pchisq((est_4[1,1]-e)^2/est_4[1,2]^2,
                                                  df=1,lower.tail=F))
fpval_4[6,]<-pval_4[4,]
fpval_4[7,]<-sapply(gridbeta0, function(e) pchisq((est_4[3,1]-e)^2/est_4[3,2]^2,
                                                  df=1,lower.tail=F))
fpval_4[8,]<-sapply(gridbeta0, function(e) pchisq((est_4[4,1]-e)^2/est_4[4,2]^2,
                                                  df=1,lower.tail=F))

plot_pval_ak(fpval_4,gridbeta0,alpha=0.1,cutoff=0.12)
