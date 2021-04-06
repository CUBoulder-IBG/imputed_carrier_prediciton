#!/usr/bin/env Rscript
n1 = as.numeric(commandArgs(trailingOnly=TRUE))[1]
n2 = as.numeric(commandArgs(trailingOnly=TRUE))[2]
setwd("/work/KellerLab/meng/acc")
#b=read.table("GLM_b0_b1.txt",head=F)[,2:3]
library(data.table)
dos=fread("ukb.vcf",head=F)
y=fread("exom.data",head=F)
ukb_sample_id=fread("head_ukb",head=F)
exom_sample_id=fread("head_exom",head=F)
info=fread("info_ukb",head=F)
index=match(as.matrix(ukb_sample_id[,1]),as.matrix(exom_sample_id[,1]))
index1=which(!is.na(index))+1
index2=na.omit(index)+1
index_snp=match(as.matrix(dos[,1]),as.matrix(y[,1]))
index1_snp=which(!is.na(index_snp))
index2_snp=na.omit(index_snp)
index_info=match(as.matrix(dos[index1_snp,1]),as.matrix(info[,1]))
info=as.matrix(info[index_info,1])

dos1=dos[index1_snp,]
dos=dos1[,..index1]
y=y[index2_snp,..index2]
index=sample(1:ncol(dos),1e5)
dos_test=dos[,..index]
dos_ref=dos[,-(..index)]
y_test=y[,..index]
y_ref=y[,-(..index)]
n_gfg=21486
b1_dos=0.9983
b0_dos=-2.395e-8

results=matrix(NA,nrow(y)*100,10)
for(i in n1:n2){
MAF=sum(na.omit(abs(as.numeric(y[i,])-2)))/(ncol(y)*2)
DATA=na.omit(cbind(as.numeric(dos_test[i,]),abs(as.numeric(y_test[i,])-2)))
for(j in 1:100){
cutoff=seq(0,0.99,0.01)[j]
p_dos=exp(DATA[,1]*b1_dos+b0_dos)/(1+exp(DATA[,1]*b1_dos+b0_dos))
#p_dos=exp(DATA[,1]*b[i,2]+b[i,1])/(1+exp(DATA[,1]*b[i,2]+b[i,1]))
data=cbind(p_dos,DATA[,2])
index1=which(data[,1]<=cutoff)
if(length(index1)>0)
data[index1,1]=0
###calculate false positive, ture positive, false negative, ture negative for each SNP
TN=length(which(apply(data,1,sum)==0))
FP_TP=length(which(data[,1]>0))
TP=length(intersect(which(data[,1]>0),which(data[,2]>0)))
FN_TP=sum(data[,2])
#Sensitivity,#of detected carrier,PPV,expected # of detected carriers in GfG
results[(i-1)*100+j,4]=TP/FN_TP
results[(i-1)*100+j,5]=TP/FP_TP
results[(i-1)*100+j,6]=(TP/FN_TP)*MAF*2*n_gfg
results[(i-1)*100+j,7]=TP
results[(i-1)*100+j,8]=TN
results[(i-1)*100+j,9]=FP_TP-TP
results[(i-1)*100+j,10]=FN_TP-TP
results[(i-1)*100+j,3]=MAF
results[(i-1)*100+j,2]=cutoff
results[(i-1)*100+j,1]=info[i]
print(i)
}
}

write.table(results,paste("ROC_",n1,"_results.txt",sep=''),col.names=F,row.names=F,quote=F)



