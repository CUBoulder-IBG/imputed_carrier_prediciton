
#Perform logistic regression on variants in UKB exome ~ UKB dosage

#Packages
library(data.table)
library(LaplacesDemon)
library(matrixStats)
#library(fastglm) #this package sucks

#Working Directory
setwd("/work/KellerLab/meng/acc")



#######################
#1 Getting data
#Read in data
dos.orig <- fread("ukb.vcf",head=F,data.table=F) #the imputed data for these 123 variants and 487410 people
exom.orig <- fread("exom.data",head=F,data.table=F) #exomic data; these don't line up
dos_sample_id <- unlist(fread("head_ukb",head=F,data.table=F))
exom_sample_id <- unlist(fread("head_exom",head=F,data.table=F))
r2.orig <- fread("info_ukb",head=F,data.table=F)
rownames(exom.orig) <- exom.orig$V1
rownames(dos.orig) <- dos.orig$V1
colnames(r2.orig) <- c("SNP","r2")

#Make exom.orig & dos.orig into matrices in order to work with them more easily
exom.mat <- as.matrix(exom.orig[,-1])
colnames(exom.mat) <- exom_sample_id
dos.mat <- as.matrix(dos.orig[,-1])
colnames(dos.mat) <- dos_sample_id
#######################




#######################
#2 Formatting data
#Put the dosage data COLUMNS (individs) in same order as the exomic data by keying the exom data
index.indivs.exom <- match(exom_sample_id,dos_sample_id)
dos1 <- dos.mat[,index.indivs.exom]

#Note there are 563 indivs in the exome data who aren't in the dosage data
summary(index.indivs.exom)

#Remove these from both exom.mat and dos1
nas.to.remove <- which(is.na(index.indivs.exom))
dos1 <- dos1[, - nas.to.remove]
exom.mat <- exom.mat[, - nas.to.remove] 

#there are now 200643-563=200080 indivs in both exom.mat and dos1
ncol(dos1)
ncol(exom.mat)

#Put the dosage data ROWS (variants) in same order as the exomic data by keying the exom data
index.snps.exom <- match(rownames(exom.mat),rownames(dos1))
dos2 <- dos1[index.snps.exom,]

#There are some values in dos2 that aren't in exom.mat. Remove them from both
rows.to.remove <- which(is.na(rownames(dos2)))
dos <- dos2[- rows.to.remove,]
exom <- 2 - exom.mat[- rows.to.remove,] #this also reverse codes so that 0 = no minor allele
exom[exom==2] <- 1 #there are 4565 2's and 63536 1's. Make the 2's 1's so we can fit logistic reg

#Now get the r2 information in the same order as above
(index.r2.exom <- match(rownames(exom),r2.orig$SNP))
r2 <- r2.orig[index.r2.exom,]

#Get MAF & MAC in exome data
r2$maf.exom <- apply(exom,1,function(x) mean(x,na.rm=T)/2)
r2$mac.exom <- apply(exom,1,sum,na.rm=TRUE)


#Remove SNPs where MAC < 3; these create fit problems below and just won't be predictable (or exist) in the GFG data anyway
maf.exom <- apply(exom,1,function(x) mean(x,na.rm=T)/2)
mac.exom <- apply(exom,1,sum,na.rm=TRUE)
(keep.snps.mac3 <- mac.exom >2)
sum(keep.snps.mac3) #92 SNPs will be retained
nrow(r2) - sum(keep.snps.mac3) #31 SNPs will be dropped
dos3 <- dos[keep.snps.mac3,]
exom3 <- exom[keep.snps.mac3,]
r23 <- r2[keep.snps.mac3,]


#######################






#######################
#3 Running logistic regression

#Vectorize the data
exom3.t <- t(exom3)
exom.out <- c(exom3.t)
dos3.t <- t(dos3)
dos.pred <- c(dos3.t)
r23.pred <- rep(r23$r2,each=ncol(exom))

#Run GLM & time it #14.6 at n=20k; 25s at n=40K; 60s at n=80K
summary(mod1 <- glm(exom.out ~ dos.pred*r23.pred,family="binomial",na.action=na.exclude))
#summary(mod2 <- glm(exom.out ~ dos.pred+r23.pred,family="binomial",na.action=na.exclude))

#Now get predicted values for each point for the best (interaction) model above
pred.alleles <- predict(mod1,type="response") 
#type="response" gives the predicted probabilities whereas type="link" give predicted log-odds of probabilities

pred.mat <- matrix(pred.alleles,ncol=ncol(exom3),nrow=nrow(exom3),byrow=TRUE)

#quick check that pred.mat & exom3 are lining up - they appear to do so
mac.pred <- apply(pred.mat,1,sum,na.rm=T)
mac.e3 <- apply(exom3,1,sum,na.rm=T)
#plot(mac.pred,mac.e3)

#######################







#######################
#4 getting PPV and Sens for each variant based on the predicted values
#get predicted for each threshold
pred.mat10 <- (pred.mat>.95)*1
pred.mat9 <- (pred.mat>.9)*1
pred.mat8 <- (pred.mat>.8)*1
pred.mat7 <- (pred.mat>.7)*1
pred.mat6 <- (pred.mat>.6)*1
pred.mat5 <- (pred.mat>.5)*1
pred.mat4 <- (pred.mat>.4)*1
pred.mat3 <- (pred.mat>.3)*1
pred.mat2 <- (pred.mat>.2)*1
pred.mat1 <- (pred.mat>.1)*1

#For each of the above, get the TN, FP, FN, TPs
#Note: rowSums2 is about 3x faster than apply
tn10 <- rowSums2(exom3==0 & pred.mat10==0,na.rm=T)
fp10 <- rowSums2(exom3==0 & pred.mat10==1,na.rm=T)
fn10 <- rowSums2(exom3==1 & pred.mat10==0,na.rm=T)
tp10 <- rowSums2(exom3==1 & pred.mat10==1,na.rm=T)

tn9 <- rowSums2(exom3==0 & pred.mat9==0,na.rm=T)
fp9 <- rowSums2(exom3==0 & pred.mat9==1,na.rm=T)
fn9 <- rowSums2(exom3==1 & pred.mat9==0,na.rm=T)
tp9 <- rowSums2(exom3==1 & pred.mat9==1,na.rm=T)

tn8 <- rowSums2(exom3==0 & pred.mat8==0,na.rm=T)
fp8 <- rowSums2(exom3==0 & pred.mat8==1,na.rm=T)
fn8 <- rowSums2(exom3==1 & pred.mat8==0,na.rm=T)
tp8 <- rowSums2(exom3==1 & pred.mat8==1,na.rm=T)

tn7 <- rowSums2(exom3==0 & pred.mat7==0,na.rm=T)
fp7 <- rowSums2(exom3==0 & pred.mat7==1,na.rm=T)
fn7 <- rowSums2(exom3==1 & pred.mat7==0,na.rm=T)
tp7 <- rowSums2(exom3==1 & pred.mat7==1,na.rm=T)

tn6 <- rowSums2(exom3==0 & pred.mat6==0,na.rm=T)
fp6 <- rowSums2(exom3==0 & pred.mat6==1,na.rm=T)
fn6 <- rowSums2(exom3==1 & pred.mat6==0,na.rm=T)
tp6 <- rowSums2(exom3==1 & pred.mat6==1,na.rm=T)

tn5 <- rowSums2(exom3==0 & pred.mat5==0,na.rm=T)
fp5 <- rowSums2(exom3==0 & pred.mat5==1,na.rm=T)
fn5 <- rowSums2(exom3==1 & pred.mat5==0,na.rm=T)
tp5 <- rowSums2(exom3==1 & pred.mat5==1,na.rm=T)

tn4 <- rowSums2(exom3==0 & pred.mat4==0,na.rm=T)
fp4 <- rowSums2(exom3==0 & pred.mat4==1,na.rm=T)
fn4 <- rowSums2(exom3==1 & pred.mat4==0,na.rm=T)
tp4 <- rowSums2(exom3==1 & pred.mat4==1,na.rm=T)

tn3 <- rowSums2(exom3==0 & pred.mat3==0,na.rm=T)
fp3 <- rowSums2(exom3==0 & pred.mat3==1,na.rm=T)
fn3 <- rowSums2(exom3==1 & pred.mat3==0,na.rm=T)
tp3 <- rowSums2(exom3==1 & pred.mat3==1,na.rm=T)

tn2 <- rowSums2(exom3==0 & pred.mat2==0,na.rm=T)
fp2 <- rowSums2(exom3==0 & pred.mat2==1,na.rm=T)
fn2 <- rowSums2(exom3==1 & pred.mat2==0,na.rm=T)
tp2 <- rowSums2(exom3==1 & pred.mat2==1,na.rm=T)

tn1 <- rowSums2(exom3==0 & pred.mat1==0,na.rm=T)
fp1 <- rowSums2(exom3==0 & pred.mat1==1,na.rm=T)
fn1 <- rowSums2(exom3==1 & pred.mat1==0,na.rm=T)
tp1 <- rowSums2(exom3==1 & pred.mat1==1,na.rm=T)

s10 <- tp10/(tp10+fn10)
p10 <- tp10/(tp10+fp10)
p10[is.nan(p10)] <- 0

s9 <- tp9/(tp9+fn9)
p9 <- tp9/(tp9+fp9)
p9[is.nan(p9)] <- 0

s8 <- tp8/(tp8+fn8)
p8 <- tp8/(tp8+fp8)
p8[is.nan(p8)] <- 0

s7 <- tp7/(tp7+fn7)
p7 <- tp7/(tp7+fp7)
p7[is.nan(p7)] <- 0

s6 <- tp6/(tp6+fn6)
p6 <- tp6/(tp6+fp6)
p6[is.nan(p6)] <- 0

s5 <- tp5/(tp5+fn5)
p5 <- tp5/(tp5+fp5)
p5[is.nan(p5)] <- 0

s4 <- tp4/(tp4+fn4)
p4 <- tp4/(tp4+fp4)
p4[is.nan(p4)] <- 0

s3 <- tp3/(tp3+fn3)
p3 <- tp3/(tp3+fp3)
p3[is.nan(p3)] <- 0

s2 <- tp2/(tp2+fn2)
p2 <- tp2/(tp2+fp2)
p2[is.nan(p2)] <- 0

s1 <- tp1/(tp1+fn1)
p1 <- tp1/(tp1+fp1)
p1[is.nan(p1)] <- 0

RES <- cbind(r23,s10,p10,s9,p9,s8,p8,s7,p7,s6,p6,s5,p5,s4,p4,s3,p3,s2,p2,s1,p1,
             tp10,fp10,fn10,
             tp9,fp9,fn9,
             tp8,fp8,fn8,
             tp7,fp7,fn7,
             tp6,fp6,fn6,
             tp5,fp5,fn5,
             tp4,fp4,fn4,
             tp3,fp3,fn3,
             tp2,fp2,fn2,
             tp1,fp1,fn1)


#Plot these results, one SNP per page
pdf("PPV.SEN.plot.pdf",width=12,height=12)
for (i in 1:nrow(RES)){
  cursnp <- RES$SNP[i]
  x <- unlist(RES[i,2:ncol(RES)])
  sen <- x[seq(4,23,2)]
  ppv <- x[seq(5,24,2)]
  plot(ppv,sen,main=paste0("Sens & PPV ",cursnp, "; MAC=",x["mac.exom"]," / 200643; r2=",round(RES$r2[i],3)),xlab="PPV",ylab="SENS",xlim=0:1,ylim=0:1)  
}
  dev.off()
  
  
  
#Get some estimates relevant to GfG 
  ng <- 30000 #estimated number in GfG
  RES$Enum.carr <- (RES$maf.exom^2 + 2*RES$maf.exom*(1-RES$maf.exom))*ng
  
  RES$Enum.called10 <-  ((RES$tp10+RES$fn10)/ncol(dos3))*ng
  RES$Enum.true.called10 <-  RES$Enum.carr * RES$s10

  RES$Enum.called9 <-  ((RES$tp9+RES$fn9)/ncol(dos3))*ng
  RES$Enum.true.called9 <-  RES$Enum.carr * RES$s9

  RES$Enum.called8 <-  ((RES$tp8+RES$fn8)/ncol(dos3))*ng
  RES$Enum.true.called8 <-  RES$Enum.carr * RES$s8

  RES$Enum.called7 <-  ((RES$tp7+RES$fn7)/ncol(dos3))*ng
  RES$Enum.true.called7 <-  RES$Enum.carr * RES$s7

  RES$Enum.called6 <-  ((RES$tp6+RES$fn6)/ncol(dos3))*ng
  RES$Enum.true.called6 <-  RES$Enum.carr * RES$s6

  RES$Enum.called5 <-  ((RES$tp5+RES$fn5)/ncol(dos3))*ng
  RES$Enum.true.called5 <-  RES$Enum.carr * RES$s5

  RES$Enum.called4 <-  ((RES$tp4+RES$fn4)/ncol(dos3))*ng
  RES$Enum.true.called4 <-  RES$Enum.carr * RES$s4

  RES$Enum.called3 <-  ((RES$tp3+RES$fn3)/ncol(dos3))*ng
  RES$Enum.true.called3 <-  RES$Enum.carr * RES$s3

  RES$Enum.called2 <-  ((RES$tp2+RES$fn2)/ncol(dos3))*ng
  RES$Enum.true.called2 <-  RES$Enum.carr * RES$s2

  RES$Enum.called1 <-  ((RES$tp1+RES$fn1)/ncol(dos3))*ng
  RES$Enum.true.called1 <-  RES$Enum.carr * RES$s1
  
  #Look
  RES[,c(1:10,55:70)]
  
#######################


 
  
#nas.dos3 <-  apply(dos3,1,function(x) sum(is.na(x)))
  
#nas.dos <-  apply(dos,1,function(x) sum(is.na(x)))
  
  
  
  
#######################
#5 getting PPV and Sens for each variant based only on dosage
  
dos.rank <- ncol(dos3) - rowRanks(dos3,ties.method="random")

  #check to make sure this worked
  summary(t(dos.rank))
  lo1 <- which(dos.rank[1,]<25)
  dos3[1:10,lo1] #yes - the first row all have relatively high dosages
  
  #get predicted for each threshold of ranked dosages
  dos.rank200 <- (dos.rank < 200)*1
  dos.rank100 <- (dos.rank < 100)*1
  dos.rank50 <- (dos.rank < 50)*1
  dos.rank25 <- (dos.rank < 25)*1
  
  #For each of the above, get the TN, FP, FN, TPs
  #Note: rowSums2 is about 3x faster than apply
  tn200 <- rowSums2(exom3==0 & dos.rank200==0,na.rm=T)
  fp200 <- rowSums2(exom3==0 & dos.rank200==1,na.rm=T)
  fn200 <- rowSums2(exom3==1 & dos.rank200==0,na.rm=T)
  tp200 <- rowSums2(exom3==1 & dos.rank200==1,na.rm=T)
  
  tn100 <- rowSums2(exom3==0 & dos.rank100==0,na.rm=T)
  fp100 <- rowSums2(exom3==0 & dos.rank100==1,na.rm=T)
  fn100 <- rowSums2(exom3==1 & dos.rank100==0,na.rm=T)
  tp100 <- rowSums2(exom3==1 & dos.rank100==1,na.rm=T)
  
  tn50 <- rowSums2(exom3==0 & dos.rank50==0,na.rm=T)
  fp50 <- rowSums2(exom3==0 & dos.rank50==1,na.rm=T)
  fn50 <- rowSums2(exom3==1 & dos.rank50==0,na.rm=T)
  tp50 <- rowSums2(exom3==1 & dos.rank50==1,na.rm=T)
  
  tn25 <- rowSums2(exom3==0 & dos.rank25==0,na.rm=T)
  fp25 <- rowSums2(exom3==0 & dos.rank25==1,na.rm=T)
  fn25 <- rowSums2(exom3==1 & dos.rank25==0,na.rm=T)
  tp25 <- rowSums2(exom3==1 & dos.rank25==1,na.rm=T)
  
  s200 <- tp200/(tp200+fn200)
  p200 <- tp200/(tp200+fp200)
  p200[is.nan(p200)] <- 0
  
  s100 <- tp100/(tp100+fn100)
  p100 <- tp100/(tp100+fp100)
  p100[is.nan(p100)] <- 0
  
  s50 <- tp50/(tp50+fn50)
  p50 <- tp50/(tp50+fp50)
  p50[is.nan(p50)] <- 0
  
  s25 <- tp25/(tp25+fn25)
  p25 <- tp25/(tp25+fp25)
  p25[is.nan(p25)] <- 0
  
  RES.rank <- cbind(r23,s200,p200,
               s100,p100,
               s50,p50,
               s25,p25,
               tp200,fp200,fn200,
               tp100,fp100,fn100,
               tp50,fp50,fn50,
               tp25,fp25,fn25)
  
  #Plot these results, one SNP per page
  pdf("PPV.SEN.plot.ranks.pdf",width=12,height=12)
  for (i in 1:nrow(RES.rank)){
    cursnp <- RES.rank$SNP[i]
    x <- unlist(RES.rank[i,2:ncol(RES.rank)])
    sen <- x[seq(4,10,2)]
    ppv <- x[seq(5,11,2)]
    
    plot(ppv,sen,main=paste0("Sens & PPV ",cursnp, "; MAC=",x["mac.exom"]," / 200643; r2=",round(RES.rank$r2[i],3)),xlab="PPV",ylab="SENS",xlim=0:1,ylim=0:1,pch=c("v","h","m","l"))
   mtext("Ranks: v=200; h=100; m=50; l=25",3,0)
  }
  dev.off()
  
#Unfortunately, a rank-based approach makes it hard to extrapolate to different samples with different sample sizes because a given rank has a different meaning depending on sample size. Being ranked 50th highest dosage out of 2M people is much higher probability of being carrier than being ranked 50th out of 100 people.
#######################
  
  
  
  
 
#######################
#6 write out results  
  
  write.table(RES,"Res.pred.prob.txt",quote=F,col.names=T,row.names=F)
  write.table(RES.rank,"Res.ranks.txt",quote=F,col.names=T,row.names=F)
  
#######################
  
  
  
  
  
  
  
  