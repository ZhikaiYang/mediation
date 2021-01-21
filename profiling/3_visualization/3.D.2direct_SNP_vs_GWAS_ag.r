setwd('/common/jyanglab/zhikaiyang/projects/mediation/largedata/dsnps_vs_gwas')
library(data.table)
f=list.files(path="/common/jyanglab/zhikaiyang/projects/mediation/largedata/gwas/output",patt=".assoc.txt")
for(ff in f)
{
out=gsub("assoc.txt","sig.region.txt",ff)
ff = paste("/common/jyanglab/zhikaiyang/projects/mediation/largedata/gwas/output/", ff, sep="")
d1=fread(ff,head=T,data.table=F)
d1=d1[d1[,13]<=1e-5,]
s=d1[,3]-50000;e=d1[,3]+50000
p=d1$p_wald
d1=cbind(d1$chr,s,e,p)
bk=1
bl=1
colnames(d1)[1]="Chr"
for(j in 2:nrow(d1))
{if(d1[j,1]==d1[j-1,1])
{
  if(d1[j,2]-d1[j-1,3]>0)
  { bk=bk+1}
}
  else{bk=bk+1} 
  bl=c(bl,bk)
}
d1=cbind(d1,bl) 
index=unique(bl)
ro=NULL
for(k in index)
{
  d2=d1[d1[,5]==k,,drop=F]
  m=min(d2[,4])
  mp=d1[which(d1[,4]==m),2][1]+50000
  interval=c(d2[1,1],d2[1,2],max(d2[,3]),m,mp,nrow(d2))
  ro=rbind(ro,interval)
}
ro=as.data.frame(ro)
colnames(ro)=c("Chr","start","end","p-value","minP_pos","N")
write.table(ro,out,col.na=T,row.na=F,quote=F,sep="\t")
}

f=list.files(path=".",patt="sig.region.txt")
r=NULL
for(i in f)
{
  ti=gsub(".sig.region.txt","",i)
 d=fread(i,head=T,data.table=F)
  date=rep(ti,nrow(d))
  d1=cbind(date,d)
  r=rbind(r,d1)
}
write.table(r,"all_traits_ag.all_res.txt",col.na=T,row.na=F,quote=F,sep="\t")







