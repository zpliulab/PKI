da <- as.matrix(read.csv("Ecoli1-trajectories100.csv",header=T,row.names = 1))
#da <- as.matrix(read.csv("GSE44183_human_expression_mat.txt",header=T,row.names = 1))
timepoint <- row.names(da)
nodename <- colnames(da)
rr<-c(rep(1,21))
tt<-as.numeric(timepoint)
da<-apply(t(da),2,as.numeric)

source('d:\\A\\work\\scRNA\\cgy\\testRcode\\testmodel\\fpca.R')

result_pca <-PCA(da,tt,rr,miss=0,ngrid=50)
time <-as.matrix(result_pca[["out21"]])
dim(time)
data1 <-as.matrix(result_pca[["yfit"]][1,])# 特征根1 拟合数据
dim(data1)
data2 <-as.matrix(result_pca[["yfit"]][2,])# 特征根1 拟合数据
data3 <-as.matrix(result_pca[["phi"]][,1])

data4 <-as.matrix(result_pca[["yfit_orig"]][1,])#原始数据估计
#fitdatainorigintimepoint <- c(data3[1],data3[100],data3[200],data3[300],data3[500],data3[600])
#data3 <-as.matrix(result_pca[["phi"]][,2])
#pdf(paste("pca_data_plot","_11.pdf",sep = ""))
p22<-plot(time,data1,col="orangered",lwd=2,lty=1,type="p",xlim = c(1,200),ylim=c(0,1),xlab = "Time",ylab = "Expression Value")+
  lines(time,data2,col="orchid",lwd=2,lty=1,type="p")+
  lines(time,result_pca[["yfit"]][3,],col="palegoldenrod",lwd=2,lty=1,type="p")+
  lines(time,result_pca[["yfit"]][4,],col="paleturquoise",lwd=2,lty=1,type="p")+
  lines(time,result_pca[["yfit"]][5,],col="palegreen",lwd=2,lty=1,type="p")+
  lines(time,data2,col="orchid",lwd=2,lty=1)+
  lines(time,result_pca[["yfit"]][3,],col="palegoldenrod",lwd=2,lty=1)+
  lines(time,result_pca[["yfit"]][4,],col="paleturquoise",lwd=2,lty=1)+
  lines(time,result_pca[["yfit"]][5,],col="palegreen",lwd=2,lty=1)+
  lines(time,data1,col="orangered",lwd=2,lty=1)+
  points(tt,da[1,],col="blue",pch=16,lty=2)+
  points(tt,data4,col="orange",lty=2,pch=17)+
  points(tt,da[2,],col="blue",pch=16,lty=2)+
  points(tt,result_pca[["yfit_orig"]][2,],col="orange",lty=2,pch=17)+
  points(tt,da[5,],col="blue",pch=16,lty=2)+
  points(tt,result_pca[["yfit_orig"]][5,],col="orange",lty=2,pch=17)+
  points(tt,da[3,],col="blue",pch=16,lty=2)+
  points(tt,result_pca[["yfit_orig"]][3,],col="orange",lty=2,pch=17)+
  points(tt,da[4,],col="blue",pch=16,lty=2)+
  points(tt,result_pca[["yfit_orig"]][4,],col="orange",lty=2,pch=17)+
  legend(150,0.8,c("G1","G2","G3","G4","G5","Estimate","Origin"),col=c("orangered","orchid","palegoldenrod","paleturquoise","palegreen","orange","blue"),
         text.col = "black",lwd=c(2,2,2,2,2,1,1), lty = c(1,1,1, 1, 1,-1,-1),pch = c(-1,-1,-1,-1,-1,16,17),
         merge = TRUE, bg = 'white')+
  abline(h = seq(0,1,0.1), col = "lightgray", lty = 3)

############new plot##########
g5<-result_pca[["yfit"]][1:5,]
G5<-t(g5)
G5df<-data.frame(G5)
colnames(G5df)<-c("G1","G2","G3","G4","G5")
G5df['Time']<-time
G5gather<-gather(G5df,Gene,Expression, -Time)

G5origin<-t(da[1:5,])
G5origindf<-data.frame(G5origin)
colnames(G5origindf)<-c("G1","G2","G3","G4","G5")
G5origindf['Time']<-tt
G5origingather<-gather(G5origindf,Gene,Expression, -Time)

G5estimate<-t(result_pca[["yfit_orig"]][1:5,])
G5estimatedf<-data.frame(G5estimate)
colnames(G5estimatedf)<-c("G1","G2","G3","G4","G5")
G5estimatedf['Time']<-tt
G5estimategather<-gather(G5estimatedf,Gene,Expression, -Time)

p_lineplot_g5<-ggplot(G5gather,aes(Time,Expression, fill=Gene))+
  geom_point(data=G5origingather,size=1.6,color='#45A7C8',shape=5,stroke=1.06)+
  geom_point(data=G5estimategather,size=1.8,color='#ffad64',shape=1,stroke=1.06)+
  geom_smooth(method="auto",aes(color=Gene),size=1,alpha=0.1)+
  theme_bw()+
  scale_fill_manual(values= pal2(5))+
  scale_color_manual(values= pal2(5))+
  facet_grid(cols = NULL)


m <- result_pca[["yfit"]]
rownames(m)<- rownames(da)
colnames(m) <- time
#write.csv(m,'fpca_E1_100.csv')

###############################################################################################
datadf<-m[1:5,]

realdf<-da[1:5,]
dim(datadf) #6 130
#row.names(datadf)<-rown
tt.grid<-colnames(datadf)
log_derivmatrix=list()
n=dim(datadf)[1]
tt.grid.d<-tt #这里修改微分取点数
#index<-seq(1,30,3)
#index<-append(index,30)
#tt.grid.d<-tt.grid
for (i in 1:n) {
  rowdata=datadf[i,]
  f <- splinefun(tt.grid,rowdata,method = "natural")
  derivrow=f(tt.grid.d, deriv = 1)
  log_derivmatrix=cbind(log_derivmatrix,derivrow)
}
row.names(log_derivmatrix)<-tt.grid.d#这里修改时间（行名）
log_derivmatrix<-t(log_derivmatrix)
row.names(log_derivmatrix)<-rown
dim(log_derivmatrix)
#dev.off() 
log_derivmatrix<-t(log_derivmatrix)
#log_derivmatrix 5 21
#datadf fpca 5 50
#realdf 5 21

#构建所需要的数据矩阵
#分组
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(ggthemes)
library(scales)
df1<-data.frame(t(datadf))
colnames(df1)<-c('G1','G2','G3','G4','G5')
df1["Time"]<-c(time)
df1 = gather(df1,key = "gene",value = "expression",-Time) 
p=ggplot(data = df1, aes(x = Time, y = expression, color =gene)) +
  geom_line()+geom_point()+theme_cleveland()

##########################################################################################

setwd("d:\\A\\work\\PKI_R")
rs2 <- as.matrix(read.csv("对比.csv",header=T))
p21<-plot(rs2[,1],rs2[,2],xlim = c(0,82),ylim=c(0,1),xlab = "The Serial Number",ylab = "R-squared")+
  points(rs2[,1],rs2[,3],col="orangered")+
  abline(h = seq(0,1,0.1), col = "lightgray", lty = 3)+
  legend(0,0.4,c("Based on Difference Value","Based on Differential Value"),col=c("black","orangered"),
       text.col = "black",pch = c(1,1),
       merge = TRUE, bg = 'white')
rs2
boxplot( rs2[,2:3], col = "grey", notch = T, horizontal = TRUE)
rs2d=data.frame(rs2)

library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(ggpubr)
gatherdata = gather(rs2d,EstimateMethod, rsquared,Based.on.Difference.Value, Based.on.Differential.Value) 
p_boxplot_Rs<-ggplot(gatherdata,aes(EstimateMethod,rsquared,fill=EstimateMethod)) + 
  geom_boxplot(width=0.5,outlier.color = cols2[3],outlier.shape = 2,outlier.size = 3) + 
  geom_jitter(shape=16, position=position_jitter(0.2),color=cols2[2])+
  stat_boxplot(geom = "errorbar",width=0.1) +
  stat_summary(fun="mean",geom="point",shape=23,
             size=4,fill="white")+
  #scale_fill_brewer(palette="Set1")+
  theme_bw()+
  scale_fill_manual(values = pal2(2))
   


##增加调色盘
cols<-c('#F6BD60','#F7EDE2','#F5CAC3','#84A59D','#F28482')
#cols<-c('#E64E00','#65B48E','#E6EB00','#E64E00')
pal<-colorRampPalette(cols)
image(x=1:20,y=1,z=as.matrix(1:20),col=pal(20))


cols2<-c('#F36993','#FFD166','#06D6A0','#68ACEA','#118AB2')
#cols2<-c('#E64E00','#65B48E','#E6EB00','#E64E00')
pal2<-colorRampPalette(cols2)
image(x=1:5,y=1,z=as.matrix(1:5),col=pal2(5))

col3<-c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3')
pal3<-colorRampPalette(cols3)
image(x=1:20,y=1,z=as.matrix(1:20),col=pal2(20))



library(patchwork)
# 第一次使用需要安装
wrap_plots(p21,p22,nrow=2, guides="collect")

