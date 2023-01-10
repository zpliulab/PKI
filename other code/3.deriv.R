setwd("d:\\A\\work\\scRNA\\2021\\Rcode")

saver2<-as.matrix(read.csv("saverdata2.csv",header=T))
row.names(saver2)<-saver2[,1]
saver2 <- saver2[,-1]
saver2<-apply(saver2,2,as.numeric)
rr1<-c(92,102,66,172,138,188)
#这组数据是经过saver之后的GSE数据，共758列

saversixpoints<-as.matrix(read.csv("saver-6points.csv",header=T))
row.names(saversixpoints)<-saversixpoints[,1]
genesymbol<-saversixpoints[,1]
saversixpoints <- saversixpoints[,-1]
saversixpoints<-apply(saversixpoints,2,as.numeric)
rr2<-c(1,1,1,1,1,1)
#这组数据是处理成了六个时间点的，取的mean

tt<-c(0,12,24,36,72,96)

#PCA对原始数据的正则化或预处理敏感（相对缩放），所以到这一步仍然未使用任何平滑化处理方法（log）

source('fpca.R')

result_pca <-PCA(saver2,tt,rr1,miss=0,ngrid=30)#修改一下 存一下30个点的

result_pca6 <-PCA(saversixpoints,tt,rr1,miss=0,ngrid=15)


time <-as.matrix(result_pca[["out21"]])

dim(time)

data1 <-as.matrix(result_pca[["yfit"]][1,])# 特征根1 拟合数据1

dim(data1)

data2 <-as.matrix(result_pca[["yfit"]][2,])# 特征根1 拟合数据2

phi <-as.matrix(result_pca[["phi"]][,1])## phi: estimated eigenfunctions at original time points tt.

data3 <-as.matrix(result_pca[["yfit_orig"]][1,])#原始数据估计
fitdatainorigintimepoint <- c(data3[1],data3[100],data3[200],data3[300],data3[500],data3[600])


plot(time,data1,col=2,lwd=2,lty=2,type="l",xlim = c(1,96),ylim=c(-5,10.5))
#lines(time,data2,col="green",lwd=2,lty=3)
points(tt,fitdatainorigintimepoint,col="blue",pch = 21,lty=1)
points(tt,saversixpoints[1,],col="orange",pch=23,lty=2)
legend(80.05,9.5,c("data1","data3","origin"),col=c("red","blue","orange"),
       text.col = "purple4",lwd=c(2,2,1), lty = c(2, 3, -1),pch = c(-1,-1,21),
       merge = TRUE, bg = 'white')



#以data1为例
f <- splinefun(time,data1,method = "natural")
f <- splinefun(x, y)
ls(envir = environment(f))
splinecoef <- get("z", envir = environment(f))
curve(f(x), 1, 96, col = "green", lwd = 1.5)

points(splinecoef, col = "purple", cex = 2)
curve(f(x, deriv = 1), 1, 96, col = 2, lwd = 1.5)
curve(f(x, deriv = 2), 1, 10, col = 2, lwd = 1.5, n = 401)
curve(f(x, deriv = 3), 1, 10, col = 2, lwd = 1.5, n = 401)
par(op)

#对yfit矩阵进行计算

curv = result_pca[['yfit']]
n=dim(curv)[1]
tt.grid = time
derivmatrix=list()

for (i in 1:n) {
  rowdata=curv[i,]
  f <- splinefun(tt.grid,rowdata,method = "natural")
  derivrow=f(tt.grid, deriv = 1)
  derivmatrix=cbind(derivmatrix,derivrow)
}
colnames(derivmatrix)<-genesymbol

write.csv(derivmatrix,"derivmatrix30.csv")
derivmatrix<-apply(derivmatrix,2,as.numeric)
write.csv(log2(derivmatrix+1),'llderivmatrix30.csv')
write.csv(curv,"fpca30.csv")

#是否对yfit进行平滑化处理？
lcurv = log2(curv+1)
derivmatrix_log=list()
for (i in 1:n) {
  rowdata=lcurv[i,]
  f <- splinefun(tt.grid,rowdata,method = "natural")
  derivrow=f(tt.grid, deriv = 1)
  derivmatrix_log=cbind(derivmatrix_log,derivrow)
}
colnames(derivmatrix)<-genesymbol
derivmatrix_log=as.matrix(derivmatrix_log)
write.csv(derivmatrix_log,"derivmatrix_log30.csv")
write.csv(lcurv,"log_fpca30.csv")

#先进行log2平滑化再计算PCA
log_saver2=log2(saver2+1)
log_saver2<-apply(log_saver2,2,as.numeric)
result_pca_log_saver <-PCA(log_saver2,tt,rr1,miss=0,ngrid=30)

log_curv = result_pca_log_saver[['yfit']]
log_derivmatrix=list()

for (i in 1:n) {
  rowdata=log_curv[i,]
  f <- splinefun(tt.grid,rowdata,method = "natural")
  derivrow=f(tt.grid, deriv = 1)
  log_derivmatrix=cbind(log_derivmatrix,derivrow)
}
colnames(log_derivmatrix)<-genesymbol
log_derivmatrix=as.matrix(log_derivmatrix)
write.csv(log_derivmatrix,"log_derivmatrix30.csv")
write.csv(log_curv,"saverlog_fpca30.csv")


#对已经经过FPCA扩增后的数据处理
getwd()
setwd("d:\\A\\work\\PKI_R\\expressiondata")
genelist<-read.csv("d:\\A\\work\\PKI_R\\150节点\\genelist.csv",row.names = 1)
timedata30<-read.csv("saverlogfpca.csv",header=T, check.names = F)

row.names(timedata30)<-make.names(timedata30[,1],TRUE)
timedata30<-timedata30[,-1]
dim(timedata30)
node<-row.names(timedata30)
datadf<-rbind()
rown<-list()
dim(timedata30)
for (i in node) {
  if(i %in% genelist[,1]) {
    print(i)
    datadf<-rbind(datadf,timedata30[i,])
    rown<-append(rown,i)
  }
  # else
  # {
  #   print('f')
  # }
}
dim(datadf) #6 130
row.names(datadf)<-rown
tt.grid<-colnames(datadf)
log_derivmatrix=list()
n=dim(datadf)[1]
tt.grid.d<-c(0,12,24,36,72,96) #这里修改微分取点数
index<-seq(1,30,3)
#index<-append(index,30)
tt.grid.d<-tt.grid
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

write.csv(log_derivmatrix,"derivmatrix30.csv")
write.csv(datadf,"saverlogfpca30.csv")