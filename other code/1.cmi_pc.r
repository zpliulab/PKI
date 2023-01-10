rm(list = ls())
setwd("d:\\A\\work\\PKI_R\\150")
getwd()
name <- "saverlog6timedata.csv"
data <- as.matrix(read.csv(name,header=T))

row.names(data)<-data[,1]
data <- data[,-1]
dim(data)

data <- t(data)
timepoint <- row.names(data)
nodename <- colnames(data)

relall <- as.matrix(read.csv("relall758.csv",header=T))
backgroundnetwork=relall[,2:3]
dim(backgroundnetwork)#417 2

node1=unique(backgroundnetwork[,1]) #41
node2=unique(backgroundnetwork[,2]) #131
node=unique(c(node1,node2)) #131 

datadf<-cbind()
rown<-list()
dim(data)
for (i in node) {
  if(i %in% nodename) {
    datadf<-cbind(datadf,data[,i])
    rown<-append(rown,i)
  }
  else
  {
    print(i)
  }
}
dim(datadf) #6 130
colnames(datadf)<-rown



ex_v <- apply(as.matrix(datadf),2,
              function(x)as.numeric(x)) 

dim(ex_v) #6 150 #6 130
pccmatrix<-cor(ex_v)
pccnode<-colnames(pccmatrix)
diag(pccmatrix)<-NA
pccmatrix<-abs(pccmatrix)
pccmatrix[!upper.tri(pccmatrix, diag = TRUE)] <- 0
pccindex=which(pccmatrix>0.4,arr.ind = T)
pcc<-list()
pccvalue<-pccmatrix[which(pccmatrix>0.4)]
pcc$num<-pccvalue
pcc$tf<-pccnode[pccindex[,1]]
pcc$target<-pccnode[pccindex[,2]]
dim(pccindex) #4253 2

pcc<-data.frame(pcc)
pccsort<-pcc[order(pcc$num,decreasing = T),]
write.csv(pccsort,"pccsort_real6logsaver_0.4.csv") 

colnames(backgroundnetwork)<-c("tf","target")
dim(backgroundnetwork)#417
backgroundnetwork_2<-rbind(backgroundnetwork,pccsort[,2:3])#535



#########high MI###########
nbins <- sqrt(NROW(ex_v))
dis_data <- discretize(ex_v,disc="equalfreq",nbins) 
colnames(dis_data) <- colnames(ex_v)
a<-mutinformation(dis_data,method="emp")
minode<-colnames(a)
diag(a)<-NA
mimatrix<-abs(a)
mimatrix[!upper.tri(mimatrix, diag = TRUE)] <- 0
miindex=which(mimatrix>0.55,arr.ind = T)
mi<-list()
mivalue<-mimatrix[which(mimatrix>0.55)]
mi$num<-mivalue
mi$tf<-minode[miindex[,1]]
mi$target<-minode[miindex[,2]]
dim(miindex) #1259
write.csv(mi,"mi_6logsaver_0.55.csv") #267 2
mi<-data.frame(mi)
misort<-mi[order(mi$num,decreasing = T),]

backgroundnetwork_3<-rbind(backgroundnetwork_2,misort[,2:3])
write.csv(backgroundnetwork_3,"bgnetwork2185.csv")
backgroundnetwork_3<-read.csv("bgnetwork2185.csv",header=T,row.names = 1)


backgroundnetwork_3<-read.csv("d:\\A\\work\\PKI_R\\150\\logsaver\\reall.csv",header=T,row.names = 1)
reN <- graph.data.frame(backgroundnetwork_3,directed = T)
reN_sim <- simplify(reN,remove.loops=T,remove.multiple = T)
#plot(reN_sim)
ed_reN_sim <- get.edgelist(reN_sim) 
dim(ed_reN_sim) #720 #390

dim(ex_v)

matPCC=rbind()
for(i in 1:dim(ed_reN_sim)[1]){
  
  PCC=cor(ex_v[,ed_reN_sim[i,1]],ex_v[,ed_reN_sim[i,2]],method = "pearson")
  if(is.na(PCC)=="FALSE"){
    if(abs(PCC)>0){
      one=c(ed_reN_sim[i,],PCC)
      matPCC=rbind(matPCC,one)
    }
  }
  
}

dim(matPCC)
write.csv(matPCC,"saver-pcc03792.csv") 


###### 0-MI  

library(infotheo)
ed_2=matPCC[,1:2]


nbins <- sqrt(NROW(ex_v))
dis_data <- discretize(ex_v,disc="equalfreq",nbins) 
colnames(dis_data) <- colnames(ex_v)

MI_0 <- rbind()
for(i in 1:dim(ed_2)[1])
{
  #i <- 1
  loc1 <- ed_2[i,1]
  loc2 <- ed_2[i,2]
  MI <- mutinformation(dis_data[,loc1],dis_data[,loc2],method="emp")
  MI_0 <- rbind(MI_0,MI)
}
median(MI_0)
mean(MI_0)
####### plot hist MI

hist(MI_0,main='MI_0_distribution',breaks =20)

MI_0_thre <-0.05663301
MI_0_left <- rbind()
for(i in 1:length(MI_0))
{
  if(MI_0[i] > MI_0_thre)
  {
    MI_0_left <- rbind(MI_0_left,c(ed_2[i,],MI_0[i]))
  }
}
dim(MI_0_left)
write.csv(MI_0_left,"MI_0_left_thre0.05663301.csv")

MI_0_net <- graph.data.frame(MI_0_left[,c(1,2)],directed = F)
pdf('./MI_0_net.pdf')
plot(MI_0_net, vertex.color="purple",vertex.size=8,
     label.font=2,label.cex=2,label.color='black',main='MI_0_net')
dev.off()

##### Calculate 1-MI if exists
library(pracma)
MI_1 <- rbind()
MI_2 <- rbind()
MI_3 <- rbind()
MI_re <- rbind()
Node_0 <- unlist(vertex_attr(MI_0_net))
for(i in 1:dim(MI_0_left)[1])
{
  #i <-1
  loc1 <- MI_0_left[i,1]
  loc2 <- MI_0_left[i,2]
  ### if they have shared one-order neighbor
  ne1 <- setdiff(Node_0[unlist(ego(MI_0_net, order = 1, MI_0_left[i,1]))],MI_0_left[i,1])
  ne2 <- setdiff(Node_0[unlist(ego(MI_0_net, order = 1, MI_0_left[i,2]))],MI_0_left[i,1])
  nn <- unique(setdiff(intersect(ne1,ne2),c(MI_0_left[i,1],MI_0_left[i,2])))
  if(isempty(nn)==F)
  {
    # 1-order
    for(j in 1:length(nn))
    {
      loc3 <- nn[j]
      con <- dis_data[,loc3]
      mi <- condinformation(dis_data[,loc1],dis_data[,loc2],con,method="emp")
      MI_1 <- rbind(MI_1,c(MI_0_left[i,c(1,2)],nn[j],1,mi))
    }
    
    
  }
  else{
    MI_re <-rbind( MI_re,c(MI_0_left[i,c(1,2)],0,'left'))
  }
}
colnames(MI_re) <- c('node1','node2','MI_order','states')
colnames(MI_1) <- c('node1','node2','neighbor1','MI_order','states')
dim(MI_1)
MI_1_value<-as.numeric(MI_1[,5])
median(MI_1_value)
mean(MI_1_value)
hist(MI_1_value,main='MI_1_distribution',breaks =50)

#### select the maximum MI and CMI
thre <-  0.174416
#### 1-order CMI
u1 <- unique(MI_1[,1])
MI_1_left_max <- rbind()
dim(MI_1)
for(i in 1:length(u1))
{
  loc <- which(MI_1[,1] %in% u1[i])
  k1 <- unique(MI_1[loc,2])
  can <- MI_1[loc,2]
  for(j in 1:length(k1))
  {
    z1 <- apply(as.matrix(MI_1[loc[which(can %in% k1[j])],dim(MI_1)[2]]),
                2,function(x) as.numeric(x))
    if(max(z1)> thre)
    {
      MI_1_left_max <- rbind(MI_1_left_max,c(u1[i],k1[j],max(z1)))
    }
  }
}
dim(MI_1_left_max) 
write.csv(MI_1_left_max,"MI_1_left_thre0.174416.csv")

MI_1_net <- graph.data.frame(MI_1_left_max[,c(1,2)],directed = F)
for(i in 1:dim(MI_1_left_max)[1])
{
  #i <-1
  loc1 <- MI_1_left_max[i,1]
  loc2 <- MI_1_left_max[i,2]
  ### if they have shared one-order neighbor
  ne1 <- setdiff(Node_0[unlist(ego(MI_1_net, order = 1, MI_1_left_max[i,1]))],MI_1_left_max[i,1])
  ne2 <- setdiff(Node_0[unlist(ego(MI_1_net, order = 1, MI_1_left_max[i,2]))],MI_1_left_max[i,1])
  nn <- unique(setdiff(intersect(ne1,ne2),c(MI_1_left_max[i,1],MI_1_left_max[i,2])))
  if(isempty(nn)==F)
  {
    
    # 2-order
    if(length(nn)>1)
    {
      for(k in 1:(length(nn)-1))
      {
        loc3 <- nn[k]
        for(b in (k+1):length(nn))
        {
          loc4 <- nn[b]
          con <- c(dis_data[,loc3],dis_data[,loc4])
          mi <- condinformation(dis_data[,loc1],dis_data[,loc2],con,method="emp")
          MI_2 <- rbind(MI_2,c(MI_0_left[i,c(1,2)],nn[k],nn[b],2,mi))
        }
      }
    }
  }
}

colnames(MI_2) <- c('node1','node2','neighbor1','neighbor2','MI_order','states')
dim(MI_2)

MI_2_value<-as.numeric(MI_2[,6])
median(MI_2_value)
mean(MI_2_value)
hist(MI_2_value,main='MI_2_distribution',breaks =50)


#### 2-order CMI
thre= 0.1591285
u2 <- unique(MI_2[,1])
MI_2_left_max <- rbind()
for(i in 1:length(u2))
{
  loc <- which(MI_2[,1] %in% u2[i])
  k1 <- unique(MI_2[loc,2])
  can <- MI_2[loc,2]
  for(j in 1:length(k1))
  {
    z1 <- apply(as.matrix(MI_2[loc[which(can %in% k1[j])],dim(MI_2)[2]]),
                2,function(x) as.numeric(x))
    if(max(z1) > thre)
    {
      MI_2_left_max <- rbind(MI_2_left_max,c(u2[i],k1[j],max(z1)))
    }
  }
}
dim(MI_2_left_max)
write.csv(MI_2_left_max,"MI_2_left_thre=0.1591285.csv")

MI_2_net <- graph.data.frame(MI_2_left_max[,c(1,2)],directed = F)
for(i in 1:dim(MI_2_left_max)[1])
{
  #i <-1
  loc1 <- MI_2_left_max[i,1]
  loc2 <- MI_2_left_max[i,2]
  ### if they have shared one-order neighbor
  ne1 <- setdiff(Node_0[unlist(ego(MI_2_net, order = 1, MI_1_left_max[i,1]))],MI_2_left_max[i,1])
  ne2 <- setdiff(Node_0[unlist(ego(MI_2_net, order = 1, MI_1_left_max[i,2]))],MI_2_left_max[i,1])
  nn <- unique(setdiff(intersect(ne1,ne2),c(MI_2_left_max[i,1],MI_2_left_max[i,2])))
  if(isempty(nn)==F)
  {
    #3 order
    if(length(nn)>2)
    {
      for(k in 1:(length(nn)-2))
      {
        loc3 <- nn[k]
        for(b in (k+1):(length(nn)-1))
        {
          loc4 <-nn[b]
          for(h in (b+1):length(nn))
          {
            loc5 <- nn[h]
            con <- c(dis_data[,loc3],dis_data[,loc4],dis_data[,loc5])
            mi <- condinformation(dis_data[,loc1],dis_data[,loc2],con,method="emp")
            MI_3 <- rbind(MI_3,c(MI_2_left_max[i,c(1,2)],nn[k],nn[b],nn[h],3,mi))
          }
        }
      }
    }
  }
}
colnames(MI_3) <- c('node1','node2','neighbor1','neighbor2','neighbor3','MI_order','states')
MI_3_value<-as.numeric(MI_3[,7])
median(MI_3_value)
mean(MI_3_value)
hist(MI_3_value,main='MI_3_distribution',breaks =50)

#### 3-order CMI
#max(apply(as.matrix(MI_3[,7]),2,function(x) as.numeric(x)))
thre3<-0.1376554
u3 <- unique(MI_3[,1])
MI_3_left_max <- rbind()
dim(MI_3)
for(i in 1:length(u3))
{
  #i<-3
  loc <- which(MI_3[,1] %in% u3[i])
  k1 <- unique(MI_3[loc,2])
  can <- MI_3[loc,2]
  for(j in 1:length(k1))
  {
    #j <-4
    z1 <- apply(as.matrix(MI_3[loc[which(can %in% k1[j])],dim(MI_3)[2]]),
                2,function(x) as.numeric(x))
    if(max(z1) > thre3)
    {
      MI_3_left_max <- rbind(MI_3_left_max,c(u2[i],k1[j],max(z1)))
    }
  }
}
dim( MI_3_left_max ) #107
dim(MI_3)
MI_3[100:105,]


all_net <- rbind(MI_1_left_max,rbind(MI_2_left_max,MI_3_left_max))
all_net <- rbind(MI_1_left_max,MI_2_left_max)
MIvalue<-as.numeric(all_net[,3])
median(MIvalue)
mean(MIvalue)
hist(MIvalue,main='MI_3_distribution',breaks =50)

thre2 <- 0.4054651
net_left <- rbind()
z1 <- unique(all_net[,1])
for(i in 1:length(z1))
{
  k1 <- which(all_net[,1] == z1[i])
  k2 <- unique(all_net[k1,2])
  can <- all_net[k1,2]
  val <- as.numeric(all_net[k1,3])
  for(j in 1:length(k2))
  {
    zz <- which(can %in% k2[j])
    va <- max(val[zz])
    if(va > thre2)
    {
      net_left <- rbind(net_left,c(z1[i],k2[j],va))
    }
  }
}

dim(unique(net_left))
unique(net_left[,2])
write.csv(MI_1,"2MI_1.csv")
write.csv(MI_2,"2MI_2.csv")
write.csv(MI_3,"2MI_3.csv")
write.csv(net_left,"netleft0.4054651.csv")