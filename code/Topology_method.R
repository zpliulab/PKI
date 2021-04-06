setwd("d:\\A\\work\\2019\\scrna\\R")

getwd()
relationship<- read.csv("rel_f1.csv",header=T,row.names = 1)
sixpointdata=data=read.csv("sixpointdata.csv",row.names = 1)

g2 <- graph_from_data_frame(relationship, directed=F)
plot(g2,layout=layout.fruchterman.reingold)

#weights   based on MI
library(KRIS)
library(entropy)
library(igraph)
library(linkcomm)
weightslist=list()
for (i in 1:length(relationship[,1])) {
  y2d = rbind( sixpointdata[relationship[i,1],], sixpointdata[relationship[i,2],])
  mi=mi.empirical(y2d, unit=c("log2"))
  print(mi)
  weightslist=append(weightslist,mi)
}

#weights based on pearson
pweightslist=list()
for (i in 1:length(relationship[,1])) {
  
  pr=cor(t(sixpointdata[relationship[i,1],]), t(sixpointdata[relationship[i,2],]))
  print(pr)
  pweightslist=append(pweightslist,abs(pr))
}


d <-  degree(g2)

cc <- transitivity(g2, type="local",weights = pweightslist)

bc <-  betweenness(g2,normalized = T,weights = pweightslist)

library(centiserve)

Lincent<-lincent(g2,weights = pweightslist) 大

Closeness_residual<-closeness.residual(g2,weights = pweightslist) 

Closeness_latora<-closeness.latora(g2,weights = pweightslist) 

Laplacian<-laplacian(g2) 

Geokpath<-geokpath(g2,weights=pweightslist)

Semilocal<-semilocal(g2)

Leverage<-leverage(g2)

Bottlenck<-bottleneck(g2)

Radiality<-radiality(g2,weights=pweightslist) 

Crossclique<- crossclique(g2)

Diffusion_degree<-diffusion.degree(g2) 

Communitycen<-communitycent(g2)

g3 <- graph_from_data_frame(relationship, directed=T)

LeaderRank<-leaderrank(g3, vids = V(g3))

Pagerank<-page.rank(g3,weights = pweightslist)

rankall=data.frame(d,cc,bc,Lincent,Closeness_latora,Closeness_residual,
                   Laplacian,Geokpath,Semilocal,Diffusion_degree,LeaderRank,
                   Pagerank$vector)
addrank=data.frame(Bottlenck,Leverage)

write.csv(addrank,"addrank.csv")

##################
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
multiplot(plot_list[[1]],plot_list[[2]],plot_list[[3]],cols = 3)


###################################

library(ggplot2)
library(reshape2)
data<- read.csv("sixpointdata.csv",header=T,row.names = 1)
data<- data[10:20,]
data$ID <- rownames(data)
data_m <- melt(data, id.vars=c("ID"))
head(data_m)
p <- ggplot(data_m, aes(x=variable,y=ID)) + 
  geom_tile(aes(fill=value))+ theme_bw()+
  scale_fill_gradient(low = "blue", high = "red")+
  scale_x_continuous(breaks = 1:6, labels = paste("t", 1:6, sep = "_"))

p <- ggplot(data_m, aes(x=variable,y=ID)) + 
  geom_tile(aes(fill = value)) +
  
  scale_fill_continuous(low = "#ED9282", high = "#30AEDE") +
  
  
  theme_bw() +
  
  theme(
    
    axis.ticks = element_blank(),
    
    panel.grid = element_blank(),
    
    panel.border = element_blank(),
    
    axis.title = element_blank(),
    
    axis.text.x = element_text(angle = 90)
    
  )

ggsave(p, filename="heatmap.pdf", width=10, height=15, units=c("cm"),colormodel="srgb")





#####boxplot
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
data=read.csv("expressiondata.csv",header = F)
data=t(data)
colnames(data)=data[1,]
data=data[-1,]
data=data.frame(data)
comb4=c('LEF1','POU5F1',	'GATA3',	'SOX2',"Timepoint")
dat=data[,comb4]
dat2 = gather(dat,key = "gene",value = "expression",-Timepoint) ###画箱线图的数据
p <- ggboxplot(dat2, x = "gene", y = "expression",
               color = "Timepoint", palette = "jco", 
               add = "jitter")


sixpointdata=data=read.csv("sixpointdata.csv",row.names = 1)
timedata=t(sixpointdata)
comb41=c('LEF1','POU5F1',	'GATA3',	'SOX2')
comb41=c('LEF1','SMAD2'	,'GATA3','SOX2')
comb41=c('LEF1'	,'POU5F1',	'NANOG','SOX2')
comb41data=timedata[,comb41]
comb41data=data.frame(comb41data)
write.csv(comb41data,"comb41data.csv")
comb41data["Time"]=c("0h","12h","24h","36h","72h","96h")
comb41data2 = gather(comb41data,key = "gene",value = "expression",-Time) 
p=ggplot(data = comb41data2, mapping = aes(x = Time, y = expression, colour = gene,group =gene)) +
  geom_line()+geom_point()


combination4=read.table("combination2.txt")
plot_list=list()
for (i in 1:8) {
  comb4=combination4[i,]
  comb4data=timedata[,as.character(comb4)]
  #print(comb4data)
  comb4data=data.frame(comb4data)
  
  comb4data["Time"]=c("0h","12h","24h","36h","72h","96h")
  comb4data2 = gather(comb4data,key = "Gene",value = "Expression",-Time) 
  
  p=ggplot(data = comb4data2, mapping = aes(x = Time, y = Expression, colour = Gene,group =Gene)) +
    geom_line()+geom_point()
  plot_list[[i]] = p
  
}


