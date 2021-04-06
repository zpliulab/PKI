# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
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



# Packages
library(ggplot2)
library(ggthemes)
#library(dplyr)

#修改路径
setwd("D:\\A\\work\\2019\\scrna\\R")
getwd()

#read data   gene-RSquared typeof(rsqdata) list
rsqdata=read.csv("tfsort2.csv",header=T,row.names=1,stringsAsFactors = FALSE)
rsqdata50=rsqdata[1:30,]

#plot1 bar 
plot1<-ggplot(data=rsqdata50,aes(x=reorder(name,num),y=num)) + 
  geom_bar(stat= 'identity',width = 0.7,aes(fill=num))+
  xlab("")+
  ylab("PKI Score")+
  theme_classic2()+
  scale_fill_gradient(high ='#4f86c6',low="#9dc3c1")+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
        panel.background=element_rect(fill="white"))+
  coord_flip()
ggsave(filename = "plot2.png",dpi = 300)



#plot2 最优子集选取过程中参数选择  
para=read.table("parameter.txt",header=T)
varx<-rep(para$α,times=3)
type<- rep(c('subset_size','change',"function"),each=9)
value<- c(para$subset_size,para$change,para$function.)
df<-data.frame(varx=varx,type=type,value=value)
plot2<-ggplot(data = df, mapping = aes(x = varx , y = value, colour = type)) + 
  geom_line(linetype="solid",size=2,)+ 
  geom_point(size=3) +
  xlab("Parameter α")+
  scale_x_continuous(limits=c(0.1,0.9),breaks = seq(0.1,0.9,0.1))+ylab('Value(%)')+
  theme_economist_white()+
  labs(colour = "Type")+
  scale_colour_discrete(h = c(100, 300) + 30, c = 50, l = 65, labels=c("R-squared value", "Function value", "Gene subset size"))


ggsave("three.png",dpi=300)



#######相关系数矩阵  list
rankdata=read.csv("allrank1.csv",header=T,stringsAsFactors = FALSE,row.names = 1)

library(vegan)
library(dplyr)
library(corrplot)
library(extrafont)
library(ggpubr)

M <- cor(rankdata[1:30,])#correlation

########散点图

rankdata=read.csv("rankall_1.csv",header=T,stringsAsFactors = FALSE,row.names = 1)
rankdata50=rankdata[1:50,]
rankdata50=scale(rankdata50)
rankdata50=data.frame(rankdata50)
p1<-ggplot(rankdata50,aes(x=PKI,y=Degree.Centrality))+
  geom_point(color="#006699",alpha=.8)+
  theme_bw()+
  xlab("PKI")+
  ylab('Degree Centrality')+
  geom_smooth(method=lm,color="#006699",alpha=.4)+ 
  stat_cor(method = "pearson")
# p2<-ggplot(rankdata50,aes(x=rankdata50$PKI,y=rankdata50$Lin.Centrality))+
#   geom_point(color="blue")+
#   theme_bw()+
#   xlab("R-squared")+
#   ylab('Lin Centrality')
p3<-ggplot(rankdata50,aes(x=PKI,y=Latora.Closeness.Centrality))+
  geom_point(color="blue",alpha=.8)+
  theme_bw()+
  xlab("R-squared")+
  ylab('Latora Closeness Centrality')+
  geom_smooth(method=lm,color="blue",alpha=.4)+ 
  stat_cor(method = "pearson")
# p3<-ggplot(rankdata50,aes(x=PKI,y=Latora.Closeness.Centrality))+
#   geom_point(color="blue")+
#   theme_bw()+
#   xlab("PKI")+
#   ylab('Latora Closeness Centrality')+
#   geom_smooth(method=lm)+ 
#   stat_cor(method = "pearson")
# p4<-ggplot(rankdata50,aes(x=PKI,y=Geodesic.K.path.Centrality))+
#   geom_point(color="blue")+
#   theme_bw()+
#   xlab("R-squared")+
#   ylab('Geodesic K-path Centrality')
# p5<-ggplot(rankdata50,aes(x=rankdata50$PKI,y=rankdata50$Semi.local.Centrality))+
#   geom_point(color="blue")+
#   theme_bw()+
#   xlab("R-squared")+
#   ylab('Semi local Centrality')+
#   geom_smooth(method=lm)+ 
#   stat_cor(method = "pearson")

p7<-ggplot(rankdata50,aes(x=PKI,y=Betweenness.Centrality))+
  geom_point(color="blue")+
  theme_bw()+
  xlab("R-squared")+
  ylab('Betweenness Centrality')+
  geom_smooth(method=lm,color="darkorange",alpha=.4)+ 
  stat_cor(method = "pearson")
p8<-ggplot(rankdata50,aes(x=rankdata50$PKI,y=rankdata50$Diffusion.Degree))+
  geom_point(color="blue")+
  theme_bw()+
  xlab("R-squared")+
  ylab('Diffusion.Degree')
p9<-ggplot(rankdata50,aes(x=PKI,y=Residual.Closeness.Centrality))+
  geom_point(color="red",alpha=.8)+
  theme_bw()+
  xlab("PKI")+
  ylab('Residual Closeness Centrality')+
  geom_smooth(method=lm,color="red",alpha=.4)+ 
  stat_cor(method = "pearson")
p10<-ggplot(rankdata50,aes(x=rankdata50$PKI,y=rankdata50$Laplacian.Centrality))+
  geom_point(color="red")+
  theme_bw()+
  xlab("R-squared")+
  ylab('Laplacian Centrality')+
  geom_smooth(method=lm,color="red",alpha=.4)+ 
  stat_cor(method = "pearson")
p11<-ggplot(rankdata50,aes(x=PKI,y=PageRank))+
  geom_point(color="#CC6699")+
  theme_bw()+
  xlab("PKI")+
  ylab('Pagerank')+
  geom_smooth(method=lm,color="#CC6699",alpha=.4)+ 
  stat_cor(method = "pearson")
p12<-ggplot(rankdata50,aes(x=PKI,y=LeaderRank))+
  geom_point(color="#CC6699",alpha=.8)+
  theme_bw()+
  xlab("PKI")+
  ylab('LeaderRank')+
  geom_smooth(method=lm,color="#CC6699",alpha=.4)+ 
  stat_cor(method = "pearson")
p13<-ggplot(rankdata50,aes(x=PKI,y=Leverage.Centrality))+
  geom_point(color="darkorange",alpha=.8)+
  theme_bw()+
  xlab("PKI")+
  ylab('Leverage Centrality')+
  geom_smooth(method=lm,color="darkorange",alpha=.4)+ 
  stat_cor(method = "pearson")
p14<-ggplot(rankdata50,aes(x=PKI,y=Bottlenck.Centrality))+
  geom_point(color="darkgreen",alpha=.8)+
  theme_bw()+
  xlab("PKI")+
  ylab('Bottlenck Centrality')+
  geom_smooth(method=lm,color="darkgreen",alpha=.4)+ 
  stat_cor(method = "pearson")

p15<-ggplot(rankdata50,aes(x=PKI,y=Cross.clique.Centrality))+
  geom_point(color="darkorange",alpha=.8)+
  theme_bw()+
  xlab("PKI")+
  ylab('Cross clique Centrality')+
  geom_smooth(method=lm,color="darkorange",alpha=.4)+ 
  stat_cor(method = "pearson")

multiplot(p1,p14,p15,p12, cols=2)

library(RColorBrewer)
####heatmap

rankdata=read.csv("allrank1.csv",header=T,row.names = 1,stringsAsFactors = FALSE)
rankdata=read.csv("rankall_1.csv",header=T,row.names = 1,stringsAsFactors = FALSE)

rankdata30=rankdata[1:50,]
M <- cor(rankdata30,method="spearman")
col3 <- colorRampPalette(c("#ff7473","#fc9d9a","white", "#0fa6ff"))
corrplot(M, method = "number",tl.col="black", cl.lim = c(0, 1))

corrplot(corr = M,method="circle",tl.pos="tp",tl.col="black",cl.lim = c(0, 1))

corrplot(corr = M, order="hclust", addrect = 3, rect.col = "black",tl.col="black")

corrplot(M,add=TRUE,type="lower", method="number", col = "black", cl.lim = c(0, 1),
         diag=FALSE,tl.pos="n",tl.col="black",cl.pos="n")

pdf("heatmap.pdf",paper='a4',family="Times")
corrplot(corr = M,order="AOE",type="upper",tl.pos="tp",tl.col="black")
corrplot(corr = M,add=TRUE, type="lower", method="number",order="AOE", col="black",
         diag=FALSE,tl.pos="n",tl.col="black",cl.pos="n")
dev.off()


M<-as.data.frame(M)
M$sum <- rowSums(M[,0:12])
M$name<-row.names(M)
sort(M$sum)
corrsum<-ggplot(data=M[,13:14],aes(x=reorder(M$name,M$sum),y=M$sum)) + theme_bw()+
  geom_bar(stat= 'identity',fill = 'snow', colour = 'black',width = 0.8)+
  xlab('Methods')+
  ylab("Sums of correlation coefficients")+
  theme(axis.text.x=element_text(size=8),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  coord_flip()


