rm(list = ls())
library(survminer)
library(tidyverse)
load('Tumor-data.rda')
load('Normal-data.rda')

source('Perturbation.R') ## from https://github.com/Marscolono/SSPGI

banet <- read.csv(file="Background.network.txt", head=F, sep="\t") #background network data
pathway <- load(file="Reactom.Rdata") #pathway data from reactom

normal_ee <- tt #expression matrix of normal samples
tumor_ee <- nn #expression matrix of cancer samples

normal_ee2 <- cbind(normal_ee, apply(normal_ee, 1, mean))

rank_normal <- rank.matrix(normal_ee2)
rank_tumor <- rank.matrix(tumor_ee)

x = cbind(rank_normal, rank_tumor)
dim(x)

n_normal = 308
n_cancer = 2167
deltarank.result <- delta.rank(net, x, n_normal, n_cancer)
save(deltarank.result, file="deltarank.result.Rdata")

### caculate the intersection-perturbation matrix
IP_normal <- IP (deltarank.result[[2]], deltarank.result[[3]])
IP_cancer <- IP (deltarank.result[[2]], deltarank.result[[4]])  

save(IP_normal, file="IP_normal.Rdata")
save(IP_cancer, file="IP_cancer.Rdata")

### select features (Sd+ kw test)
data <- t(cbind(IP_cancer, IP_normal))
group <- as.factor(c(rep("cancer", 2167), rep("normal", 308)))

kw.test <- function(x){
  data<- as.data.frame(cbind(x, group))
  colnames(data) <- c("value", "group")
  p<- kruskal.test(value ~ group, data =data)$p.value
  return(p)
}
P.value <- apply(data, 2, kw.test)
f1 <- order(P.value)[1:30000]

sd = apply(EPm_cancer, 1, sd)   
f2 =rev(order(sd))[1:30000]

fea.loc <- intersect(f1, f2) #location of feature
feature <- EPm_cancer[fea.loc,] #feature matrix used for clustering
save(feature,file = 'feature.rda')

####################################################consensus cluster
library(ConsensusClusterPlus)
dir.create('ConsensusCluster/')
results = ConsensusClusterPlus(d = feature,
                               maxK=10,
                               reps=1000,
                               pItem=0.8,
                               pFeature=1,
                               title='ConsensusCluster/',
                               innerLinkage="complete",
                               finalLinkage="complete",
                               clusterAlg="pam",
                               distance='spearman', 
                               seed=123456,
                               plot="pdf")

icl <- calcICL(results,title = 'ConsensusCluster/',plot = 'pdf')

Kvec = 2:10
x1 = 0.1; x2 = 0.9
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") 
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
optK

PAC <- as.data.frame(PAC)
PAC$K <- 2:10
PAC$tt <- 'B'
PAC$tt[5] <- 'A'
library(ggplot2)
library(ggthemr)
ggthemr('flat dark')
ggplot(PAC,aes(factor(K),PAC,group=1))+
  geom_line(size=0.8,color='white')+
  geom_point(size=4,shape=21,stroke = 1.5,aes(fill=tt),color='white')+
  scale_fill_npg()+
  labs(y='Proportion of ambiguous clustering',x='Cluster number K')+
  annotate('segment',x=3.4,xend = 4.6,y=0.21,yend = 0.2,
           size=1.5,arrow=arrow(),alpha=0.8,color='#EEDF7C')+
  theme(axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 14,colour = 'black'),
        plot.title = element_text(size=15,hjust=0.5),
        legend.position = 'none')
ggthemr_reset()

clusterNum=6
cluster=results[[clusterNum]][["consensusClass"]]

sub <- data.frame(Sample=names(cluster),Cluster=cluster)
sub$Cluster <- paste0('GINS',sub$Cluster)
table(sub$Cluster)

head(sub)
cc <- sub$Cluster
names(cc) <- sub$Sample
cc2 <- sort(cc)

my <- results[[6]][["ml"]]
rownames(my) <- sub$Sample
colnames(my) <- sub$Sample
my2 <- my[names(cc2),names(cc2)]

library(pheatmap)
library(ggsci)
col <- c('#B0997F','#93a8ac','#ffc857','#61a5c2','#119da4','#FF6666')
pheatmap(1-my2,show_colnames = F,show_rownames = F,
         #treeheight_row = 20,treeheight_col = 20,
         cluster_rows = F,cluster_cols = F,
         clustering_method = 'complete',
         color = colorRampPalette(c("#C75D30","white"))(50),
         annotation_names_row = F,annotation_names_col = F,
         annotation_row = data.frame(Subtype=cc2),
         annotation_col = data.frame(Subtype=cc2),
         annotation_colors = list(Subtype=c('GINS1'=pal_npg()(10)[1],'GINS2'=pal_npg()(10)[2],
                                            'GINS3'=pal_npg()(10)[3],'GINS4'=pal_npg()(10)[4],
                                            'GINS5'=pal_npg()(10)[7],'GINS6'=pal_npg()(10)[6])))
graph2pdf(file='pheat-cluster.pdf',width=5.5,height=4.5)

ss <- merge(sub,tt,by.x=1,by.y=2)
table(ss$Cluster)
table(ss$Cluster,ss$GEO)

library(survival)
fit <- survfit(Surv(OS.time,OS)~Cluster,ss)
ggsurvplot(fit,pval = T)
fit

fit2 <- survfit(Surv(RFS.time,RFS)~Cluster,ss)
ggsurvplot(fit2,pval = T)
fit2
table(ss$Cluster)
table(ss$Cluster,ss$GEO)
sub3 <- ss[,1:2]

customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}
set_theme <- function() {
  theme_set(theme_bw(base_rect_size = 1.5)+
              theme(panel.background = element_rect(fill = "#f0f8ff", color = NA),
                    panel.grid.minor  = element_blank(),
                    panel.grid.major.x = element_blank(),
                    panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
                    axis.line = element_line(color = "#606F7B"),
                    legend.position = 'none'))
}

ggthemr::ggthemr_reset()
surv_pvalue(fit)
pp <- ggsurvplot(fit,
           palette = col,
           conf.int = FALSE,
           size =1.3,#线条粗细
           pval = T,
           pval.method = T,
           #test.for.trend=T,
           #pval.coord = c(0.8,0.8), #numeric vector, of length 2, specifying the x and y coordinates of the p-value. Default values are NULL.
           #pval.method.coord = (0.8,0.8), #the same as pval.coord but for displaying log.rank.weights name
           legend.labs=paste0('GINS',1:6), 
           legend.title="",
           legend='none',
           font.main = 15,
           xlab= "Time in years",
           ylab=NULL,
           risk.table=TRUE,
           risk.table.pos = 'out',
           tables.col = "strata",
           break.time.by = 2,
           risk.table.title="Number at risk",
           risk.table.height=.32,
           risk.table.y.text.col = T, ## risk table左侧是否用对应分层变量的颜色注释
           risk.table.y.text = T,
           risk.table.y.title = F,
           ncensor.plot = F,
           ggtheme = set_theme())
pp
dev.off()
pp$plot <- customize_labels(pp$plot,
                            font.x        = c(14, "bold", "darkred"), 
                            font.y        = c(14, "bold", "darkred"),
                            font.xtickslab = c(12, "plain", "black"),
                            font.ytickslab = c(12, "plain",'black'))
pp$plot <- pp$plot+
  labs(title = 'Overall survival')+
  theme(plot.title = element_text(face = "bold",colour = "darkred",size = 18,hjust = 0.5))
pp$table <- customize_labels(pp$table,
                             font.title  = c(14, "bold", "darkgreen"), 
                             font.x        = c(14, "bold", "darkred"), 
                             font.y        = c(14, "bold", "darkred"),
                             font.xtickslab = c(12, "plain", "black"),
                             font.ytickslab = c(12, "bold"))
pp
library(export)
graph2pdf(file='Subtype-OS-KM.pdf',width=7,height=6.8)


pp <- ggsurvplot(fit2,
                 palette = col,
                 conf.int = FALSE,
                 size =1.3,#线条粗细
                 pval = T,
                 pval.method = T,
                 #pval.coord = c(0.8,0.8), #numeric vector, of length 2, specifying the x and y coordinates of the p-value. Default values are NULL.
                 #pval.method.coord = (0.8,0.8), #the same as pval.coord but for displaying log.rank.weights name
                 legend.labs=paste0('GINS',1:6), 
                 legend='none',
                 legend.title="",
                 xlab="Time in years",
                 ylab=NULL,
                 risk.table=TRUE,
                 risk.table.pos = 'out',
                 tables.col = "strata",
                 break.time.by = 2,
                 risk.table.title="Number at risk",
                 risk.table.height=.32,
                 risk.table.y.text.col = T, ## risk table左侧是否用对应分层变量的颜色注释
                 risk.table.y.text = T,
                 risk.table.y.title = F,
                 ncensor.plot = F,
                 ggtheme = set_theme())
pp
dev.off()
pp$plot <- customize_labels(pp$plot,
                            font.x        = c(14, "bold", "darkred"), 
                            font.y        = c(14, "bold", "darkred"),
                            font.xtickslab = c(12, "plain", "black"),
                            font.ytickslab = c(12, "plain",'black'))
pp$plot <- pp$plot+
  labs(title = 'Relapse-free survival')+
  theme(plot.title = element_text(face = "bold",colour = "darkred",size = 18,hjust = 0.5))

pp$table <- customize_labels(pp$table,
                             font.title  = c(14, "bold", "darkgreen"), 
                             font.x        = c(14, "bold", "darkred"), 
                             font.y        = c(14, "bold", "darkred"),
                             font.xtickslab = c(12, "plain", "black"),
                             font.ytickslab = c(12, "bold"))
pp
graph2pdf(file='Subtype-RFS-KM.pdf',width=7,height=6.8)

# -------------------------------------------------------------------------

ibrary(umap)
umapdata <- umap(t(feature), n_components = 2, random_state = 15,k=100)

kk <- as.data.frame(umapdata$layout)
kk$ID <- colnames(feature)
kk <- merge(kk,ss[,1:2],by.x=3,by.y=1)

ggplot(kk,aes(V1,V2,color=Subtype,fill=Subtype))+
  scale_fill_manual(values = col)+
  #scale_color_manual(values = col)+
  geom_point(size=2.5,shape=21,color='grey20')+
  labs(x='UMAP1',y='UMAP2')+
  theme_bw(base_rect_size = 1.5)+
  theme(axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 14,colour = 'darkred',face='bold'),
        plot.title = element_text(size=15,hjust=0.5),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.position = 'none')+
  annotate('text',x = -2.4,y = 2.05,label='GINS1',fontface='bold.italic')+
  annotate('text',x = 2.6,y = -1.1,label='GINS2',fontface='bold.italic')+
  annotate('text',x = 3,y = 2,label='GINS5',fontface='bold.italic')+
  annotate('text',x = 2.1,y = -2.3,label='GINS4',fontface='bold.italic')+
  annotate('text',x = -2.2,y = -1.3,label='GINS3',fontface='bold.italic')+
  annotate('text',x = -1.4,y = -2.7,label='GINS6',fontface='bold.italic')
ggsave(filename = 'UMAP.pdf',width = 5,height = 4.9)

# -------------------------------------------------------------------------

library(citccmst)
load(list.files(system.file("extdata", package="citccmst"), full.names=TRUE))
citvalid.exp.annot <- data.frame(id=rownames(citvalid.exp.norm), stringsAsFactors=FALSE, row.names = rownames(citvalid.exp.norm))
ann <- readxl::read_xls('pmed.1001453.s013.xls',skip = 1)[,2:3]
ann$`Probe Set ID`%in%rownames(citvalid.exp.norm)

citvalid.citccmst <- cit.assignCcmst(data=tumor_ee,
                                     data.annot=citvalid.exp.annot,
                                     data.colId="id",
                                     data.colMap="id" ,
                                     citccmst.colMap="Probe.Set.ID",
                                     dist.method="dqda",
                                     plot=T)
table(citvalid.citccmst$citccmst)


















