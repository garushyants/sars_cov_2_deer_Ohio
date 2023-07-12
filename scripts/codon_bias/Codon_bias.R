library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(stringr)
library(reshape2)
library(gridExtra)
library("rstudioapi")

path<-getSourceEditorContext()$path
setwd(dirname(path))
setwd("../../data/")

folder = "../figures/codon_bias"

########################################################
#Read annotated vcfs with cluster mutations
AnnotatedVCF<-read.csv("Delta/Delta.clusters.snpeff.corrected.vcf", skip=8, header =F, sep='\t')
#AnnotatedVCFAlpha<-read.csv("Alpha/Alpha.clusters.snpeff.corrected.vcf", skip=8, header =F, sep='\t')

#parse vcf
AnnotatedVCFINFOA<-separate(data = AnnotatedVCF, col = V8,
                           into = c("Cluster","INFOelse"),
                           sep = ',')
AnnotatedVCFINFOB<-separate(data = AnnotatedVCFINFOA, col = INFOelse,
                           into = c("Node","INFOelse"),
                           sep = ';')
AnnotatedVCFINFOC<-separate(data = AnnotatedVCFINFOB, col = INFOelse,
                            into = c("A1","Effect","Codon","ProtChange","A2","gene","A3","Type","A4","A5","A6"),
                            sep ='[|]')

CodonDf<-separate(data = AnnotatedVCFINFOC, col = Codon, into =c("CodonBefore","CodonAfter"), sep = "/")

AnnotatedVCFINFO<-CodonDf[,c(1:9,11:13,14,16,18)]
AnnotatedVCFINFO$Cluster <- gsub('Cluster=', '', AnnotatedVCFINFO$Cluster)
AnnotatedVCFINFO$Node <- gsub('Node=', '', AnnotatedVCFINFO$Node)
AnnotatedVCFINFO$Change<-paste(AnnotatedVCFINFO$V4,">",AnnotatedVCFINFO$V5,
                               sep="")

#############

Clusterlevel<-data.frame(Cluster=unique(AnnotatedVCFINFO$Cluster), level=seq(1,length(unique(AnnotatedVCFINFO$Cluster))))
AnnotatedVCFINFOLevNoPr<-merge(AnnotatedVCFINFO,Clusterlevel, by="Cluster")

#########################################
#Select synonymous only
SynOnly<-subset(AnnotatedVCFINFOLevNoPr,
                    AnnotatedVCFINFOLevNoPr$Effect == "SILENT" &
                  AnnotatedVCFINFOLevNoPr$Type == "CODING")

SynOnly$CodonAUp<-toupper(SynOnly$CodonAfter)
SynOnly$CodonBUp<-toupper(SynOnly$CodonBefore)

CodonUsageSyn<-merge(as.data.frame(table(SynOnly$CodonAUp)),
      as.data.frame(table(SynOnly$CodonBUp)),
      by="Var1",all=T)
names(CodonUsageSyn)<-c("Codon","After","Before")
#########################################

######################################################
#Upload deer and human data
CodonUsageDeerHuman<-read.csv("codon_bias/Codon_usage_human_deer.mod.csv", header = T)

CodonUsageDeerHumanRD<-CodonUsageDeerHuman%>% group_by(AA) %>% 
  mutate(RSCU_Deer=Odocoileus.virginianus.texanus*length(Codon)/sum(Odocoileus.virginianus.texanus)) %>%
  ungroup()
CodonUsageDeerHumanR<-CodonUsageDeerHumanRD%>% group_by(AA) %>% 
  mutate(RSCU_Human=Homo.sapiens*length(Codon)/sum(Homo.sapiens)) %>%
  ungroup()

LargeDf<-merge(subset(CodonUsageDeerHumanR, CodonUsageDeerHumanR$AA !='STOP'),
               CodonUsageSyn, by="Codon", all=T)


RSCUDf<-melt(LargeDf[,c('Codon','RSCU_Human','RSCU_Deer')])
P1<-ggplot(data=RSCUDf,
       aes(y=Codon, x= variable,fill=value))+
  geom_tile(size=0.2, color ="white")+
  scale_fill_distiller(name = "RSCU", palette = "RdPu") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90),
        axis.title.x = element_blank())
P1
###
RSCUChangeA<-merge(SynOnly, CodonUsageDeerHumanR, by.x="CodonAUp", by.y ="Codon")[,c(1,22,23)]
names(RSCUChangeA)<-c("Codon","AA","RSCU_Deer")
RSCUChangeA$Dataset<-rep("after",length(RSCUChangeA$Codon))
RSCUChangeB<-merge(SynOnly, CodonUsageDeerHumanR, by.x="CodonBUp", by.y ="Codon")[,c(1,22,23)]
names(RSCUChangeB)<-c("Codon","AA","RSCU_Deer")
RSCUChangeB$Dataset<-rep("before",length(RSCUChangeB$Codon))
RSCUChangeToPlot<-rbind(RSCUChangeA,RSCUChangeB)
RSCUChangeToPlot$Dataset<-factor(RSCUChangeToPlot$Dataset, levels = c("before","after"))


P2<-ggplot(data=RSCUChangeToPlot,
       aes(x=Dataset, y=RSCU_Deer, fill = Dataset))+
  geom_boxplot() +
  geom_jitter(color="black", size=1.2, alpha=0.5,
              width=0.2) +
  theme_minimal()+
  stat_compare_means(label.x=1.2, label.y =2.2, size =3)
P2
P12<-ggarrange(P1,P2, ncol=2,
               labels = c("a","b"),
          widths = c(0.4,0.6))
P12
ggsave("Codon_usage.png", plot=P12, path=folder,
       width =20,height=20, dpi =400, units="cm")



sum(LargeDf$RSCU_Deer*LargeDf$After, na.rm = T)
sum(LargeDf$RSCU_Deer*LargeDf$Before, na.rm = T)
















