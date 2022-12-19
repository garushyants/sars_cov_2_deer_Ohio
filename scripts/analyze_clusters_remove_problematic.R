library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(stringr)
library(reshape2)
library(gridExtra)


#set path
setwd("/panfs/pan1.be-md.ncbi.nlm.nih.gov/virbac/sars_cov_2_deers/20220921_final_vcfs")

folder = "figures_20220929"
###########Read SARS-CoV-2 genomic data from gff
#Read annotation and get data
GFFDf <- read.csv("../vep_data_20220506/Sars_cov_2.ASM985889v3.101.gff3", skip = 8, header = F, sep="\t")
geneDf<-subset(GFFDf,GFFDf$V3 == "gene")
geneDfA<-separate(data = geneDf, col = V9, into = c("A1","gene","A3","A4","A5","A6","A7"), sep = ";")
geneDfA$gene <- gsub('Name=', '', geneDfA$gene)

ymin<-c(0,-1)
ymax<-c(1,0)
geneDfA$textpos<-rep(c(1.1,-1.1),6,length.out=length(geneDfA$V1))
geneDfA$hjust<-rep(c(0.1,1),6,length.out=length(geneDfA$V1))
geneDfA$ymin<-rep(ymin, 6,length.out = length(geneDfA$V1))
geneDfA$ymax<-rep(ymax, 6,length.out = length(geneDfA$V1))


#adding info about parts of ORF1ab
ORF1ab<-read.csv("../vep_data_20220506/NC_045512.2_ORF1ab.tsv", header = F, sep="\t")
ORF1abINFO<-separate(data = ORF1ab, col = V9, into = c("A1","A2","A3","A4","gene","A6"), sep = ";")
ORF1abINFO$gene <- gsub('product=', '', ORF1abINFO$gene)
ORF1abINFO$ymin<-rep(0,length(ORF1abINFO$V1))
ORF1abINFO$ymax<-rep(1,length(ORF1abINFO$V1))
ORF1abINFO$textpos<-rep(1.1,length(ORF1abINFO$V1))
ORF1abINFO$hjust<-rep(0,length(ORF1abINFO$V1))

#merge initial gff and ORF1ab info
colsToKeep<-c("V4","V5","ymin","ymax","gene","textpos","hjust")

SARSCoV2Genes<-rbind(geneDfA[geneDfA$gene != "ORF1ab",colsToKeep],ORF1abINFO[,colsToKeep])
RdRpline<-data.frame(V4=c(13442),V5=c(16236),ymin=c(0),ymax =c(1),gene=c("RdRp"),textpos=c(1.1),hjust=c(0.1))
SARSCoV2GenesCor<-rbind(subset(SARSCoV2Genes, SARSCoV2Genes$gene != "RNA-dependent RNA polymerase"),RdRpline)
SARSCoV2GenesCor$gene<-factor(SARSCoV2GenesCor$gene, levels =SARSCoV2GenesCor[order(SARSCoV2GenesCor$V4),]$gene)
##################################################
##Draw plain SARS-CoV-2 genome
SARSCoVGenome <- ggplot(data = SARSCoV2GenesCor)+
  geom_rect(aes(xmin = V4, xmax = V5, ymin = ymin, ymax = ymax,
                fill = gene),
            color = "black",
            size =0.2)+
  geom_rect(xmin=0,xmax=29903,ymin=-0.05,ymax=0.05)+
  geom_text(aes(x=V4+(V5-V4)/2, 
                y=textpos,
                label = gene,
                hjust =hjust), 
            size =2.5,
            angle = 45)+
  ylim(c(-2.6,+2.7))+
  scale_fill_manual(values =c("#8dd3c7","#8dd3c7","#8dd3c7","#8dd3c7","#8dd3c7",
                              "#8dd3c7","#8dd3c7","#8dd3c7","#8dd3c7","#8dd3c7",
                              "#8dd3c7","#8dd3c7","#8dd3c7","#8dd3c7",
                              "#8dd3c7", "#ffffb3", "#bebada", "#fb8072",
                              "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
                              "#d9d9d9", "#bc80bd", "#ccebc5"))+
  scale_x_continuous(expand = c(0.01, 0), limits = c(0,30000))+
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(color = "white", size =12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_text(color = "white", size =10),
        axis.ticks.y = element_line(color = "white"),
        plot.margin=unit(c(0.2,1,0,1), "cm"))
SARSCoVGenome

########################################################
########################################################
#Read annotated vcfs with cluster mutations
setwd("/panfs/pan1.be-md.ncbi.nlm.nih.gov/virbac/sars_cov_2_deers/20220926_redefined_delta_clusters")

#load cluster names
ClusterNames<-read.csv("Clusters_rename_rule.txt", sep ="\t", header =F)
colnames(ClusterNames)<-c("Cluster", "ClusterNode")

AnnotatedVCF<-read.csv("Delta.clusters.snpeff.redefined_clusters.vcf", skip=8, header =F, sep='\t')
AnnotatedVCFINFOA<-separate(data = AnnotatedVCF, col = V8,
                           into = c("ClusterNode","INFOelse"),
                           sep = ',')
AnnotatedVCFINFOB<-separate(data = AnnotatedVCFINFOA, col = INFOelse,
                           into = c("Node","INFOelse"),
                           sep = ';')
AnnotatedVCFINFOC<-separate(data = AnnotatedVCFINFOB, col = INFOelse,
                            into = c("MutType","A1"),
                            sep ='[(]')
AnnotatedVCFINFOD<-separate(data = AnnotatedVCFINFOC, col = A1,
                            into = c("A1","A2","Codon","Prot_long","A3","gene","Type","A4","A5",
                                     "A6","A7"),
                            sep ='[|]')
AnnotatedVCFINFOE<-separate(data = AnnotatedVCFINFOD, col = Prot_long,
                            into = c("ProtChange","A"),
                            sep ='[/]')

AnnotatedVCFINFO<-AnnotatedVCFINFOE[,c(1:10,13,14,17,18)]
AnnotatedVCFINFO$ClusterNode <- gsub('Cluster=', '', AnnotatedVCFINFO$ClusterNode)
AnnotatedVCFINFO$Node <- gsub('Node=', '', AnnotatedVCFINFO$Node)
AnnotatedVCFINFO$ProtChange <- gsub('p.', '', AnnotatedVCFINFO$ProtChange)
AnnotatedVCFINFO$MutType <- gsub('EFF=', '', AnnotatedVCFINFO$MutType)
AnnotatedVCFINFO$Change<-paste(AnnotatedVCFINFO$V4,">",AnnotatedVCFINFO$V5,
                               sep="")
AnnotatedVCFINFOProperNames<-merge(AnnotatedVCFINFO, ClusterNames, by ="ClusterNode")


Clusterlevel<-data.frame(Cluster=unique(AnnotatedVCFINFOProperNames$Cluster), 
                         level=seq(1,length(unique(AnnotatedVCFINFOProperNames$Cluster))))
AnnotatedVCFINFOLev<-merge(AnnotatedVCFINFOProperNames,Clusterlevel, by="Cluster")



###remove problematic sites
#Check this additional list of Ambiquous sites from here: https://github.com/W-L/ProblematicSites_SARS-CoV2

ProbLematicSitesDf<-read.csv("problematic_sites_sarsCov2.vcf", comment.char = "#", header = F, sep = "\t")
ProbLematicSitesDfMask<-subset(ProbLematicSitesDf, ProbLematicSitesDf$V7 == "mask")

#I am masking problematic sites
AnnotatedVCFINFOLevNoPr<-subset(AnnotatedVCFINFOLev,
                                    !(AnnotatedVCFINFOLev$V2 %in% ProbLematicSitesDfMask$V2))

#########################################
#Save all to file for supplementary table
ToSave<-AnnotatedVCFINFOLevNoPr[,c(1,3:4,6,7,10:14,16)]

write.csv(ToSave,file="Deltas_noproblematic_suppltableN3.csv",row.names = F, quote=F)

#########################################
##
#Check the number of mutations per cluster per gene
MutsPerGene<-AnnotatedVCFINFOLevNoPr %>% group_by(Cluster,gene) %>% count()

#Check selection in S
subset(MutsPerGene, MutsPerGene$gene=='S' & MutsPerGene$n > 5)

# #####I can try to remove sites within ARTIC primers
# Artic<-read.csv("ARTICv3_nCoV-2019.primer.no_alt.bed", header = F, sep = "\t")
# 
# mIntervalsdf<-Artic[,2:3]
# mIntervals<-mIntervalsdf[order(mIntervalsdf$V2),]
# 
# #Subset only positions outside ARTIC v3 primers 
# library(intervals)
# starts = findInterval(AnnotatedVCFINFOLevNoPr$V2, mIntervals$V2)
# ends = findInterval(AnnotatedVCFINFOLevNoPr$V2, mIntervals$V3)
# 
# PosOutsideArtic<-AnnotatedVCFINFOLevNoPr$V2[starts != (ends + 1L)]
# 
# AnnotatedVCFINFOLevNoPrNoArtic<-subset(AnnotatedVCFINFOLevNoPr,
#                                        AnnotatedVCFINFOLevNoPr$V2 %in% PosOutsideArtic)
##

AnnotatedVCFINFOLevNoPrNoSYN<-subset(AnnotatedVCFINFOLevNoPr,AnnotatedVCFINFOLevNoPr$MutType == "missense_variant" |
                                       AnnotatedVCFINFOLevNoPr$MutType == "stop_gained")

###
#Find common variants
MutCounts<-AnnotatedVCFINFOLevNoPrNoSYN %>% group_by(V2,V4,V5) %>%
  summarise(count = length(Node))
MutCountsRecurrent<-subset(MutCounts,MutCounts$count>1)

#################
#plot variants on genome

AlongGenomePlot<-ggplot()+
  geom_vline(data = MutCountsRecurrent,aes(xintercept = V2),
             color = "#9ebcda",
             size=1.5,
             alpha=0.5)+
  geom_jitter(data=AnnotatedVCFINFOLevNoPrNoSYN,
             aes(x = V2,
                  y = level,
                  fill = MutType),
              color = "black",
              shape = 25, 
             height= 0.1,
              size =3,
              alpha =0.6)+
  scale_fill_manual(values = c("#8dd3c7","#ffffb3","#bebada","#fb8072",
                               "#80b1d3","#fdb462","#b3de69","#fccde5",
                               "#d9d9d9","#bc80bd","#ccebc5","#ffed6f"),
                    name = "Mutation type")+
  scale_y_continuous(name ="Cluster",breaks = Clusterlevel$level,labels = Clusterlevel$Cluster)+
  scale_x_continuous(expand= c(0.01,0),limits = c(0,30000),breaks=seq(0,30000,2000))+
  theme_minimal()+
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(size =10))
        #plot.margin=unit(c(0.3,1,0.3,1), "cm"))
AlongGenomePlot

###
#Combine with genome
Figure1<-ggarrange(AlongGenomePlot,
          SARSCoVGenome,
          ncol=1,
          align = "v",
          heights = c(1,0.5))
Figure1

ggsave(plot = Figure1, path=folder,file = "Figure_SNVs_in_clusters_along_genome_no_problematic_missense.svg",
       width = 22, height = 25, dpi = 300, units = "cm")

#########################
#Draw synonymous and missense and see if there are difference

AlongGenomePlotMisSyn<-ggplot(data=AnnotatedVCFINFOLevNoPr)+
  geom_vline(data = MutCountsRecurrent,aes(xintercept = V2),
             color = "#d9d9d9",
             alpha=0.8)+
  geom_point(aes(x = V2,
                 y = level,
                 fill = MutType),
             color = "black",
             shape = 25, 
             size =4,
             alpha =0.6)+
  scale_fill_manual(values = c("#8dd3c7","#ffffb3","#bebada","#fb8072",
                               "#80b1d3","#fdb462","#b3de69","#fccde5",
                               "#d9d9d9","#bc80bd","#ccebc5","#ffed6f"),
                    name = "Mutation type")+
  scale_y_continuous(name ="Cluster",breaks = Clusterlevel$level,labels = Clusterlevel$Cluster)+
  scale_x_continuous(expand= c(0.01,0),limits = c(0,30000),breaks=seq(0,30000,2000))+
  facet_wrap(~MutType, ncol =1)+
  theme_classic()+
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(size =10))
#plot.margin=unit(c(0.3,1,0.3,1), "cm"))
AlongGenomePlotMisSyn

Figure2<-ggarrange(AlongGenomePlotMisSyn,
                   SARSCoVGenome,
                   ncol=1,
                   align = "v",
                   heights = c(1,0.4))
Figure2

ggsave(plot = Figure2, path=folder,file = "Figure_SNVs_in_clusters_along_genome_no_problematic_by_type.svg",
       width = 32, height = 32, dpi = 300, units = "cm")












