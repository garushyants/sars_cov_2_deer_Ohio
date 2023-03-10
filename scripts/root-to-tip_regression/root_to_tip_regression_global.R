library(ggplot2)
library(ggpubr)
library(stringr)
library("rstudioapi")
library(DescTools)
path<-getSourceEditorContext()$path
setwd(SplitPath(path)$dirname)


setwd('../../data/global_dataset')
folder<-"../../figures/root-to-tip"
###Read output of root_to_tips_regression.py
RootToTip<-read.csv("Combined.root_to_tips_counts.tsv", header = T, sep='\t',
                    stringsAsFactors = F)

###Read dates
Dates<-read.csv("Combined.all.names_and_dates.csv", header = T,sep = ',')

###get metadata for the samples

Metadata<-read.csv("Combined.all.host.lineage", sep = '\t', header=F)

##get lineages definitions

Lineages<-read.csv("Lineage_renaming.csv",header=F)

MetadataLineages<-merge(x=Metadata,y=Lineages,by.x='V3', by.y ='V1',all.x=T)
#MetadataLineages$V2.y[is.na(MetadataLineages$V2.y)]<-"Other"
MetadataLineages$coloringTag<-paste(MetadataLineages$V2.x,MetadataLineages$V2.y,sep="_")

#replacing NA
rep_str=c("_NA"="")
MetadataLineages$coloringTag <- str_replace_all(MetadataLineages$coloringTag, rep_str)

###Merge dataframes

RootToTipWithDates<-merge(x=RootToTip,y=Dates, by.x='Leaf', by.y='name')

RootToTipWithDatesMeta<-merge(x=RootToTipWithDates, y=MetadataLineages, by.x='Leaf', by.y='V1')


#####
#Save final dataset
ToSave<-RootToTipWithDatesMeta

names(ToSave)<-c("strain","All","Syn","Mis","Ssyn","Smis","date","clade","host","VOC","coloringTag")

write.table(ToSave, file="Combined.all.root_to_tip_from_R.tsv",
          quote=F, sep="\t",
          row.names = F)

#####

excluded_for_regression<-c("deer","deer_Alpha","deer_Delta","deer_Omicron","human_Omicron")
ForRegression<-subset(RootToTipWithDatesMeta, !(RootToTipWithDatesMeta$coloringTag %in% excluded_for_regression))

##colors
mypallete<-c("#4393c3","#08519c","#6a3d9a","#cab2d6",
  "#e3e69e","#f2db96","#fdae61","#f46d43",
  "#e08214")

lmcolor<-"#bdbdbd"


#Plot1 Plot All mutations
lm_eqn <- function(df,x,y){
  m <- lm(y ~ x, df)
  eq <- substitute(paste(a," subst/site/year",sep=""),
                   list(a = formatC(unname(coef(m)[2]/29903), format = "e", digits = 2)))
  as.character(as.expression(eq))
}

PlotAll<-ggplot()+
  geom_point(data=RootToTipWithDatesMeta,
             aes(x=date, y=All, color = coloringTag),
             shape=19)+
  geom_smooth(data=ForRegression,
              aes(x=date, y=All),
              color = lmcolor,
              method='lm')+
  annotate("text",
           label = lm_eqn(ForRegression,ForRegression$date,ForRegression$All),
           x = 2022,
           y = 2,
           color = lmcolor,
           parse=T)+
  scale_x_continuous(name="years", breaks= seq(2020,2023,0.5),limits = c(2019.5,2023.0))+
  scale_y_continuous(name="# of mutations", breaks = seq(0,100,5))+
  scale_color_manual(name ="",values=mypallete)+
  theme_classic2()
PlotAll

#Plot2 Synonymous
PlotSyn<-ggplot()+
  geom_point(data=RootToTipWithDatesMeta,
             aes(x=date, y=Syn, color = coloringTag),
             shape=19)+
  geom_smooth(data=ForRegression,
              aes(x=date, y=Syn),
              color = lmcolor,
              method='lm')+
  scale_x_continuous(name="years", breaks= seq(2020,2023,0.5),limits = c(2019.5,2023.0))+
  scale_y_continuous(name="# of synonymous mutations", breaks = seq(0,100,5))+
  scale_color_manual(name ="",values=mypallete)+
  theme_classic2()
PlotSyn

##Plot missense
PlotMis<-ggplot()+
  geom_point(data=RootToTipWithDatesMeta,
             aes(x=date, y=Mis, color = coloringTag),
             shape=19)+
  geom_smooth(data=ForRegression,
              aes(x=date, y=Mis),
              color = lmcolor,
              method='lm')+
  scale_x_continuous(name="years", breaks= seq(2020,2023,0.5),limits = c(2019.5,2023.0))+
  scale_y_continuous(name="# of missense mutations", breaks = seq(0,100,5))+
  scale_color_manual(name ="",values=mypallete)+
  theme_classic2()
PlotMis

#Plot S syn
PlotSsyn<-ggplot()+
  geom_point(data=RootToTipWithDatesMeta,
             aes(x=date, y=Ssyn, color = coloringTag),
             shape=19)+
  geom_smooth(data=ForRegression,
              aes(x=date, y=Ssyn),
              color = lmcolor,
              method='lm')+
  scale_x_continuous(name="years", breaks= seq(2020,2023,0.5),limits = c(2019.5,2023.0))+
  scale_y_continuous(name="# of synonymous mutations in Spike", breaks = seq(0,100,1), limits = c(0,10))+
  scale_color_manual(name ="",values=mypallete)+
  theme_classic2()
PlotSsyn

#Plot S missense
#Plot S syn
PlotSmis<-ggplot()+
  geom_point(data=RootToTipWithDatesMeta,
             aes(x=date, y=Smis, color = coloringTag),
             shape=19)+
  geom_smooth(data=ForRegression,
              aes(x=date, y=Smis),
              color = lmcolor,
              method='lm')+
  scale_x_continuous(name="years", breaks= seq(2020,2023,0.5),limits = c(2019.5,2023.0))+
  scale_y_continuous(name="# of missense mutations in Spike", breaks = seq(0,100,5))+
  scale_color_manual(name ="",values=mypallete)+
  theme_classic2()
PlotSmis

#Combine together

top<-ggarrange(PlotAll, labels = c("a"))
middle<-ggarrange(PlotSyn,PlotMis,
          PlotSsyn,PlotSmis,
          ncol=2,
          nrow=2,
          labels=c('b','c','d','e'),
          legend = "none")

FinalPlot<-ggarrange(top,
          bottom,
          ncol =2,
          widths = c(0.8,1))
FinalPlot

###Save figure

# ggsave("Figure_root_to_tip_regression.png",path=folder,
#       width=45,height=20,units='cm',
#       dpi=300, plot=FinalPlot)
# 
# ggsave("Figure_root_to_tip_regression.svg",path=folder,
#        width=45,height=20,units='cm',
#        dpi=300, plot=FinalPlot)


###figure for main text
MainPlot<-ggarrange(PlotSyn,PlotMis,
                  ncol=2,
                  labels=c('a','b'),
                  common.legend = T,
                  legend = "right")

ggsave("Figure_Main_root_to_tip_regression.svg",path=folder,
       width=30,height=11,units='cm',
       dpi=300, plot=MainPlot)

###############################
###Do calculations for all genes

M<-read.csv("Combined_dataset.root_to_tips_counts.M.tsv", header = T, sep='\t',
                    stringsAsFactors = F)
N<-read.csv("Combined_dataset.root_to_tips_counts.N.tsv", header = T, sep='\t',
            stringsAsFactors = F)
ORF1ab<-read.csv("Combined_dataset.root_to_tips_counts.ORF1ab.tsv", header = T, sep='\t',
            stringsAsFactors = F)
ORF10<-read.csv("Combined_dataset.root_to_tips_counts.ORF10.tsv", header = T, sep='\t',
            stringsAsFactors = F)
ORF3a<-read.csv("Combined_dataset.root_to_tips_counts.ORF3a.tsv", header = T, sep='\t',
            stringsAsFactors = F)
ORF6<-read.csv("Combined_dataset.root_to_tips_counts.ORF6.tsv", header = T, sep='\t',
            stringsAsFactors = F)
ORF7a<-read.csv("Combined_dataset.root_to_tips_counts.ORF7a.tsv", header = T, sep='\t',
            stringsAsFactors = F)
ORF7b<-read.csv("Combined_dataset.root_to_tips_counts.ORF7b.tsv", header = T, sep='\t',
            stringsAsFactors = F)
ORF8<-read.csv("Combined_dataset.root_to_tips_counts.ORF8.tsv", header = T, sep='\t',
                stringsAsFactors = F)
#Merge them
keeps<-c(1,5,6)
MN<-merge(M[,keeps], N[,keeps], by = "Leaf")
MN1ab<-merge(MN, ORF1ab[,keeps], by = "Leaf")
MN1ab10<-merge(MN1ab, ORF10[,keeps], by = "Leaf")
MN1ab103<-merge(MN1ab10, ORF3a[,keeps], by = "Leaf")
MN1ab1036<-merge(MN1ab103, ORF6[,keeps], by = "Leaf")
MN1ab10367a<-merge(MN1ab1036, ORF7a[,keeps], by = "Leaf")
MN1ab10367<-merge(MN1ab10367a, ORF7b[,keeps], by = "Leaf")
MN1ab103678<-merge(MN1ab10367, ORF8[,keeps], by = "Leaf")
#Add dates and colors
GenesDates<-merge(x=MN1ab103678,y=Dates, by.x='Leaf', by.y='name')
GenesDatesMeta<-merge(x=GenesDates, y=MetadataLineages, by.x='Leaf', by.y='V1')

#For regression
GenesForRegression<-subset(GenesDatesMeta, !(GenesDatesMeta$coloringTag %in% excluded_for_regression))

#########
###plot

#functions
PlotByGene<-function(gene,syn,mis,laba="a",labb="b"){
  #Plot gene syn
  # gene<-"N"
  # syn<"Nsyn"
  # mis<-"Nmis"
  PlotGenesyn<-ggplot()+
    geom_point(data=GenesDatesMeta,
               aes_string(x="date", y=syn, color = "coloringTag"),
               shape=19)+
    geom_smooth(data=GenesForRegression,
                aes_string(x="date", y=syn),
                color = lmcolor,
                method='lm')+
    scale_x_continuous(name="years", breaks= seq(2020,2023,0.5),limits = c(2019.5,2023.0))+
    scale_y_continuous(name=paste("# of synonymous mutations in ",gene,sep=""), breaks = seq(0,100,1))+
    scale_color_manual(name ="",values=mypallete)+
    theme_classic2()
  #PlotGenesyn
  
  #Plot gene missense
  PlotGenemis<-ggplot()+
    geom_point(data=GenesDatesMeta,
               aes_string(x="date", y=mis, color = "coloringTag"),
               shape=19)+
    geom_smooth(data=GenesForRegression,
                aes_string(x="date", y=mis),
                color = lmcolor,
                method='lm')+
    scale_x_continuous(name="years", breaks= seq(2020,2023,0.5),limits = c(2019.5,2023.0))+
    scale_y_continuous(name=paste("# of missense mutations in ",gene,sep=""), breaks = seq(0,500,2))+
    scale_color_manual(name ="",values=mypallete)+
    theme_classic2()
  #PlotGenemis
  
  GenePlot<-ggarrange(PlotGenesyn,
                      PlotGenemis,
            ncol =2,
            labels = c(laba,labb),
            common.legend = T,
            legend = "right")

  return(GenePlot)
}

#Final plots
NPlot<-PlotByGene("N","Nsyn","Nmis")
NPlot
ggsave("N_mutations_global_dataset.png",path=folder,
       plot=NPlot,
       width=30,height=15,units='cm',
       dpi=300)
MPlot<-PlotByGene("M","Msyn","Mmis")
MPlot
ggsave("M_mutations_global_dataset.png",path=folder,
       plot=MPlot,
       width=30,height=15,units='cm',
       dpi=300)
ORF10Plot<-PlotByGene("ORF10","ORF10syn","ORF10mis")
ORF10Plot
ORF1abPlot<-PlotByGene("ORF1ab","ORF1absyn","ORF1abmis","f","g")
ORF1abPlot

FinalPlotORF1ab<-ggarrange(FinalPlot,
                           ORF1abPlot,
                           heights=c(1,0.7),
                           ncol=1)

ggsave("Figure_root_to_tip_regression_withORF1ab.png",path=folder,
       width=45,height=35,units='cm',
       dpi=300, plot=FinalPlotORF1ab)

ggsave("Figure_root_to_tip_regression_withORF1ab.svg",path=folder,
       width=45,height=35,units='cm',
       dpi=300, plot=FinalPlotORF1ab)
ORF3aPlot<-PlotByGene("ORF3a","ORF3asyn","ORF3amis")
ORF3aPlot
ggsave("ORF3a_mutations_global_dataset.png",path=folder,
       plot=ORF3aPlot,
       width=30,height=15,units='cm',
       dpi=300)
ORF6Plot<-PlotByGene("ORF6","ORF6syn","ORF6mis")
ORF6Plot
ggsave("ORF6_mutations_global_dataset.png",path=folder,
       plot=ORF6Plot,
       width=30,height=15,units='cm',
       dpi=300)
ORF7aPlot<-PlotByGene("ORF7a","ORF7asyn","ORF7amis")
ORF7aPlot
ggsave("ORF7a_mutations_global_dataset.png",path=folder,
       plot=ORF7aPlot,
       width=30,height=15,units='cm',
       dpi=300)
ORF7bPlot<-PlotByGene("ORF7b","ORF7bsyn","ORF7bmis")
ORF7bPlot
ggsave("ORF7b_mutations_global_dataset.png",path=folder,
       plot=ORF7bPlot,
       width=30,height=15,units='cm',
       dpi=300)
ORF8Plot<-PlotByGene("ORF8","ORF8syn","ORF8mis")
ORF8Plot
ggsave("ORF8_mutations_global_dataset.png",path=folder,
       plot=ORF8Plot,
       width=30,height=15,units='cm',
       dpi=300)








