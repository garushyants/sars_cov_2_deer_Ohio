library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(gridExtra)
library(stringr)
library(reshape2)
library(gridExtra)
library(zoo)
#library(grid)
library(DescTools)

path<-getSourceEditorContext()$path
setwd(dirname(path))
setwd("../../data/lineages_frequencies")

########################################################
#Rule to substitute AA
aacodefile<-read.csv("../genetic_code.txt", header = F, sep =" ")
aacode<-unique(aacodefile[,c(2,3)])

###Function to parse vcf and find recurrent events in clusters
parseVCF<-function(filename)
{
  #filename<-"Delta.clusters.snpeff.corrected.vcf"
  AnnotatedVCF<-read.csv(filename, header =F, skip=8, sep='\t')
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
  
  #AnnotatedVCFINFO<-AnnotatedVCFINFOE[,c(1:10,13,14,17,18)]
  AnnotatedVCFINFOE$ClusterNode <- gsub('Cluster=', '', AnnotatedVCFINFOE$ClusterNode)
  AnnotatedVCFINFOE$Node <- gsub('Node=', '', AnnotatedVCFINFOE$Node)
  AnnotatedVCFINFOE$ProtChangeLong <- gsub('p.', '', AnnotatedVCFINFOE$ProtChange, fixed=T)
  #generating short values
  REFAA<-substring(AnnotatedVCFINFOE$ProtChangeLong,1,3)
  REFAAshort<-as.character(aacode$V2[match(REFAA, aacode$V3)])
  ALTAA<-substring(AnnotatedVCFINFOE$ProtChangeLong,nchar(AnnotatedVCFINFOE$ProtChangeLong)-2,nchar(AnnotatedVCFINFOE$ProtChangeLong))
  ALTAAshort<-as.character(aacode$V2[match(ALTAA, aacode$V3)])
  POS<-substring(AnnotatedVCFINFOE$ProtChangeLong,4,nchar(AnnotatedVCFINFOE$ProtChangeLong)-3)
  #merge
  AnnotatedVCFINFOE$ProtChange<-paste(AnnotatedVCFINFOE$gene,":",REFAAshort,POS,ALTAAshort,sep="")
  
  AnnotatedVCFINFOE$MutType <- gsub('EFF=', '', AnnotatedVCFINFOE$MutType)
  INFOLong<-AnnotatedVCFINFOE[,c("V2","V4","V5","MutType","ProtChangeLong","ProtChange", "gene", "Type", "A2")]
  INFOShort<-INFOLong %>% 
    group_by(V2,V4,V5,MutType,ProtChangeLong,ProtChange,gene,Type) %>% summarise(events=length(A2))
  INFOShortMissense<-subset(INFOShort,(INFOShort$MutType != "synonymous_variant" &
                                         INFOShort$Type == "protein_coding"))
  return(INFOShortMissense)
}
##############
#Get tables to analyze for Alpha and Delta clusters
DeltaVCF<-parseVCF("Delta.clusters.snpeff.corrected.vcf")
AlphaVCF<-parseVCF("Alpha.clusters.snpeff.corrected.vcf")

##############
#Read data on mutations counts from alignment (only for clusters)
getVariantsFromAli<-function(filename)
{
  FromAli<-read.csv(filename, header = F, sep = " ")
  #aggregate data
  SampleIDsAggr<-aggregate(V4 ~ V1+V2+V3, data = FromAli, paste, collapse = ",")
  LineagesAggr <- FromAli%>%
    group_by(V1,V2,V3) %>%
    summarise(lineage = paste(unique(V5), collapse = ","))
  FromAliCounts <- FromAli%>%
    group_by(V1,V2,V3) %>%
    summarise(count=length(V4))
  #merge aggregated
  df_list <- list(SampleIDsAggr, FromAliCounts, LineagesAggr)
  FromAliAggregated<-Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  return(FromAliAggregated)
}

DeltaFromAliAggregated<-getVariantsFromAli("Delta.variants_from_ali.strict.only_clusters.tsv")
AlphaFromAliAggregated<-getVariantsFromAli("Alpha.variants_from_ali.strict.only_clusters.tsv")

#I read data about singletons because I will require it later on
DeltaFromAliSingletons<-getVariantsFromAli("Delta.variants_from_ali.strict.only_singletons.tsv")

######################
#merge data from alignment with VCF data
DeltaVariantsInfo<-merge(DeltaVCF, DeltaFromAliAggregated, by.x = c("V2","V5"), by.y = c("V1","V3"))
#The final number is lower because all non-missense (3) and reversions are removed
AlphaVariantsInfo<-merge(AlphaVCF, AlphaFromAliAggregated, by.x = c("V2","V5"), by.y = c("V1","V3"))
#In case of Alpha lineage is not provided for NY samples that is why I substitute it with B.1.1.7 for all
AlphaVariantsInfo$lineage<-rep("B.1.1.7", length(AlphaVariantsInfo$V2))

##add mutation in outbreak.info format

DeltaVariantsInfo$mutationup<-str_replace_all(DeltaVariantsInfo$ProtChange,"ORF1ab","ORF1a")
DeltaVariantsInfo$mutation<-tolower(DeltaVariantsInfo$mutationup)

AlphaVariantsInfo$mutationup<-str_replace_all(AlphaVariantsInfo$ProtChange,"ORF1ab","ORF1a")
AlphaVariantsInfo$mutation<-tolower(AlphaVariantsInfo$mutationup)
########################################
#Going to outbreak.info to select the varaints that can be potentially of interest
#Select lineages of interest
#DeltaLineagesToCheck<-c("AY.3", "AY.25", "AY.44", "AY.75","AY.103","B.1.617.2")
AlphaLineagesToCheck<-c("B.1.1.7")


#########################################
####The idea of this one is just to get the list of mutations and check them
MutsSpike<-c("S:L18F","S:T29I","S:H69Y","S:G75V","S:Y145D","S:G476S","S:N501Y","S:P681H")
DeltaLineagesToCheck<-c("AY.3", "AY.25", "AY.44", "AY.75", "AY.103","B.1.617.2")

# DeltaVCF$mutationup<-str_replace_all(DeltaVCF$V13,"ORF1ab","ORF1a")
# DeltaVCF$mutation<-tolower(DeltaVCF$mutationup)
# PreDeltaVCFToCheck<-DeltaVCF[DeltaVCF$V9 >1,]
# DeltaVCFToCheck<-unique(rbind(PreDeltaVCFToCheck,
#                           subset(DeltaVCF,DeltaVCF$mutationup %in% MutsSpike)))
DeltaMutsToCheck<-as.vector(DeltaVariantsInfo$mutationup)
DeltaMutsToChecklow<-tolower(DeltaMutsToCheck)

#devtools::install_github("outbreak-info/R-outbreak-info")
library(outbreakinfo)

authenticateUser()

folder<-"figures_lineages_freqs"
###
#Plot the prevalence of this lineages in Ohio
Del_Oh = getPrevalence(pangolin_lineage = c("AY.3", "AY.25", "AY.44", "AY.75", "AY.103","B.1.617.2"), 
                       location="Ohio")
#subset the time frame of interest
#select only interval when deer samples were collected
DeltaLineagesOhio<-plotPrevalenceOverTime(Del_Oh, title = "")+
  scale_x_date(limits =as.Date(c("2021-07-01","2022-12-01")),
               breaks = "1 month",
               date_labels = "%b %Y")+
  facet_wrap(~lineage, ncol=1)+
  labs(caption="")+
  theme(legend.position = "none")
DeltaLineagesOhio

ggsave("Delta_lineages_Ohio.svg",
       path = folder,
       plot = DeltaLineagesOhio,
       width = 16,
       height = 21,
       dpi =400,
       units= "cm")

#Plot the prevalence of mutations of interest in Ohio in Spike

#Create df to plot
#mutations =  c(MutsSpike[1]), 
Del_US_Delta_lineages = getPrevalence(location = "United States", pangolin_lineage = DeltaLineagesToCheck)
Del_US_muts=getPrevalence(mutations = c(MutsSpike[1]),location = "United States", pangolin_lineage = DeltaLineagesToCheck)
Del_US_muts$lineage = str_replace_all(str_split(Del_US_muts$lineage," AND ",simplify =T)[,1],"[()]","")
Del_US_muts$Mut<-rep(MutsSpike[1],length(Del_US_muts$lineage_count))
###Merged DFs to calculate the rates correctly
DfMergedCorrectly<-merge(Del_US_muts[,c(1,3,9,10,11)],Del_US_Delta_lineages[,c(1,3,9,10)], by=c("date","lineage"))
###
for (Mut in MutsSpike[2:length(MutsSpike)]) {
  loopdf = getPrevalence(mutations =  c(Mut), location = "United States", pangolin_lineage = DeltaLineagesToCheck)
  loopdf$lineage = str_replace_all(str_split(loopdf$lineage," AND ",simplify =T)[,1],"[()]","")
  loopdf$Mut<-rep(Mut,length(loopdf$lineage_count))
  loopdfMergedCorrectly<-merge(loopdf[,c(1,3,9,10,11)],Del_US_Delta_lineages[,c(1,3,9,10)], by=c("date","lineage"))
  Del_US_muts<-rbind(Del_US_muts,loopdf)
  DfMergedCorrectly<-rbind(DfMergedCorrectly,loopdfMergedCorrectly)
  ###The sleep line below is required because outbreak.info limits API requests sometimes for 1 per minute, 
  #of course this in the end takes long, but at least it works without breaking everything
  Sys.sleep(60)
}

#set mutations in correct order
Del_US_muts$Mut<-factor(Del_US_muts$Mut, levels = MutsSpike)
DfMergedCorrectly$Mut<-factor(DfMergedCorrectly$Mut, levels = MutsSpike)
DfMergedCorrectly$Rate<-DfMergedCorrectly$lineage_count.x/DfMergedCorrectly$lineage_count.y
############
##This is the old version of plot that use direct information from outbreak.info directly, so it calculate frequencies 
##on the wrong background
DeltaLineagesUSMut<-ggplot(Del_US_muts, 
                           aes(x = date, y = proportion, colour = linadj,
                               fill = linadj, group = linadj)) +
  geom_ribbon(aes(ymin = proportion_ci_lower, ymax = proportion_ci_upper), alpha = 0.25, size = 0) +
  geom_line(size = 0.7) +
  scale_y_continuous(labels = scales::percent, limits = c(0,0.005)) +
  #plotPrevalenceOverTime(Del_US_muts, title = "")+
  scale_x_date(name ="", limits =as.Date(c("2021-07-01","2022-12-01")),
               breaks = "1 month",
               date_labels = "%b %Y")+
  facet_grid(cols =vars(linadj),rows =vars(Mut))+
  labs(caption="")+
  theme_minimal()+
  theme(legend.position = "none",
        plot.title=element_text(hjust=0.5),
        axis.text.x = element_text(angle=90),
        panel.spacing = unit(1, "lines"))
DeltaLineagesUSMut
#save
ggsave("Delta_Spike_recurrent_muts_AF_humans_US.svg",path = folder,plot=DeltaLineagesUSMut, 
       width = 25,
       height = 25,
       dpi =400,
       units ="cm")

#######
#This is a new less fancy plot that calculates the frequencies of mutations over time on the lineage background
#Not among all samples collected on particular date
DeltaMutsInLineagesNew<-ggplot(DfMergedCorrectly, 
                           aes(x = date, y = Rate, colour = lineage,
                               fill = lineage, group = lineage))+
  geom_line(size = 0.7)+
  facet_grid(cols =vars(lineage),rows =vars(Mut))+
  labs(caption="")+
  theme_minimal()+
  scale_y_continuous(labels = scales::percent)+
  scale_x_date(name ="", limits =as.Date(c("2021-12-01","2022-06-01")),
               breaks = "1 month",
               date_labels = "%b %Y")+
  theme(legend.position = "none",
        plot.title=element_text(hjust=0.5),
        axis.text.x = element_text(angle=90),
        panel.spacing = unit(1, "lines"))
DeltaMutsInLineagesNew

###And separate plot for S:L18F only
DfMergedL18FOnly<-DfMergedCorrectly[which(DfMergedCorrectly$Mut == "S:L18F"),]
# DfMergedL18FOnly$WeekMutCounts<-rollsumr(DfMergedL18FOnly$lineage_count.x, k = 7, fill = NA, align="left")
# DfMergedL18FOnly$WeekLinCounts<-rollsumr(DfMergedL18FOnly$lineage_count.y, k = 7, fill = NA, align="left")
# DfMergedL18FOnly$WeekRate<-DfMergedL18FOnly$WeekMutCounts/DfMergedL18FOnly$WeekLinCounts
#working on the second axis
ylim.prim <- c(min(DfMergedL18FOnly$Rate,na.rm=T), max(DfMergedL18FOnly$Rate,na.rm=T))
ylim.sec <- c(min(DfMergedL18FOnly$lineage_count.y), max(DfMergedL18FOnly$lineage_count.y))

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]
#
L18FInLineagesNew<-ggplot(DfMergedL18FOnly, 
                               aes(x = date, y = Rate, group = lineage))+
  geom_vline(xintercept = as.Date("2021-11-01"), size=0.5, color ="#4292c6")+
  geom_vline(xintercept = as.Date("2021-12-01"), size=0.5, color ="#4292c6")+
  geom_line(aes(y = a + lineage_count.y*b), color = "light grey",
            size = 0.5) +
  geom_line(size = 0.7, color = "#bd0026")+
  scale_y_continuous("AF", sec.axis = sec_axis(~ (. - a)/b,
                                               name = "Sample counts"))+
  coord_cartesian(ylim=c(0, 0.5))+
  facet_wrap(~lineage)+#,scales="free_y")+
  labs(caption="")+
  theme_minimal()+
  scale_x_date(name ="", limits =as.Date(c("2021-07-01","2022-06-01")),
               breaks = "1 month",
               date_labels = "%b %Y")+
  theme(legend.position = "none",
        plot.title=element_text(hjust=0.5),
        axis.text.x = element_text(angle=90),
        panel.spacing = unit(1, "lines"))
L18FInLineagesNew
################
ggsave("Spike_L18F_by_lineage.png",path = folder,plot=L18FInLineagesNew, 
       width = 35,
       height = 20,
       dpi =300,
       units ="cm")
ggsave("Spike_L18F_by_lineage.svg",path = folder,plot=L18FInLineagesNew, 
       width = 35,
       height = 20,
       dpi =300,
       units ="cm")
###############
#Plot AF of selected mutations in those lineages globally, not only in US
#retrieve data
#Common Delta mutations globally
Delta_mutsFreq =getMutationsByLineage(pangolin_lineage = DeltaLineagesToCheck,frequency = 0.5)
#get AF for mutatations in region (modified outbreak.info function)
get_mutations<-function(pangolin_lineage, frequency=0.75, logInfo=TRUE, location)
{
  df <- purrr::map_df(pangolin_lineage, function(lineage) getGenomicData(query_url="lineage-mutations", 
                                                                  pangolin_lineage = lineage, 
                                                                  frequency = 0, logInfo = logInfo,
                                                                  location= location))
  
  if(!is.null(df) && nrow(df) != 0){
    mutations = df %>%
      filter(prevalence >= frequency) %>%
      pull(mutation) %>%
      unique()
    
    df <- df %>%
      filter(mutation %in% mutations)
  }
  return(df)
}
#apply function
char_muts = get_mutations(pangolin_lineage = DeltaLineagesToCheck,frequency = 0.0000001, 
                          location ="United States")
char_muts_to_plot=subset(char_muts,char_muts$mutation %in% c(DeltaMutsToChecklow,tolower(MutsSpike)))

Delta_mutsFreq_to_plot<-subset(Delta_mutsFreq, Delta_mutsFreq$gene == "S")
###plot common Delta mutations
DFtoPlot<-char_muts_to_plot
#create labels and order
DFtoPlot$mutlabel<-as.vector(str_replace(toupper(DFtoPlot$mutation),"A:","a:"))
DFtoPlot$mutlabel<-factor(DFtoPlot$mutlabel,
                          levels = as.vector(DeltaVCFToCheck[order(-DeltaVCFToCheck$V1),]$mutationup))#,toupper(unique(Delta_mutsFreq_to_plot$mutation))))
####Save this data to file
DeltaPrevalenceToSave<-reshape(DFtoPlot[,c("mutlabel","lineage","prevalence")], idvar = "mutlabel", timevar = "lineage", direction ="wide")
write.csv(DeltaPrevalenceToSave, file ="DeltaPrevalenceUS.csv", quote =F, row.names =F)
###
#plotMutationHeatmap(DFtoPlot, title = "", lightBorders = FALSE)
borderColour = "#FFFFFF"
MUTATIONPALETTE = c('#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a')
#add blanks
blank = crossing(unique(DFtoPlot$query_key), unique(DFtoPlot$mutlabel))
names(blank) <-c("query_key","mutlabel")

#Plot
GISAIDAFplot<-ggplot(data=DFtoPlot, aes(y=mutlabel,x=query_key,fill=prevalence))+
  geom_tile(colour = borderColour, fill = "#dedede", data = blank,size=1) +
  geom_tile(color = borderColour,size=1)+
  theme_minimal()+
  coord_fixed() +
  scale_fill_gradientn(colours = MUTATIONPALETTE, limits = c(0,0.23), labels = scales::percent)+
  labs(y ="",x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        panel.grid = element_blank(),
        legend.position = "left")
GISAIDAFplot

############
#Plot common Delta S mutations
Delta_mutsFreq_to_plot$mutation<-factor(toupper(Delta_mutsFreq_to_plot$mutation), 
                                        levels = as.vector(toupper(unique(Delta_mutsFreq_to_plot[order(-Delta_mutsFreq_to_plot$codon_num),]$mutation))))

DeltaCommonPlot<-ggplot(data=Delta_mutsFreq_to_plot, aes(y=mutation,x=query_key,fill=prevalence))+
  #geom_tile(colour = borderColour, fill = "#dedede", data = blank,size=1) +
  geom_tile(color = borderColour,size=1)+
  theme_minimal()+
  coord_fixed() +
  scale_fill_gradientn(colours = c("#fcfbfd","#efedf5","#dadaeb","#bcbddc","#9e9ac8","#807dba","#6a51a3","#54278f","#3f007d"),
                       limits = c(0,1), labels = scales::percent)+
  labs(y ="",x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        panel.grid = element_blank(),
        legend.position = "left")
DeltaCommonPlot

#save
ggsave("Delta_Muts_Human_US_GISAID_CommonS.svg",path = folder,plot=DeltaCommonPlot, 
       width = 10,
       height = 10,
       dpi =400,
       units ="cm")

#
DeltaVCFToCheck$mutationup<-factor(DeltaVCFToCheck$mutationup, 
                                   levels = as.vector(DeltaVCFToCheck[order(-DeltaVCFToCheck$V1),]$mutationup))
HomoplasyPlot<-ggplot(data=DeltaVCFToCheck, aes(y=mutationup,x=V6,fill=V9))+
  geom_tile(colour = borderColour)+
  geom_text(aes(label =V9 ))+
  scale_fill_gradientn(colors = c("#bdd7e7","#6baed6","#3182bd","#08519c"))+
  theme_minimal()+
  coord_fixed() +
  labs(x="",y="")+
  theme(legend.position = "none",
        axis.text = element_blank())
HomoplasyPlot

#AFplot
AFPlot<-ggplot(data=DeltaVCFToCheck, aes(y=mutationup,x=V6,fill=V15/60))+
  geom_tile(colour = borderColour)+
  #geom_text(aes(label =V15 ))+
  scale_fill_gradientn(colours = c('#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a'),
                       limits = c(0,0.23), labels = scales::percent)+
  scale_y_discrete(labels=DeltaVCFToCheck[order(-DeltaVCFToCheck$V1),]$V14,
                   position="right")+
  theme_minimal()+
  coord_fixed() +
  labs(x="",y="")+
  theme(legend.position = "none",
        axis.text.x = element_blank())
AFPlot

GISAIDAFplot+HomoplasyPlot+AFPlot

ggsave("Delta_Deer_Muts_GISAID_Freqs.svg",path = folder,plot=GISAIDAFplot+HomoplasyPlot+AFPlot, 
       width = 20,
       height = 20,
       dpi =400,
       units ="cm")


#########################################
#################################
#Plot the same for Alpha
#########################################
AlphaMutsSpike<-c("S:V3G","S:S12F","S:T29I","S:R246I","S:S247G","S:Y248C")
AlphaLineagesToCheck<-c("B.1.1.7")

AlphaVCF$mutationup<-str_replace_all(AlphaVCF$V13,"ORF1ab","ORF1a")
AlphaVCF$mutation<-tolower(AlphaVCF$mutationup)
PreAlphaVCFToCheck<-AlphaVCF[AlphaVCF$V9 >1,]
AlphaVCFToCheck<-unique(rbind(PreAlphaVCFToCheck,
                              subset(AlphaVCF,AlphaVCF$mutationup %in% AlphaMutsSpike)))
AlphaMutsToCheck<-as.vector(AlphaVCFToCheck$mutationup)
AlphaMutsToChecklow<-tolower(AlphaMutsToCheck)

#Plot the prevalence of this lineages in Ohio
Alpha_Oh = getPrevalence(pangolin_lineage = AlphaLineagesToCheck, location="Ohio")
#subset the time frame of interest
#select only interval when deer samples were collected
AlphaLineagesOhio<-plotPrevalenceOverTime(Alpha_Oh, title = "")+
  scale_x_date(limits =as.Date(c("2021-07-01","2022-02-01")),
               breaks = "1 month",
               date_labels = "%b %Y")+
  facet_wrap(~lineage, ncol=1)+
  labs(caption="")+
  theme(legend.position = "none")
AlphaLineagesOhio

ggsave("Alpha_lineage_Ohio.svg",
       path = folder,
       plot = AlphaLineagesOhio,
       width = 16,
       height = 10,
       dpi =400,
       units= "cm")

#####
#mutations
Alpha_mutsFreq =getMutationsByLineage(pangolin_lineage = AlphaLineagesToCheck,frequency = 0.5)
char_mutsA = get_mutations(pangolin_lineage = AlphaLineagesToCheck,frequency = 0.0000001,
                           location="United States")
char_muts_to_plotA=subset(char_mutsA,char_mutsA$mutation %in% AlphaMutsToChecklow)

Alpha_mutsFreq_to_plot<-subset(Alpha_mutsFreq, Alpha_mutsFreq$gene == "S")

#DFtoPlot<-rbind(char_muts_to_plot,Delta_mutsFreq_to_plot)
DFtoPlotA<-char_muts_to_plotA
#create labels and order

DFtoPlotA$mutlabel<-factor(as.vector(str_replace(toupper(DFtoPlotA$mutation),"A:","a:")),
                          levels = c(as.vector(AlphaVCFToCheck[order(-AlphaVCFToCheck$V1),]$mutationup)))#,toupper(unique(Delta_mutsFreq_to_plot$mutation))))
###
AlphaPrevalenceToSave<-reshape(DFtoPlotA[,c("mutlabel","lineage","prevalence")], idvar = "mutlabel", timevar = "lineage", direction ="wide")
write.csv(AlphaPrevalenceToSave, file ="AlphaPrevalenceUS.csv", quote =F, row.names =F)

#plotMutationHeatmap(DFtoPlot, title = "", lightBorders = FALSE)
borderColour = "#FFFFFF"
MUTATIONPALETTEALPHA = c("#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026")
#add blanks
blank = crossing(unique(DFtoPlotA$query_key), unique(DFtoPlotA$mutlabel))
names(blank) <-c("query_key","mutlabel")

#Plot
AlphaAFplot<-ggplot(data=DFtoPlotA, aes(y=mutlabel,x=query_key,fill=prevalence))+
  geom_tile(colour = borderColour, fill = "#dedede", data = blank,size=1) +
  geom_tile(color = borderColour,size=1)+
  theme_minimal()+
  coord_fixed() +
  scale_fill_gradientn(colours = MUTATIONPALETTEALPHA, limits = c(0,0.23), labels = scales::percent)+
  labs(y ="",x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        panel.grid = element_blank(),
        legend.position = "left")
AlphaAFplot

#common
Alpha_mutsFreq_to_plot$mutation<-factor(toupper(Alpha_mutsFreq_to_plot$mutation), 
                                        levels = as.vector(toupper(unique(Alpha_mutsFreq_to_plot[order(-Alpha_mutsFreq_to_plot$codon_num),]$mutation))))

AlphaCommonPlot<-ggplot(data=Alpha_mutsFreq_to_plot, aes(y=mutation,x=query_key,fill=prevalence))+
  #geom_tile(colour = borderColour, fill = "#dedede", data = blank,size=1) +
  geom_tile(color = borderColour,size=1)+
  theme_minimal()+
  coord_fixed() +
  scale_fill_gradientn(colours = c("#fcfbfd","#efedf5","#dadaeb","#bcbddc","#9e9ac8","#807dba","#6a51a3","#54278f","#3f007d"),
                       limits = c(0,1), labels = scales::percent)+
  labs(y ="",x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        panel.grid = element_blank(),
        legend.position = "left")
AlphaCommonPlot

ggsave("Alpha_Muts_Human_US_GISAID_CommonS.svg",path = folder,plot=AlphaCommonPlot, 
       width = 10,
       height = 10,
       dpi =400,
       units ="cm")

#homoplasy
AlphaVCFToCheck$mutationup<-factor(AlphaVCFToCheck$mutationup, levels = as.vector(AlphaVCFToCheck[order(-AlphaVCFToCheck$V1),]$mutationup))
AlphaHomoplasyPlot<-ggplot(data=AlphaVCFToCheck, aes(y=mutationup,x=V6,fill=V9))+
  geom_tile(colour = borderColour)+
  geom_text(aes(label =V9 ))+
  scale_fill_gradientn(colors = c("#bdd7e7","#6baed6","#3182bd","#08519c"))+
  theme_minimal()+
  coord_fixed() +
  labs(x="",y="")+
  theme(legend.position = "none",
        axis.text = element_blank())
AlphaHomoplasyPlot

#AFplot
AAFPlot<-ggplot(data=AlphaVCFToCheck, aes(y=mutationup,x=V6,fill=V15/33))+
  geom_tile(colour = borderColour)+
  #geom_text(aes(label =V15 ))+
  scale_fill_gradientn(colours = MUTATIONPALETTEALPHA,
                       limits = c(0,0.3), labels = scales::percent)+
  scale_y_discrete(labels=AlphaVCFToCheck[order(-AlphaVCFToCheck$V1),]$V14,
                   position="right")+
  theme_minimal()+
  coord_fixed() +
  labs(x="",y="")+
  theme(legend.position = "none",
        axis.text.x = element_blank())
AAFPlot

ggsave("Alpha_Deer_Muts_GISAID_Freqs.svg",path = folder,plot=AlphaAFplot+AlphaHomoplasyPlot+AAFPlot, 
       width = 16,
       height = 20,
       dpi =400,
       units ="cm")

####Plot AF over time for individual mutations
#Create df to plot
Alpha_US_muts = getPrevalence(mutations =  c(AlphaMutsSpike[1]), location = "United States", pangolin_lineage = AlphaLineagesToCheck)
Alpha_US_muts$linadj = str_replace_all(str_split(Alpha_US_muts$lineage," AND ",simplify =T)[,1],"[()]","")
Alpha_US_muts$Mut<-rep(MutsSpike[1],length(Alpha_US_muts$lineage_count))
for (Mut in AlphaMutsSpike[2:length(AlphaMutsSpike)]) {
  loopdf = getPrevalence(mutations =  c(Mut), location = "United States", pangolin_lineage = AlphaLineagesToCheck)
  loopdf$linadj = str_replace_all(str_split(loopdf$lineage," AND ",simplify =T)[,1],"[()]","")
  loopdf$Mut<-rep(Mut,length(loopdf$lineage_count))
  Alpha_US_muts<-rbind(Alpha_US_muts,loopdf)
}

#set mutations in correct order
Alpha_US_muts$Mut<-factor(Alpha_US_muts$Mut, levels = AlphaMutsSpike)
##plot
AlphaLineagesUSMut<-ggplot(Alpha_US_muts, 
                           aes(x = date, y = proportion, colour = linadj,
                               fill = linadj, group = linadj)) +
  geom_ribbon(aes(ymin = proportion_ci_lower, ymax = proportion_ci_upper), alpha = 0.25, size = 0) +
  geom_line(size = 0.7) +
  scale_y_continuous(labels = scales::percent, limits = c(0,0.005)) +
  #plotPrevalenceOverTime(Del_US_muts, title = "")+
  scale_x_date(name ="", limits =as.Date(c("2021-07-01","2022-02-28")),
               breaks = "1 month",
               date_labels = "%b %Y")+
  facet_grid(cols =vars(linadj),rows =vars(Mut))+
  labs(caption="")+
  theme_minimal()+
  theme(legend.position = "none",
        plot.title=element_text(hjust=0.5),
        axis.text.x = element_text(angle=90),
        panel.spacing = unit(1, "lines"))
AlphaLineagesUSMut

ggsave("Alpha_Spike_recurrent_muts_AF_humans_US.svg",
       path = folder,
       plot=AlphaLineagesUSMut, 
       width = 12,
       height = 18,
       dpi =400,
       units ="cm")