library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(stringr)
library(reshape2)
library(cowplot)
#library(zoo)
library("rstudioapi")

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
DeltaLineagesToCheck<-c("AY.3", "AY.25", "AY.44", "AY.75", "AY.103","B.1.617.2")
AlphaLineagesToCheck<-c("B.1.1.7")


###Authentificate in outbreak.info
#devtools::install_github("outbreak-info/R-outbreak-info")
library(outbreakinfo)

authenticateUser()

folder<-"../../figures/lineages_frequencies"

####Get the list of mutations in Spike
AllDeltaMutsSpike<-DeltaVariantsInfo[DeltaVariantsInfo$gene=="S",]$ProtChange
AllAlphaMutsSpike<-AlphaVariantsInfo[AlphaVariantsInfo$gene=="S",]$ProtChange

####Get mutations characteristic for lineages (frequent globally)
Delta_mutsFreq =getMutationsByLineage(pangolin_lineage = DeltaLineagesToCheck,frequency = 0.5)
Alpha_mutsFreq =getMutationsByLineage(pangolin_lineage = AlphaLineagesToCheck,frequency = 0.5)

#############
#Now I can select interesting mutations in Spike

#First let's select the recurrent ones or common for Alpha and Delta datasets
DeltaMutsSpikeDf<-subset(DeltaVariantsInfo,
                       DeltaVariantsInfo$gene=="S" & 
                         (DeltaVariantsInfo$mutation %in% Alpha_mutsFreq$mutation | 
                            DeltaVariantsInfo$mutation %in% AlphaVariantsInfo$mutation |
                            DeltaVariantsInfo$events > 1))

AlphaMutsSpikeDf<-subset(AlphaVariantsInfo,
                         AlphaVariantsInfo$gene=="S" & 
                           (AlphaVariantsInfo$mutation %in% Delta_mutsFreq$mutation | 
                              AlphaVariantsInfo$mutation %in% DeltaVariantsInfo$mutation |
                              AlphaVariantsInfo$events > 1))

#Then let's add mutations with at least some mentions in the literature
DeltaKnown<-c("S:G75V","S:Y145D","S:K150I","S:257D","S:G476S",
              "S:P527L", "S:A681H","S:T859I")
AlphaKnown<-c("S:V3G","S:S12F","S:R246I","S:S247G","S:Y248C",
              "S:P384S","S:V445A","S:483A")
DeltaMutsSpikeDfKnown<-DeltaVariantsInfo[DeltaVariantsInfo$ProtChange %in% DeltaKnown,]
AlphaMutsSpikeDfKnown<-AlphaVariantsInfo[AlphaVariantsInfo$ProtChange %in% AlphaKnown,]

#Combining
DeltaMutsSpike<-rbind(DeltaMutsSpikeDf,DeltaMutsSpikeDfKnown)
AlphaMutsSpike<-rbind(AlphaMutsSpikeDf,AlphaMutsSpikeDfKnown)

############
#In other genes I am only selecting recurrent mutations
DeltaMutsOther<-DeltaVariantsInfo[(DeltaVariantsInfo$events >1 & DeltaVariantsInfo$gene !="S"),]
AlphaMutsOther<-AlphaVariantsInfo[(AlphaVariantsInfo$events >1 & AlphaVariantsInfo$gene !="S"),]

###########
#Get frequencies of variants in US
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
Delta_char_muts = get_mutations(pangolin_lineage = DeltaLineagesToCheck,frequency = 0.0000001, 
                          location ="United States")
Alpha_char_muts = get_mutations(pangolin_lineage = AlphaLineagesToCheck,frequency = 0.0000001, 
                                location ="United States")



#################################
#####Generating table with frequencies for all mutations of interest for supplimentary data
DeltaDfToSaveMutsFreqLong<-Delta_char_muts[(Delta_char_muts$mutation %in% c(DeltaMutsSpike$mutation,DeltaMutsOther$mutation)),]
AlphaDfToSaveMutsFreqLong<-Alpha_char_muts[(Alpha_char_muts$mutation %in% c(AlphaMutsSpike$mutation,AlphaMutsOther$mutation)),]
#reshape
DeltaDfToSaveMutsFreqWide<-reshape(DeltaDfToSaveMutsFreqLong[,c("mutation","lineage","prevalence")], idvar = "mutation", timevar = "lineage", direction ="wide")
DeltaToSaveNoSingl<-merge(rbind(DeltaMutsSpike,DeltaMutsOther),DeltaDfToSaveMutsFreqWide, by="mutation", all=T)[,c(2:4,7:8,10,12:14,16:21)]
DeltaToSaveSingl<-merge(DeltaToSaveNoSingl,DeltaFromAliSingletons[,c(1,3:6)],
                           by.x = c("V2","V5"), by.y = c("V1","V3"))

names(DeltaToSaveSingl)[names(DeltaToSaveSingl) == 'count.x'] <- 'count'
names(DeltaToSaveSingl)[names(DeltaToSaveSingl) == 'lineage.x'] <- 'lineage'
DeltaToSaveFinal<-dplyr::bind_rows(DeltaToSaveNoSingl[!(DeltaToSaveNoSingl$ProtChange %in% DeltaToSaveSingl$ProtChange),],DeltaToSaveSingl)
DeltaToSaveFinal$dataset<-rep("Delta",length(DeltaToSaveFinal$V2))

#
AlphaDfToSaveMutsFreqWide<-reshape(AlphaDfToSaveMutsFreqLong[,c("mutation","lineage","prevalence")], idvar = "mutation", timevar = "lineage", direction ="wide")
AlphaToSaveNoSingl<-merge(rbind(AlphaMutsSpike,AlphaMutsOther),AlphaDfToSaveMutsFreqWide, by="mutation", all=T)[,c(2:4,7:8,10,12:14,16)]
AlphaToSaveNoSingl$dataset<-rep("Alpha",length(AlphaToSaveNoSingl$V2))

###
MergeBothDatasets<-dplyr::bind_rows(DeltaToSaveFinal,AlphaToSaveNoSingl)
names(MergeBothDatasets)<-c("POS","ALT","REF","Protein change","gene","EventsInClusters","SamplesWithMutation","NumberOfSamplesInClustersWithMutation",
                            "Lineages","AY.3","AY.25","AY.44","AY.75","AY.103","B.1.617.2","SingletonsWithMutation","NumberOfSingletonsWithMutation",
                            "LineagesSingletons","Dataset","B.1.1.7")
MergeBothDatasets[is.na(MergeBothDatasets)] <- ""
write.table(MergeBothDatasets, file ="Table_S8.tsv", quote =F, row.names =F, sep="\t")
##################################

##################################
####Ploting variants
#Get required data
GetPrevalenceOfMutationsOverTime<-function(mutations,lineages)
{
  # mutations<-DeltaMutsSpike$mutation
  # lineages<-DeltaLineagesToCheck
  Del_US_Delta_lineages = getPrevalence(location = "United States", pangolin_lineage = lineages)
  Del_US_muts=getPrevalence(mutations = c(mutations[1]),location = "United States", pangolin_lineage = lineages)
  Del_US_muts$lineage = str_replace_all(str_split(Del_US_muts$lineage," AND ",simplify =T)[,1],"[()]","")
  Del_US_muts$Mut<-rep(DeltaMutsSpike$mutation[1],length(Del_US_muts$lineage_count))
  ###Merged DFs to calculate the rates correctly
  DfMergedCorrectly<-merge(Del_US_muts[,c(1,3,9,10,11)],Del_US_Delta_lineages[,c(1,3,9,10)], by=c("date","lineage"))
  ###
  for (Mut in mutations[2:length(mutations)]) {
    loopdf = getPrevalence(mutations =  c(Mut), location = "United States", pangolin_lineage = lineages)
    if (nrow(loopdf) >0) #here I am only selecting the ones for which data is available
    {
      loopdf$lineage = str_replace_all(str_split(loopdf$lineage," AND ",simplify =T)[,1],"[()]","")
      loopdf$Mut<-rep(Mut,length(loopdf$lineage_count))
      loopdfMergedCorrectly<-merge(loopdf[,c(1,3,9,10,11)],Del_US_Delta_lineages[,c(1,3,9,10)], by=c("date","lineage"))
      Del_US_muts<-rbind(Del_US_muts,loopdf)
      DfMergedCorrectly<-rbind(DfMergedCorrectly,loopdfMergedCorrectly)
      ###The sleep line below is required because outbreak.info limits API requests sometimes for 1 per minute, 
      #of course this in the end takes long, but at least it works without breaking everything
      Sys.sleep(15)
    }
  }
  
  #set mutations in correct order
  Del_US_muts$Mut<-factor(Del_US_muts$Mut, levels = mutations)
  DfMergedCorrectly$Mut<-factor(DfMergedCorrectly$Mut, levels = mutations)
  DfMergedCorrectly$Rate<-DfMergedCorrectly$lineage_count.x/DfMergedCorrectly$lineage_count.y
  
  return(DfMergedCorrectly)
}

DeltaDfMergedCorrectly<-GetPrevalenceOfMutationsOverTime(DeltaMutsSpike$mutation,DeltaLineagesToCheck)

#Plot the prevalence of mutations of interest in United States in Spike over time
##Delta
DeltaMutsInLineages<-ggplot(DeltaDfMergedCorrectly, 
                               aes(x = date, y = Rate, colour = lineage,
                                   fill = lineage, group = lineage))+
  geom_line(linewidth = 0.7)+
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
DeltaMutsInLineages

ggsave("Delta_spike_muts_over_time_in_humans.png",path = folder,plot=DeltaMutsInLineages, 
       width = 25,
       height = 40,
       dpi =300,
       units ="cm")

###And separate plot for S:L18F only
DfMergedL18FOnly<-DeltaDfMergedCorrectly[which(DeltaDfMergedCorrectly$Mut == "s:l18f"),]

ylim.prim <- c(min(DfMergedL18FOnly$Rate,na.rm=T), max(DfMergedL18FOnly$Rate,na.rm=T))
ylim.sec <- c(min(DfMergedL18FOnly$lineage_count.y), max(DfMergedL18FOnly$lineage_count.y))

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]
#
L18FInLineages<-ggplot(DfMergedL18FOnly, 
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
L18FInLineages
################
ggsave("Spike_L18F_by_lineage.png",path = folder,plot=L18FInLineages, 
       width = 35,
       height = 20,
       dpi =300,
       units ="cm")
ggsave("Spike_L18F_by_lineage.svg",path = folder,plot=L18FInLineages, 
       width = 35,
       height = 20,
       dpi =300,
       units ="cm")

#Alpha
AlphaDfMergedCorrectly<-GetPrevalenceOfMutationsOverTime(AlphaMutsSpike$mutation,AlphaLineagesToCheck)

#Plot the prevalence of mutations of interest in United States in Spike over time
##Delta
AlphaMutsInLineages<-ggplot(AlphaDfMergedCorrectly, 
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
AlphaMutsInLineages

#############################################################
#Plots for Figure 4

#Selecting colors
borderColour = "#FFFFFF"
MUTATIONPALETTE = c('#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a')
MUTATIONPALETTEALPHA = c("#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026")

plotOutbreakInfoAFPerLineageOutbreak<-function(borderColor,pallete,df,allDeerSpikeMutsToPlot)
{
  df$mutlabel<-toupper(df$mutation)
  Muts<-vector()
  
  if(missing(allDeerSpikeMutsToPlot)) {
    Muts<-factor(df$mutlabel,levels = as.vector(unique(df[order(-df$codon_num),]$mutlabel)))
    #df$mutlabel<-factor(df$mutlabel,levels = as.vector(unique(df[order(-df$codon_num),]$mutlabel)))
  } else {
    Muts<-factor(allDeerSpikeMutsToPlot,levels=allDeerSpikeMutsToPlot)
    #df$mutlabel<-factor(df$mutlabel,levels = as.vector(allDeerSpikeMutsToPlot))
  }
  
  blank = crossing(unique(df$query_key), Muts)
  names(blank) <-c("query_key","mutlabel")
  
 AFplot<-ggplot(data=df, aes(y=mutlabel,x=query_key,fill=prevalence))+
    geom_tile(colour = borderColour, fill = "#dedede", data = blank,size=1) +
    geom_tile(color = borderColour,size=1)+
    theme_minimal()+
    coord_fixed() +
    scale_fill_gradientn(colours = pallete, limits = c(0.0,1.0), labels = scales::percent)+
    labs(y ="",x="")+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
          panel.grid = element_blank(),
          legend.position = "left")
  AFplot
  return(AFplot)
  
}

#plots
#Delta
Delta_mutsFreq_to_plot<-subset(Delta_mutsFreq, Delta_mutsFreq$gene == "S")
Delta_spikemutsDeer_toPlot<-DeltaDfToSaveMutsFreqLong[DeltaDfToSaveMutsFreqLong$gene == "S",]
DeltaCommonMuts<-plotOutbreakInfoAFPerLineageOutbreak(borderColour,MUTATIONPALETTE,Delta_mutsFreq_to_plot)
DeltaCommonMuts
DeltaDeerMuts<-plotOutbreakInfoAFPerLineageOutbreak(borderColour,MUTATIONPALETTE,Delta_spikemutsDeer_toPlot,DeltaMutsSpike[order(-DeltaMutsSpike$V2),]$ProtChange)
DeltaDeerMuts
#Alpha
Alpha_mutsFreq_to_plot<-subset(Alpha_mutsFreq, Alpha_mutsFreq$gene == "S")
Alpha_spikemutsDeer_toPlot<-AlphaDfToSaveMutsFreqLong[AlphaDfToSaveMutsFreqLong$gene == "S",]
AlphaCommonMuts<-plotOutbreakInfoAFPerLineageOutbreak(borderColour,MUTATIONPALETTEALPHA,Alpha_mutsFreq_to_plot)
AlphaCommonMuts
AlphaDeerMuts<-plotOutbreakInfoAFPerLineageOutbreak(borderColour,MUTATIONPALETTEALPHA,Alpha_spikemutsDeer_toPlot,AlphaMutsSpike[order(-AlphaMutsSpike$V2),]$ProtChange)
AlphaDeerMuts

###Plot AF and number of events in deer
HomoplasyPlot<-function(df,cltotal)
{
  DeltaDeerToPlot<-df[df$gene == "S",]
  DeltaDeerToPlot$ProtChange<-factor(DeltaDeerToPlot$ProtChange, levels=as.factor(DeltaDeerToPlot[order(-DeltaDeerToPlot$V2),]$ProtChange))
  DeltaDeerToPlot$x<-rep("# of events", length(DeltaDeerToPlot$V2))
  #
  plot<-ggplot(data=DeltaDeerToPlot, aes(y=ProtChange,x=x,fill=events))+
    geom_tile(colour = borderColour)+
    geom_text(aes(label =events ))+
    scale_fill_gradientn(colors = c("#bdd7e7","#6baed6","#3182bd","#08519c"))+
    theme_minimal()+
    coord_fixed() +
    labs(x="",y="") +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))
  return(plot)
}
  
AFPlot<-function(df,cltotal,pallete)
{
  DeltaDeerToPlot<-df[df$gene == "S",]
  DeltaDeerToPlot$ProtChange<-factor(DeltaDeerToPlot$ProtChange, levels=as.factor(DeltaDeerToPlot[order(-DeltaDeerToPlot$V2),]$ProtChange))
  DeltaDeerToPlot$AFinDeer<-DeltaDeerToPlot$count/cltotal
  DeltaDeerToPlot$x<-rep("AF in deer", length(DeltaDeerToPlot$V2))
  #
  plot<-ggplot(data=DeltaDeerToPlot, aes(y=ProtChange,x=x,fill=AFinDeer))+
    geom_tile(colour = borderColour)+
    scale_fill_gradientn(colors = pallete)+
    theme_minimal()+
    coord_fixed() +
    labs(x="",y="") +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))
  return(plot)
}



DclusterTotal<-47
AclusterTotal<-33

DeltaAFPlot<-AFPlot(DeltaToSaveFinal,DclusterTotal,MUTATIONPALETTE)
DeltaHomoplasyPlot<-HomoplasyPlot(DeltaToSaveFinal)

DeltaAll<-plot_grid(DeltaCommonMuts,NULL, NULL,
                    DeltaDeerMuts,DeltaAFPlot, DeltaHomoplasyPlot,
                    ncol=3,
                    nrow = 2,
                    byrow = TRUE,
                    align ="v",axis="tb",
                    rel_heights = c(0.88,1),
                    rel_widths = c(6,1,1))


AlphaAFPlot<-AFPlot(AlphaToSaveNoSingl,AclusterTotal,MUTATIONPALETTEALPHA)
AlphaHomoplasyPlot<-HomoplasyPlot(AlphaToSaveNoSingl)

AlphaAll<-plot_grid(AlphaCommonMuts,NULL, NULL,
                    AlphaDeerMuts,AlphaAFPlot, AlphaHomoplasyPlot,
                    ncol=3,
                    nrow = 2,
                    byrow = TRUE,
                    align ="v",axis="tb",
                    rel_heights = c(1,1),
                    rel_widths = c(6,1,1))

###Saving final plots
ggsave("Delta_AF_all.svg",path = folder,plot=DeltaAll, 
       width = 20,
       height = 20,
       dpi =300,
       units ="cm")

ggsave("Alpha_AF_all.svg",path = folder,plot=AlphaAll, 
       width = 20,
       height = 20,
       dpi =300,
       units ="cm")


