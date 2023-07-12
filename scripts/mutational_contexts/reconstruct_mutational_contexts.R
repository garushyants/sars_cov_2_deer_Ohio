library(MutationalPatterns)
library(BSgenome.SARScov2.NCBI.ASM985889v2)
library(ggplot2)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(stringr)
library("rstudioapi")
library(DescTools)

path<-getSourceEditorContext()$path
setwd(SplitPath(path)$dirname)
setwd("../../data/mutational_contexts")
folder<-"../../figures/mutational_contexts"
##Load vcfs
sample_names <-c("Deer Alpha","Deer Delta","Human")
vcf_files<-c("Alpha.clusters.snpeff.nopr.forMutationalSignatures.vcf",
             "Delta.clusters.snpeff.nopr.forMutationalSignatures.vcf",
             "UShER_all.vcf")

sample_namesS <-c("Deer Alpha","Deer Delta","Human")
vcf_filesS<-c("Alpha.clusters.snpeff.nopr.forMutationalSignatures.synonymous.vcf",
             "Delta.clusters.snpeff.nopr.forMutationalSignatures.synonymous.vcf",
             "UShER_synonymous.vcf")

###genome
ref_genome<-"BSgenome.SARScov2.NCBI.ASM985889v2"
genome<-getBSgenome(ref_genome)

vcfsA <- read_vcfs_as_granges(vcf_files, sample_names, genome = ref_genome)
summary(vcfsA)
vcfsS <- read_vcfs_as_granges(vcf_filesS, sample_namesS, genome = ref_genome)
summary(vcfsS)

###
#Plot
type_occurrencesA <- mut_type_occurrences(vcfsA, ref_genome)
type_occurrencesS <- mut_type_occurrences(vcfsS, ref_genome)

#plot muts by type all
SingleMutA<-plot_spectrum(type_occurrencesA,by=rownames(type_occurrencesA),legend = TRUE) +
  theme_minimal()+
  theme(text = element_text(size=10),
        legend.title = element_text(size=8))
SingleMutA

#plot muts by type synonymous
SingleMutS<-plot_spectrum(type_occurrencesS,by=rownames(type_occurrencesS),legend = TRUE) +
  theme_minimal()+
  theme(text = element_text(size=10),
        legend.title = element_text(size=8))
SingleMutS

#plot contexts all
mut_matA <- mut_matrix(vcf_list = vcfsA, ref_genome = ref_genome)
head(mut_matA)
Mut_96_plotA<-plot_96_profile(mut_matA) +
  theme_minimal()+
  theme(text = element_text(size=10),
        axis.text.x = element_text(size=6, angle=90, hjust = 0.5))
Mut_96_plotA

#plot contexts synonymous
mut_matS <- mut_matrix(vcf_list = vcfsS, ref_genome = ref_genome)
head(mut_matS)
Mut_96_plotS<-plot_96_profile(mut_matS) +
  theme_minimal()+
  theme(text = element_text(size=10),
        axis.text.x = element_text(size=6, angle=90, hjust = 0.5))
Mut_96_plotS


#Save figures
PlotToSaveA<-ggarrange(SingleMutA,Mut_96_plotA,
          ncol = 1,
          heights = c(0.8,1),
          labels = c("a","b"),
          font.label = list(size = 10, face = "plain"))
PlotToSaveA
ggsave("SARS-CoV-2_deer_clusters_human_all_mutcontexts.svg", plot = PlotToSaveA,path=folder,
       width =20, height =20, units='cm', dpi =400)

PlotToSaveS<-ggarrange(SingleMutS,Mut_96_plotS,
                       ncol = 1,
                       heights = c(0.8,1),
                       labels = c("a","b"),
                       font.label = list(size = 10, face = "plain"))
PlotToSaveS
ggsave("SARS-CoV-2_deer_clusters_human_synonymous_mutcontexts.svg", path=folder,plot = PlotToSaveS,
       width =20, height =20, units='cm', dpi =400)

###############
#Draw subsamples for Human all sample
setwd("./subsamples_human_435_10/")
subsample_names <-as.character(seq(1,10))
vcf_files_subsample<-c("all1.header.vcf","all2.header.vcf","all3.header.vcf","all4.header.vcf","all5.header.vcf",
             "all6.header.vcf","all7.header.vcf","all8.header.vcf","all9.header.vcf","all10.header.vcf")

vcfsSub <- read_vcfs_as_granges(vcf_files_subsample, subsample_names, genome = ref_genome, 
                                predefined_dbs_mbs = T)#,
                                #remove_duplicate_variants = F)

type_occurrencesSubsample <- mut_type_occurrences(vcfsSub, ref_genome)

#plot muts by type all
SingleMutSubsample<-plot_spectrum(type_occurrencesSubsample,legend = TRUE,
                                  indv_points = T,
                                  error_bars = "stdev") +
  theme_minimal()+
  theme(text = element_text(size=10),
        legend.title = element_text(size=8))+
  scale_y_continuous(limits = c(0,0.75))
SingleMutSubsample


# mut_matSub <- mut_matrix(vcf_list = vcfsSub, ref_genome = ref_genome)
# Mut_96_plotSub<-plot_96_profile(mut_matSub) +
#   theme_minimal()+
#   theme(text = element_text(size=10),
#         axis.text.x = element_text(size=6, angle=90, hjust = 0.5))
# Mut_96_plotSub
setwd("./../")
ggsave("SARS-CoV-2_human_subsamples_spectrum_all.svg", path=folder,plot = SingleMutSubsample,
       width =15, height =20, units='cm', dpi =300)

