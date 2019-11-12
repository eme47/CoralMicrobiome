#install binaries from Bioconductor and packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
BiocManager::install("decontam")

source("http://bioconductor.org/biocLite.R")
biocLite(supressUpdates=FALSE)

library(decontam)
packageVersion("decontam")

library("ggplot2")
packageVersion("ggplot2")


#install DADA2
source("http://bioconductor.org/biocLite.R")
biocLite(supressUpdates=FALSE)
biocLite("Shortread", suppressUpdates = FALSE)
biocLite("devtools")
library("devtools")
packageVersion("devtools")

devtools::install_github("benjjneb/dada2")
packageVersion("dada2")
library("dada2")

version


path="/Users/eeggleston/Desktop/CoralMicrobiome"
list.files(path)

#bring in forward reads
Fs=sort(list.files(path,pattern="_R1_001.fastq", full.names = TRUE))
sample.names=sapply(strsplit(basename(Fs),"_"),'[',1)
sample.names
Fs

#inspect read quality
plotQualityProfile(Fs[1:2])

#filter and trim
filterFs=file.path(path, "filtered_Fs", paste0(sample.names, "_F_filt.fastq.gz"))

outF=filterAndTrim(Fs, filterFs, truncLen = c(240), trimLeft=10, maxN=0, maxEE=c(2), truncQ = 2,rm.phix = TRUE,compress = TRUE,multithread = TRUE)
head(outF)

errFs=learnErrors(filterFs,multithread = TRUE)  
plotErrors(errFs, nominalQ = TRUE)

#dereplicate reads with high abundance and give them sample names
derepsFs=derepFastq(filterFs, verbose=TRUE)
names(derepsFs) <- sample.names
head(derepsFs)

dadaF=dada(derepsFs,err=errFs, multithread = TRUE)
dadaF[[1]]
dadaF[[7]]
dadaF[[13]]

seqtabF=makeSequenceTable(dadaF)
dim(seqtabF)

#inspect distribution of sequence lengths
table(nchar(getSequences(seqtabF)))

#remove chimeras
seqtabF.nochim=removeBimeraDenovo(seqtabF,method="consensus", multithread=TRUE,verbose=TRUE)
dim(seqtabF.nochim)
class(seqtabF.nochim)
head(seqtabF.nochim)
sum(seqtabF.nochim)/sum(seqtabF)
#About 16% chimeric reads

getNF=function(x) sum(getUniques(x))
getNF

summary_tabF=data.frame(row.names=sample.names, dada2_input=outF[,1],
                       filtered=outF[,2],dadaF=sapply(dadaFs, getNF),
                       nonchim=rowSums(seqtabF.nochim),
                       final_perc_reads_retained=round(rowSums(seqtabF.nochim)/outF[,1]*100, 1))
summary_tabF

getwd()
list.files(path)
#taxaF=assignTaxonomy(seqtabF.nochim, "/Users/eeggleston/Desktop/CoralMicrobiome/SILVA_132_SSURef_Nr99_tax_silva.fasta", multithread = FALSE)
taxaF=assignTaxonomy(seqtabF.nochim, "/Users/eeggleston/Desktop/CoralMicrobiome/silva_nr_v132_train_set.fa.gz", multithread = TRUE)

taxa.printF=taxaF
rownames(taxa.printF)=NULL
head(taxa.printF)
tail(taxa.printF)

#giving our seq headers more manageable names (ASV_1, ASV_2...)
asvF_seqs=colnames(seqtabF.nochim)
asvF_seqs
asvF_headers= vector(dim(seqtabF.nochim)[2], mode="character")
asvF_headers
for (i in 1:dim(seqtabF.nochim)[2]){
  asvF_headers[i]=paste(">ASV", sep="_")
}
asvF_headers
head(asvF_seqs)
head(asvF_headers)

#making and writing out a fasta of our final ASV seqs
asvF_fasta=c(rbind(asvF_headers,asvF_seqs))
write(asvF_fasta, "ASVsFs.fa")
head(asvF_fasta)

#count table
asvF_tab=t(seqtabF.nochim)
row.names(asvF_tab)= sub(">", "", asvF_headers)
head(asvF_tab)
write.table(asvF_tab, "ASVsFs_counts.tsv", sep = "\t", quote=F, col.names = NA)

#tax table
asvF_tax=taxaF
row.names(asvF_tax)= sub(">", "", asvF_headers)
write.table(asvF_tab, "ASVsFs_taxonomy.tsv", sep="\t", quote=F, col.names = NA)

#phyloseq analysis
library(phyloseq)
packageVersion("phyloseq")
sd <- read.csv("coralMD.csv", header = TRUE, sep = ",")
row.names(sd)<- sd$sample
samples= sample_data(sd)
psF <- phyloseq(otu_table(seqtabF.nochim,taxa_are_rows = FALSE),sample_data(sd),tax_table(taxaF))
psF

str(psF)
stF <- subset_taxa(psF,
                       Family !="Mitochondria" &
                         Class !="Chloroplast")

stF

stF2 <- subset_taxa(psF,
                   Family !="Mitochondria" &
                     Class !="Chloroplast" &
                     Order !="Chloroplast")

stF2
###PHYLOSEQ

plot_richness(psF, measures=c("Shannon", "Simpson"))
plot_richness(stF, measures=c("Shannon", "Simpson"))

plot_richness(stF, measures="Shannon", color = "collection_site")

estimate_richness(stF, measures=c("Observed","Shannon", "InvSimpson"))

#Includes Mitochondria and Chloroplast ASvs
psF.prop <- transform_sample_counts(psF, function(otu) otu/sum(otu))
ord.nmds.brayF <- ordinate(psF.prop, method="NMDS", distance="bray")


plot_ordination(psF.prop, ord.nmds.brayF, title="Bray NMDS",color="collection_site")

#subset taxa Forward reads (minus Mitochondria and Chloroplast)
stF.prop <- transform_sample_counts(stF, function(otu) otu/sum(otu))
ord.nmds.braystF <- ordinate(stF.prop, method="NMDS", distance="bray")

###FIGURE OF NMDS FOR MANUSCRIPT
plot_ordination(stF.prop, ord.nmds.braystF, title="Bray NMDS",color="collection_site")

###OTHER PHYLOSEQ ANALYSES
#heatmap
ps_abundstF= filter_taxa(stF, function(x) sum (x > total*0.20)>0, TRUE)
ps_abundstF



plot_heatmap(stF, method= "NMDS", distance = "bray")

plot_heatmap(ps_abundF, method= "NMDS", distance = "bray", 
             taxa.label = "Genus", taxa.order = "Order", 
             low = "beige", high = "red", na.value = "beige")


totalF=median(sample_sums(psF))
standfF= function(x, t=total) round (t *(x/sum(x)))
ps_transformF=transform_sample_counts(psF,standfF)
ps_transformF

# raw bar chart
top20stF <- names(sort(taxa_sums(stF), decreasing=TRUE))[1:20]
ps.top20stF <- transform_sample_counts(stF, function(OTU) OTU/sum(OTU))
ps.top20stF <- prune_taxa(top20stF, ps.top20stF)

plot_bar(ps.top20stF, fill="Genus")
plot_bar(ps.top20stF, fill="Family")
plot_bar(ps.top20stF, fill="Phylum")
plot_bar(ps.top20stF, fill="Kingdom")
plot_bar(ps.top20stF, fill="Family", x= "collection_site")

##breakaway
devtools::install_github("adw96/breakaway")
library(breakaway)
devtools::install_github("adw96/DivNet")
library(DivNet)
packageVersion("breakaway")
packageVersion("DivNet")
install.packages("tidyverse")
library(tidyverse)
library(phyloseq)
data("GlobalPatterns")
GlobalPatterns %>% sample_data

psF %>% sample_data
stF %>% sample_data

reef_F= stF %>%
  subset_samples(collection_site %in% c("mermaid_reef",
                                        "sandy_cay")) %>%
  tax_glom("Order")
reef_F

reef_F2= stF %>%
  subset_samples(collection_site %in% c("mermaid_reef",
                                        "sandy_cay")) %>%
  tax_glom("Family")
reef_F2


# water = GlobalPatterns %>%
#   subset_samples(SampleType %in% c("Freshwater", 
#                                    "Freshwater (creek)",
#                                    "Ocean",
#                                    "Sediment (estuary)")) %>%
#   tax_glom("Order")
# 
# water
observed_f=sample_richness(reef_F)
summary(observed_f)
plot(observed_f, reef_F, color="sample")
data.frame("observed_richness"= (observed_f %>% summary)$estimate,
           "depth"= phyloseq::sample_sums(reef_F),
           "type" = reef_F %>% sample_data %>% get_variable("collection_site")) %>%
  ggplot(aes(x=depth, y=observed_richness, color=type))+
  geom_point()

observed_f2=sample_richness(reef_F2)
summary(observed_f2)
plot(observed_f2, reef_F2, color="sample")
data.frame("observed_richness"= (observed_f2 %>% summary)$estimate,
           "depth"= phyloseq::sample_sums(reef_F2),
           "type" = reef_F2 %>% sample_data %>% get_variable("collection_site")) %>%
  ggplot(aes(x=depth, y=observed_richness, color=type))+
  geom_point()

# observed_c=sample_richness(water)  
# summary(observed_c)
# plot(observed_c, water, color="SampleType")
# 
# data.frame("observed_richness" = (observed_c %>% summary)$estimate,
#            "depth" = phyloseq::sample_sums(water),
#            "type" = water %>% sample_data %>% get_variable("SampleType")) %>%
#   ggplot(aes(x=depth, y=observed_richness, color=type))+
#   geom_point()
fa=breakaway(reef_F)
fa
plot(fa, reef_F, color="collection_site")
plot(fa, reef_F, color="sample")

fa2=breakaway(reef_F2)
fa2
plot(fa2, reef_F2, color="collection_site")
plot(fa2, reef_F2, color="sample")


# ba=breakaway(water)
# ba
# plot(ba,water, color="SampleType")
# summary(ba) %>%
#   add_column("SampleNames"= water %>% otu_table %>% sample_names)

summary(fa)
summary(fa) %>%
  add_column("SampleNames"= reef_F %>% otu_table %>% sample_names)

summary(fa2)
summary(fa2) %>%
  add_column("SampleNames"= reef_F2 %>% otu_table %>% sample_names)


reef_F2 %>%
  chao_bunge %>%
  plot(reef_F2)
#error with alpha diversity estimation, therefore used observed alpha diversity from Phyloseq
##in lieu of the breakaway estimates

water %>%
  chao_bunge %>%
  plot(water, color="SampleType")

bt=betta(summary(ba)$estimate,
         summary(ba)$error,
         make_design_matrix(water, "SampleType"))
bt$table
summary(fa2)$estimate
summary(fa2)$error
bF <- betta(summary(fa2)$estimate,
            summary(fa2)$error,
            make_design_matrix(reef_F2,"collection_site"))
###DIVERSITY DIFFERENCES BETWEEN SANDY CAY and MERMAID###
bF$table
write.table(bF$table, "betaFsOnly.tsv", sep="\t", quote = F, col.names=NA)

library(DivNet)
#level of Order
dv_reefF <- divnet(reef_F, ncores = 1, tuning="test")
dv_reefF %>% names
dv_reefF$shannon %>%
  summary %>%
  add_column("SampleNames"= reef_F %>% otu_table %>% sample_names)

plot(reef_F %>% sample_shannon,
     reef_F,
     col="collection_site") +ylim (0,3.5)

dv_reefF_st <- reef_F %>%
  divnet(X="collection_site", ncore =8)

plot(dv_reefF_st$shannon,
     reef_F,
     col="collection_site")

testDiversity(dv_reefF_st, "shannon")

plot(dv_reefF_st$simpson, 
     reef_F, 
     col = "collection_site")
testDiversity(dv_reefF_st, "simpson")

bcF<- dv_reefF_st$'bray-curtis'
bcF
bcF %>% unique

simplifyBeta(dv_reefF_st, reef_F, "bray-curtis", "collection_site")
simplifyBeta(dv_reefF_st, reef_F, "euclidean", "collection_site")

simplifyBeta(dv_reefF_st, reef_F, "bray-curtis", "collection_site") %>%
  ggplot(aes(x= interaction(Covar1,Covar2),
             y= beta_est))+
  geom_point()+
  geom_linerange(aes(ymin=lower, ymax=upper)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab ("Estimates of Bray-Curtis distance")



##OTHER ANALYSES FOR PHYLOGENY##
#make tree
install.packages("phangorn")
packageVersion("phangorn")
library(DECIPHER)
packageVersion("DECIPHER")
seqs= getSequences(seqtabF.nochim)
seqs2=getSequences(sfF2)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
library(phangorn)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
fit= pml(treeNJ,data=phang.align)

bs <- bootstrap.pml(fit,bs=100, optNni=TRUE, multicore = TRUE)
fitGTR <- update(fit, k=4, inv= 0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic", control = 
                      pml.control(trace=0))

psFtree <- phyloseq(otu_table(seqtabF.nochim,taxa_are_rows = FALSE),
                    sample_data(sd),tax_table(taxaF),
                    phy_tree(fitGTR$tree))
psFbsTree <- phyloseq(otu_table(seqtabF.nochim,taxa_are_rows = FALSE),
                    sample_data(sd),tax_table(taxaF),
                    phy_tree(fit$tree,bs))


psFtree
psFbsTree
library(phyloseq)
library(ggplot2)
packageVersion("ggplot2")
head(phy_tree(psFtree)$node.label)
phy_tree(psFtree)$node.label = substr(phy_tree(psFtree)$node.label, 1, 4)
ntaxa(psFtree)
physeq= prune_taxa(taxa_names(psFtree)[1:50],psFtree)

plot_tree(physeq)

psFbsTree

head(phy_tree(psFbsTree)$node.label)
head(phy_tree(psFbsTree)$node)

ntaxa(psFbsTree)
physeq2= prune_taxa(taxa_names(psFbsTree)[1:50],psFbsTree)


plot_tree(physeq)
plot_tree(physeq2)
plot_tree(physeq2, "treeonly", nodeplotblank)
plot_tree(physeq2, nodelabf=nodeplotblank, label.tips="taxa_names", ladderize="left")
plot_tree(physeq, nodelabf=nodeplotboot(), ladderize="left", color="collection_site")
plot_tree(physeq2, nodelabf=nodeplotboot(), ladderize="left", color="collection_site")
plot_tree(physeq2, nodelabf=nodeplotboot(), ladderize="left",
          color="collection_site", size= "abundance",
          base.spacing = 0.07)
plot_tree(physeq, nodelabf=nodeplotboot(), ladderize="left",
          color="collection_site", size= "abundance",
          base.spacing = 0.07,label.tips = "Class",
          shape = "Kingdom")
plot_tree(physeq, nodelabf=nodeplotboot(), ladderize="left",
          color="collection_site", size= "abundance",
          base.spacing = 0.07,
          shape = "Family")
plot_tree(physeq, nodelabf=nodeplotboot(), ladderize="left",
          color="collection_site", size= "abundance",
          base.spacing = 0.07,
          shape = "Genus")
plot_tree(physeq, nodelabf=nodeplotboot(), ladderize="left",
          color="collection_site", size= "abundance",
          base.spacing = 0.07,
          shape = "Class")
plot_tree(physeq, nodelabf=nodeplotboot(), ladderize="left",
          color="collection_site", size= "abundance",
          base.spacing = 0.07,
          shape = "Phylum")

endo <- subset_taxa(psFbsTree, Family =="Endozoicomonadaceae")
plot_tree(endo, color = "collection_site", 
          shape = "Class", label.tips = "Class",
          size = "abundance", plot.margin = 0.7)
alpha <- subset_taxa(psFbsTree, Class =="Alphaproteobacteria")
plot_tree(alpha, color = "collection_site", 
          shape = "Family", label.tips = "Family",
          size = "abundance", plot.margin = 0.7)
nomito <- subset_taxa(psFbsTree, !Class =="Vibrionaceae")
gammaprot <- subset_taxa(psFbsTree, Class =="Gammaproteobacteria")
plot_tree(nomito, color = "collection_site", 
          shape = "Family", label.tips = "Class",
          size = "abundance", base.spacing = 0.1,
          plot.margin = 1)
proteos <- subset_taxa(psFbsTree, Phylum =="Proteobacteria")
plot_tree(proteos, color = "collection_site", 
          shape = "Class", label.tips = "Family",
          size = "abundance", base.spacing = 0.05)

cyanos <- subset_taxa(psFbsTree, Phylum =="Cyanobacteria", )
plot_tree(cyanos, color = "collection_site", 
          shape = "Order", label.tips = "Order",
          size = "abundance", base.spacing = 0.1,
          plot.margin = 1)

