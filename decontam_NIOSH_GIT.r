
#Load Packages
library(DESeq2)
library(edgeR)
library(limma)
library(BiocManager)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(scales)
library(data.table)
library(pathfindR)
library(fBasics)
library(forcats)
library(maptools)
library(phyloseq)
library(vegan)
library(themetagenomics)
library(tidyr)
library(tableone)
library(decontam)
library(dplyr)
library(grid)
library(cowplot)
library(colorspace)

#Set Theme
theme<-
    theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), 
    legend.position="none")

#Choose Alpha/FDR
alpha = 0.01


##Load the files needed
file = "niosh.open.filter.biom"
map = "Map.master.R.filter.lns.a3.txt"

# Load the abundace table and mapping table 
abundance.table = import_biom(file, taxaPrefix=F)
mapping.table=sample_data(read.table(map, header=T, sep="\t", row.names=1))

#Merge abundance and meta data into a phyloseq object
lung.physeq=phyloseq(otu_table(abundance.table),tax_table(abundance.table), mapping.table)

#Make the columnnames of the pyloseq object of the phylogenetic tree
colnames(tax_table(lung.physeq))=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Load the tree file (use the unannotated.tree)
treefile = "97_otus.tree"
tree.obj = import_qiime(treefilename = treefile) 

#Now merge the three separate phyloseq objects into a single object
otu.table = merge_phyloseq(lung.physeq, mapping.table, tree.obj)

# Remove taxa with 0 abundance
otu.table = subset_taxa(otu.table, rowSums(otu_table(otu.table)) != 0)

#Keep Samples for Paper
otu.table.paper = subset_samples(otu.table, Paper_Active_BKG==1)

#Save File
save.image(file="NIOSH.decontam.RData")


#Create a TRUE/FALSE Variable for Subject and Control
sample_data(otu.table.paper)$is.neg <- sample_data(otu.table.paper)$Sample_Type == "BKG"
contamdf.prev <- isContaminant(otu.table.paper, method="prevalence", neg="is.neg",threshold=0.5)
table(contamdf.prev$contaminant)

#Normalize Count Table
normalizeSample = function(x) {
    x/sum(x)
}

#Relative Abundance Method
ps.pa <- transform_sample_counts(otu.table.paper, normalizeSample)

#Prevelance Method
#ps.pa <- transform_sample_counts(otu.table.paper, function(abund) 1*(abund>0))

#Create table of Controls and Samples
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_Type == "BKG", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_Type!= "BKG", ps.pa)

#Mean
df.pa <- data.frame(pa.pos=rowMeans(otu_table(ps.pa.pos)), pa.neg=rowMeans(otu_table(ps.pa.neg)),
                      contaminant=contamdf.prev$contaminant, row.names=rownames(contamdf.prev))

#Prevalance
#df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
#                      contaminant=contamdf.prev$contaminant, row.names=rownames(contamdf.prev))


#ADD Phylogenetic Tree (OTUID) to the table and rank order by Frequency and Prevelance 
df.pa = cbind(as(df.pa, "data.frame"), as(tax_table(otu.table.paper)[rownames(df.pa), ], "matrix"))
#Replace OTU with Taxa
df.pa$row2 <- paste(df.pa$Domain,df.pa$Phylum,df.pa$Class,df.pa$Order,df.pa$Family,df.pa$Genus,df.pa$OTU)
#Replace Spaces with .
df.pa$row2 <- gsub('\\s+', '.', df.pa$row2)
#Just Highest Level and OTU
df.pa[df.pa=="g__"]<-NA
df.pa[df.pa=="f__"]<-NA
df.pa[df.pa=="o__"]<-NA
df.pa[df.pa=="c__"]<-NA

#Fix Taxa Annotation
df.pa$gs <- ifelse(is.na(df.pa$Genus),  paste(df.pa$Family,rownames(df.pa),sep="_"), paste(df.pa$Genus,rownames(df.pa),sep="_"))
df.pa$gs <- ifelse(is.na(df.pa$Family), paste(df.pa$Order,rownames(df.pa),sep="_"),df.pa$gs)
df.pa$gs <- ifelse(is.na(df.pa$Order),  paste(df.pa$Class,rownames(df.pa),sep="_"),df.pa$gs)
df.pa$gs <- ifelse(is.na(df.pa$Class),  paste(df.pa$Phylum,rownames(df.pa),sep="_"),df.pa$gs)


#For Mean
df.pa <- df.pa[df.pa$pa.pos > 10^-9,]
df.pa <- df.pa[df.pa$pa.neg > 10^-8,]

#Mean
pdf(file="NIOSH_Mean_Rel_Abundance_of_Contam_or_NoContam_PREV_edit.pdf", width=5, height=5)
    ggplot(data=df.pa, aes(y=pa.neg, x=pa.pos, color=contaminant)) + 
    #ggplot(data=df.pa, aes(y=neg, x=pos, color=contaminant)) + 
    #geom_point(aes(size =size), alpha=0.65,shape=19) + #Chose Colors and size for dots
    geom_point(alpha=0.65,shape=19) + #Chose Colors and size for dots
    #geom_text_repel(aes(label=ifelse(df.pa$pa.neg > 25 | df.pa$pa.pos > 150, as.character(df.pa$gs),''),size=abundance),force=25,segment.colour="grey",segment.alpha=0.5) +
    #geom_text_repel(aes(label=ifelse(df.pa$pa.neg > 10^-2 | df.pa$pa.pos > 10^-2, as.character(df.pa$gs),'')),size=2,force=25,segment.colour="grey",segment.alpha=0.5) +
    #geom_text_repel(aes(label=ifelse(df.pa$contaminant=="TRUE" & df.pa$pa.neg > 10^-3, as.character(df.pa$gs),'')),size=2,force=25,segment.colour="grey",segment.alpha=0.5) +
    geom_text_repel(aes(label=ifelse(rownames(df.pa)==813945, as.character(df.pa$gs),'')),size=2,force=25,segment.colour="grey",segment.alpha=0.5) +
    #ylab("Rel Abundance (Negative Controls)") + xlab("Rel Abundance (True Samples)")+
    scale_color_manual(values=c("red","#3B9E55"))+
    scale_y_continuous(name="Rel Abundance (Negative Controls)",trans="log10", breaks = trans_breaks('log10', function(x) 10^x),limits=c(10^-6,10^0), labels = trans_format('log10', math_format(10^.x))) +
    scale_x_continuous(name="Rel Abundance (True Samples)",trans="log10", breaks = trans_breaks('log10', function(x) 10^x),limits=c(10^-9,10^0), labels = trans_format('log10', math_format(10^.x))) +
    #scale_y_continuous(name="Prevalence (Negative Controls)",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
    #scale_x_continuous(name="Prevalence (True Samples)",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
    #scale_size(name="Relative Abundance",
    #    breaks=fivenum(df.pa$size),
    #    labels=format((fivenum(df.pa$abundance)),scientific=TRUE))+
    #theme(legend.key=element_blank())
    theme(axis.text=element_text(size=15),panel.grid.major = element_blank(),legend.position="none", 
    panel.grid.minor = element_blank(), panel.background= element_blank(), 
    axis.line = element_line(colour = "black"), strip.text.y = element_blank(),strip.text = element_text(size = 20,face="bold"), 
    strip.background = element_blank(), axis.title=element_text(size=15,face="bold"))
dev.off()

#Prevalance ---> You need to go above and recalculate prevalance before running this figure code
pdf(file="NIOSH_Prevalance_of_Contam_or_NoContam_PREV_log.pdf", width=5, height=5)
    ggplot(data=df.pa, aes(y=pa.neg, x=pa.pos, color=contaminant)) + 
    #ggplot(data=df.pa, aes(y=neg, x=pos, color=contaminant)) + 
    #geom_point(aes(size =size), alpha=0.65,shape=19) + #Chose Colors and size for dots
    geom_point(alpha=0.65,shape=19) + #Chose Colors and size for dots
    #geom_text_repel(aes(label=ifelse(df.pa$pa.neg > 25 | df.pa$pa.pos > 150, as.character(df.pa$gs),''),size=abundance),force=25,segment.colour="grey",segment.alpha=0.5) +
    #geom_text_repel(aes(label=ifelse(df.pa$pa.neg > 10^-2 | df.pa$pa.pos > 10^-2, as.character(df.pa$gs),'')),size=2,force=25,segment.colour="grey",segment.alpha=0.5) +
    #geom_text_repel(aes(label=ifelse(df.pa$contaminant=="TRUE" & df.pa$pa.neg > 10^-3, as.character(df.pa$gs),'')),size=2,force=25,segment.colour="grey",segment.alpha=0.5) +
    geom_text_repel(aes(label=ifelse(rownames(df.pa)==813945, as.character(df.pa$gs),'')),size=2,force=25,segment.colour="grey",segment.alpha=0.5) +
    #ylab("Rel Abundance (Negative Controls)") + xlab("Rel Abundance (True Samples)")+
    scale_color_manual(values=c("red","#3B9E55"))+
    scale_y_continuous(name="Prevalence (Negative Controls)",trans="log10", breaks = trans_breaks('log10', function(x) 10^x),limits=c(10^-1,10^3), labels = trans_format('log10', math_format(10^.x))) +
    scale_x_continuous(name="Prevalence (True Samples)",trans="log10", breaks = trans_breaks('log10', function(x) 10^x),limits=c(10^-1,10^3), labels = trans_format('log10', math_format(10^.x))) +
    #scale_y_continuous(name="Prevalence (Negative Controls)",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
    #scale_x_continuous(name="Prevalence (True Samples)",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
    #xlab("Prevalence (True Samples)")+
    #ylab("Prevalence (Negative Controls)")+
    #scale_size(name="Relative Abundance",
    #    breaks=fivenum(df.pa$size),
    #    labels=format((fivenum(df.pa$abundance)),scientific=TRUE))+
    #theme(legend.key=element_blank())
    theme(axis.text=element_text(size=15),panel.grid.major = element_blank(),legend.position="none", 
    panel.grid.minor = element_blank(), panel.background= element_blank(), 
    axis.line = element_line(colour = "black"), strip.text.y = element_blank(),strip.text = element_text(size = 20,face="bold"), 
    strip.background = element_blank(), axis.title=element_text(size=15,face="bold"))
dev.off()


