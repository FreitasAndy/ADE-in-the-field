#### 
#### Analysis from field experiment in the Amazon Rainforest
#### Cecropia, Schizolobium, Handroanthus
#### Control x ADE
####


library(readxl)
library(ggplot2)
library(tidyverse)
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(dplyr)
library(magrittr)
library(microbiome)
library(microeco)
library(ggpubr)
library(corrplot)
library(ALDEx2)
library(data.table)
library(knitr)
library(file2meco)
theme_set(theme_bw())


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


#####################
# DADA2 for ITS
path <- "./ITS"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#plotQualityProfile(fnFs[1:6])
#plotQualityProfile(fnRs[1:6])
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out.fun <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen = 50,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=F)
head(out.fun)
derepFs <- derepFastq(filtFs, verbose=F)
derepRS <- derepFastq(filtRs, verbose = F)
errF <- learnErrors(derepFs, multithread=TRUE)
errR <- learnErrors(derepRS, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
 dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[4]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out.fun, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
taxa <- assignTaxonomy(seqtab.nochim,
                       "C:/Users/andersonfreitas/Downloads/sh_general_release_dynamic_04.04.2024.fasta",
                       multithread=TRUE, tryRC = TRUE)
# save.image(file = "./fungi_Dada.RData")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

samdf <- read_excel("./Planilha Sequenciamento.xlsx",
                    sheet = "Planilha2")
samdf <- samdf[order(samdf$SampleID, decreasing = F),]
samdf <- as.data.frame(samdf)
row.names(samdf) <- samdf$SampleID
otu <- as.data.frame(t(seqtab.nochim))
write.csv(otu, file = "./otu_fungi.csv")
write.csv(taxa, file = "./tax_fungi.csv")
ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE), 
               sample_data(samdf), 
               tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps.fun <- merge_phyloseq(ps, dna)
taxa_names(ps.fun) <- paste0("ASV", seq(ntaxa(ps)))
ps.fun



save.image(file = "./fungi_Dada_good.RData")

######################################
# DADA2 for 16S
path <- "./16S"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnRs[1:6])
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(20,20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=F)
save.image(file = "./bacteria_ate_filterandtrim.RData")
head(out)
derepFs <- derepFastq(filtFs, verbose=F)
derepRS <- derepFastq(filtRs, verbose = F)
errF <- learnErrors(filtFs, multithread=TRUE, USE_QUALS = FALSE)
errR <- learnErrors(filtRs, multithread=TRUE, USE_QUALS = FALSE) #thank you Guillermouceda on GitHub
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[4]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
taxa <- assignTaxonomy(seqtab.nochim,
                       "C:/Users/Dinos/Downloads/silva_nr99_v138.1_train_set.fa.gz",
                       multithread=TRUE, tryRC = TRUE)
taxa <- addSpecies(taxa, "C:/Users/Dinos/Downloads/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

samdf <- read_excel("./Planilha Sequenciamento.xlsx",
                    sheet = "Sheet1")
samdf <- samdf[order(samdf$SampleID, decreasing = F),]
samdf <- as.data.frame(samdf)
row.names(samdf) <- samdf$SampleID
otu <- as.data.frame(t(seqtab.nochim))
ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE), 
               sample_data(samdf), 
               tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps.bac <- merge_phyloseq(ps, dna)
taxa_names(ps.bac) <- paste0("ASV", seq(ntaxa(ps)))
ps.bac

#save.image("F:/OneDrive/Anderson-BackUp/FAPEAM/Paper5_FieldExperiment/Fungi and Bac done.RData")



##### Analysis ####

df.fun <- phyloseq2meco(ps.fun)
df.fun$tidy_dataset()
df.fun$tax_table %<>% tidy_taxonomy

df.bac <- phyloseq2meco(ps.bac)
df.bac$tidy_dataset()
df.bac$tax_table %<>% tidy_taxonomy




### Phyla Abundance

t1 <- trans_abund$new(dataset = df.bac, taxrank = "Phylum", ntaxa = 10)
t2 <- trans_abund$new(dataset = df.fun, taxrank = "Phylum", ntaxa = 10)

t1 <- trans_abund$new(dataset = df.fun, taxrank = "Phylum", ntaxa = 10, groupmean = "Treatment")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
pb.plot <- g1 + theme_classic() + theme(axis.title.y = element_text(size = 18)) +
  labs(title = "Prokaryotes", tag = "A")

t2 <- trans_abund$new(dataset = df.fun, taxrank = "Phylum", ntaxa = 10, groupmean = "Treatment")
g2 <- t2$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
pf.plot <- g2 + theme_classic() + theme(axis.title.y = element_text(size = 18)) +
  labs(title = "Fungi", tag = "B")

plot.hey <- ggarrange(pb.plot, pf.plot, ncol = 1, nrow = 2, common.legend = F, legend = "bottom")

ggsave(filename = "./Figures/phyla_abundance.svg", plot = plot.hey, device = "svg",
       dpi = 1200,  width = 17,  height = 24, units = "cm")



### Microbial Biomass

biomass2 <- read_excel("Mapping File.xlsx")
biomass  <- filter(biomass2, Plant != "Acacia" & Plant != "Cecropia")
ggplot(data = biomass, aes(x = Treatment, y = `DNA Count`, fill = Substrate)) +
  geom_boxplot() +
  theme_classic()
biomass$Treatment <- factor(biomass$Treatment,
                            levels = c("Cecropia Control", "Cecropia ADE",
                                       "Schizolobium Control", "Schizolobium ADE",
                                       "Handroanthus Control", "Handroanthus ADE")
)
df1 <- data_summary(biomass[,c(5:9)], varname="DNA Count",
                    groupnames=c("Depth", "Substrate", "Plant", "Treatment"))
df1$Plant <- factor(df1$Plant, levels = c("Schizolobium", "Handroanthus"))
df1$Substrate <- factor(df1$Substrate, levels = c("Control", "ADE"))
mb.plot = 
  ggplot(data = df1, aes(x=Plant, y=`DNA Count`, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=`DNA Count`, ymax=`DNA Count`+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Tree", y= "ng of DNA") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Estimated microbial biomass") + guides(color = "none") +
  facet_grid(.~Depth)
mb.plot 
mb.k <- rstatix::dunn_test(biomass, `DNA Count` ~ Treatment, p.adjust.method = "fdr")
mb.k

########


# Testing Normality - Handroanthus
ipe.df <- filter(data.ex, Tree == "Handroanthus")

# height
set.seed(123) # Para reprodutibilidade
dados_exp <- data.frame(
  Block = ipe.df$Block,
  Treatment = ipe.df$Treatment,
  Response = ipe.df$Height
  )
print(dados_exp)

dados_tratamento_A <- subset(dados_exp, Treatment == "Control")$Response
dados_tratamento_B <- subset(dados_exp, Treatment == "ADE")$Response

shapiro.test(dados_tratamento_A)
shapiro.test(dados_tratamento_B)

modelo_anova <- aov(Response ~ Treatment + Block, data = dados_exp)
summary(modelo_anova)
plot(modelo_anova, which = 2) #homogeneidade dos resíduos


# stem
dados_exp <- data.frame(
  Block = ipe.df$Block,
  Treatment = ipe.df$Treatment,
  Response = ipe.df$Stem
)
print(dados_exp)

dados_tratamento_A <- subset(dados_exp, Treatment == "Control")$Response
dados_tratamento_B <- subset(dados_exp, Treatment == "ADE")$Response

shapiro.test(dados_tratamento_A)
shapiro.test(dados_tratamento_B)

modelo_anova <- aov(Response ~ Treatment + Block, data = dados_exp)
summary(modelo_anova)
plot(modelo_anova, which = 2)




## Testing normality - Schizolobium

par.df <- filter(data.ex, Tree == "Schizolobium")

# height
set.seed(123) # Para reprodutibilidade
dados_exp <- data.frame(
  Block = par.df$Block,
  Treatment = par.df$Treatment,
  Response = par.df$Height
)
print(dados_exp)

dados_tratamento_A <- subset(dados_exp, Treatment == "Control")$Response
dados_tratamento_B <- subset(dados_exp, Treatment == "ADE")$Response

shapiro.test(dados_tratamento_A)
shapiro.test(dados_tratamento_B)

modelo_anova <- aov(Response ~ Treatment + Block, data = dados_exp)
summary(modelo_anova)
plot(modelo_anova, which = 2) #homogeneidade dos resíduos


# stem
dados_exp <- data.frame(
  Block = par.df$Block,
  Treatment = par.df$Treatment,
  Response = par.df$Stem
)
print(dados_exp)

dados_tratamento_A <- subset(dados_exp, Treatment == "Control")$Response
dados_tratamento_B <- subset(dados_exp, Treatment == "ADE")$Response

shapiro.test(dados_tratamento_A)
shapiro.test(dados_tratamento_B)

modelo_anova <- aov(Response ~ Treatment + Block, data = dados_exp)
summary(modelo_anova)
plot(modelo_anova, which = 2)



## Plot - Height
df1 <- data_summary(data.ex[,c(5,7,8)], varname="Height",
                    groupnames=c("Tree", "Treatment"))
df1$Treatment <- factor(df1$Treatment, levels = c("Control", "ADE"))
height.plot = 
  ggplot(data = df1, aes(x=Tree, y=Height, fill = Treatment)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Height, ymax=Height+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Tree", y= "Height (cm)", tag = "A") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Plant Height") + guides(color = "none")
height.plot 


## Plot - Stem
df1 <- data_summary(data.ex[,c(6,7,8)], varname="Stem",
                    groupnames=c("Tree", "Treatment"))
df1$Treatment <- factor(df1$Treatment, levels = c("Control", "ADE"))
stem.plot = 
  ggplot(data = df1, aes(x=Tree, y=Stem, fill = Treatment)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Stem, ymax=Stem+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Tree", y= "Stem Diameter (cm)", tag = "B") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Stem Diameter") + guides(color = "none")
stem.plot 

dual.plot <- ggpubr::ggarrange(height.plot, stem.plot, ncol = 2, 
                               common.legend = TRUE, legend = "bottom")
dual.plot

ggsave(filename = "./Figures/Height_and_Stem.svg", plot = dual.plot, device = "svg",
       dpi = 1200,  width = 22,  height = 12, units = "cm")



### Alpha diversity
#ps.bac
#ps.fun


#Bacteria - Observed
set.seed(27)
input.bac.R <- rarefy_even_depth(ps.bac)
richness <- estimate_richness(input.bac.R)
head(richness)
head(sample_data(ps.bac))
alpha.all = cbind(richness, as.data.frame(sample_data(input.bac.R)))
div.stat <- alpha.all %>%
  group_by(Plant) %>%
  dunn_test(Observed ~ Substrate)
div.stat #no diff in Observed

div.stat <- alpha.all %>%
  group_by(Plant) %>%
  dunn_test(Simpson ~ Substrate)
div.stat #no diff in Observed

df5 <- data_summary(alpha.all, varname="Observed",
                    groupnames=c("Substrate", "Plant", "Treatment"))
head(df5)
df5$Substrate <- factor(df5$Substrate, levels = c("Control", "ADE"))
alphabac.plot = 
  ggplot(data = df5, aes(x=Plant, y=Observed, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Observed, ymax=Observed+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "", y= "Number of Observed Taxa", tag = "A") +
  theme(axis.text.x = element_text(angle = -15, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Observed Diversity - Bacteria") + guides(color = "none") +
  theme_classic()
alphabac.plot


#Bacteria InvSimpson
df5 <- data_summary(alpha.all, varname="InvSimpson",
                    groupnames=c("Substrate", "Plant", "Treatment"))
head(df5)
df5$Substrate <- factor(df5$Substrate, levels = c("Control", "ADE"))
alphabac2.plot = 
  ggplot(data = df5, aes(x=Plant, y=InvSimpson, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=InvSimpson, ymax=InvSimpson+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "", y= "Inverse Simpson Index", tag = "B") +
  theme(axis.text.x = element_text(angle = -15, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Inverse Simpson Diversity - Bacteria") + guides(color = "none") +
  theme_classic()
alphabac2.plot



#Fungi Observed
set.seed(27)
input.fun.R <- rarefy_even_depth(ps.fun)
richness <- estimate_richness(input.fun.R)
head(richness)
head(sample_data(ps.fun))
alpha.all = cbind(richness, as.data.frame(sample_data(input.fun.R)))
div.stat <- alpha.all %>%
  group_by(Plant) %>%
  dunn_test(Observed ~ Substrate)
div.stat #no diff in Observed

div.stat <- alpha.all %>%
  group_by(Plant) %>%
  dunn_test(InvSimpson ~ Substrate)
div.stat #no diff in Observed

df5 <- data_summary(alpha.all, varname="Observed",
                    groupnames=c("Substrate", "Plant", "Treatment"))
head(df5)
df5$Substrate <- factor(df5$Substrate, levels = c("Control", "ADE"))
alphafun.plot = 
  ggplot(data = df5, aes(x=Plant, y=Observed, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Observed, ymax=Observed+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "", y= "Number of Observed Taxa", tag = "A") +
  theme(axis.text.x = element_text(angle = -15, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Observed Diversity - Fungi") + guides(color = "none") +
  theme_classic()
alphafun.plot


#Fungi InvSimpson
df5 <- data_summary(alpha.all, varname="InvSimpson",
                    groupnames=c("Substrate", "Plant", "Treatment"))
head(df5)
df5$Substrate <- factor(df5$Substrate, levels = c("Control", "ADE"))
alphafun2.plot = 
  ggplot(data = df5, aes(x=Plant, y=InvSimpson, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=InvSimpson, ymax=InvSimpson+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "", y= "Inverse Simpson Index", tag = "B") +
  theme(axis.text.x = element_text(angle = -15, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Inverse Simpson Diversity - Fungi") + guides(color = "none") +
  theme_classic()
alphafun2.plot


alpha.t <- ggpubr::ggarrange(alphabac.plot, alphafun.plot, alphabac2.plot, alphafun2.plot,
                             ncol = 2, nrow = 2, common.legend = TRUE, align = "hv", legend = "bottom")
alpha.t

ggsave(filename = "./Figures/AlphaDiversity.svg", plot = alpha.t, device = "svg",
       dpi = 1200,  width = 20,  height = 20, units = "cm")



### Beta Diversity - Bacteria

inputbac.clr = microbiome::transform(ps.bac, "clr")
df        = as(sample_data(inputbac.clr), "data.frame")
ds        = phyloseq::distance(inputbac.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Plant*Substrate, data = df, permutations = 999)
permanova
#And the plot
ze=c("#79b643", "#bb59bd", "#52ac7a", "#cf426b", "#4dadcf","#cf5234",
     "#7475cb", "#cea13d", "#be6c91", "#748037", "#be774e")

input_ord = ordinate(inputbac.clr, "NMDS" , "euclidean") 
p4 = plot_ordination(inputbac.clr, input_ord, color = "Treatment", shape = "Plant")
p1.ctrl.ade = p4 + geom_point(aes(shape = Plant), size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = -45, y = -45, hjust = 0.2 , 
           label = bquote(''~R^2~'= 0.14  |  p = 0.001'), size = 4)+
  annotate("text", x = -45, y = -48, hjust = 0.2 , 
           label = bquote('Stress = 0.156'), size = 3)+
  scale_fill_manual(values = ze) +
  stat_ellipse() +
  theme(axis.text.x = element_text(size = 20))+
  theme_bw() 
p1.ctrl.ade

ggplot2::ggsave(filename = "./Figures/Beta_euclidean_NMDS.svg", plot = p1.ctrl.ade,
                device = "svg", dpi = 1200,  width = 20,  height = 13,  units = "cm")




### Beta Diversity - Fungi

inputfun.clr = microbiome::transform(ps.fun, "clr")
df        = as(sample_data(inputfun.clr), "data.frame")
ds        = phyloseq::distance(inputfun.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Plant*Substrate, data = df, permutations = 999)
permanova
#And the plot
ze=c("#79b643", "#bb59bd", "#52ac7a", "#cf426b", "#4dadcf","#cf5234",
     "#7475cb", "#cea13d", "#be6c91", "#748037", "#be774e")

input_ord = ordinate(inputfun.clr, "NMDS" , "euclidean") 
p4 = plot_ordination(inputfun.clr, input_ord, color = "Treatment", shape = "Plant")
p1.ctrl.ade = p4 + geom_point(aes(shape = Plant), size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = 10, y = 77, hjust = 0.2 , 
           label = bquote(''~R^2~'= 0.35  |  p = 0.001'), size = 4)+
  annotate("text", x = 10, y = 70, hjust = 0.2 , 
           label = bquote('Stress = 0.085'), size = 3)+
  scale_fill_manual(values = ze) +
  stat_ellipse() +
  theme(axis.text.x = element_text(size = 20))+
  theme_bw() 
p1.ctrl.ade

ggplot2::ggsave(filename = "./Figures/Beta_Fungi_euclidean_NMDS.svg", plot = p1.ctrl.ade,
                device = "svg", dpi = 1200,  width = 20,  height = 13,  units = "cm")




input.bac.agg <- microbiome::aggregate_rare(ps.bac, level = "Family", prevalence = 1/1000, detection = 1/1000)
input.fun.agg <- microbiome::aggregate_rare(ps.fun, level = "Family", prevalence = 1/1000, detection = 1/1000)


### ALDEx2 - Bacteria - Handroanthus
bac.ipe = subset_samples(input.bac.agg, Plant == "Handroanthus")
mi         = as.data.frame((otu_table(bac.ipe)))
var        = sample_data(bac.ipe)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
bac.ipe.out= x[(x$we.eBH<="0.05"),]
#View(bac.ipe.out)
#colnames(bac.ipe.out)
filt       = cbind(as.data.frame(row.names(bac.ipe.out)),bac.ipe.out[c(2,3,21,22,24)])
colnames(filt)
colnames(filt)[which(names(filt) == "row.names(bac.ipe.out)")] <- "Family"
write.csv(x = filt, file = "./Figures/Sheets/ALDEx2_Family_Handroanthus_Bacteria.csv")
#tax_table(ps.fun)[c(1539,1672),]


### ALDEx2 - Bacteria - Schizolobium
bac.par = subset_samples(input.bac.agg, Plant == "Schizolobium")
mi         = as.data.frame((otu_table(bac.par)))
var        = sample_data(bac.par)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
bac.par.out= x[(x$we.eBH<="0.05"),]
#View(bac.par.out)
#colnames(bac.par.out)
filt       = cbind(as.data.frame(row.names(bac.par.out)),bac.par.out[c(2,3,21,22,24)])
colnames(filt)
colnames(filt)[which(names(filt) == "row.names(bac.par.out)")] <- "Family"
write.csv(x = filt, file = "./Figures/Sheets/ALDEx2_Family_Handroanthus_Bacteria.csv")
#tax_table(ps.bac)[c(1539,1672),] 
#not significant



### ALDEx2 - Fungi - Schizolobium
fun.par = subset_samples(input.fun.agg, Plant == "Schizolobium")
mi         = as.data.frame((otu_table(fun.par)))
var        = sample_data(fun.par)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
fun.par.out= x[(x$we.eBH<="0.05"),]
#View(fun.par.out)
#colnames(fun.par.out)
filt       = cbind(as.data.frame(row.names(fun.par.out)),fun.par.out[c(2,3,23,24,26)])
colnames(filt)
colnames(filt)[which(names(filt) == "row.names(fun.par.out)")] <- "Family"
write.csv(x = filt, file = "./Figures/Sheets/ALDEx2_Family_Schizolobium_Fungi.csv")
#tax_table(ps.fun)[c(4,13,36,41,51,80,120,137,196,218,253,311,364,544,708),] 



### ALDEx2 - Fungi - Handroanthus
fun.ipe = subset_samples(input.fun.agg, Plant == "Handroanthus")
mi         = as.data.frame((otu_table(fun.ipe)))
var        = sample_data(fun.ipe)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
fun.ipe.out= x[(x$we.eBH<="0.05"),]
#View(fun.ipe.out)
#colnames(fun.ipe.out)
filt       = cbind(as.data.frame(row.names(fun.ipe.out)),fun.ipe.out[c(2,3,21,22,24)])
colnames(filt)
colnames(filt)[which(names(filt) == "row.names(fun.ipe.out)")] <- "Family"
write.csv(x = filt, file = "./Figures/Sheets/ALDEx2_Family_Handroanthus_Fungi.csv")
#write.csv(x = tax_table(ps.fun)[c(4,6,7,9,12,14,24,27,29,30,33,  
#37,38,41,42,44,45,52,54,55,56,57,  
#59,60,76,80,104,118,124,131,145,148,176,184,189, 
#194,227,256,309,312,316,319,321,363,372,377, 
#382,397,406,409,419,430,445,453,462,470,488, 
#491,538,592,631,650,658,676,711,722,746,822, 
#835,847,877,907,949,950,1364),], file = "./Figures/Sheets/Tax_Ipe_Fungi.csv")




##                ##
## ALDEx2 - Genus ##
##                ##


bacgen.agg <- microbiome::aggregate_rare(ps.bac, level = "Genus", prevalence = 1/1000, detection = 1/1000)
bacfun.agg <- microbiome::aggregate_rare(ps.fun, level = "Genus", prevalence = 1/1000, detection = 1/1000)


### ALDEx2 - Bacteria - Handroanthus
bac.ipe = subset_samples(bacgen.agg, Plant == "Handroanthus")
mi         = as.data.frame((otu_table(bac.ipe)))
var        = sample_data(bac.ipe)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
bac.ipe.out= x[(x$we.eBH<="0.05"),]
#View(bac.ipe.out)
#colnames(bac.ipe.out)
filt       = cbind(as.data.frame(row.names(bac.ipe.out)),bac.ipe.out[c(2,3,21,22,24)])
colnames(filt)
colnames(filt)[which(names(filt) == "row.names(bac.ipe.out)")] <- "Genus"
write.csv(x = filt, file = "./Figures/Sheets/ALDEx2_Genus_Handroanthus_Bacteria.csv")



### ALDEx2 - Bacteria - Schizolobium
bac.par = subset_samples(bacgen.agg, Plant == "Schizolobium")
mi         = as.data.frame((otu_table(bac.par)))
var        = sample_data(bac.par)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
bac.par.out= x[(x$we.eBH<="0.05"),]
#View(bac.par.out)
#colnames(bac.par.out)
filt       = cbind(as.data.frame(row.names(bac.par.out)),bac.par.out[c(2,3,21,22,24)])
colnames(filt)
colnames(filt)[which(names(filt) == "row.names(bac.par.out)")] <- "Genus"
write.csv(x = filt, file = "./Figures/Sheets/ALDEx2_Genus_Schizolobium_Bacteria.csv")



### ALDEx2 - Fungi - Schizolobium
fun.par = subset_samples(bacfun.agg, Plant == "Schizolobium")
mi         = as.data.frame((otu_table(fun.par)))
var        = sample_data(fun.par)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
fun.par.out= x[(x$we.eBH<="0.05"),]
#View(fun.par.out)
#colnames(fun.par.out)
filt       = cbind(as.data.frame(row.names(fun.par.out)),fun.par.out[c(2,3,23,24,26)])
colnames(filt)
colnames(filt)[which(names(filt) == "row.names(fun.par.out)")] <- "Genus"
write.csv(x = filt, file = "./Figures/Sheets/ALDEx2_Genus_Schizolobium_Fungi.csv")


### ALDEx2 - Fungi - Handroanthus
fun.ipe = subset_samples(bacfun.agg, Plant == "Handroanthus")
mi         = as.data.frame((otu_table(fun.ipe)))
var        = sample_data(fun.ipe)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
fun.ipe.out= x[(x$we.eBH<="0.05"),]
#View(fun.ipe.out)
#colnames(fun.ipe.out)
filt       = cbind(as.data.frame(row.names(fun.ipe.out)),fun.ipe.out[c(2,3,21,22,24)])
colnames(filt)
colnames(filt)[which(names(filt) == "row.names(fun.ipe.out)")] <- "Family"
write.csv(x = filt, file = "./Figures/Sheets/ALDEx2_Genus_Handroanthus_Fungi.csv")


ps.bac
ps.fun


### FAPROTAX

final.r6 <- phyloseq2meco(ps.bac)
final.r6
t1 <- trans_func$new(final.r6)
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = TRUE)
genes = cbind(as.data.frame(t1$res_spe_func_perc),sample_data(ps.bac))

genes[is.na(genes)] <- 0
colnames(genes)
my.variables <- genes[c(1:48)]
for(i in 1:length(my.variables)) {
  aaa = kruskal.test(x = as.matrix(my.variables[,i]), g = genes$Treatment)
  if(aaa$p.value <= 0.05) {
    print(colnames(my.variables[i]))
    print(aaa)
  } else {
    print(colnames(my.variables[i]))
    print("No difference")
  }
}

#Writing table
write.csv(genes, file = "./Figures/Sheets/genesbac.csv")


## FUNGUILD

final.r6 <- phyloseq2meco(ps.fun)
final.r6
t1 <- trans_func$new(final.r6)
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = TRUE)
genes = cbind(as.data.frame(t1$res_spe_func_perc),sample_data(ps.fun))

genes[is.na(genes)] <- 0
colnames(genes)
my.variables <- genes[c(1:24)]
for(i in 1:length(my.variables)) {
  aaa = kruskal.test(x = as.matrix(my.variables[,i]), g = genes$Treatment)
  if(aaa$p.value <= 0.05) {
    print(colnames(my.variables[i]))
    print(aaa)
  } else {
    print(colnames(my.variables[i]))
    print("No difference")
  }
}

#Writing table
write.csv(genes, file = "./Figures/Sheets/genesfun.csv")



### NETWORKS

#Joining ASV tables
taxa_names(ps.fun) <- paste("Fun_", taxa_names(ps.fun), sep="")

otu1 <- otu_table(ps.bac)
otu2 <- otu_table(ps.fun)
otu1 <- otu1[, order(sample_names(ps.bac))]
otu2 <- otu2[, order(sample_names(ps.fun))]

otu_combined <- rbind(otu1, otu2)

taxa1 <- tax_table(ps.bac)
taxa2 <- tax_table(ps.fun)

tax_combined <- rbind(taxa1, taxa2)

sample_data_combined <- sample_data(ps.fun)  # same as psf

ps_combined <- phyloseq(otu_table(otu_combined, taxa_are_rows=TRUE),
                        tax_table(tax_combined),
                        sample_data_combined)
ps.gen <- aggregate_rare(ps_combined, level = "Genus", prevalence = 1/1000, detection = 1/1000)



#Network for Schizolobium Control
Schi.ctrl <- subset_samples(ps.gen, Treatment == "Schizolobium Control")
Schi.ctrl <- phyloseq2meco(Schi.ctrl)
Schi.ctrl$tax_table %<>% tidy_taxonomy
Schi.ctrl$cal_abund()
t1 <- trans_network$new(dataset = Schi.ctrl, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.001)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
Schi.ctrl.csv = t1$res_network_attr
write.csv(Schi.ctrl.csv, "./Figures/Schi_Ctrl_atributes.csv")
t1$save_network(filepath = "./Figures/Schi_Control_network.gexf")



#Network for Schizolobium ADE
Schi.ade <- subset_samples(ps.gen, Treatment == "Schizolobium ADE")
Schi.ade <- phyloseq2meco(Schi.ade)
Schi.ade$tax_table %<>% tidy_taxonomy
Schi.ade$cal_abund()
t1 <- trans_network$new(dataset = Schi.ade, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.001)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
Schi.ade.csv = t1$res_network_attr
write.csv(Schi.ade.csv, "./Figures/Schi_ade_atributes.csv")
t1$save_network(filepath = "./Figures/Schi_ADE_network.gexf")



#Network for Handroanthus Control
Ipe.ctrl <- subset_samples(ps.gen, Treatment == "Handroanthus Control")
Ipe.ctrl <- phyloseq2meco(Ipe.ctrl)
Ipe.ctrl$tax_table %<>% tidy_taxonomy
Ipe.ctrl$cal_abund()
t1 <- trans_network$new(dataset = Ipe.ctrl, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.001)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
Ipe.ctrl.csv = t1$res_network_attr
write.csv(Ipe.ctrl.csv, "./Figures/Ipe_Ctrl_atributes.csv")
t1$save_network(filepath = "./Figures/Ipe_Control_network.gexf")



#Network for Handroanthus ADE
Ipe.ade <- subset_samples(ps.gen, Treatment == "Handroanthus ADE")
Ipe.ade <- phyloseq2meco(Ipe.ade)
Ipe.ade$tax_table %<>% tidy_taxonomy
Ipe.ade$cal_abund()
t1 <- trans_network$new(dataset = Ipe.ade, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.001)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
Ipe.ade.csv = t1$res_network_attr
write.csv(Ipe.ade.csv, "./Figures/Ipe_ade_atributes.csv")
t1$save_network(filepath = "./Figures/Ipe_ADE_network.gexf")


##
#Kruskal Chemistry
chem <- read_excel("C:/Users/andersonfreitas/OneDrive/Anderson-BackUp/FAPEAM/Paper2_FieldPlaces/Chemical_Composition.xlsx", 
                   sheet = "Planilha2")
chem[is.na(chem)] <- 0
colnames(chem)
my.variables <- chem[c(4:20)]
for(i in 1:length(my.variables)) {
  aaa = kruskal.test(x = as.matrix(my.variables[,i]), g = chem$Source)
  if(aaa$p.value <= 0.05) {
    print(colnames(my.variables[i]))
    print(aaa)
  } else {
    print(colnames(my.variables[i]))
    print("No difference")
  }
}
