###
### Final analysis/visualization script for Jarvis et al., 2017, JNCI-CS
###

### ANALYSIS: Mutational Signatures
require(deconstructSigs)
require(ggplot2)

require(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38

#1. Read in and start formatting the dataframe
cosmic_mut_sort <- with(cosmic_mut, cosmic_mut[order(cosmic_mut[,"sample"]),])
rownames(cosmic_mut_sort) <- NULL
cosmic_mut_sort$sample <- as.factor(cosmic_mut_sort$sample)

deconstructSigs_input <- cosmic_mut_all_sort[,c(1:2,6,16)]
deconstructSigs_input$ref <- substr(deconstructSigs_input$mut, 1, 1) 
deconstructSigs_input$alt <- substr(deconstructSigs_input$mut, 3, 3) 
deconstructSigs_input <- subset(deconstructSigs_input, select = c("chr", "pos", "ref", "alt", "sample"))

#1a. Due to the size of the dataframe, we must split it into 3 discrete sections so the matricies can be created
# efficiently. Be sure not to split up mutations within a cell line between different files. 
# Beyond that requirement, breakpoints for files are arbitrary. 
deconstructSigs_input_1 <- deconstructSigs_input[c(1:224680),]
deconstructSigs_input_2 <- deconstructSigs_input[c(224681:443115),]
deconstructSigs_input_3 <- deconstructSigs_input[c(443116:663075),]

#2. Create nx96 matrix, mapping number of trinucleotide muts to each sample (cell line)
mut.counts_1 <- mut.to.sigs.input(mut.ref = deconstructSigs_input_1, sample.id = "sample", chr = "chr", pos = "pos", ref = "ref", alt = "alt", bsg = hg38)
mut.counts_2 <- mut.to.sigs.input(mut.ref = deconstructSigs_input_2, sample.id = "sample", chr = "chr", pos = "pos", ref = "ref", alt = "alt", bsg = hg38)
mut.counts_3 <- mut.to.sigs.input(mut.ref = deconstructSigs_input_3, sample.id = "sample", chr = "chr", pos = "pos", ref = "ref", alt = "alt", bsg = hg38)

mut.counts_12 <- rbind(mut.counts_1, mut.counts_2)
mut.counts <- rbind(mut.counts_12, mut.counts_3)

#3. Match sample mutations to known signature mutational profiles:
#3a. Load the reference signature file first
signatures.nature2013 <- load("~/Desktop/RCRH Sequence Analysis/signatures.nature2013.rda") #If this doesn't work, load from the 'Files' tab in the view panel (if file is in the wd)

#3b. Get sigature context for file
context <- getTriContextFraction(mut.counts.ref = mut.counts, trimer.counts.method = "default")
context$tca_tct <- context[,"T[C>T]A"] + context[,"T[C>T]T"] + context[,"T[C>G]A"] + context[,"T[C>G]T"]
context$sample <- rownames(context)
tca_tct <- subset(context, select = c("sample", "tca_tct"))
rownames(tca_tct) <- NULL
context$sample <- NULL
context$tca_tct <- NULL

#4. Create a function to write output.sigs for every cell line into one table
#4a. Split the dataframe into individual file
rm(output.sigs.final)
output.sigs.final <- as.data.frame(whichSignatures(context,
                                                   sample.id = "ZR-75-30",
                                                   signatures.cosmic,
                                                   contexts.needed = F))
for(i in (1:nrow(context))) {
  output.sigs <- as.data.frame(whichSignatures(context,
                                               sample.id = rownames(context[i,]),
                                               signatures.cosmic,
                                               contexts.needed = F))
  output.sigs.final <- rbind(output.sigs.final, output.sigs)
}


output.sigs.final <- output.sigs.final[-c(1021),]
output.sigs.final$zAPOBEC.Sig <- output.sigs.final$weights.Signature.2 + output.sigs.final$weights.Signature.13

output.sigs.final <- output.sigs.final[,c(1:30,319,320)]
output.sigs.final$sample <- rownames(output.sigs.final)
rownames(output.sigs.final) <- NULL

sigs_tissues <- merge(output.sigs.final, cell_line_mutload, by = "sample") 
sigs_tissues <- sigs_tissues[,-c(3,14,34)]

#5. Each mutation plot was created separatly by replaceing the tissue variable name in line 79, 
# and running through the ggplot command (through line 134). Note that a lettering scheme is applied to individuals 
# signatures to allow for correct sorting.
sigs_individual <- subset(sigs_tissues, tissue == "large_intestine")
sigs_individual <- sigs_individual[,-c(32)]

sigs_melt <- melt(sigs_individual, id = "sample")
colnames(sigs_melt) <- c("sample", "sig", "value")
sigs_melt[,"sig"] <- gsub("weights.", "", sigs_melt[,"sig"])

sigs_melt[,"sig"] <- gsub("Signature.10", "I", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.11", "J", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.12", "K", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.14", "L", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.15", "M", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.16", "N", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.17", "O", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.18", "P", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.19", "Q", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.20", "R", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.21", "S", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.22", "T", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.23", "U", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.24", "V", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.25", "W", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.26", "X", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.27", "Y", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.28", "Z", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.29", "ZZ", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.30", "ZZZ", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("unknown", "ZZZZ", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("zAPOBEC.Sig", "ZZZZZ", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.1", "A", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.3", "B", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.4", "C", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.5", "D", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.6", "E", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.7", "F", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.8", "G", sigs_melt[,"sig"])
sigs_melt[,"sig"] <- gsub("Signature.9", "H", sigs_melt[,"sig"])

list <- sigs_individual[order(sigs_individual$zAPOBEC.Sig),] 
list1 <- as.vector(list[,"sample"])

# Bar plots for mutational signature proportion (Used in FIGURE_2)
#
ggplot(sigs_melt, aes(sample, value, fill = sig)) +
  geom_col() +
  #scale_fill_brewer() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  xlab("Cancer Cell Line") +
  ylab("Mutational Signature Proportion") +
  scale_x_discrete(limits = list1) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

### Enrichment Score Calculations (uses file created by the "count_trinuc_muts.pl" script)
#
library(stringr)
library(plyr)

cosmic_mut_all_sort <- read.table(file = "cosmic_mut_all_sort.txt", header = T, sep = "\t", stringsAsFactors = T)
colnames(test) <- c("chr", "pos", "5_tetnuc", "3_tetnuc", "trinuc", "mut", "trinuc_mut", "strand", "context", "C_count", "TC_count", "TCA_count", "TCT_count", "YTCA_count", "RTCA_count", "sample")

enrich_pre <- rbind(enrich1, enrich2)
enrich_post <- rbind(enrich3, enrich4)
enrich_tot <- rbind(enrich_pre, enrich_post)
enrich_tot <- unique(enrich_tot)

enrich_tot$Mut_TCW <- "0"
enrich_tot$Mut_C <- "0"
enrich_tot$Con_TCW <- "0"
enrich_tot$Con_C <- "0"

mut <- data.frame(do.call('rbind', strsplit(as.character(enrich_tot$mut),'>',fixed=T)))
enrich_tot$mut_ref <- mut[,1]

enrich_C <- subset(enrich_tot, mut_ref == "C")
enrich_CtoK <- subset(enrich_C, mut != "C>A") # Remove C>A mutations!

# Mut_C
enrich_CtoK[which(enrich_CtoK$mut_ref == "C"),"Mut_C"] <- "1"

# Mut_TCW
enrich_CtoK[which(enrich_CtoK$trinuc_mut == "T[C>G]A"),"Mut_TCW"] <- "1"
enrich_CtoK[which(enrich_CtoK$trinuc_mut == "T[C>G]T"),"Mut_TCW"] <- "1"
enrich_CtoK[which(enrich_CtoK$trinuc_mut == "T[C>T]A"),"Mut_TCW"] <- "1"
enrich_CtoK[which(enrich_CtoK$trinuc_mut == "T[C>T]T"),"Mut_TCW"] <- "1"

# Con_C
enrich_CtoK$Con_C <- str_count(enrich_CtoK$context, "C") + str_count(enrich_CtoK$context, "G")

# Con_TCW 
enrich_CtoK$Con_TCW <- str_count(enrich_CtoK$context, "TCA") + str_count(enrich_CtoK$context, "TCT") + str_count(enrich_CtoK$context, "TGA") + str_count(enrich_CtoK$context, "TGT")

# Aggregate
enrich_final <- enrich_CtoK[,16:20]
enrich_final$Mut_TCW <- as.integer(enrich_final$Mut_TCW)
enrich_final$Mut_C <- as.integer(enrich_final$Mut_C)

enrich_final <- ddply(enrich_final, "sample", numcolwise(sum))

rownames(enrich_final) <- enrich_final$sample
enrich_final$sample <- NULL
enrich_final$enrich_score <- (enrich_final$Mut_TCW / enrich_final$Con_TCW) / (enrich_final$Mut_C / enrich_final$Con_C)

enrich_matrix <- as.data.frame(enrich_final$Mut_TCW) 
enrich_matrix$Mut_Denom <- enrich_final$Mut_C - enrich_final$Mut_TCW
enrich_matrix$Con_TCW <- enrich_final$Con_TCW
enrich_matrix$Con_Denom <- enrich_final$Con_C - enrich_final$Con_TCW
rownames(enrich_matrix) <- rownames(enrich_final)
colnames(enrich_matrix) <- c("Mut_TCW", "Mut_Denom", "Con_TCW", "Con_Denom")

enrich_matrix <- as.matrix(enrich_matrix)

exe_fisher <- function(x) {
  m <- matrix(unlist(x), ncol = 2, nrow = 2, byrow = T)
  f <- fisher.test(m)
  return(as.data.frame(f$p.value))
}

fishers <- t(as.data.frame(apply(enrich_matrix, 1, exe_fisher)))
fishers <- as.data.frame(fishers)

enrich_final$fisher_pval <- fishers$V1
enrich_final$bh_adj_qval <- p.adjust(enrich_final$fisher_pval, method = "BH")

enrich_final$Mut_Ratio <- enrich_final$Mut_TCW/(enrich_final$Mut_C - enrich_final$Mut_TCW)
enrich_final$Con_Ratio <- enrich_final$Con_TCW/(enrich_final$Con_C - enrich_final$Con_TCW)
enrich_final[which(enrich_final$Mut_Ratio < enrich_final$Con_Ratio), "bh_adj_qval"] <- 1
enrich_final$Mut_Ratio <- NULL
enrich_final$Con_Ratio <- NULL
enrich_final$sample <- rownames(enrich_final)
rownames(enrich_final) <- NULL


# FIGURE_1: Median mutations in each cell line (uses columns 5 and 8 from the CosmicCLP_MutantExport.tsv file)
# 

# Format the tissue type info
cosmic_tissue_type <- read.table(file = "cosmic_tissue_type.txt", header = T, stringsAsFactors = F, fill = T)
cosmic_tissue_type <- cosmic_tissue_type[,c(1:2)]
colnames(cosmic_tissue_type) <- c("sample", "tissue")
cosmic_tissue_type <- unique(cosmic_tissue_type)

# Combine tissue and mutation information
cosmic_mut_tissue <- merge(cosmic_mut_all_sort, cosmic_tissue_type, by = "sample", all.x = T)
cell_line_mutload <- as.data.frame(table(cosmic_mut_tissue$sample))
colnames(cell_line_mutload) <- c("sample", "mut_tot")
cell_line_mutload <- merge(cell_line_mutload, cosmic_tissue_type, by = "sample", all.x = T)
cell_line_mutload[775,3] <- "upper_aerodigestive_tract"

# Replace the tissue variable in line 218 with each tissue type sequentially and iteratively run through line 220
# to generate a complete quantile table
mut_sub <- subset(cell_line_mutload, tissue == "endometrium")
x <- as.data.frame(t(quantile(mut_sub$mut_tot)))
mut_med_quantiles <- rbind(mut_med_quantiles, x)

rownames(mut_med_quantiles) <- a
mut_med_quantiles$tissue <- rownames(mut_med_quantiles)
colnames(mut_med_quantiles) <- c("low", "first", "med", "third", "high", "tissue")

ggplot(mut_med_quantiles, aes(tissue, med)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1), legend.position = "none") +
  geom_errorbar(aes(ymin=first, ymax=third), width=.3) +
  scale_y_continuous(limits = c(0,3501),
                     breaks = c(0,1200,2400,3600)) +
  xlab("Tissue Type") +
  ylab("Median Number of Mutations") + 
  #geom_text(aes(label = freq), vjust = -0.2) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_x_discrete(limits = c("pleura","bone","kidney","prostate","pancreas",
                              "central_nervous_system","vulva","small_intestine","soft_tissue","upper_aerodigestive_tract",
                              "autonomic_ganglia","urinary_tract","breast","thyroid","ovary",
                              "testis","NS","haematopoietic_and_lymphoid_tissue","salivary_gland","liver",
                              "biliary_tract","adrenal_gland","cervix", "stomach","large_intestine",
                              "lung","oesophagus","skin","placenta","endometrium")) +
  ggtitle("Median Mutations by Tissue Type")

# FIGURE_1: Dots for cell numbers
# 
number <- as.data.frame(table(cosmic_tissue_type$tissue))
colnames(number) <- c("tissue", "freq")
ggplot(number, aes(tissue, freq)) +
  geom_point(size = 4) +
  ylim(0,200) +
  scale_x_discrete(limits = c("pleura","bone","kidney","prostate","pancreas",
                              "central_nervous_system","vulva","small_intestine","soft_tissue","upper_aerodigestive_tract",
                              "autonomic_ganglia","urinary_tract","breast","thyroid","ovary",
                              "testis","NS","haematopoietic_and_lymphoid_tissue","salivary_gland","liver",
                              "biliary_tract","adrenal_gland","cervix", "stomach","large_intestine",
                              "lung","oesophagus","skin","placenta","endometrium"))


# FIGURE_2: Plotting mutload vs cell line order (model lines with shaded intervals)
#
sigs_tissues_individual <- subset(sigs_tissues, tissue == "large_intestine")

sigs_tissues_individual_1 <- sigs_tissues_individual[order(sigs_tissues_individual$zAPOBEC.Sig),]
rownames(sigs_tissues_individual_1) <- c(1:nrow(sigs_tissues_individual_1))
sigs_tissues_individual_1[,"order"] <- rownames(sigs_tissues_individual_1)

ggplot(sigs_tissues_individual_1, aes(as.numeric(order), mut_tot)) +
  geom_point(shape = 18, size = 4) +
  geom_smooth(span = 0.75) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  xlab("Breast Cancer Cell Line") +
  ylab("Mut Burden") +
  ylim(0,1600) +
theme_bw() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

# FIGURE_3: TCW, enrichment, and APOBEC sig correlation plots
# 
library(gridExtra)
sigs_enrich <- merge(sigs_tissues, enrich_final, by = "sample")
sigs_enrich_tcw <- merge(sigs_enrich, tca_tct, by = "sample")

# Scatterplots
x <- ggplot(sigs_enrich_tcw, aes(enrich_score, tca_tct)) +
  geom_point() +
  xlim(0,5) +
  ylim(0,0.5)
y <- ggplot(sigs_enrich_tcw, aes(zAPOBEC.Sig, tca_tct)) +
  geom_point() +
  xlim(0,0.6) +
  ylim(0,0.5)
grid.arrange(x,y, ncol = 2, nrow = 1)

# FIGURE_4: Images were created using a separate script that used values from the Supplementary Table 1 
# (counts at all trinucleotide contexts).