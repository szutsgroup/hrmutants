{\rtf1\ansi\ansicpg1252\cocoartf1404\cocoasubrtf470
{\fonttbl\f0\fnil\fcharset0 Menlo-Regular;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww16060\viewh17860\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 library(tidyverse)\
library(superheat)\
library(MutationalPatterns)\
library(BSgenome.Ggallus.UCSC.galGal4)\
library(BSgenome.Hsapiens.NCBI.GRCh38)\
\
\
#reading in postprocessed SNV data\
snv <- read.delim("all_SNVs_postprocessed.isomut", stringsAsFactors = FALSE)\
\
\
#creating a decoy vcf object\
vcf_list <- vector(mode = "list", length = length(unique(snv$genotype)))\
names(vcf_list) <- unique(snv$genotype)\
\
for (i in unique(snv$genotype)) \{\
  t <- filter(snv, genotype == i)\
  \
  vcf_list[[i]] <- with(t, GRanges(seqnames = paste0("chr", chr), ranges = IRanges(start = pos, end = pos), seqinfo = seqinfo(Ggallus)))\
  names(vcf_list[[i]] ) <- paste0("chr", t$chr, ":", t$pos, "_", t$ref, "/", t$mut)\
  mcols(vcf_list[[i]] )$paramRangeID <- factor(NA)\
  mcols(vcf_list[[i]] )$REF <- DNAStringSet(t$ref)\
\
  mcols(vcf_list[[i]] )$ALT <- DNAStringSetList(unname(sapply(t$mut, DNAStringSet)))\
  mcols(vcf_list[[i]] )$QUAL <- as.numeric(NA)\
  mcols(vcf_list[[i]])$FILTER <- "PASS"\
\}\
\
#estimating a suitable component number\
mut_mat <- mut_matrix(vcf_list = vcf_list, ref_genome = ref_genome)\
mut_mat <- mut_mat + 0.0001\
estimate <- nmf(mut_mat, rank=1:5, method="brunet", nrun=10, seed=123456)\
plot(estimate)\
ggsave(file = "SNV_plots/denovo_NMF/NMF_component_number_test.pdf", width = 8, height = 5)\
\
#2 components seem to be OK\
#extracting 2 components\
nmf_res_2 <- extract_signatures(mut_mat, rank = 2, nrun = 100)\
colnames(nmf_res_2$signatures) <- c("Signature A", "Signature B")\
rownames(nmf_res_2$contribution) <- c("Signature A", "Signature B")\
\
#saving some plots and numeric data\
plot_96_profile(nmf_res_2$signatures, condensed = TRUE, ymax = 0.075)\
ggsave("SNV_plots/denovo_NMF/NMF_denovo_2comp.pdf",  width = 8, height = 5)\
\
plot_contribution(nmf_res_2$contribution, nmf_res_2$signatures, mode = "absolute")\
ggsave("SNV_plots/denovo_NMF/NMF_denovo_comparison.pdf", width = 8, height = 5)\
\
write.table(mut_mat, file = "triplet_mutation_matrix.txt", sep = "\\t", quote = FALSE, row.names = FALSE)\
\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 write.table(nmf_res_2$contribution, file = "extracted_signatures.txt", sep = "\\t", quote = FALSE, row.names = FALSE)\
\
#"humanizing" the denovo components (correcting for differences between human and chicken triplet frequencies)\
humanized <- nmf_res_2$signatures / (triplets96$Gallus.gallus * sum(as.numeric(seqlengths(Ggallus)))) * (triplets96$Homo.sapiens * sum(as.numeric(seqlengths(Hsapiens))))\
\
\
#cosine similarity and Spearman correlation\
cos_sim_denovo_signatures_hum = cos_sim_matrix(humanized, cancer_signatures)\
cancer_signatures <- read.delim("https://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt", stringsAsFactors = FALSE)[,1:33]\
cos_sim_denovo_signatures_hum %>% \
    superheat( bottom.label.text.angle = 90, bottom.label.text.size = 4, title = "Cosinus similarity", bottom.label.size = .6, left.label.size = .3, left.label.text.size = 4, legend.breaks = seq(0, 1, by = .2), heat.lim = c(0, 1))\
\
cor(humanized, cancer_signatures, method = "spearman") %>% \
    superheat( bottom.label.text.angle = 90, bottom.label.text.size = 4, title = "Spearman correlation", bottom.label.size = .6, left.label.size = .3, left.label.text.size = 4)\
\
\
\
}
