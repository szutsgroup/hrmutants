{\rtf1\ansi\ansicpg1252\cocoartf1404\cocoasubrtf470
{\fonttbl\f0\fnil\fcharset0 Menlo-Regular;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww19860\viewh16600\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Isomut postprocessing\
\
# Aim: find a Fischer score tresholds where at most 5 SNVs, at most 1 deletions and at most 1 insertions are left in any of the starting clones\
\
library(tidyverse)\
library(BSgenome.Ggallus.UCSC.galGal4)\
\
\
#treatments\
treat_list <- list(starting_clone = c("HR5", "HR21", "HR1", "HR9", "HR29", \
                                   "HR33", "HR13", "HR17", "HR37", "HR25", \
                                   "HR41"),\
                             mock = c("HR6", "HR7", "HR8", "HR22", "HR23", \
                                      "HR24", "HR2", "HR2", "HR3", "HR10", \
                                      "HR11", "HR12", "HR30", "HR31", "HR32", \
                                      "HR34", "HR35", "HR36", "HR14", "HR15", \
                                      "HR16", "HR18", "HR19", "HR20", "HR38", \
                                      "HR39", "HR40", "HR26", "HR27", "HR28", \
                                      "HR42", "HR43", "HR44"))\
#genotypes\
geno_list <- list(WT = c("HR1", "HR2", "HR3", "HR4"),\
                BRCA1 = c("HR5", "HR6", "HR7", "HR8"),\
                BRCA2 = c("HR21", "HR22", "HR23", "HR24"),\
              RAD51C = c("HR9", "HR10", "HR11", "HR12"),\
               RAD52 = c("HR29", "HR30", "HR31", "HR32"),\
               RAD54 = c("HR33", "HR34", "HR35", "HR36"),\
               XRCC2 = c("HR13", "HR14", "HR15", "HR16"),\
               XRCC3 = c("HR17", "HR18", "HR19", "HR20"),\
                 ATM = c("HR37", "HR38", "HR39", "HR40"),\
                CHK2 = c("HR41", "HR42", "HR43", "HR44"),\
               PALB2 = c("HR25", "HR26", "HR27", "HR28"))\
\
\
#reading in isomut SNV output file\
snv <- read.delim("isomut_output/all_SNVs.isomut")\
snv$sample <- substr(snv$X.sample_name, 0, 5)\
\
#cumulative plot: SNV counts vs. Fischer score\
csnv <- matrix(0, ncol = length(table(snv$sample)), nrow = 301)\
for (i in 1:301) \{ \
    csnv[i,] <- as.vector(tabulate(factor(filter(snv, score > (i-1)/10)$X.sample_name, levels = names(table(snv$X.sample_name))), nbins = length(table(snv$sample)))) \
\}\
csnv <- data.frame(csnv, stringsAsFactors = FALSE)\
names(csnv) <- names(table(snv$sample))\
csnv$treshold <- seq(0, 30, by = 0.1)\
csnvm <- gather(csnv, )\
csnvm$status <- ifelse(csnvm$variable %in% treat_list$starting_clone, "starting clone", "other")\
ggplot(data = csnvm, aes(x = treshold, y = log(value, base = 10), group = variable, color = status)) + geom_line() + ylab("log10 (Mutation number) ") + xlab("Score treshold") + geom_hline(yintercept = log10(5), linetype = 2)\
\
#calculate suitable Fischer treshold\
nvScore <- csnv[ ,c(treat_list$starting_clone, "treshold")] %>% \
  apply(1, function(x) ifelse(sum(x[1:length(treat_list$starting_clone)] < 5) == length(treat_list$starting_clone), x[length(treat_list$starting_clone) + 1], NA)) %>% \
  min(na.rm = TRUE)\
\
#filtered snv set\
snv_filt <- filter(snv, score > snvScore, chr != "MT")\
snv_filt$genotype <- NA\
for (i in 1:length(geno_list)) \{\
  snv_filt[which(snv_filt$sample %in% geno_list[[i]]), "genotype"] <- names(geno_list)[i]\
\}\
snv_filt$treatment <- NA\
for (i in 1:length(treat_list)) \{\
  snv_filt[which(snv_filt$sample %in% treat_list[[i]]), "treatment"] <- names(treat_list)[i]\
\}\
\
\
\
#reading in isomut indel output file\
indel <- read.delim(file.path(outDir, "all_indels.isomut"))\
indel$sample <- substr(indel$X.sample_name, 0, 5)\
\
#Fischer scores separately for insertions and deletions\
#first insertions\
ins <- filter(indel, type == "INS")\
cins <- matrix(0, ncol = samplenum, nrow = 201)\
for (i in 1:201) \{ \
    cins[i,] <- as.vector(tabulate(factor(filter(ins, score > (i-1)/10)$sample, levels = names(table(indel$sample))), nbins = samplenum)) \
\}\
cins <- data.frame(cins, stringsAsFactors = FALSE)\
names(cins) <- names(table(indel$sample))\
cins$treshold <- seq(0, 20, by = 0.1)\
cinsm <- melt(cins, id = "treshold")\
cinsm$status <- ifelse(cinsm$variable %in% setdiff(treat_list$starting_clone, problematic), "starting clone", "other")\
ggplot(data = cinsm, aes(x = treshold, y = log(value, base = 10), group = variable, color = status)) + geom_line() + ylab("log10 (Insertion number)") + xlab("Score treshold") + geom_hline(yintercept = log10(1), linetype = 2)\
\
#deletions\
del <- filter(indel, type == "DEL")\
cdel <- matrix(0, ncol = samplenum, nrow = 201)\
for (i in 1:201) \{ \
    cdel[i,] <- as.vector(tabulate(factor(filter(del, score > (i-1)/10)$sample, levels = names(table(indel$sample))), nbins = samplenum)) \
\}\
cdel <- data.frame(cdel, stringsAsFactors = FALSE)\
names(cdel) <- names(table(indel$sample))\
cdel$treshold <- seq(0, 20, by = 0.1)\
cdelm <- melt(cdel, id = "treshold")\
cdelm$status <- ifelse(cdelm$variable %in% setdiff(treat_list$starting_clone, problematic), "starting clone", "other")\
ggplot(data = cdelm, aes(x = treshold, y = log(value, base = 10), group = variable, color = status)) + geom_line() + ylab("log10 (Deletion number)") + xlab("Score treshold") + geom_hline(yintercept = log10(1), linetype = 2)\
\
\
write.table(snv_filt, file = "all_SNVs_postprocessed.isomut", sep = "\\t", row.names = FALSE, quote = FALSE)\
write.table(del_filt, file = "all_DELs_postprocessed.isomut", sep = "\\t", row.names = FALSE, quote = FALSE)\
write.table(ins_filt, file = "all_INSs_postprocessed.isomut", sep = "\\t", row.names = FALSE, quote = FALSE)\
\
}
