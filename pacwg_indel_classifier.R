library(tidyverse)
library(pracma)

# preformatting
rbind.data.frame(del_filt, ins_filt) -> x
x$genotype <- toupper(x$genotype)
x <- filter(x, treatment != "starting_clone")
x$label <- apply(x, 1, function(x) paste0(x[12], "_", x[13]))

# creating output table
# 83 is the number of PACWG indel categories
restab <- matrix(0, ncol = length(table(x$label)), nrow = 83)
colnames(restab) <- names(table(x$label))


# for each genotype/treatment combination separately
# calculate number of indels in each subcategory, and fill in the respective column of the output table
for (gg in unique(x$label)) {
  
  filter(x, type == "DEL", ref %in% c("C", "G"), nchar(ref) == 1, !is.na(sample), label == gg) -> a
  if (nrow(a) == 0) {
    a <- data.frame(hpl = -1)
  } else {
    with(a, GRanges(paste0("chr", chr), ranges = IRanges(start = pos+1, end = pos+10))) %>% getSeq(Ggallus, .) -> seqsa
    
    a$nei <- as.character(seqsa)
    pattern_a <- paste0("^", a$ref, "*", sep = "")
    
    
    a$nei2 <- NA
    for (i in seq_along(pattern_a)) {
      a$nei2[i] <- gsub(pattern_a[i], "", a$nei[i])
    }
    a$hpl <- (nchar(a$nei) - nchar(a$nei2))
    a[which(a$hpl > 5), "hpl"] <- "6+"
  }
  a$hpl <- factor(a$hpl, levels = c(1:5, "6+"))
  
  #===============
  
  filter(x, type == "DEL", ref %in% c("A", "T"), nchar(ref) == 1, !is.na(sample),label == gg) -> b
  if (nrow(b) == 0) {
    b <- data.frame(hpl = -1)
  } else {
    with(b, GRanges(paste0("chr", chr), ranges = IRanges(start = pos+1, end = pos+10))) %>% getSeq(Ggallus, .) -> seqsb
    
    b$nei <- as.character(seqsb)
    pattern_b <- paste0("^", b$ref, "*", sep = "")
    
    
    b$nei2 <- NA
    for (i in seq_along(pattern_b)) {
      b$nei2[i] <- gsub(pattern_b[i], "", b$nei[i])
    }
    b$hpl <- (nchar(b$nei) - nchar(b$nei2))
    b[which(b$hpl > 5), "hpl"] <- "6+"
  }
  b$hpl <- factor(b$hpl, levels = c(1:5, "6+"))
  
  #===============
  
  filter(x, type == "INS", mut %in% c("C", "G"), nchar(mut) == 1, !is.na(sample), label == gg) -> c
  if (nrow(c) == 0) {
    c <- data.frame(hpl = -1)
  } else {
    with(c, GRanges(paste0("chr", chr), ranges = IRanges(start = pos +1, end = pos+10))) %>% getSeq(Ggallus, .) -> seqsc
    
    c$nei <- as.character(seqsc)
    pattern_c <- paste0("^", c$mut, "*", sep = "")
    
    
    c$nei2 <- NA
    for (i in seq_along(pattern_c)) {
      c$nei2[i] <- gsub(pattern_c[i], "", c$nei[i])
    }
    c$hpl <- (nchar(c$nei) - nchar(c$nei2))
    c[which(c$hpl > 4), "hpl"] <- "5+"
  }
  c$hpl <- factor(c$hpl, levels = c(0:4, "5+"))
  
  #===============
  
  filter(x, type == "INS", mut %in% c("A", "T"), nchar(mut) == 1, !is.na(sample), label == gg) -> d
  if (nrow(d) == 0) {
    d <- data.frame(hpl = -1)
  } else {
    with(d, GRanges(paste0("chr", chr), ranges = IRanges(start = pos +1, end = pos+10))) %>% getSeq(Ggallus, .) -> seqsd
    
    d$nei <- as.character(seqsd)
    pattern_d <- paste0("^", d$mut, "*", sep = "")
    
    
    d$nei2 <- NA
    for (i in seq_along(pattern_d)) {
      d$nei2[i] <- gsub(pattern_d[i], "", d$nei[i])
    }
    d$hpl <- (nchar(d$nei) - nchar(d$nei2))
    d[which(d$hpl > 4), "hpl"] <- "5+"
  }
  d$hpl <- factor(d$hpl, levels = c(0:4, "5+"))
  #===============
  checkHomLen <- function(deleted, next50) {
    ret <- 0
    for (i in 1:nchar(deleted)) {
      if (substr(deleted, 1, i) == substr(next50, 1, i)) ret <- i
    }
    return(ret)
  }
  
  del_context$mhl <- apply(del_context, 1, function(x) checkHomLen(x[8], x[9]))
  
  y <- del_context
  y$Genotype <- toupper(y$Genotype)
  y$label <- apply(y, 1, function(x) paste0(x[2], "_", x[3]))
  
  filter(y, Category %in% c("repeat", "no homology"), nchar(Deleted) > 1, !is.na(Sample), label == gg) -> e
  if (nrow(e) == 0) {
    e <- list(ea = data.frame(nrep = -1) %>% mutate(nrep = factor(nrep, levels = c(1:5, "6+"))),
              eb = data.frame(nrep = -1) %>% mutate(nrep = factor(nrep, levels = c(1:5, "6+"))),
              ec = data.frame(nrep = -1) %>% mutate(nrep = factor(nrep, levels = c(1:5, "6+"))),
              ed = data.frame(nrep = -1) %>% mutate(nrep = factor(nrep, levels = c(1:5, "6+"))))
  } else {
    pattern_e <- paste0("^(", e$Deleted, ")*")
    
    e$nrep <- NA
    for (i in seq_along(pattern_e)) {
      tmp <- gsub(pattern_e[i], "x", e$Next50[i])
      e$nrep[i] <- (nchar(e$Next50[i]) - (nchar(tmp) - 1)) / nchar(e$Deleted[i]) + 1
    }
    e[which(e$nrep > 5), "nrep"] <- "6+"
    
    e <- list(ea = filter(e, nchar(Deleted) == 2) %>% mutate(nrep = factor(nrep, levels = c(1:5, "6+"))),
              eb = filter(e, nchar(Deleted) == 3) %>% mutate(nrep = factor(nrep, levels = c(1:5, "6+"))),
              ec = filter(e, nchar(Deleted) == 4) %>% mutate(nrep = factor(nrep, levels = c(1:5, "6+"))),
              ed = filter(e, nchar(Deleted) >= 5) %>% mutate(nrep = factor(nrep, levels = c(1:5, "6+"))))
  }
  
  #===============
  z <- ins_context
  z$label <- apply(z, 1, function(x) paste0(x[8], "_", x[9]))
  
  filter(z,  nchar(Inserted) > 1, !is.na(Sample), label == gg ) -> f
  if (nrow(f) == 0) {
    f <- list(fa = data.frame(nrep = -1) %>% mutate(nrep = factor(nrep, levels = c(0:4, "5+"))),
              fb = data.frame(nrep = -1) %>% mutate(nrep = factor(nrep, levels = c(0:4, "5+"))),
              fc = data.frame(nrep = -1) %>% mutate(nrep = factor(nrep, levels = c(0:4, "5+"))),
              fd = data.frame(nrep = -1) %>% mutate(nrep = factor(nrep, levels = c(0:4, "5+"))))
  } else {
    pattern_f <- paste0("^(", f$Inserted, ")*")
    
    f$nrep <- NA
    for (i in seq_along(pattern_f)) {
      tmp <- gsub(pattern_f[i], "x", f$Next50[i])
      f$nrep[i] <- (nchar(f$Next50[i]) - (nchar(tmp) - 1)) / nchar(f$Inserted[i])
    }
    f[which(f$nrep > 4), "nrep"] <- "5+"

    
    f <- list(fa = filter(f, nchar(Inserted) == 2) %>% mutate(nrep = factor(nrep, levels = c(0:4, "5+"))),
              fb = filter(f, nchar(Inserted) == 3) %>% mutate(nrep = factor(nrep, levels = c(0:4, "5+"))),
              fc = filter(f, nchar(Inserted) == 4) %>% mutate(nrep = factor(nrep, levels = c(0:4, "5+"))),
              fd = filter(f, nchar(Inserted) >= 5) %>% mutate(nrep = factor(nrep, levels = c(0:4, "5+"))))
  }
  #===============
  
  filter(y, Category == "homology", nchar(Deleted) > 1, !is.na(Sample), label == gg) -> g
  
  if (nrow(g) == 0) {
    g <- list(ga = data.frame(mhl = -1) %>% mutate(mhl = factor(mhl, levels = c(1))),
              gb = data.frame(mhl = -1) %>% mutate(mhl = factor(mhl, levels = c(1:2))),
              gc = data.frame(mhl = -1) %>% mutate(mhl = factor(mhl, levels = c(1:3))),
              gd = data.frame(mhl = -1) %>% mutate(mhl = ifelse(mhl > 4, "5+", mhl), mhl = factor(mhl, levels = c(1:4, "5+"))))
  } else {
        g <- list(ga = filter(g, nchar(Deleted) == 2) %>% mutate(mhl = factor(mhl, levels = c(1))),
              gb = filter(g, nchar(Deleted) == 3) %>% mutate(mhl = factor(mhl, levels = c(1:2))),
              gc = filter(g, nchar(Deleted) == 4) %>% mutate(mhl = factor(mhl, levels = c(1:3))),
              gd = filter(g, nchar(Deleted) >= 5) %>% mutate(mhl = ifelse(mhl > 4, "5+", mhl), mhl = factor(mhl, levels = c(1:4, "5+"))))
  }
  
  ret <- c(table(a$hpl), table(b$hpl), table(c$hpl), table(d$hpl), table(e$ea$nrep), table(e$eb$nrep), table(e$ec$nrep), table(e$ed$nrep), table(f$fa$nrep), table(f$fb$nrep), table(f$fc$nrep), table(f$fd$nrep), table(g$ga$mhl), table(g$gb$mhl), table(g$gc$mhl), table(g$gd$mhl))
  
  restab[,gg] <- ret
  
  pdf(paste0("_", gg, "_id_spectrum.pdf"), width = 20, height = 3)
  barplot(ret, cex.names = .7, ylim = c(0, 20))
  dev.off()
print(gg)
}

# load PACWG signatures, fomat table
id_sig <- read.delim("~/Downloads/sigProfiler_ID_signatures.csv", sep = ",")
rownames(id_sig) <- id_sig$Mutation.Type
id_sig$Mutation.Type <- NULL
id_sig <- as.matrix(id_sig)

# normalize to 1 for all genotypes
apply(restab, 2, function(x) x/sum(x)) -> restab
rownames(restab) <- rownames(id_sig)

# subsetting indel signatures
id_select <-paste0("ID", c(1,3,5,6,8))
id_select0 <- paste0("ID", 1:17)

# helper function for consistent coloring
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colz <- gg_color_hue(length(id_select))

# calculating non-negative least squares coefficients
nnls_tab <- matrix(0, ncol = ncol(restab), nrow = ncol(id_sig[,id_select]))
colnames(nnls_tab) <- colnames(restab)
rownames(nnls_tab) <- colnames(id_sig[,id_select])
for (n in colnames(nnls_tab)) {
  nnls_tab[,n] <- lsqnonneg(id_sig[,id_select], restab[,n])$x
}

# calculating RMSDs
rmsd_tab <- data.frame(Genotype = colSums(restab), RMSD = NA, Mutnum = as.numeric(table(x$label)))
for (i in 1:11) {
  sweep(id_sig[,id_select], 2, nnls_tab[,i], "*") %>% rowSums() -> a
  restab[,i] -> b
  sqrt(mean((a-b)**2)) -> rmsd_tab$RMSD[i]
}

# attaching unassigned indels
rbind(nnls_tab, table(x$label) - colSums(nnls_tab)) -> nnls_tab
rownames(nnls_tab)[length(id_select) + 1] <- "Unassigned"

# creating contribution plot object
as.data.frame(nnls_tab) %>%
  mutate(Signature = factor(rownames(.), levels = c("Unassigned", id_select))) %>%
  gather(Genotype, Contribution, -Signature) %>%
  mutate(Genotype = gsub("_mock", "", Genotype),
         Genotype = factor(Genotype, levels = c("WT", "BRCA1", "RAD51C", "XRCC2", "XRCC3", "BRCA2", "PALB2", "RAD52", "RAD54", "ATM", "CHK2")),
         Mutnum = rep(as.numeric(table(x$label)), each = length(id_select) + 1)) %>%
  ggplot(aes(x = Genotype, y = Contribution, fill = Signature)) + geom_col(color = "black") + theme_classic() + theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
                                                                                           axis.title.y = element_text(size = 16),
                                                                                           axis.text = element_text(size = 14),
                                                                                           axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
   xlab("") + ylab("Absolute contribution") + scale_fill_manual(values = c("grey70", colz)) -> g1

# creating RMSD plot object
mutate(rmsd_tab, Genotype = gsub("_mock", "", rownames(rmsd_tab)),
       Genotype = factor(Genotype, levels = c("WT", "BRCA1", "RAD51C", "XRCC2", "XRCC3", "BRCA2", "PALB2", "RAD52", "RAD54", "ATM", "CHK2"))) %>%
  ggplot(aes(x = Genotype, y = RMSD/Mutnum)) + geom_col() + theme_classic() + theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
                                                                                    axis.title.y = element_text(size = 16),
                                                                                    axis.text = element_text(size = 14),
                                                                                    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("") + ggtitle("RMSD / Indel number") + scale_x_discrete(labels = rep("", 11)) + theme(plot.title = element_text(hjust = 0, size = 12))  -> g2

# plotting together
library(gridExtra)
gz <- list(ggplotGrob(g2), ggplotGrob(g1))
gz[[1]]$widths <- gz[[2]]$widths
arrangeGrob(grobs = gz, heights = c(1, 2.5))
