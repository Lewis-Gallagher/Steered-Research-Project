# Read in and remove first column, not sure why its there
setwd("/home/lewisg/BS7120/Bowtie/Results/Rn4")
library(limma)
mic_raw1 <- read.table("SRR1178015_rsem_map.genes.results", header =TRUE)
mic_raw2 <- read.table("SRR1178022_rsem_map.genes.results", header =TRUE)
mic_raw3 <- read.table("SRR1178048_rsem_map.genes.results", header =TRUE)
con_raw1 <- read.table("SRR1178019_rsem_map.genes.results", header =TRUE)
con_raw2 <- read.table("SRR1178024_rsem_map.genes.results", header =TRUE)
colnames(mic_raw1) <- c("gene_id", "expected_count", "tauvalue", "transcript_id")
colnames(mic_raw2) <- c("gene_id", "expected_count", "tauvalue", "transcript_id")
colnames(mic_raw3) <- c("gene_id", "expected_count", "tauvalue", "transcript_id")
colnames(con_raw1) <- c("gene_id", "expected_count", "tauvalue", "transcript_id")
colnames(con_raw2) <- c("gene_id", "expected_count", "tauvalue", "transcript_id")


#Make matrix of count data
mic_counts <- data.frame(con_raw1$expected_count,
                         con_raw2$expected_count,
                         mic_raw1$expected_count,
                         mic_raw2$expected_count,
                         mic_raw3$expected_count)
MIC_matrix <- data.matrix(mic_counts)
colnames(MIC_matrix) <- c("con1_ExpCounts","con2_ExpCounts",
                          "mic1_ExpCounts","mic2_ExpCounts",
                          "mic3_ExpCounts")
rownames(MIC_matrix) <- (con_raw1$gene_id)

# Normalisation and Filtering (from edgeR userguide)
library(edgeR)
group <- c(1,1,2,2,2)
dge <- DGEList(counts=MIC_matrix, genes = rownames(MIC_matrix), group=group)

# Filter
keep <- rowSums(cpm(dge)>1) >= 2
dge <- dge[keep,,keep.lib.size=FALSE]
dge <- calcNormFactors(dge)

# Make the experiment matrix model
design <- model.matrix(~group, data = dge$samples)
colnames(design) <- (c("CON","MIC"))

# Estimate dispersion using Bayes method
dge <- estimateDisp(dge, design, robust=TRUE)
dge$common.dispersion

# Fit
fit <- glmQLFit(dge, design, robust = TRUE)
plotQLDisp(fit)
tr <- glmTreat(fit, coef=2, lfc = 1.5)
trTags <- topTags(tr, n = Inf, adjust = "BH", sort.by = "PValue", p.value=0.05)


# Create summary tables
sumTR <- data.frame(summary(de <- decideTestsDGE(tr, adjust.method = "BH", p.value = 0.05, lfc = 1.5)))
colnames(sumTR) <- c("signal", "chemical", "freq")
sumTR

trTagsTable <- data.frame(trTags$table)
colnames(trTagsTable) <- c("ensembl_gene_id", "logFC", "unshrunk_logFC", "logCPM", "PValue", "FDR")


# Some plots
plotBCV(dge)
plotMDS.DGEList(dge)
barplot(dge$samples$lib.size*1e-6, ylab="Library size (millions)")
plotMD(tr)
abline(h=c(-1, 1), col="blue")


# Retrieve gene symbols with BiomaRt
library(biomaRt)

geneid <- c(trTagsTable$ensembl_gene_id)

# Rn6
# ensembl <- listMarts()
# ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                    dataset = "rnorvegicus_gene_ensembl")
# 
# genes <- getBM(attributes=c("ensembl_gene_id",
#                             "external_gene_name"),
#                filters="ensembl_gene_id",
#                values=geneid, mart=ensembl,
#                uniqueRows = TRUE)

# # Rn4
ensembl <- listMarts(host = "may2012.archive.ensembl.org")
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "rnorvegicus_gene_ensembl",
                   host = "may2012.archive.ensembl.org")

genes <- getBM(attributes=c("ensembl_gene_id",
                            "external_gene_id"),
               filters="ensembl_gene_id",
               values=geneid, mart=ensembl,
               uniqueRows = TRUE)

trTagsTable <- merge(x = trTagsTable, y = genes, by = "ensembl_gene_id", all = TRUE)
trTagsTable <- trTagsTable[,c(1,7,2,3,4,5,6)]
colnames(trTagsTable) <- c("ensembl_gene_id", "symbol", "logFC", "unshrunk_logFC", "logCPM", "PValue", "FDR")
trTagsTable <- trTagsTable[order(trTagsTable$PValue),]



write.csv(trTagsTable, "rn4_tr_isoformsAll.csv")
