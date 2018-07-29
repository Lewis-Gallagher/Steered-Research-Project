# Performs differential gene expression analysis on RSEM input files *.genes.results
# Gene symbols are dependant upon the genome build being used
# Outputs .csv file containing DGE results sorted by p-value and .csv file containting
# normalised read counts


library(edgeR)
library(biomaRt)
library(limma)
setwd("/home/lewisg/BS7120/Bowtie/bowtie-0.12.7/Results/Rn4")

genesResults <- dir(path = ".", pattern = "*.genes.results")
genesResults
geneList <- lapply(genesResults, read.table)

# Know which SRR numbers are controls - the order matters
names(geneList) <- c("mic_raw1", "con_raw1", "mic_raw2", "con_raw2", "mic_raw3")


# Make matrix of count data i.e. second columm of each table
mic_counts <- data.frame(lapply(geneList, "[",2))

MIC_matrix <- data.matrix(mic_counts)
rm(mic_counts)
MIC_matrix <- MIC_matrix[,c(2,4,1,3,5)] # rearragne so controls are first columns
colnames(MIC_matrix) <- c("con1_ExpCounts","con2_ExpCounts",
                          "mic1_ExpCounts","mic2_ExpCounts",
                          "mic3_ExpCounts")
geneid <- geneList$con_raw1$V1
rownames(MIC_matrix) <- geneid

# Normalisation and Filtering (from edgeR userguide)

group <- c(1,1,2,2,2)
dge <- DGEList(counts = MIC_matrix, genes = geneid, group = group)

# Filter https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
keep <- rowSums(cpm(dge)>1) >= 2
dge <- dge[keep,,keep.lib.size = FALSE]
nc <- cpm(dge, prior.count = 2, log = TRUE, normalized.lib.sizes = FALSE)
dge <- calcNormFactors(dge, method = "TMM")

# Make the experiment matrix model
design <- model.matrix(~group, data = dge$samples)
colnames(design) <- (c("CON","MIC"))

# Estimate dispersion using Bayes methhod
dge <- estimateDisp(dge, design, robust=TRUE)
dge$common.dispersion

# Fit
fit <- glmQLFit(dge, design, robust = TRUE)
tr <- glmTreat(fit, coef = 2, lfc = 1.2) # filter lfc

# Create summary tables
trTags <- topTags(tr, n = Inf, adjust = "BH", sort.by = "PValue", p.value = 0.05) # filter pvalue

sumTR <- data.frame(summary(de <- decideTestsDGE(
    tr, adjust.method = "BH", lfc = 1.2, p.value = 0.05))) # filter lfc and pvalue
colnames(sumTR) <- c("signal", "chemical", "freq")
sumTR

trTagsTable <- data.frame(trTags$table)
colnames(trTagsTable) <- c("ensembl_gene_id","logFC",
                           "unshrunk_logFC", "logCPM",
                           "PValue", "FDR")


# Some plots
plotQLDisp(fit)
plotBCV(dge)
plotMDS.DGEList(dge)
barplot(dge$samples$lib.size*1e-6, ylab="Library size (millions)",
        xlab = "Library", main = "Library sizes",
        names.arg = c("Con1", "Con2", "MIC1", "MIC2", "MIC3"))
plotMD(tr)
abline(h=c(-1, 1), col="yellow")



# Retrieve gene symbols with BiomaRt
# The biomart changes depending on genome build.  Use the relevant one

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

# Rn4
ensembl <- listMarts(host = "may2012.archive.ensembl.org")
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "rnorvegicus_gene_ensembl",
                   host = "may2012.archive.ensembl.org")

genes <- getBM(attributes=c("ensembl_gene_id",
                            "external_gene_id"),
               filters="ensembl_gene_id",
               values=geneid, mart=ensembl,
               uniqueRows = TRUE)


# Merge the found genes table with our top hits to retieve gene symbols
trTagsTable <- merge(x = trTagsTable, y = genes, by = "ensembl_gene_id", all = FALSE)
trTagsTable <- trTagsTable[,c(1,7,2,3,4,5,6)]
colnames(trTagsTable) <- c("ensembl_gene_id", "symbol", "logFC", "unshrunk_logFC", "logCPM", "PValue", "FDR")
trTagsTable <- trTagsTable[order(trTagsTable$PValue),]
head(trTagsTable)

# Create table of normalised counts
nc <- data.frame(nc)
id <- rownames(nc)
nc <- cbind(ensembl_gene_id = id, nc) 
nc <- merge(x = nc, y = genes, by = "ensembl_gene_id", all = FALSE)
nc <- nc[,c(1,7,2,3,4,5,6)]


write.csv(trTagsTable, "/home/lewisg/BS7120/Bowtie/Results/Bowtie1/Rn6/b1_rn6_genes_Unfiltered.csv")
write.csv(nc, "/home/lewisg/BS7120/Bowtie/Results/Bowtie1/Rn6/rn6_genes_normCounts.csv")
