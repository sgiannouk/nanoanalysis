### Stavros Giannoukakos ###
### ONT ANALYSIS FOR DIFFERENTIAL TRANSCRIPT USAGE (DTU) 


args <- commandArgs(TRUE)
if (length(args) == 7) {
  # Input filtered matrix (output of step 8)
  matrix <- args[1]
  # CSV file used for running TALON
  input_groups <- args[2]
  # Output direcotry where all stats will be saves
  main_outdir <- args[3]
  # Minimum number of total mapped sequence reads 
  # for a gene isoform to be considered
  minFeatureExpr <- as.numeric(args[4])
  # Adjusted p-value threshold for differential expression analysis
  adjPValueThreshold <- as.numeric(args[5])
  # Minimum required log2 fold change for differential expression analysis
  lfcThreshold <- as.numeric(args[6])
  # Number of cores to be used for the differential expression analysis
  threads <- as.numeric(args[7])
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

# matrix <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_2/filt_talon_abundance.csv"
# input_groups <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_2/talon_input.csv"
# main_outdir <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_2/diffExpr_analysis"
# minFeatureExpr <- 10  # Minimum feature counts
# adjPValueThreshold <- 0.01
# lfcThreshold <- 1
# threads <- 2


library("tidyr")
library("stageR")
library("DEXSeq")
library("DRIMSeq")
library("ggplot2")
library("BiocParallel")
BPPARAM = MulticoreParam(workers=threads)



##### DIFFERENTIAL TRANSCRIPT USAGE (DTU) ANALYSIS USING DRIMSEQ/DEXSEQ/STAGER #####
print("RUNNING DIFFERENTIAL TRANSCRIPT USAGE (DTU) ANALYSIS USING DRIMSEQ/DEXSEQ/STAGER")

outdir <- file.path(main_outdir, "DiffTranscriptUsage")
dir.create(outdir, showWarnings = FALSE)
setwd(outdir)

# Input the filtered expression matrix
expr_file <- read.csv(matrix)
print(paste("Total number of unique genes:", length(unique(expr_file$annot_gene_id)), sep=" " ))
print(paste("Total number of unique trancripts:", length(unique(expr_file$annot_transcript_id)), sep=" " ))

# Obtain number of input samples
n <- length(expr_file[ ,10:length(expr_file)])
print(paste("Total number of input samples:", n, sep=" " ))

# Reading input csv file containing the groups
group_samples <- read.csv(input_groups, header=F)[ ,1:3]
group_samples <- group_samples[order(group_samples$V2), ] # Order data frame

# Preparing a dataframe containing information about the samples 
group_samples <- data.frame(sample_id = group_samples$V1, condition = group_samples$V2, batch=group_samples$V3)

sampletypevalues <- factor(unique(group_samples$condition)) # Obtaining the sample groups
print(paste("Total number of input groups: ", length(sampletypevalues)," (", sampletypevalues[1], " and ", sampletypevalues[2],")", sep="" ))


# Obtaining the counts table
matfile <- expr_file[ ,c(2,1,10:length(expr_file))]  # Obtaining only the annot_gene_id and the counts
colnames(matfile)[1:2] <- c("feature_id","gene_id")
# Create a dmDSdata object that is the starting point of DRIMSeq
drimseq <- dmDSdata(counts = matfile, samples = group_samples)
print("General input stats:")
print(drimseq)   # General stats


# Check what is the minimal number of replicates per condition
print("Number of replicates per condition")
print(table((DRIMSeq::samples(drimseq))$condition))

# Filtering of lowly expressed transcript
# The parameters are the suggested by ONT
drimseq <- dmFilter(drimseq, 
                    min_samps_feature_expr = 1,
                    min_feature_expr = minFeatureExpr)

# Obtaining the counts after dmFiltering
filtered_counts <- counts(drimseq)

# Printing stats
print(paste("In total ", length(unique(filtered_counts$gene_id)), 
            " out of ", length(unique(matfile$gene_id))," genes passed the min. thresholds.",sep=""))
print(paste("In total ", length(unique(filtered_counts$feature_id)), 
            " out of ", length(unique(matfile$feature_id))," genes passed the min. thresholds.",sep=""))

# Creating with design matrix
if (length(unique(group_samples$batch) == 1)) {
  design <- model.matrix( ~ condition, data = DRIMSeq::samples(drimseq))
} else {
  design <- model.matrix( ~ condition + batch, data = DRIMSeq::samples(drimseq)) }
rm(n, threads, main_outdir, minFeatureExpr, matfile)


### Differential transcript usage using DEXSeq ###
sample.data<-DRIMSeq::samples(drimseq)  # Reformating and obtaining the metadata
row.names(sample.data) <- sample.data$sample_id; sample.data$sample_id <- NULL
sample.data$condition <- factor(sample.data$condition)
count.data <- data.frame(counts(drimseq)[,-c(1:2)])  # Obtaining the read counts of the filtered data frame

# Constructing the DEXSeqDataSet object that will store our data
dxd <- DEXSeqDataSet(countData=count.data, 
                     sampleData=sample.data, 
                     design= ~sample + exon + condition:exon, 
                     featureID=filtered_counts$feature_id, 
                     groupID=filtered_counts$gene_id)

# Normalisation  - DEXSeq uses the same method as DESeq and DESeq2
dxd <- estimateSizeFactors(dxd)

# Estimating the variability of the data, in order to test for differential exon usage
dxd <- estimateDispersions(dxd, BPPARAM=BPPARAM)
rm(sample.data, count.data)

# Having the dispersion estimates and the size factors, we can now test for differential exon usage. 
# For each gene, DEXSeq fits a generalized linear model with the formula ~sample + exon + condition:exon 
# and compare it to the smaller model (the null model) ~ sample + exon
dxd <- testForDEU(dxd, reducedModel=~sample + exon, BPPARAM=BPPARAM)

# Estimating relative exon usage fold changes
dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=BPPARAM)

# Obtaining the intermediate and final results
dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)

# # MA-plot - Drawing the expression levels over the exons to highlight differential exon usage
# png(paste(outdir,"/DEXSeq_MAPlot.png",sep=""), units='px', height=900, width=1600, res=90)
# plotMA(dxr, cex=0.8, alpha=0.05)
# title(main = "MA-plot")
# dev.off()

# Optimising the DEXSeq results
dxr_res <- data.frame(dxr)
# Renaming the log2FC
colnames(dxr_res)[10] <- "logFC"


# # MA-plot - Drawing the expression levels over the exons to highlight differential exon usage
# logUp <- which(dxr_res$logFC >= lfcThreshold)
# logDown <- which(dxr_res$logFC <= -lfcThreshold)
# withStat <- which(dxr_res$padj <= adjPValueThreshold)
# colours <- c(noDifference="dimgray", upRegulated="mediumseagreen", downRegulated="indianred3")
# gene <- rep("noDifference", nrow(dxr_res))
# gene[logUp[logUp %in% withStat]] <- "upRegulated"
# gene[logDown[logDown %in% withStat]] <- "downRegulated"
# ggplot(data.frame(dxr_res), aes(y=logFC, x=exonBaseMean)) + 
#       geom_point(size=1.2) + 
#       geom_hline(yintercept = -lfcThreshold, color="indianred3") + 
#       geom_hline(yintercept = lfcThreshold, color="mediumseagreen") +
#       theme_bw() +
#       aes(colour=gene) + 
#       scale_colour_manual(name="Genes", values=colours) +
#       xlab("log(CountsPerMillion)") +
#       ylab("log(FoldChange)") +
#       ggtitle("MA plot - log(FC) vs. log(CPM) on gene level data")
# ggsave(file=paste(outdir,"/DEXSeq_MAplot.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)


### StageR analysis on the DEXSeq differential transcript usage results ###
pConfirmation <- matrix(dxr$pvalue, ncol=1)
dimnames(pConfirmation) <- list(c(dxr$featureID), c("transcript"))
pScreen <- perGeneQValue(dxr)
tx2gene <- data.frame(row.names = dxr$featureID , transcript = dxr$featureID, gene = dxr$groupID)
           
stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)
# Transcript-level adjusted p-values for genes not passing the screening stage are set to NA by default
stageRObj <- stageWiseAdjustment(object=stageRObj, method="dtu", alpha=adjPValueThreshold)
padj <- getAdjustedPValues(stageRObj, order=TRUE, onlySignificantGenes=FALSE)

# Creating the output matrix with the stats from DEXSeq and stageR
dexseq_results <- data.frame(dxr)
# Choosing the necessary columns to maintain
dexseq_results <- dexseq_results[ ,c(1,2,10,6,7)]; row.names(dexseq_results) <- NULL
# Renaming the columns
colnames(dexseq_results)[1:3] <- c("gene_id", "transcript_id", "log2FoldChange")
# Merging the two data frames
dexseq_results <- merge(dexseq_results, padj, by.x="transcript_id", by.y="txID")
# Rearranging 
dexseq_results$geneID <- NULL; dexseq_results <- dexseq_results[ ,c(2,1,3:7)]
# Renaming
colnames(dexseq_results)[6:7] <- c("padj_gene", "padj_transcript")
# Incorporating the gene and transcript names in the data frame
colnames(expr_file)[1:2] <- c("gene_id", "transcript_id")
dexseq_results <- merge(dexseq_results, expr_file[ ,c(2:4)], by="transcript_id", all.x=T)
# Rearranging
dexseq_results <- dexseq_results[ ,c(2,1,8,9,3:7)]
# Renaming the new columns
colnames(dexseq_results)[3:4] <- c("gene_name", "transcript_name")
# Incorporating the sample data in the data frame
dexseq_results <- merge(counts(drimseq)[ ,c(1:length(counts(drimseq)))], dexseq_results, by.x=c("gene_id", "feature_id"), by.y=c("gene_id","transcript_id"))
# Rearranging
num <- length(group_samples$sample_id)
idx <- length(dexseq_results)
dexseq_results <- dexseq_results[ ,c(1, 2, num+3, num+4, 3:(2+num), (idx-4):idx)]
# Renaming
colnames(dexseq_results)[2] <- "transcript_id"

dexseq_results <- dexseq_results[order(dexseq_results$padj_gene), ]
# Exporting the normalised results table containing all features along with the output stats  from DEXSeq and stageR
write.table(dexseq_results, file=paste(outdir,"/",sampletypevalues[1],"VS",sampletypevalues[2],"_DEXSeqStageR_allDTU.csv", sep=""), sep="\t", row.names = F, quote=FALSE)

# Filtering out transcripts that are not showing DTU
candidates <- union(which(dexseq_results$padj_gene <= adjPValueThreshold), which(dexseq_results$padj_transcript < adjPValueThreshold))
dexseq_results_filt <- dexseq_results[candidates, ]

# Ordering by padj_transcript
dexseq_results_filt <- dexseq_results_filt[order(dexseq_results_filt$padj_gene), ]

# Exporting the normalised results table containing the selected features with adjusted p value lower than the user-input value
write.table(dexseq_results_filt, file=paste(outdir,"/",sampletypevalues[1],"VS",sampletypevalues[2],"_DEXSeqStageRBelow", gsub("[.]", "", adjPValueThreshold), ".csv", sep=""), sep="\t", row.names = F, quote=FALSE)
rm(filtered_counts, stageRObj)

# Data frame manipulation towards melting the sample counts
dexseq_results_filt_forplot <- dexseq_results_filt[ ,!(names(dexseq_results_filt) %in% c("log2FoldChange","pvalue","padj"))]
# Melting
dexseq_results_filt_forplot <- dexseq_results_filt_forplot %>% gather(key='sample', value='norm_count', -gene_id, -transcript_id, -gene_name, -transcript_name, -padj_gene, -padj_transcript)
# Adding group identity
dexseq_results_filt_forplot$Groups <- group_samples[match(dexseq_results_filt_forplot$sample, group_samples$sample_id), ]$condition


# Plotting genes showing DTU 
for(gene in unique(dexseq_results_filt$gene_name)){
    gdf <- dexseq_results_filt_forplot[which(dexseq_results_filt_forplot$gene_name==gene), ]
    
    ggplot(gdf, aes(x=transcript_name, y=norm_count)) + 
           geom_boxplot(aes(fill=Groups), position="dodge") + 
           geom_dotplot(binaxis="y", stackdir="center", dotsize=0.6, aes(fill=Groups), position="dodge") + 
           theme_bw() +
           scale_fill_brewer(palette="Paired") +
           theme(axis.text.x = element_text(angle = 30,hjust=1)) +
           labs(title=paste("Boxplots showing transcript expression\nlevels across conditions for gene", gene), x="DTUs", y="Transcript read count")
ggsave(file=paste(outdir,"/DEXSeqStageR_DTU_", gene, ".png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)}
