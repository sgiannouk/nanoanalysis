### Stavros Giannoukakos ###
### ONT ANALYSIS FOR DIFFERENTIAL GENE EXPRESSION (DGE) 


args <- commandArgs(TRUE)
if (length(args) == 6) {
  # Input filtered matrix (output of step 8)
  matrix <- args[1]
  # CSV file used for running TALON
  input_groups <- args[2]
  # Output direcotry where all stats will be saves
  main_outdir <- args[3]
  # Minimum gene counts
  minGeneExpr <- as.numeric(args[4])
  # Adjusted p-value threshold for differential expression analysis
  adjPValueThreshold <- as.numeric(args[5])
  # Minimum required log2 fold change for differential expression analysis
  lfcThreshold <- as.numeric(args[6])
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

# matrix <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_2/prefilt_talon_abundance.tsv"
# input_groups <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_2/talon_input.csv"
# main_outdir <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_2/diffExpr_analysis"
# minGeneExpr <- 10  # Minimum gene counts
# adjPValueThreshold <- 0.01
# lfcThreshold <- 1


library("dplyr")
library("edgeR")
library("DRIMSeq")
library("ggplot2")
library("reshape")
options(scipen = 999)



##### DIFFERENTIAL GENE EXPRESSION (DGE) ANALYSIS USING DRIMSEQ/EDGER #####
print("RUNNING DIFFERENTIAL GENE EXPRESSION (DGE) ANALYSIS USING DRIMSEQ/EDGER")

outdir <- file.path(main_outdir, "DiffGeneExpr")
dir.create(outdir, showWarnings = FALSE)
setwd(outdir)

# Input the filtered expression matrix
expr_file <- read.table(matrix, sep="\t", header=T)
expr_file <- expr_file[ ,3:length(expr_file)]
# Ommiting the genomic transcripts 
expr_file <- expr_file[which(expr_file$transcript_novelty != "Genomic"), ]

# Exporting the unfiltered table
# write.table(expr_file, file=paste(outdir,"/edgeR_prefiltered_table.csv", sep=""), sep="\t", row.names = F, quote=FALSE)

print(paste("Total number of unique genes:", length(unique(expr_file$annot_gene_id)), sep=" " ))

# Obtain number of input samples
n <- length(expr_file[ ,10:length(expr_file)])
print(paste("Total number of input samples:", n, sep=" " ))
# Reading input csv file containing the groups
group_samples <- read.csv(input_groups, header=F)[ ,1:3]
group_samples <- group_samples[order(group_samples$V2), ] # Order data frame
# Preparing a dataframe containing information about the samples 
group_samples <- data.frame(sample_id = group_samples$V1, condition = group_samples$V2, batch=group_samples$V3)

sampletypevalues <- factor(unique(group_samples$condition)) # Obtaining the sample groups
print(paste("Total number of input groups: ", 
            length(sampletypevalues)," (", sampletypevalues[1], " and ", sampletypevalues[2],")", sep="" ))


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

# Filtering of lowly expressed genes
# The parameters are the suggested by ONT
drimseq <- dmFilter(drimseq, 
                    min_samps_gene_expr = 1,
                    min_samps_feature_expr = 1,
                    min_gene_expr = minGeneExpr,
                    min_feature_expr = 2)

# Obtaining the counts after dmFiltering
filtered_counts <- counts(drimseq)

# Printing stats
print(paste("In total ", length(unique(filtered_counts$gene_id)), 
            " out of ", length(unique(matfile$gene_id))," genes passed the min. thresholds.",sep=""))

# Creating with design matrix
if (length(unique(group_samples$batch) == 1)) {
  design <- model.matrix( ~ condition, data = DRIMSeq::samples(drimseq))
} else {
  design <- model.matrix( ~ condition + batch, data = DRIMSeq::samples(drimseq)) }

# Summarising all reads per gene name and removing the duplicate  rows
gene_filt_counts <- filtered_counts; gene_filt_counts$feature_id <- NULL;
gene_filt_counts <- data.frame(dplyr::group_by(gene_filt_counts, gene_id) %>% dplyr::summarise_all(sum))

# Making the annot_transcript_id as row names and removing it as from 1st column
row.names(gene_filt_counts) <- gene_filt_counts$gene_id; gene_filt_counts$gene_id <- NULL

# Storing the raw read counts table in a simple list-based data object called a DGEList.
edgeR_table <- DGEList(gene_filt_counts, group=group_samples$condition)

# Normalisation for RNA composition by finding a set of scaling factors for the library sizes 
# that minimize the log-fold changes between the samples for most genes. The default method
# for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples.
edgeR_table <- calcNormFactors(edgeR_table)

# Estimating common dispersion and tagwise dispersions 
edgeR_table <- estimateDisp(edgeR_table, design)

# Performing quasi-likelihood F-tests
fit <- glmQLFit(edgeR_table, design)

qlf <- glmQLFTest(fit)

# Obtaining the results table and sorting by p-value
edger_res_overall <- topTags(qlf, n=Inf, sort.by="PValue")[[1]]
edger_res <- edger_res_overall[,c(1,4,5)];  colnames(edger_res) <- c("log2FoldChange",  "pval", "padj")

# Output basic stats
print(paste("Results indicating the gene set with log2FoldChange value > |", lfcThreshold,"| and with adjasted p-value < ", adjPValueThreshold, sep=""))
is.de <- decideTests(qlf, p.value=adjPValueThreshold, lfc=lfcThreshold)
print(summary(is.de))


# Selecting the samples
selected_samples <- colnames(gene_filt_counts[ ,(which(group_samples$condition==sampletypevalues[1] | group_samples$condition==sampletypevalues[2]))])
# Merging cpm values and stats from EdgeR
total_results <- data.frame(merge(cpm(edgeR_table)[ ,selected_samples], edger_res, by=0))
# Importing the gene names and remove duplicates
total_results <- merge(expr_file[ ,c("annot_gene_id", "annot_gene_name")], total_results, by.y="Row.names", by.x="annot_gene_id", all.y=TRUE, no.dups=TRUE)
total_results <- total_results[!duplicated(total_results$annot_gene_id),  ]
# Renaming the first two columns
colnames(total_results)[1] <- "gene_id"; colnames(total_results)[2] <- "gene_name"
# Sorting the table by pvalue and adjasted p-value
total_results <- total_results[with(total_results, order(pval, padj)), ]
# Exporting the normalised results table containing the selected features with log2FoldChange greater than 1 and adjusted p value lower than the user-input value
write.table(total_results, file=paste(outdir,"/",sampletypevalues[1],"VS",sampletypevalues[2],"_edgeR_allGenes.csv", sep=""), sep="\t", row.names = F, quote=FALSE)

# Output only results where logFC
results_selected <- total_results[(abs(total_results$log2FoldChange)>=lfcThreshold) & (total_results$padj<=adjPValueThreshold), ]
# Sorting the table by pvalue and adjasted p-value
results_selected <- results_selected[with(results_selected, order(pval, padj)), ]
# Exporting the normalised results table containing the selected features with log2FoldChange greater than 1 and adjusted p value lower than the user-input value
write.table(results_selected, file=paste(outdir,"/",sampletypevalues[1],"VS",sampletypevalues[2],"_edgeR_topGenesBelow", gsub("[.]", "", adjPValueThreshold), "LFC", gsub("[.]", "", lfcThreshold), ".csv", sep=""), sep="\t", row.names = F, quote=FALSE)
rm(design, drimseq, fit, qlf)


# MA-plot
logUp <- which(edger_res_overall$logFC >= lfcThreshold)
logDown <- which(edger_res_overall$logFC <= -lfcThreshold)
withStat <- which(edger_res_overall$FDR <= adjPValueThreshold)
colours <- c(noDifference="dimgray", upRegulated="mediumseagreen", downRegulated="indianred3")
gene <- rep("noDifference", nrow(edger_res_overall))
gene[logUp[logUp %in% withStat]] <- "upRegulated"
gene[logDown[logDown %in% withStat]] <- "downRegulated"
ggplot(data.frame(edger_res_overall), aes(y=logFC, x=logCPM)) +
        geom_point(size=1.2) +
        geom_hline(yintercept = -lfcThreshold, color="indianred3") +
        geom_hline(yintercept = lfcThreshold, color="mediumseagreen") +
        theme_bw() +
        aes(colour=gene) +
        scale_colour_manual(name="Genes", values=colours) +
        labs(title="MA plot - log(FC) vs. log(CPM) on gene level data", x="log(CountsPerMillion)", y="log(FoldChange)")
ggsave(file=paste(outdir,"/edgeR_MAplot.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
rm(is.de, logUp, logDown, withStat, colours, gene)


# MA-plot
withNovel <- expr_file[grepl("^TALONT[^*]+$", expr_file$annot_transcript_id), ]
withNovel <- as.vector(withNovel$annot_gene_id)
withNovel <- which(rownames(edger_res_overall) %in% withNovel)
edger_res_overall$novelty <- 16
edger_res_overall[withNovel, length(edger_res_overall)] <- 18
edger_res_overall$novelty = as.factor(edger_res_overall$novelty)
logUp <- which(edger_res_overall$logFC >= lfcThreshold)
logDown <- which(edger_res_overall$logFC <= -lfcThreshold)
withStat <- which(edger_res_overall$FDR <= adjPValueThreshold)
colours <- c(noDifference="dimgray", upRegulated="mediumseagreen", downRegulated="indianred3")
gene <- rep("noDifference", nrow(edger_res_overall))
gene[logUp[logUp %in% withStat]] <- "upRegulated"
gene[logDown[logDown %in% withStat]] <- "downRegulated"
ggplot(edger_res_overall, aes(x=logCPM, y=logFC)) +
       geom_point(aes(shape=factor(novelty,  labels=c("Known", "Novel")))) +
       geom_hline(yintercept = -lfcThreshold, color="indianred3") +
       geom_hline(yintercept = lfcThreshold, color="mediumseagreen") +
       theme_bw() +
       aes(colour=gene) +
       scale_colour_manual(name="Genes", values=colours) +
       labs(title="MA plot - log(FC) vs. log(CPM) on gene level data with Known/Novel-transcript indication", x="log(CountsPerMillion)", y="log(FoldChange)", shape="Known/Novel")
ggsave(file=paste(outdir,"/edgeR_MAplot_withTranscriptIndications.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
edger_res_overall$novelty <- NULL
rm(logUp, logDown, withStat, colours, gene,withNovel)


# Volcano plot
# Compute significance, with a maximum of 350 for the p-values set to 0 due to limitation of computation precision
edger_res_overall <- merge(expr_file[ ,c("annot_gene_id", "annot_gene_name")], edger_res_overall, by.x="annot_gene_id", by.y = 0, all.y=T)
edger_res_overall <- edger_res_overall[!duplicated(edger_res_overall$annot_gene_id),  ]
row.names(edger_res_overall) <- edger_res_overall$annot_gene_id; edger_res_overall$annot_gene_id <- NULL; colnames(edger_res_overall)[1] <- "gene_name"
edger_res_overall$significant <- (-log10(edger_res_overall$FDR))
edger_res_overall[is.infinite(edger_res_overall$significant), "significant"] <- 350
# Select transcripts with a defined p-value (DESeq assigns NA to some transcripts)
genes.to.plot <- !is.na(edger_res_overall$PValue)
## Volcano plot of adjusted p-values
cols <- densCols(edger_res_overall$logFC, edger_res_overall$significant)
cols[edger_res_overall$PValue==0] <- "purple"
edger_res_overall$pch <- 19
edger_res_overall$pch[edger_res_overall$PValue == 0] <- 6
print("Generating the volcano plot...")
png(paste(outdir,"/edgeR_volcanoPlot.png",sep=""), units='px', height=900, width=1600, res=100)
plot(edger_res_overall$logFC,
     edger_res_overall$significant,
     col=cols,
     panel.first=grid(),
     main="Volcano plot",
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     pch=edger_res_overall$pch, cex=0.4)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(adjPValueThreshold), col="brown")
## Plot the names of a reasonable number of transcripts, by selecting those begin not only significant but also having a strong effect size
gn.selected <- (abs(edger_res_overall$logFC)>=lfcThreshold & edger_res_overall$FDR<=adjPValueThreshold)
text(edger_res_overall$logFC[gn.selected], edger_res_overall$significant[gn.selected],lab=edger_res_overall$gene_name[gn.selected], cex=0.6)
dev.off()
rm(genes.to.plot, cols, gn.selected)
edger_res_overall$significant <- NULL; edger_res_overall$pch  <- NULL; #edger_res_overall$gene_name <- NULL
rm(gene_filt_counts)


### TOP 30 DE GENES
# Plotting the normalized count values for the top 30 differentially expressed genes (by padj values)
## Order results by by pvalue and adjasted p-value values and get the first 30 genes
total_results <- total_results[with(total_results, order(pval, padj)), ]
# Obtaining the top 30 most significant genes and the per sample normalised (CPM) data
top30_sigNorm <- total_results[1:30, !names(total_results) %in% c("log2FoldChange", "pval", "padj")]
# Gathering the columns to have normalized counts to a single column
melt_top30_sigNorm <- melt(top30_sigNorm)
melt_top30_sigNorm$group <- gsub("_[^_]+$", "", melt_top30_sigNorm$variable)
melt_top30_sigNorm$value[melt_top30_sigNorm$value == 0] <- 0.1

ggplot(melt_top30_sigNorm, aes(x = gene_name, y = value, color = group)) +
       geom_point(size=2.5) +
       scale_y_log10(labels = function(x) format(x, scientific = F)) +
       xlab("Genes") +
       ylab("log10(CPM)") +
       ggtitle("Top 30 Significant DE Genes") +
       theme_bw() +
       scale_colour_manual(name="", values=c("#66CC99", "#877598")) +
       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
       theme(plot.title = element_text(hjust = 0.5)) +
       theme(legend.position="bottom")
ggsave(file=paste(outdir,"/edgeR_top30MostSingificantGenes.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
rm(melt_top30_sigNorm, top30_sigNorm)
rm(matrix, filtered_counts, edgeR_table, edger_res, results_selected)


### Checking number of  transcripts per condition
# Obtaining the counts table
matfile <- expr_file[ ,c(2,1,10:length(expr_file))]  # Obtaining only the annot_gene_id and the counts
colnames(matfile)[1:2] <- c("feature_id","gene_id")

# Create a dmDSdata object that is the starting point of DRIMSeq
drimseq <- dmDSdata(counts = matfile, samples = group_samples)

# Filtering of lowly expressed transcript
# The parameters are the suggested by ONT
drimseq <- dmFilter(drimseq, 
                    min_samps_gene_expr = 1,
                    min_samps_feature_expr = 1,
                    min_gene_expr = 2,
                    min_feature_expr = 2)


# Obtaining the counts after dmFiltering
filtered_counts <- counts(drimseq)

# Selecting the samples from the first group
selected_group1 <- cbind(filtered_counts[ ,c(1,2)], filtered_counts[ ,which(colnames(filtered_counts) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[1]])])
# Removing unexpressed genes in this group (genes with 0 in all samples)
selected_group1 <- selected_group1[rowSums(selected_group1[ ,which(colnames(selected_group1) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[1]])]) > 0, ] 
# Transform to count table
selected_group1 <- data.frame(table(selected_group1$gene_id))
# Collapsing and counting transcripts
selected_group1 <- data.frame(table(selected_group1$Freq))
# Normalising for total sequencing depth
selected_group1 <- cbind(selected_group1[ ,1], data.frame(cpm(selected_group1[ ,2])))
colnames(selected_group1) <- c("Var1", "Freq")
selected_group1$group <- as.character(sampletypevalues[1])  # Renaming all genes to samplegroup


# Selecting the samples from the second group
selected_group2 <- cbind(filtered_counts[,c(1,2)], filtered_counts[ ,which(colnames(filtered_counts) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[2]])])
# Removing unexpressed genes in this group (genes with 0 in all samples)
selected_group2 <- selected_group2[rowSums(selected_group2[ ,which(colnames(selected_group2) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[2]])]) > 0, ]
# Transform to count table
selected_group2 <- data.frame(table(selected_group2$gene_id))
# Collapsing and counting transcripts
selected_group2 <- data.frame(table(selected_group2$Freq))
# Normalising for total sequencing depth
selected_group2 <- cbind(selected_group2[ ,1], data.frame(cpm(selected_group2[ ,2])))
colnames(selected_group2) <- c("Var1", "Freq")
selected_group2$group <- as.character(sampletypevalues[2])  # Renaming all genes to samplegroup

# Merging the melted data frames
merged_melted_tables <- rbind(selected_group1[1:40, ], selected_group2[1:40, ])
# write.table(merged_melted_tables, file=paste(outdir,"/table_of_NumofTrPerCondition.csv", sep=""), sep="\t", row.names = F, quote=FALSE)
rm(selected_group1, selected_group2)

ggplot(merged_melted_tables, aes(x = Var1, y = Freq, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual("Groups", values = c("#66CC99", "#877598")) +
  theme_bw() + 
  theme(legend.position="bottom") +
  theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title="Number of isoform usage per condition (subset of up to 40 transcripts)", x="Number of transcprits", y="Relative frequency") +
  ggsave(file=paste(outdir, "/DGE_TrNumPerCondition.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)


# ########### NUMBER OF ISOFORMS PER FOLD CHANGE ###########
# num_of_isoforms_plot <- function(mat, who){
#             # Selecting the samples from the first group
#             selected_group1 <- cbind(gene_id=mat[ ,1], mat[ ,which(colnames(mat) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[1]])])
#             # Removing unexpressed genes in this group (genes with 0 in all samples)
#             selected_group1 <- selected_group1[rowSums(selected_group1[ ,which(colnames(selected_group1) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[1]])]) > 0, ]
#             # Normalising for sequencing depth
#             selected_group1 <- cbind(selected_group1[ ,1:2], data.frame(cpm(selected_group1[ ,3:length(selected_group1)])))
#             # Transform to count table
#             selected_group1 <- data.frame(plyr::count(selected_group1, "gene_id"))
#             # Collapsing and counting transcripts
#             selected_group1 <- data.frame(table(Var1=selected_group1$freq))
#             selected_group1$group <- as.character(sampletypevalues[1])  # Renaming all genes to samplegroup
#     
#             # Selecting the samples from the second group
#             selected_group2 <- cbind(gene_id=mat[ ,1], mat[ ,which(colnames(mat) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[2]])])
#             # Removing unexpressed genes in this group (genes with 0 in all samples)
#             selected_group2 <- selected_group2[rowSums(selected_group2[ ,which(colnames(selected_group2) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[2]])]) > 0, ]
#             # Normalising for sequencing depth
#             selected_group2 <- cbind(selected_group2[ ,1:2], data.frame(cpm(selected_group2[ ,3:length(selected_group2)])))
#             # Transform to count table
#             selected_group2 <- data.frame(plyr::count(selected_group2, "gene_id"))
#             # Collapsing and counting transcripts
#             selected_group2 <- data.frame(table(Var1=selected_group2$freq))
#             selected_group2$group <- as.character(sampletypevalues[2])  # Renaming all genes to samplegroup
#             # Merging the melted data frames
#             merged_melted_tables <- rbind(selected_group1[1:55, ], selected_group2[1:55, ])
#             write.table(merged_melted_tables, file=paste(outdir,"/table_of_NumofTrPerCondition_", who,".csv", sep=""), sep="\t", row.names = F, quote=FALSE)
#     
#             ggplot(merged_melted_tables, aes(x = Var1, y = Freq, fill = group)) +
#                    geom_bar(stat = "identity", position = "dodge") +
#                    scale_fill_manual("Groups", values = c("#66CC99", "#877598")) +
#                    theme_bw() +
#                    theme(legend.position="bottom") +
#                    theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#                    labs(title=paste(who,": Number of isoform usage per condition ((subset of up to 55 transcripts))",sep=""), x="Number of isoforms", y="Relative Frequency") +
#              ggsave(file=paste(outdir, "/TrNumPerCondition_", who,".png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
#              rm(selected_group1, selected_group2, merged_melted_tables)
# }
# 
# logUp <- rownames(edger_res_overall)[which(edger_res_overall$logFC >= lfcThreshold & edger_res_overall$FDR <= adjPValueThreshold)]
# logUp <- data.frame(subset(expr_file, annot_gene_id %in% logUp))
# num_of_isoforms_plot(logUp, "LogUp")
# 
# logDown <- rownames(edger_res_overall)[which(edger_res_overall$logFC <= -lfcThreshold & edger_res_overall$FDR <= adjPValueThreshold)]
# logDown <- subset(expr_file, annot_gene_id %in% logDown)
# num_of_isoforms_plot(logDown, "LogDown")
# 
# logNonSign <- rownames(edger_res_overall)[which(edger_res_overall$FDR > adjPValueThreshold)]
# logNonSign <- subset(expr_file, annot_gene_id %in% logNonSign)
# num_of_isoforms_plot(logNonSign, "LogNonSign")
# rm(logUp, logDown, logNonSign, num_of_isoforms_plot)
