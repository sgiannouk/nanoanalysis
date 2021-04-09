### Stavros Giannoukakos ###
### ONT ANALYSIS FOR DIFFERENTIAL GENE EXPRESSION (DGE) 


args <- commandArgs(TRUE)
if (length(args) == 5) {
  # Input filtered matrix (output of step 8)
  matrix <- args[1]
  # CSV file used for running TALON
  input_groups <- args[2]
  # Output direcotry where all stats will be saves
  main_outdir <- args[3]
  # Adjusted p-value threshold for differential expression analysis
  adjPValueThreshold <- as.numeric(args[4])
  # Minimum required log2 fold change for differential expression analysis
  lfcThreshold <- as.numeric(args[5])
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

# matrix <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_3/prefilt_talon_abundance.tsv"
# input_groups <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_3/talon_input.csv"
# main_outdir <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_3/diffExpr_analysis"
# adjPValueThreshold <- 0.05
# lfcThreshold <- 2


suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("DRIMSeq"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape"))
options(ggrepel.max.overlaps = Inf)
options(scipen = 999)



##### DIFFERENTIAL GENE EXPRESSION (DGE) ANALYSIS USING EDGER #####
print("RUNNING DIFFERENTIAL GENE EXPRESSION (DGE) ANALYSIS USING EDGER")
outdir <- file.path(main_outdir, "diffExpr_DGE")
dir.create(outdir, showWarnings = FALSE)

# Input the filtered expression matrix
expr_file <- read.table(matrix, sep="\t", header=T)
expr_file <- expr_file[ ,3:length(expr_file)]

# Omitting the genomic transcripts 
expr_file <- expr_file[which(expr_file$transcript_novelty != "Genomic"), ]

# Obtain number of input samples and unique genes
print(paste("Total number of input samples:", length(expr_file[ ,10:length(expr_file)]), sep=" " ))
print(paste("Total number of unique genes:", length(unique(expr_file$annot_gene_id)), sep=" " ))

# Reading input csv file containing the groups
group_samples <- read.csv(input_groups, header=F)[ ,1:3]
group_samples <- group_samples[order(group_samples$V2), ] # Order data frame
colnames(group_samples) <- c("sample_id", "condition", "batch")

sampletypevalues <- factor(unique(group_samples$condition)) # Obtaining the sample groups
print(paste("Total number of input groups: ", length(sampletypevalues)," (", sampletypevalues[1], "/", sampletypevalues[2],")", sep="" ))


# Obtaining the counts table
matfile <- expr_file[ ,c(1,10:length(expr_file))]  # Obtaining only the annot_gene_id and the counts
colnames(matfile)[1] <- "gene_id"

# Summarising all reads per gene name and removing the duplicate  rows
matfile <- data.frame(dplyr::group_by(matfile, gene_id) %>% dplyr::summarise_all(sum))
row.names(matfile) <- matfile$gene_id; matfile$gene_id <- NULL

# Creating with design matrix
design <- model.matrix( ~ condition, data = group_samples)

# Storing the raw read counts table in a simple list-based data object called a DGEList.
edgeR_table <- DGEList(matfile, group=group_samples$condition)

# Filtering out lowly expressed genes
keep <- filterByExpr(edgeR_table, group=group_samples$condition)
edgeR_table <- edgeR_table[keep, , keep.lib.sizes=FALSE]

# Obtain the filtered counts
filtered_counts <- edgeR_table$counts

# Normalisation for RNA composition by finding a set of scaling factors for the library sizes 
# that minimize the log-fold changes between the samples for most genes. The default method
# for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples.
edgeR_table <- calcNormFactors(edgeR_table, method = "TMM")

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
print(summary(decideTests(qlf, p.value=adjPValueThreshold, lfc=lfcThreshold)))


# Selecting the samples
selected_samples <- colnames(filtered_counts[ ,(which(group_samples$condition==sampletypevalues[1] | group_samples$condition==sampletypevalues[2]))])
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
rm(edgeR_table, results_selected, fit, qlf, keep, selected_samples, edger_res, input_groups, matrix)




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
        geom_hline(yintercept = -lfcThreshold, color="indianred3", size=.3) +
        geom_hline(yintercept = lfcThreshold, color="mediumseagreen", size=.3) +
        theme_minimal() +
        aes(colour=gene) +
        scale_colour_manual(name="Genes", values=colours) +
        labs(title="MA plot - log(FC) vs. log(CPM) on gene level data", x="log(CountsPerMillion)", y="log(FoldChange)")
ggsave(file=paste(outdir,"/edgeR_MAplot.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
rm(logUp, logDown, withStat, colours, gene)


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
       geom_hline(yintercept = -lfcThreshold, color="indianred3", size=.3) +
       geom_hline(yintercept = lfcThreshold, color="mediumseagreen", size=.3) +
       theme_minimal() +
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
abline(v=c(-lfcThreshold,lfcThreshold), col="brown")
abline(h=-log10(adjPValueThreshold), col="brown")
## Plot the names of a reasonable number of transcripts, by selecting those begin not only significant but also having a strong effect size
gn.selected <- (abs(edger_res_overall$logFC)>=lfcThreshold & edger_res_overall$FDR<=adjPValueThreshold)
text(edger_res_overall$logFC[gn.selected], edger_res_overall$significant[gn.selected],lab=edger_res_overall$gene_name[gn.selected], cex=0.6)
dev.off()
rm(genes.to.plot, cols, gn.selected)
edger_res_overall$significant <- NULL; edger_res_overall$pch  <- NULL; #edger_res_overall$gene_name <- NULL



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
       geom_point(size=2.4, shape=15) +
       scale_y_log10(labels = function(x) format(x, scientific = F)) +
       theme_bw() +
       scale_colour_manual(name="", values=c("#7ebdb4", "#f4a548")) +
       theme(axis.text.x = element_text(angle = 45, hjust = 1),
             plot.title = element_text(hjust = 0.5),
             legend.position="bottom") + 
       labs(title="Top 30 Significant DE Genes", x="Genes", y="log10(CPM)")
ggsave(file=paste(outdir,"/edgeR_top30MostSingificantGenes.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
rm(melt_top30_sigNorm, top30_sigNorm)




### PLOTTING THE PERCENTAGE OF ISOFORMS PER NUMBER OF TRANSCRIPTS ###
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
                    min_gene_expr = 5,
                    min_feature_expr = 1)


# Obtaining the counts after dmFiltering
filtered_counts_new <- counts(drimseq)
# Selecting the samples from the first group
selected_group1 <- cbind(filtered_counts_new[ ,c(1,2)], filtered_counts_new[ ,which(colnames(filtered_counts_new) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[1]])])
# Removing unexpressed genes in this group (genes with 0 in all samples)
selected_group1 <- selected_group1[rowSums(selected_group1[ ,which(colnames(selected_group1) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[1]])]) > 0, ] 
# Transform to count table
selected_group1 <- data.frame(table(selected_group1$gene_id))
# Collapsing and counting transcripts
selected_group1 <- data.frame(table(selected_group1$Freq))
# Normalising for total sequencing depth
selected_group1 <- data.frame(cbind(selected_group1[ ,1], round(selected_group1[ ,2]/sum(selected_group1[ ,2])*100,3)))
# selected_group1 <- cbind(selected_group1[ ,1], data.frame(cpm(selected_group1[ ,2])))
colnames(selected_group1) <- c("Var1", "Freq")
selected_group1$group <- as.character("Healthy samples")  # Renaming all genes to samplegroup
selected_group1 <- selected_group1[order(selected_group1$Var1), ]


# Selecting the samples from the second group
selected_group2 <- cbind(filtered_counts_new[,c(1,2)], filtered_counts_new[ ,which(colnames(filtered_counts_new) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[2]])])
# Removing unexpressed genes in this group (genes with 0 in all samples)
selected_group2 <- selected_group2[rowSums(selected_group2[ ,which(colnames(selected_group2) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[2]])]) > 0, ]
# Transform to count table
selected_group2 <- data.frame(table(selected_group2$gene_id))
# Collapsing and counting transcripts
selected_group2 <- data.frame(table(selected_group2$Freq))
# Normalising for total sequencing depth
selected_group2 <- data.frame(cbind(selected_group2[ ,1], round(selected_group2[ ,2]/sum(selected_group2[ ,2])*100,3)))
# selected_group2 <- cbind(selected_group2[ ,1], data.frame(cpm(selected_group2[ ,2],log=T)))
colnames(selected_group2) <- c("Var1", "Freq")
selected_group2$group <- as.character("Cancer samples")  # Renaming all genes to samplegroup
selected_group2 <- selected_group2[order(selected_group2$Var1), ]

# Merging the melted data frames
merged_melted_tables <- rbind(selected_group1[1:101,], selected_group2[1:101,])
merged_melted_tables$group <- factor(merged_melted_tables$group, levels=c("Healthy samples","Cancer samples"))
# rm(selected_group1, selected_group2)

ggplot(merged_melted_tables, aes(x = Var1, y = Freq, fill = group)) + 
       geom_bar(stat = "identity", position="fill", width=1.1, alpha=.9) +
       scale_fill_manual("", values = c("#78c4d4","#8f4f4f")) +
       theme_minimal() +
       scale_x_continuous(limits= c(0, 101), breaks=seq(0,100,5)) +
       # scale_y_continuous(labels = function(x) paste0(x, "%")) +
       theme(panel.border = element_blank(),
             panel.grid.minor = element_blank(),
             panel.grid.major.x = element_blank(),
             panel.grid.major.y = element_line(size=.1, color="gray78"),
             axis.title.x = element_text(vjust=-2),
             legend.title = element_blank(),
             legend.text=element_text(size=11),
             legend.position = "bottom") +
       labs(title="Number of isoform usage per condition",
            x="Number of isoforms",
            y="Relative frequency (percentage)")
ggsave(file=paste(outdir, "/DGE_TrPrcPerCondition.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
rm(list=setdiff(ls(), c("main_outdir", "filtered_counts", "group_samples", "design", "adjPValueThreshold", "lfcThreshold")))




################################################################
################ FUNCTIONAL ENRICHMENT ANALYSIS ################
suppressPackageStartupMessages(library("SPIA"))
suppressPackageStartupMessages(library("DOSE"))
suppressPackageStartupMessages(library("ReactomePA"))
suppressPackageStartupMessages(library("enrichplot"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("clusterProfiler"))
options(scipen = 0)


print("RUNNING FUNCTIONAL ENRICHMENT ANALYSIS")
outdir <- file.path(main_outdir, "diffExpr_DGE/funcEnrichment_analysis")
dir.create(outdir, showWarnings = FALSE)

row.names(filtered_counts) <- gsub("\\..*","",row.names(filtered_counts))
# Matching the ENSEMBL gene_id with the ENTREZ IDs and Symbols
gene.annot <- bitr(row.names(filtered_counts),
                   fromType = "ENSEMBL",
                   toType = c("ENTREZID", "SYMBOL"),
                   OrgDb = "org.Hs.eg.db")

# Removing duplicates
gene.annot <- gene.annot[!duplicated(gene.annot$ENSEMBL), ]
# Merging the count matrix with the obtained annotation
filtered_counts <- merge(filtered_counts, gene.annot, by.x=0, by.y="ENSEMBL")
filtered_counts <- filtered_counts[ ,c(1, length(filtered_counts), (length(filtered_counts)-1),(4:length(filtered_counts)-2))]  # Rearranging the columns
colnames(filtered_counts)[1] <- "ENSEMBL"  # Renaming the 1st column

# Running edgeR again
edgeR_table <- DGEList(filtered_counts[ ,4:length(filtered_counts)], genes=filtered_counts[ ,1:3], group=group_samples$condition)
edgeR_table <- calcNormFactors(edgeR_table)
edgeR_table <- estimateDisp(edgeR_table, design)
fit <- glmQLFit(edgeR_table, design)
qlf <- glmQLFTest(fit)
print(summary(decideTests(qlf, p.value=adjPValueThreshold, lfc=lfcThreshold)))

# Obtaining the most significant genes
res <- topTags(qlf, n=Inf, sort.by="PValue", p.value = adjPValueThreshold)[[1]]
res <- res[abs(res$logFC) >= lfcThreshold, ]
res <- res[,c(3,4)]  # Choosing only ENTREZ and logFC
res <- setNames(res$logFC, res$ENTREZID)
res <- sort(res, decreasing = T)
de.entrez <- names(res)

ego <- enrichGO(de.entrez, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
# Investigate how the significant GO terms are distributed over the GO graph (GO DAG graph)
goplot(ego, showCategory = 6)
ggsave(file=paste(outdir, "/1.enrichment.DAGgraph.png",sep=""), width = 10, height = 8, units = "in", dpi = 900)


# Visualize enriched terms
barplot(ego, showCategory=20)
ggsave(file=paste(outdir, "/2.enrichment.enrichplot-bar.png",sep=""), width = 20, height = 8, units = "in", dpi = 900)

# Dot plot is similar to bar plot with the capability to encode another score as dot size
dotplot(ego, showCategory=30)
ggsave(file=paste(outdir, "/3.enrichment.enrichplot-dot.png",sep=""), width = 20, height = 10, units = "in", dpi = 900)

go <- enrichGO(de.entrez, OrgDb = "org.Hs.eg.db", ont="all")
dotplot(go, split="ONTOLOGY") + 
        facet_grid(ONTOLOGY~., scale="free")
ggsave(file=paste(outdir, "/4.enrichment.enrichplot-gocategories.png",sep=""), width = 20, height = 12, units = "in", dpi = 900)


# Gene-Concept Network
ego2 <- simplify(ego)  # remove redundant GO terms
enrichplot::cnetplot(ego2, showCategory = 5, foldChange=res)
ggsave(file=paste(outdir, "/5.enrichment.geneconcept-network.png",sep=""), width = 10, height = 6, units = "in", dpi = 900)

# UpSet Plot
# The upsetplot is an alternative to cnetplot for visualizing the complex association between 
# genes and gene sets. It emphasizes the gene overlapping among different gene sets.
png(paste(outdir, "/6.enrichment.upset.png",sep=""), units='px', width=1600, height=900, res=100)
upsetplot(ego, n=10)
dev.off()

# Heatmap-like functional classification
heatplot(ego2, showCategory = 15, foldChange=res)
ggsave(file=paste(outdir, "/7.enrichment.heatmap.png",sep=""), width = 18, height = 3, units = "in", dpi = 900)

# The SPIA (Signaling Pathway Impact Analysis) tool can be used to integrate the lists of 
# differentially expressed genes, their fold changes, and pathway topology to identify affected pathways
spia_result <- spia(de=res, all=as.vector(filtered_counts$ENTREZID), organism="hsa", nB=4000)

# Output of the SPIA significant results
spia_result.sign <- spia_result[spia_result$pGFdr<=adjPValueThreshold, ]
write.table(spia_result.sign, file=paste(outdir,"/spia.significant.csv", sep=""), sep="\t", row.names = F, quote=FALSE)

png(paste(outdir, "/8.enrichment.spia-results.png",sep=""), units='px', width=1600, height=900, res=100)
plotP(spia_result, threshold=0.05)
dev.off()
