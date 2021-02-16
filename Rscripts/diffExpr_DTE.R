### Stavros Giannoukakos ###
### ONT ANALYSIS FOR DIFFERENTIAL TRANSCRIPT EXPRESSION (DTE)


args <- commandArgs(TRUE)
if (length(args) == 6) {
  # Input filtered matrix (output of step 8)
  matrix <- args[1]
  # CSV file used for running TALON
  input_groups <- args[2]
  # Output direcotry where all stats will be saves
  main_outdir <- args[3]
  # Minimum transcript counts
  minFeatureExpr <- as.numeric(args[4])
  # Adjusted p-value threshold for differential expression analysis
  adjPValueThreshold <- as.numeric(args[5])
  # Minimum required log2 fold change for differential expression analysis
  lfcThreshold <- as.numeric(args[6])
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

matrix <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_3/filt_talon_abundance.csv"
input_groups <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_3/talon_input.csv"
main_outdir <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_3/diffExpr_analysis"
minFeatureExpr <- 10  # Minimum transcript counts
adjPValueThreshold <- 0.05
lfcThreshold <- 2


library("readr")
library("dplyr")
library("edgeR")
library("ggpubr")
library("DRIMSeq")
library("ggplot2")
library("reshape")
library("patchwork")
library("ggchicklet")
options(scipen = 999)



##### DIFFERENTIAL TRANSCRIPT EXPRESSION (DTE) ANALYSIS USING DRIMSEQ/EDGER #####
print("RUNNING DIFFERENTIAL TRANSCRIPT EXPRESSION (DTE) ANALYSIS USING DRIMSEQ/EDGER")

outdir <- file.path(main_outdir, "diffExpr_DTE")
dir.create(outdir, showWarnings = FALSE)
setwd(outdir)

# Input the filtered expression matrix
expr_file <- read.csv(matrix, header=T)
# Exporting the unfiltered table
write.table(expr_file, file=paste(outdir,"/edgeR_filtered_table.tsv", sep=""), sep="\t", row.names = F, quote=FALSE)

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
print(paste("Total number of input groups: ", 
            length(sampletypevalues)," (", sampletypevalues[1], " and ", sampletypevalues[2],")", sep="" ))


# Obtaining the counts table
matfile <- expr_file[ ,c(2,1,10:length(expr_file))]  # Obtaining only the annot_transcript_id and the counts
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
print(paste("In total ", length(unique(filtered_counts$feature_id)), "/", length(unique(matfile$feature_id)),
            " transcripts passed the min. threshold of ", minFeatureExpr," counts.",sep=""))

# Creating with design matrix
if (length(unique(group_samples$batch) == 1)) {
  design <- model.matrix( ~ condition, data = DRIMSeq::samples(drimseq))
} else {
  design <- model.matrix( ~ condition + batch, data = DRIMSeq::samples(drimseq)) }


# Removing gene_id column
transcript_filt_counts <- filtered_counts; transcript_filt_counts$gene_id <- NULL

# Making the annot_transcript_id as row names and removing it as from 1st column
row.names(transcript_filt_counts) <- transcript_filt_counts$feature_id; transcript_filt_counts$feature_id <- NULL

# Storing the raw read counts table in a simple list-based data object called a DGEList.
edgeR_table <- DGEList(transcript_filt_counts, group=group_samples$condition)

# Normalisation for RNA composition by finding a set of scaling factors for the library sizes 
# that minimize the log-fold changes between the samples for most transcripts. The default method
# for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples.
edgeR_table <- calcNormFactors(edgeR_table)

# Estimating common dispersion and tagwise dispersions 
edgeR_table <- estimateDisp(edgeR_table, design)

# # Plotting the per-transcript dispersion estimates
# png(paste(outdir,"/edgeR_DispersionPlot.png",sep=""), units='px', height=900, width=1600, res=90)
# plotBCV(edgeR_table)
# title(main = "Dispersion Estimates")
# dev.off()

# Performing quasi-likelihood F-tests
fit <- glmQLFit(edgeR_table, design)

# # Plotting the quasi-likelihood dispersion  estimates
# png(paste(outdir,"/edgeR_QLDispersionPlot.png",sep=""), units='px', height=900, width=1600, res=90)
# plotQLDisp(fit)
# title(main = "QL Dispersion Estimates")
# dev.off()

qlf <- glmQLFTest(fit)

# Obtaining the results table and sorting by p-value
edger_res_overall <- topTags(qlf, n=Inf, sort.by="PValue")[[1]]
edger_res <- edger_res_overall[,c(1,4,5)];  colnames(edger_res) <- c("log2FoldChange",  "pval", "padj")

# Output basic stats
print(paste("Results indicating the transcript set with log2FoldChange value > |", lfcThreshold,"| and with adjasted p-value < ", adjPValueThreshold, sep=""))
is.de <- decideTests(qlf, p.value=adjPValueThreshold, lfc=lfcThreshold)
print(summary(is.de))



# Selecting the samples
selected_samples <- colnames(transcript_filt_counts[ ,(which(group_samples$condition==sampletypevalues[1] | group_samples$condition==sampletypevalues[2]))])
# Merging cpm values and stats from EdgeR
total_results <- data.frame(merge(cpm(edgeR_table)[ ,selected_samples], edger_res, by=0))
# Importing the gene names and remove duplicates
total_results <- merge(expr_file[ ,c("annot_transcript_id", "annot_transcript_name", "annot_gene_id", "annot_gene_name")], total_results, by.y="Row.names", by.x="annot_transcript_id", all.y=TRUE, no.dups=TRUE)
total_results <- total_results[!duplicated(total_results$annot_transcript_id),  ]
# Renaming the first two columns
colnames(total_results)[1:4] <- c("transcript_id", "transcript_name", "gene_id", "gene_name")
# Sorting the table by pvalue and adjasted p-value
total_results <- total_results[with(total_results, order(pval, padj)), ]
# Exporting the normalised results table containing the selected features with log2FoldChange greater than 1 and adjusted p value lower than the user-input value
write.table(total_results, file=paste(outdir,"/",sampletypevalues[1],"VS",sampletypevalues[2],"_edgeR_allTranscripts.csv", sep=""), sep="\t", row.names = F, quote=FALSE)

# Output only results where logFC 
results_selected <- total_results[(abs(total_results$log2FoldChange)>=lfcThreshold) & (total_results$padj<=adjPValueThreshold), ]
# Sorting the table by pvalue and adjasted p-value
results_selected <- results_selected[with(results_selected, order(pval, padj)), ]
# Exporting the normalised results table containing the selected features with log2FoldChange greater than 1 and adjusted p value lower than the user-input value
write.table(results_selected, file=paste(outdir,"/",sampletypevalues[1],"VS",sampletypevalues[2],"_edgeR_topTranscriptsBelow", gsub("[.]", "", adjPValueThreshold), "LFC", gsub("[.]", "", lfcThreshold), ".csv", sep=""), sep="\t", row.names = F, quote=FALSE)

# Export a list of all novel differentially expressed isoforms to proceed for Functional Annotation
list_of_de_transcripts <- results_selected[,1]  # List of all transcripts that were differentially expressed
novel_de_transcripts <- grep('TALONT', list_of_de_transcripts, value=TRUE)  # Obtaining only the novel transcripts that were DEd
write.table(novel_de_transcripts, file=paste(outdir,"/novel_de_transcripts_for_funcAnnotation.tsv", sep=""), sep="\t", row.names = F, col.names = F, quote=FALSE)



# MA-plot - Drawing the expression levels over the exons to highlight differential exon usage
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
        scale_colour_manual(name="Transcripts", values=colours) +
        labs(title="MA plot - log(FC) vs. log(CPM) on transcript level data", x="log(CountsPerMillion)", y="log(FoldChange)")
ggsave(file=paste(outdir,"/edgeR_MAplot.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
rm(logUp, logDown, withStat, colours, gene)

# MA-plot
withNovel <- expr_file[grepl("^TALON[^*]+$", expr_file$annot_transcript_id), ]
withNovel <- as.vector(withNovel$annot_transcript_id)
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
p1 <- ggplot(edger_res_overall, aes(x=logCPM, y=logFC)) + 
             geom_point(aes(shape=factor(novelty,  labels=c("Known", "Novel")))) +
             geom_hline(yintercept = -lfcThreshold, color="indianred3") +
             geom_hline(yintercept = lfcThreshold, color="mediumseagreen") +
             theme_minimal() +
             aes(colour=gene) +
             scale_colour_manual(name="Transcripts", values=colours) +
             labs(title="MA plot - log(FC) vs. log(CPM) on transcript level data with Known/Novel-transcript indication", 
                  x="log(CountsPerMillion)", y="log(FoldChange)", shape="Known/Novel") +
             theme(legend.position = "right")
ggsave(file=paste(outdir,"/edgeR_MAplot_withTranscriptIndications.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)


# Obtaining the novelty information of the significant transcript
significant_overview  <- expr_file[expr_file$annot_transcript_id %in% results_selected$transcript_id,c(8:9)]
significant_overview$ISM_subtype <- sub("None", "", significant_overview$ISM_subtype)
significant_overview <- data.frame(paste(significant_overview$transcript_novelty, significant_overview$ISM_subtype, sep=" "))
# Transform to count table
significant_overview <- data.frame(table(significant_overview[,1]))
# Rename ISM to ISM None as it should be
significant_overview$Var1 <- gsub("ISM.$", "ISM None", significant_overview$Var1)
# Get the percentage of each group and sort the groups
significant_overview$perc <-round(significant_overview[ ,2]/sum(significant_overview[ ,2])*100,3)
significant_overview <- significant_overview[order(-significant_overview$Freq), ]
significant_overview$Var1 <- factor(significant_overview$Var1, levels=rev(significant_overview$Var1))
# Choosing colours for each group
legend_colors <- c("#ae4544", "#d8cb98", "#a4ad6f", "#cc7c3a", "#436f82", "#7c5981", "#cccccc", "#ffcda3", "#ef4f4f")

p2 <- ggplot(significant_overview, aes(x=Var1, y=perc, fill=Var1)) +
            geom_chicklet(width = .45) +
             theme_minimal() +
             coord_flip() +
             scale_y_continuous(limits= c(0, round(max(significant_overview$perc), digits = -1)), 
                                breaks=seq(0, round(max(significant_overview$perc), digits = -1),5)) +
             theme(plot.title = element_text(colour = "#a6a6a4", size=11),
                   panel.grid.major.y = element_blank(),
                   panel.border = element_blank(),
                   axis.ticks.y = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.grid.major.x = element_line(size=.1, color="gray78"),
                   legend.position = "None") +
             labs(x="", y="Percentage of transcripts") +
             scale_fill_manual("", values = legend_colors[1:length(significant_overview$Var1)])

figure <- p1 + p2 + plot_layout(ncol = 1, heights=c(1,.25))
ggsave(file=paste(outdir,"/edgeR_MAplot_withTranscriptIndicationsNgroupSun.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
rm(logUp, logDown, withStat, colours, gene,withNovel, p1, p2)

# Volcano plot
# Compute significance, with a maximum of 350 for the p-values set to 0 due to limitation of computation precision
edger_res_overall <- merge(expr_file[ ,c("annot_transcript_id", "annot_transcript_name")], edger_res_overall, by.x="annot_transcript_id", by.y = 0, all.y=T)
edger_res_overall <- edger_res_overall[!duplicated(edger_res_overall$annot_transcript_id),  ]
row.names(edger_res_overall) <- edger_res_overall$annot_transcript_id; edger_res_overall$annot_transcript_id <- NULL; colnames(edger_res_overall)[1] <- "transcript_name"
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
     xlab="Effect size: log2(Fold-Change)",
     ylab="-log10(adjusted p-value)",
     pch=edger_res_overall$pch, cex=0.4)
abline(v=0)
abline(v=c(-lfcThreshold,lfcThreshold), col="brown")
abline(h=-log10(adjPValueThreshold), col="brown")
## Plot the names of a reasonable number of transcripts, by selecting those begin not only significant but also having a strong effect size
gn.selected <- (abs(edger_res_overall$logFC)>=lfcThreshold & edger_res_overall$FDR<=adjPValueThreshold)
text(edger_res_overall$logFC[gn.selected], edger_res_overall$significant[gn.selected],lab=edger_res_overall$transcript_name[gn.selected], cex=0.6)
dev.off()
rm(genes.to.plot, cols, gn.selected)
edger_res_overall$significant <- NULL; edger_res_overall$pch  <- NULL; edger_res_overall$transcript_name <- NULL
rm(is.de,  fit, qlf, transcript_filt_counts, edger_res, edger_res_overall, edgeR_table, design)


### TOP 30 DE TRANSCRIPTS
# Plotting the normalized count values for the top 30 differentially expressed transcripts (by padj values)
# Plotting the normalized count values for the top 30 differentially expressed genes (by padj values)
## Order results by by pvalue and adjasted p-value values and get the first 30 genes
total_results <- total_results[with(total_results, order(pval, padj)), ]
# Obtaining the top 30 most significant genes and the per sample normalised (CPM) data
top30_sigNorm <- total_results[1:30, !names(total_results) %in% c("gene_name", "gene_id", "log2FoldChange", "pval", "padj")]
# Gathering the columns to have normalized counts to a single column
melt_top30_sigNorm <- melt(top30_sigNorm)
melt_top30_sigNorm$group <- gsub("_[^_]+$", "", melt_top30_sigNorm$variable)
melt_top30_sigNorm$value[melt_top30_sigNorm$value == 0] <- 0.1

ggplot(melt_top30_sigNorm, aes(x = transcript_name, y = value, color = group)) +
  geom_point(size=2.5) +
  scale_y_log10(labels = function(x) format(x, scientific = F)) +
  theme_bw() +
  scale_colour_manual(name="", values=c("#66CC99", "#877598")) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom") +
  labs(title="Top 30 Significant DE Transcripts", x="Transcripts", y="log10(CPM)")
ggsave(file=paste(outdir,"/edgeR_top30MostSingificantTranscripts.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
rm(melt_top30_sigNorm, top30_sigNorm, n, input_groups, legend_colors, minFeatureExpr, 
   significant_overview, filtered_counts,drimseq, figure, matfile, matrix)
