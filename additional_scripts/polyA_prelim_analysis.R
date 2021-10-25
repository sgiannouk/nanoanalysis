### Stavros Giannoukakos ###
### ONT POLYA-TAIL ANALYSIS

args <- commandArgs(TRUE)
if (length(args) == 2) {
  # Input/Output direcotry
  polyAdir <- args[1]
  # Number of cores to use
  cores.to.use = args[2]
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

# polyAdir <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_3/polyA_estimation_new"
# cores.to.use <- 2  # Number of cores to use



suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("doParallel"))
registerDoParallel(cores=cores.to.use)


### POLYA-TAIL ANALYSIS ###
print("RUNNING POLYA-TAIL PRELIMINARY ANALYSIS")

outdir <- file.path(polyAdir, "dpa_results")
dir.create(outdir, showWarnings = F)
setwd(outdir)


# Input of all polyA estimations from all samples into a single data frame
polyAfiles <- list.files(path = polyAdir, pattern = "*.transcripts.pass.tsv", all.files = T, full.names = T, recursive = T)
dataset <- ldply(polyAfiles, read.table, sep = "\t", header=T, .parallel=T)
rm(polyAfiles, cores.to.use)

# Pro-filtering QC
data_table = data.frame(table(dataset$qc_tag))
data_table$Var1 <- factor(data_table$Var1, levels = c("SUFFCLIP", "NOREGION", "ADAPTER", "READ_FAILED_LOAD", "PASS"))
ggplot(data_table, aes(x=Var1, y=Freq, fill=Var1)) + 
       geom_bar(stat="identity", width=.3) +
       theme_minimal() +
       scale_fill_brewer(palette="Set3") +
       theme(panel.border = element_blank(),
             panel.grid.minor = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.major.x = element_line(size=.1, color="gray78"),
             axis.title.x = element_text(vjust=-2),
             legend.title = element_blank(),
             legend.text=element_text(size=11),
             legend.position = "None") +
       labs(title="Overall polyA QC tags (prior filtering)", x="", y="Frequency") +
       coord_flip()
ggsave(file=paste(outdir, "/qc_prefilter.png",sep=""), width = 10, height = 3, units = "in", dpi = 900)
rm(data_table)


# Barplot showing qc_tag content
gc_tag_mat <- plyr::count(dataset, c("sample", "qc_tag"))
gc_tag_mat$qc_tag <- factor(gc_tag_mat$qc_tag, levels = c("SUFFCLIP", "NOREGION", "READ_FAILED_LOAD", "ADAPTER", "PASS"))
ggplot(gc_tag_mat, aes(fill=qc_tag, y=freq, x=sample)) + 
       geom_bar(stat="identity", position="fill") +
       theme_minimal() +
       scale_y_continuous(labels = scales::percent) +
       scale_fill_brewer(palette="Set3") +
       labs(title="PolyA QC tags across samples (filtered)", x="Sample names", y="Frequency (%)", fill = "QC tags") +
       theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file=paste(outdir, "/qc_tagsplot_filtered.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
rm(gc_tag_mat)

# Keeping only "PASS" transcripts
dataset <- dataset[dataset$qc_tag=="PASS", ]
dataset <- dataset[with(dataset, order(sample, group, contig)), ]
write.table(dataset, file=paste(outdir,"/all_tails.tsv", sep=""), sep="\t", row.names = F, quote=FALSE)


# Plotting Violin chart
ggplot(dataset, aes(x=sample, y=polya_length, fill=group)) +
       geom_violin(trim=FALSE) +
       geom_boxplot(width=0.1) +
       theme_minimal() +
       scale_fill_brewer(palette="Set1") +
       theme(panel.border = element_blank(),
             panel.grid.minor = element_blank(),
             panel.grid.major.x = element_blank(),
             panel.grid.major.y = element_line(size=.1, color="gray78"),
             axis.title.x = element_text(vjust=-2),
             legend.title = element_blank(),
             legend.text=element_text(size=11),
             legend.position = "None") +
        labs(title="Violin plot of polyA lenghts", x="Samples", y="PolyA length")
ggsave(file=paste(outdir, "/qc_polyAdistrib-violin_filtered.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)


# Plotting distribution
# ggplot(distrib_table, aes(x=polya_length, fill=group)) + 
#        geom_histogram(binwidth=.5, alpha=.5, position="identity")

