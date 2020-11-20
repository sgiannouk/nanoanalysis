### Stavros Giannoukakos ###
### ONT ANALYSIS transcript type analysis 

args <- commandArgs(TRUE)
if (length(args) == 2) {
  matrix = args[1]
  outdir <- args[2]
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

# outdir <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_2"
# matrix <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_2/perTranscript_expression_matrix.csv"


library("dplyr")
library("ggplot2")
library("reshape2")
library("RColorBrewer")
setwd(outdir)

# Input the edited expression matrix
expr_file <- read.csv(matrix, header=TRUE, row.names=1)
# Remove transcript names
row.names(expr_file) <- NULL
# Convert sample columns into numeric
expr_file[ ,2:ncol(expr_file)] <- sapply(expr_file[ ,2:ncol(expr_file)], as.numeric)
# Collapsing  all reads per categories
expr_file <- aggregate(expr_file[ ,2:ncol(expr_file)], by=list(Category=expr_file[,1]), FUN=sum)
# Output the matrix
write.table(expr_file, file=paste(outdir,"/transcript-type_summarisation.tsv", sep=""), sep="\t", row.names = F, quote=FALSE)
# Remove lowly expressed categories 
expr_file <- expr_file[rowSums(expr_file[2:ncol(expr_file)]) > 20, ]
# Order categories based on sum of rows
expr_file <- expr_file[order(rowSums(expr_file[2:ncol(expr_file)])), ]
# Divide each column by its sum (getting percentages)
expr_file[ ,2:ncol(expr_file)] <- as.data.frame(lapply(expr_file[ ,2:ncol(expr_file)], function(x) x/sum(x)))


# Prepare the matrix for plotting
mtx <- melt(expr_file, id.vars="Category")
# Keep the correct order
mtx$variable <- factor(mtx$variable, levels=colnames(expr_file)[2:length(colnames(expr_file))])

num <- length(unique(expr_file$Category))

# ggplot 
ggplot(mtx, aes(x = variable  , y = value, fill = factor(Category, levels=expr_file$Category))) +
geom_bar(stat = "identity", position = "stack", width = 0.5) +
theme_bw() +
scale_fill_manual("", values = brewer.pal(num, "Spectral")) +
scale_y_continuous(limits= c(0, 1), breaks=seq(0,1,.1), labels = scales::percent_format(accuracy = 5L)) +
theme(text=element_text(family="Source Sans Pro"),
      plot.title = element_text(colour = "#a6a6a4", size=13),
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(size=.1, color="gray78"),
      axis.title.x = element_text(vjust=-2),
      legend.title = element_blank(),
      legend.text=element_text(size=11),
      legend.position = "bottom") +
labs(title="Transcript type summarisation (ref. genome/Talon)",
     x="",
     y="Percentage of reads (%)")
ggsave(file="transcript_type_summarisation.png", width = 10, height = 6, units = "in", dpi = 1200)
