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

# outdir <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/R_analysis/R_script"
# matrix <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/R_analysis/R_script/perGene_expression_matrix.csv"

library("pals")
library("dplyr")
library("plotly")
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
# Divide each column by its sum (getting percentages)
expr_file[ ,2:ncol(expr_file)] <- as.data.frame(lapply(expr_file[ ,2:ncol(expr_file)], function(x) x/sum(x)))


# Prepare the matrix for plotting
mtx <- melt(expr_file, id.vars="Category")
# Keep the correct order
mtx$variable <- factor(mtx$variable, levels=sort(colnames(expr_file)[2:length(colnames(expr_file))]))

# Colouring 
num <- length(unique(expr_file$Category))
colors  <- colorRampPalette(brewer.pal(8, "Dark2"))(num)

# ggplot 
p <- ggplot(mtx, aes(x = variable  , y = value, fill = Category)) +
     geom_bar(stat = "identity", position = "stack", width = 0.3) +
     theme_bw() +
     scale_fill_manual("", values = colors) +
     scale_y_continuous(labels = scales::percent) +
     theme(legend.position = "bottom") +
     ylab("Percentage of reads (%)") +
     xlab("") +
     ggtitle("Transcript type summarisation  (ref. genome)")
ggsave(file="transcript_type_summarisation.png", width = 10, height = 6, units = "in", dpi = 1200)

ptly <- ggplotly(p, originalData = T, dynamicTicks = T) %>% 
        layout(yaxis = list(tickformat = "%"))

htmlwidgets::saveWidget(ptly, "transcript_type_summarisation.html")
