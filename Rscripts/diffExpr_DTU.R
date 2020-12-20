### Stavros Giannoukakos ###
### ONT ANALYSIS FOR DIFFERENTIAL TRANSCRIPT USAGE (DEU) AND DIFFERENTIAL EXON USAGE (DEU)


args <- commandArgs(TRUE)
if (length(args) == 5) {
  # Input filtered matrix (output of step 8)
  matrix <- args[1]
  # CSV file used for running TALON
  input_groups <- args[2]
  # Output directory where all stats will be saved
  main_outdir <- args[3]
  # Fasta file with spliced exons for each transcript
  ref_transcriptome <- args[4]
  # Transcriptome annotation from the TALON database
  ref_annotation <- args[5]
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

# matrix <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_2/new_filter/filt_talon_abundance.csv"
# input_groups <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_2/talon_input.csv"
# main_outdir <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_2/diffExpr_analysis"
# ref_transcriptome <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_2/new_filter/reference_transcriptome.fasta"
# ref_annotation <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation_2/new_filter/database_talon.gtf"

library("IsoformSwitchAnalyzeR")
options(scipen = 999)


##### DIFFERENTIAL EXON USAGE (DEU) ANALYSIS USING IsoformSwitchAnalyzeR #####
print("RUNNING DIFFERENTIAL EXON USAGE (DEU) ANALYSIS USING IsoformSwitchAnalyzeR")

outdir <- file.path(main_outdir, "DiffTranscriptUsage_NEWWWW")
dir.create(outdir, showWarnings = FALSE)
setwd(outdir)


# Reading input csv file containing the groups
group_samples <- data.frame(readr::read_csv(file = input_groups, col_names=F)[ ,1:2])
colnames(group_samples) <- c("sampleID","condition")
group_samples <- group_samples[order(group_samples$condition), ] # Order data frame
# Preparing a dataframe containing information about the samples 
sampletypevalues <- factor(unique(group_samples$condition)) # Obtaining the sample groups


# Importing the talon filtered matrix
matfile <- data.frame(readr::read_csv(file = matrix))
rownames(matfile) <- matfile$annot_transcript_id
tallonCols <- c("annot_gene_id",
                "annot_transcript_id",
                "annot_gene_name",
                "annot_transcript_name",
                "n_exons",
                "length",
                "gene_novelty",
                "transcript_novelty",
                "ISM_subtype")

# Removing the annotation columns
expression <- matfile[ ,setdiff(colnames(matfile), tallonCols)]

### Create the switchAnalyzeRlist with all necessary files
talonSwitch <- importRdata(isoformCountMatrix   = expression,
                           designMatrix         = group_samples,
                           isoformExonAnnoation = ref_annotation,
                           isoformNtFasta       = ref_transcriptome)
summary(talonSwitch)
rm(expression, matfile, tallonCols, ref_annotation, ref_transcriptome, input_groups, matrix, main_outdir)

# Pre-filtering step. Remove single isoform genes or non-expressed isoforms
talonSwitch <- preFilter(switchAnalyzeRlist       = talonSwitch,
                         # geneExpressionCutoff     = 5,
                         removeSingleIsoformGenes = TRUE)

talonSwitch <- isoformSwitchAnalysisPart1(switchAnalyzeRlist   = talonSwitch,
                                          switchTestMethod     = 'DRIMSeq',
                                          orfMethod            = "longest",
                                          alpha                = 0.05,
                                          dIFcutoff            = 0.1,
                                          pathToOutput         = outdir,
                                          outputSequences      = TRUE, 
                                          prepareForWebServers = FALSE)


print(extractSwitchSummary(talonSwitch))
# Exporting the IsoformSwitchAnalyzeR test results
write.table(talonSwitch$isoformFeatures, file=paste(outdir,"/isoformSwitchAnalysis_results.csv", sep=""), sep=",", row.names=F, quote=F)

# Creating the external files that need to be filled by

print(paste("Pfam -- run:$ pfam_scan.pl -fasta ",outdir,"/isoformSwitchAnalyzeR_isoform_AA.fasta -dir /home/stavros/playground/progs/PfamScan/databases", sep=""))
file.create(paste(outdir,"/results_pfam.txt", sep=""))


print("CPAT -- nt | http://lilab.research.bcm.edu/cpat/")
file.create(paste(outdir,"/results_cpat.txt", sep=""))

print("IUPred2A -- AA | https://iupred2a.elte.hu/")
file.create(paste(outdir,"/results_IUPred2A.txt", sep=""))

print("SignalP -- AA | http://www.cbs.dtu.dk/services/SignalP/")
file.create(paste(outdir,"/results_SignalP.txt", sep=""))


readline(prompt="Press [enter] to continue")

talonSwitch <- isoformSwitchAnalysisPart2(switchAnalyzeRlist       = talonSwitch,
                                          alpha                    = 0.05,
                                          dIFcutoff                = 0.1,
                                          n                        = Inf,
                                          removeNoncodinORFs       = TRUE,
                                          codingCutoff             = 0.5,  # For CPC2: The cutoff suggested is 0.5 for all species
                                          pathToPFAMresultFile     = paste(outdir,"/results_pfam.txt", sep=""),
                                          # pathToCPC2resultFile = paste(outdir,"/results_cpc2.csv", sep=""),
                                          pathToCPATresultFile     = paste(outdir,"/results_cpat.txt", sep=""),
                                          pathToIUPred2AresultFile = paste(outdir,"/results_IUPred2A.txt", sep=""),
                                          pathToSignalPresultFile  = paste(outdir,"/results_SignalP.txt", sep=""),
                                          consequencesToAnalyze    = c('intron_retention',
                                                                       'coding_potential',
                                                                       'ORF_seq_similarity',
                                                                       'NMD_status',
                                                                       'domains_identified',
                                                                       'IDR_identified',
                                                                       'IDR_type',
                                                                       'signal_peptide_identified'),
                                          pathToOutput            = outdir,
                                          fileType                = 'png',
                                          asFractionTotal         = FALSE,
                                          outputPlots             = FALSE,
                                          quiet                   = FALSE)

# Extracting the (top) switching genes/isoforms (with functional consequences).
extractTopSwitches(talonSwitch, filterForConsequences = TRUE, n=Inf)
isoswitch_genes <- extractTopSwitches(talonSwitch, filterForConsequences = TRUE, n=Inf)$gene_name
# The number of isoform switches with predicted functional consequences 
extractSwitchSummary(talonSwitch, filterForConsequences = TRUE)

for (gene_name in isoswitch_genes){
     print(gene_name)
     png(paste(outdir,"/switchPlot_", gene_name,".png",sep=""), units='px', height=1400, width=2900, res=290)
     switchPlot(talonSwitch, gene=gene_name)
     dev.off()
}

png(paste(outdir,"/common_switch_consequences.png",sep=""), units='px', height=800, width=1800, res=100)
extractConsequenceSummary(talonSwitch,
                          consequencesToAnalyze='all',
                          plotGenes = FALSE,
                          asFractionTotal = FALSE,
                          removeEmptyConsequences = TRUE,
                          localTheme=theme_bw())
dev.off()

# png(paste(outdir,"/switchPlot_", gene_name,".png",sep=""), units='px', height=900, width=1600, res=100)
# extractConsequenceEnrichment(talonSwitch,
#                              consequencesToAnalyze='all',
#                              analysisOppositeConsequence = TRUE,
#                              returnResult = FALSE)

# 
# extractSplicingEnrichment(talonSwitch,
#                           returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
# )
# 
# ggplot(data=talonSwitch$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
#   geom_point(aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), size=1 )+
#   geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
#   geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
#   facet_wrap( ~ condition_2) +
#   # facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
#   labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
#   theme_bw()
