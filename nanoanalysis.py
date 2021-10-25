###Stavros Giannoukakos### 
#Version of the program
__version__ = "v0.2.2"

import time
import argparse
import subprocess
# from Bio import SeqIO
from pathlib import Path
from datetime import datetime
from collections import Counter
import shutil, glob, sys, os


ont_data =  "/home/stavros/playground/ont_basecalling/guppy_v3_basecalling"
annotation = {"NonTransf":"control", "Tumour":"treatment"}
chosen_samples = ("NonTransf_1",  "NonTransf_2",  "NonTransf_3", "Tumour_1",  "Tumour_2",  "Tumour_4")
# chosen_samples = ("Sample_3") #"Sample_2" , 
# chosen_samples = ("NonTransf_1",  "NonTransf_2",  "NonTransf_3", "Tumour_1",  "Tumour_2",  "Tumour_4", "Sample_3", "VU001PL")


### References
refGenomeGRCh38 = "/home/stavros/references/reference_genome/GRCh38_GencodeV31_primAssembly/GRCh38.primary_assembly.genome.fa"
refGenomeGRCh38_traclean = "/home/stavros/references/reference_genome/GRCh38_GencodeV31_primAssembly/GRCh38.primary_assembly.genome.edited.fa"
refAnnot = "/home/stavros/references/reference_annotation/GRCh38_gencode.v35.primary_assembly.annotation.gtf"
refAnnot_bed = "/home/stavros/references/reference_annotation/GRCh38_gencode.v35.primary_assembly.annotation.bed"
refAnnot_bed12 = "/home/stavros/references/reference_annotation/bed12/gencode.v35.primary_assembly.annotation.bed12"
ensembl_db = "/home/stavros/.cache/pyensembl/GRCh38/ensembl95/agfusion.homo_sapiens.95.db"
### R Scripts
rscripts = "{0}/Rscripts".format(os.path.dirname(os.path.realpath(__file__)))
### JAFFA directory
jaffa_pipeline = "/home/stavros/playground/progs/JAFFA_v2.1"
### Talon analysis
transcriptclean = "python3 /home/stavros/playground/progs/TranscriptClean/TranscriptClean.py"  ### TranscriptClean
sj_ref = "/home/stavros/playground/progs/TranscriptClean/GRCh38_SJs.ref"  ### TranscriptClean ref. SJ
map_antisense = "python3 /home/stavros/playground/progs/TALON/src/talon/post/map_antisense_genes_to_sense.py"
### Trinotate analysis
trinity_sqlite = "/home/stavros/references/trinity_db/Trinotate.sqlite"
pfam_dir = "/home/stavros/playground/progs/PfamScan/databases"
pfam_db = "/home/stavros/references/trinity_db/Pfam-A.hmm"
uniprot = "/home/stavros/references/trinity_db/uniprot_sprot.pep"
tmhmm_exe = "/home/stavros/playground/progs/tmhmm-2.0c/bin/tmhmm"
rnammer_exe = "/home/stavros/playground/progs/rnammer_v1.2/rnammer"
rnammerTransc = "/home/stavros/playground/progs/Trinotate/util/rnammer_support/RnammerTranscriptome.pl"



usage = "nanodRNA_analysis [options]"
epilog = " -- June 2019 | Stavros Giannoukakos -- "
description = "DESCRIPTION"

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage, description=description, epilog=epilog)
# Number of threads/CPUs to be used
parser.add_argument('-t', '--threads', dest='threads', default=str(60), metavar='', 
                	help="Number of threads to be used in the analysis")
# Adjusted p-value threshold for differential expression analysis
parser.add_argument('-adjpval', '--adjPValueThreshold', dest='adjPValueThreshold', default=str(0.05), metavar='', 
                	help="Adjusted p-value threshold for differential\nexpression analysis")
# Minimum required log2 fold change for differential expression analysis
parser.add_argument('-lfc', '--lfcThreshold', dest='lfcThreshold', default=str(2), metavar='', 
                	help="Minimum required log2 fold change for diffe-\nrential expression analysis")
# Minimum counts for the ISM group
parser.add_argument('-ith', '--ismTheshold', dest='ismTheshold', default=str(20), metavar='', 
                	help="Min. counts filtering threshold for the ISM\nSuffix transcripts (either group)")
# Minimum polyA counts
parser.add_argument('-mpa', '--minPolyA', dest='minPolyA', default=str(10), metavar='', 
                	help="Min. number of transcripts containing a polyA estimation")
# Top N genes to be used for the heatmap
parser.add_argument('-n', '--n_top', dest='n_top', default=str(60), metavar='', 
                	help="Top N genes to be used for the heatmap")
# Display the version of the pipeline 
parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}\n{epilog}')
# Get the options and return them
args = parser.parse_args()

script_dir = os.path.dirname(os.path.realpath(__file__))
startTime = datetime.now()

# Main folder hosting the analysis
analysis_dir = os.path.join(script_dir, "batch3_analysis_v3.6.1")
# analysis_dir = os.path.join(script_dir, "batch1_analysis_v3.6.1")
# analysis_dir = os.path.join(script_dir, "batch1Nbatch3_analysis_v3.6.1")

prepr_dir = os.path.join(analysis_dir, "preprocessed_data")
alignments_dir = os.path.join(analysis_dir, "alignments")
reports_dir = os.path.join(analysis_dir, "reports")
# Reporting directories
initial_qc_reports = os.path.join(analysis_dir, "reports/QC_initial")
postAlignment_reports = os.path.join(analysis_dir, "reports/QC_post")
coverage_reports = os.path.join(analysis_dir, "reports/coverage")
pipeline_reports = os.path.join(analysis_dir, "reports/pipeline_reports")
# Downstream analysis directories
expression_analysis_dir = os.path.join(analysis_dir, "expression_analysis")
fusions_events_jaffa = os.path.join(analysis_dir, "jaffa_fusion_analysis")
antisense_dir = os.path.join(analysis_dir, "expression_analysis/antisense_genes_analysis")
differential_expression = os.path.join(analysis_dir, "expression_analysis/differential_expression")
polyA_analysis_dir = os.path.join(analysis_dir, "polyA_estimation")
# Subdirectories
methylation_dir = os.path.join(analysis_dir, "methylation_analysis")
eligos_methylation_dir = os.path.join(analysis_dir, "eligos2_meth_analysis_m6a")
dte = os.path.join(analysis_dir, "expression_analysis/differential_expression/diffExpr_DTE")
predict_product = os.path.join(differential_expression, "diffExpr_DTE/predict_isoform_productivity")
functional_annot = os.path.join(differential_expression, "diffExpr_DTE/functional_annotation")
mdt_dir = os.path.join(differential_expression, "most_dominant_transcript")
temp_fasta = os.path.join(methylation_dir, "temp")
if not os.path.exists(pipeline_reports): os.makedirs(pipeline_reports)



def quality_control(seq_summary_file, sample_id, raw_data_dir):
	# Producing preliminary QC reports 
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} QUALITY CONTROL OF THE INPUT SAMPLES')


	if not os.path.exists(initial_qc_reports): os.makedirs(initial_qc_reports)

	# Using NanoPlot (github.com/wdecoster/NanoPlot) to extract basic information regarding the sequencing
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  NanoPlot - Quality Control of the input {sample_id} raw data: in progress ..')
	nanoPlot = " ".join([
	"NanoPlot",  # Call NanoPlot
	"--threads", args.threads,  # Number of threads to be used by the script
	"--summary", seq_summary_file,  # Input of <sequencing_summary> generated by Albacore1.0.0 / Guppy 2.1.3+
	"--color saddlebrown",  # Color for the plots
	"--prefix", os.path.join(initial_qc_reports, f"{sample_id}_"),  # Create the report in the reports directory
	"--tsv_stats",  # Output the stats file as a properly formatted TSV
	"--format png",  # Output format of the plots
	"--dpi 900"]) # Set the dpi for saving images in high resolution
	subprocess.run(nanoPlot, shell=True)

	# Move log files to the 'pipeline_reports dir'
	subprocess.run(f"mv {initial_qc_reports}/*.log {pipeline_reports}", shell=True)
	return

def alignment_against_ref(fastq_pass, sample_id, raw_data_dir, seq_summary_file):
	""" Using Minimap2 to align the raw data against the reference genome """
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} ALIGNING AGAINSTE THE REF. GENOME')

	
	if not os.path.exists(alignments_dir): os.makedirs(alignments_dir)

	# Remove reads with length smaller than 47 nucleotides
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Filtering out reads with lenght less than 47nt: in progress ..')
	filtered_fastq = f'{alignments_dir}/{sample_id}.filt.fastq'
	# subprocess.run(f'seqtk seq -L 47 {fastq_pass} > {filtered_fastq}', shell=True)

	### ALIGN THE RAW READS AGAINST THE REFERENCE GENOME
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Minimap2 - Mapping {sample_id} against the reference genome: in progress ..')
	bam_file = f"{alignments_dir}/{sample_id}.genome.bam"
	minimap2_genome = " ".join([
	"minimap2",  # Call minimap2 (v2.17-r941)
	"-t", args.threads,  # Number of threads to use
	"-ax splice",   # Long-read spliced alignment mode and output in SAM format (-a)
	"-k 13",  # k-mer size
	"-uf",  # Find canonical splicing sites GT-AG - f: transcript strand
	"--secondary=no",  # Do not report any secondary alignments
	"--MD",  # output the MD tag
	refGenomeGRCh38,  # Inputting the reference genome
	filtered_fastq,  # Input .fastq.gz file
	"|", "samtools sort",  # Calling 'samtools sort' to sort the output alignment file
	"--threads", args.threads,  # Number of threads to be used by 'samtools sort'
	"--output-fmt BAM",  # Specify output format
	"-",  # Input from standard output
	"-o", bam_file,  # Sorted output  BAM file
	"2>>", os.path.join(pipeline_reports, "alignment_minimap2-report.txt")])  # Directory where all reports reside
	# subprocess.run(minimap2_genome, shell=True)
	# subprocess.run(f'samtools index -@ {args.threads} {bam_file}', shell=True)
	
	mapping_qc(sample_id, seq_summary_file, bam_file)
	return

def mapping_qc(sample_id, seq_summary_file, bam_file):
	""" Outputting multiple alignment statistics """
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} GENERATING ALIGNMENT STATS')
	

	if not os.path.exists(postAlignment_reports): os.makedirs(postAlignment_reports)
	
	# ### EXPORTING ALIGNMENT STATS
	# print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  PycoQC - Generating post-alignment stats of {sample_id}: in progress ..')
	# pycoQC = " ".join([
	# "pycoQC",  # Call pycoQC
	# "--quiet",  # Reduce verbosity
	# "--summary_file", seq_summary_file, 
	# "--bam_file", bam_file,  # Input of bam files from Minimap2
	# "--html_outfile", os.path.join(postAlignment_reports, f"{sample_id}.pycoQC-report.html"),  # Create the report in the reports directory
	# "--report_title", "\"Post-alignment quality control report\"",  # A title to be used in the html report
	# "2>>", os.path.join(pipeline_reports, "postalignment_pycoQC-report.txt")])
	# subprocess.run(pycoQC, shell=True)

	# # Picard CollectAlignmentSummaryMetrics
	# print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Picard - Collecting alignment summary metrics of {sample_id}: in progress ..')
	# CollectAlignmentSummaryMetrics = ' '.join([
	# "PicardCommandLine CollectAlignmentSummaryMetrics",  # Call picard-tools CollectAlignmentSummaryMetrics
	# "VERBOSITY= ERROR",  # Control verbosity of logging
	# "QUIET= true",  # Whether to suppress job-summary info on System.err
	# f"INPUT= {bam_file}",  # Input BAM file
	# f"OUTPUT= {postAlignment_reports}/{sample_id}.alignment_metrics.txt",  # Output
	# f"REFERENCE_SEQUENCE= {refGenomeGRCh38}",  # Reference sequence file
	# "2>>", os.path.join(pipeline_reports, "postalignment_collectAlignmentSummaryMetrics-report.txt")])
	# subprocess.run(CollectAlignmentSummaryMetrics, shell=True) 

	# # Check duplicate reads
	# print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Picard - Extracting read duplication stats of {sample_id}: in progress ..')
	# duplicate_reads = ' '.join([
	# "PicardCommandLine MarkDuplicates",  # Call MarkDuplicates
	# f"INPUT= {bam_file}",  # Input BAM file
	# f"OUTPUT= {alignments_dir}/{sample_id}.genome.dedup.bam",
	# f"METRICS_FILE= {postAlignment_reports}/{sample_id}.mark_duplicates.txt",  # Output file
	# "2>>", os.path.join(pipeline_reports, "postalignment_picard_markDuplicate-report.txt")])
	# subprocess.run(duplicate_reads, shell=True)
	# subprocess.run(f'rm {alignments_dir}/{sample_id}.genome.dedup.bam', shell=True)

	# # BAM stats
	# print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  BAM stats - Generating post-alignment stats of {sample_id}: in progress ..')
	# bam_stat = ' '.join([
	# "bam_stat.py",
	# "-i", bam_file,  # Input BAM file
	# f"> {postAlignment_reports}/{sample_id}.bamstat.txt",  # Output file
	# "2>>", os.path.join(pipeline_reports, "postalignment_bamstats-report.txt")])
	# subprocess.run(bam_stat, shell=True)

	# # BAM read distribution
	# print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  RSeQC - Generating read distribution stats of {sample_id}: in progress ..')
	# read_distribution = ' '.join([
	# "read_distribution.py", 
	# "-i", bam_file,  # Input BAM file
	# "-r", refAnnot_bed12,  # Reference in bed12
	# f"> {postAlignment_reports}/{sample_id}.fragSize",  # Output file
	# "2>>", os.path.join(pipeline_reports, "postalignment_read_distribution-report.txt")])
	# subprocess.run(read_distribution, shell=True)

	# # Check the strandness of the reads
	# print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  RSeQC - Generating read strandness stats of {sample_id}: in progress ..')
	# strandness = ' '.join([
	# "infer_experiment.py",  # Call samtools infer_experiment
	# "-i", bam_file,  # Input BAM file
	# "-r", refAnnot_bed,  # Reference gene model in bed format
	# f"> {postAlignment_reports}/{sample_id}.strandness.txt",  # Output file
	# "2>>", os.path.join(pipeline_reports, "postalignment_strandness-report.txt")])
	# subprocess.run(strandness, shell=True)

	# # Gene body coverage
	# print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  RSeQC - Generating gene body coverage of {sample_id}: in progress ..')
	# gene_coverage = ' '.join([
	# "geneBody_coverage.py",  # Call samtools geneBody_coverage
	# "-i", bam_file,  # Input BAM file
	# "-r", refAnnot_bed12,  # Reference gene model in bed format
	# "-f png",  # Output file format
	# "-o", f"{postAlignment_reports}/gene_coverage.{sample_id}",  # Output file
	# "2>>", os.path.join(pipeline_reports, "postalignment_gene_coverage-report.txt")])
	# subprocess.run(gene_coverage, shell=True)

	# # Number of reads mapped to each chromosome
	# print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  RSeQC - Generating mapping stats of {sample_id}: in progress ..')
	# mapping_pos = ' '.join([
	# "samtools idxstats",  # Call samtools idxstats
	# bam_file,  # Input BAM file
	# f"> {postAlignment_reports}/{sample_id}.samtools_idxstats.txt",  # Output file
	# "2>>", os.path.join(pipeline_reports, "postalignment_samtools_idxstats-report.txt")])
	# subprocess.run(mapping_pos, shell=True)

	# QualiMap QC
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  QualiMap rnaseq - Performing QC control on {sample_id}: in progress ..')
	qualimp = ' '.join([
	"qualimap rnaseq",  # Call qualimap rnaseq
	"--java-mem-size=15G",  # Set Java memory heap size
	"-bam", bam_file,  # Input BAM file
	"-gtf", refAnnot,  # Annotations file in Ensembl GTF format
	"-outdir", f'{coverage_reports}/{sample_id}/',  # Output folder for HTML report and raw data
	# "2>>", os.path.join(pipeline_reports, "postalignment_qualimap_rnaseq-report.txt")
	])
	subprocess.run(qualimp, shell=True)

	# QualiMap QC
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  QualiMap bamqc - Performing QC control on {sample_id}: in progress ..')
	qualimp = ' '.join([
	"qualimap bamqc",  # Call qualimap bamqc
	"-nt", args.threads,  #  Number of thread
	"--java-mem-size=15G",  # Set Java memory heap size
	"--genome-gc-distr hg19",  # Species to compare with genome GC distribution
	"-bam", bam_file,  # Input BAM file
	"--feature-file", refAnnot,  # Annotations file in Ensembl GTF format
	"-outdir", f'{coverage_reports}/{sample_id}/',  # Output folder for HTML report and raw data
	# "2>>", os.path.join(pipeline_reports, "postalignment_qualimap_bamqc-report.txt")
	])
	subprocess.run(qualimp, shell=True)

	# # Wub Alignment based QC plots
	# print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  WUB - Alignment based QC plots of {sample_id}: in progress ..')
	# alignment_qc = ' '.join([
	# "bam_alignment_qc.py",
	# "-f", refGenomeGRCh38,  # Input reference file
	# "-x -Q",  # Do not plot per-reference information/ Keep qiet
	# f"-r", "{postAlignment_reports}/{sample_id}.bam_alignment_qc.pdf",  # Output pdf file
	# f"-p", "{postAlignment_reports}/{sample_id}.bam_alignment_qc.pk",  # Output pk file
	# bam_file,
	# "2>>", os.path.join(pipeline_reports, "postalignment_wub-report.txt")])
	# subprocess.run(alignment_qc, shell=True)

	# # AlignQC 
	# print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  AlignQC - Generating post-alignment stats of {sample_id}: in progress ..')
	# alignqc = ' '.join([
	# "alignqc analyze",  # Calling AlignQC analyze
	# "--threads", str(int(int(args.threads)/2)),
	# "--genome", refGenomeGRCh38,  # Reference in .fasta
	# f"--gtf", "{refAnnot}.gz",  # Input reference annotation in gzipped form
	# f"--output", "{postAlignment_reports}/{sample_id}.alignqc.xhtml",  # Output pdf file
	# bam_file,
	# "2>>", os.path.join(pipeline_reports, "postalignment_alignqc-report.txt")])
	# subprocess.run(alignqc, shell=True)
	return

def polyA_estimation(sample_id, sum_file, fastq_pass, raw_data_dir):
	""" PolyA length estimation using Nanopolish polyA """
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} POLYA LENGTH ESTIMATION')


	polyA_analysis_dir_idv = os.path.join(polyA_analysis_dir, sample_id)
	if not os.path.exists(polyA_analysis_dir_idv): os.makedirs(polyA_analysis_dir_idv)
	filtered_fastq = f'{alignments_dir}/{sample_id}.filt.fastq'
	

	# Nanopolish index
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Nanopolish index - Indexing the output of the guppy basecaller of sample {sample_id}: in progress ..')
	indexing = ' '.join([
	"nanopolish index",  # Calling Nanopolish index
	"--sequencing-summary", sum_file,  # the sequencing summary file from Guppy
	"--directory", f"{raw_data_dir}/workspace/fast5_pass",  # Input BAM file
	filtered_fastq,
	"2>>", os.path.join(pipeline_reports, "nanopolish_index-report.txt")])
	subprocess.run(indexing, shell=True)

	# Nanopolish polyA
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Nanopolish polyA -  Estimate the polyadenylated tail lengths of {sample_id}: in progress ..')
	polyA_est = ' '.join([
	"nanopolish polya",  # Calling Nanopolish polyA
	"--threads", args.threads,  # Number of threads to use
	"--genome", refGenomeGRCh38,  # The reference genome assembly that  was used
	"--reads", filtered_fastq,  # The raw 1D ONT direct RNA reads in fastq
	"--bam", f"{alignments_dir}/{sample_id}.genome.bam",  # The reads aligned to the genome assembly in BAM format
	f"> {polyA_analysis_dir_idv}/{sample_id}.polya_results.tsv",  # Output file
	"2>>", os.path.join(pipeline_reports, "nanopolish_polyA-report.txt")])
	subprocess.run(polyA_est, shell=True)
	return

def methylation_detection_xpore(sample_id, fastq_pass, raw_data_dir):
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} METHYLATION DETECTION')
	
	methylation_dir_ind = f"{methylation_dir}/{sample_id}"
	if not os.path.exists(methylation_dir_ind): os.makedirs(methylation_dir_ind)


	# 1. Align nanopore events to reference k-mers
	eventalign = " ".join([
	"nanopolish eventalign",  # Calling nanopolish eventalign
	"--threads", args.threads,  # Use NUM threads
	"--scale-events",  # Scale events to the model, rather than vice-versa
	"--signal-index",  # write the raw signal start and end index values for the event to the tsv output
	"--min-mapping-quality 0",  # Only use reads with mapping quality at least 0
	"--genome", refGenomeGRCh38,  # The genome we are computing a consensus for is in FILE
	"--reads", f"{alignments_dir}/{sample_id}.filt.fastq",  # The 2D ONT reads are in fasta FILE
	"--bam", f"{alignments_dir}/{sample_id}.genome.bam",  # The reads aligned to the genome assembly are in bam
	"--summary", f"{methylation_dir_ind}/{sample_id}.eventalig-summary.txt",  # Summarize the alignment of each read/strand in FILE
	">", f"{methylation_dir_ind}/{sample_id}.eventalign.txt",  # Output file
	"2>>", os.path.join(pipeline_reports, "methylation_eventalign-report.txt")])
	# subprocess.run(eventalign, shell=True)


	# 2. Calling xpore-dataprep to prepare the data for analysis
	xpore_dataprep = " ".join([
	"xpore dataprep",  # Calling xpore-dataprep
	"--n_processes", args.threads,  # Number of processes to run
	# "--gtf_path_or_url", f"{expression_analysis_dir}/database_ALL_talon.gtf",  # GTF of the reference
	# "--transcript_fasta_paths_or_urls", f"{expression_analysis_dir}/reference_transcriptome.fasta",  # Reference transcriptome
	"--eventalign", f"{methylation_dir_ind}/{sample_id}.eventalign.txt",  # Eventalign filepath
	# "--summary", f"{methylation_dir_ind}/{sample_id}.eventalig-summary.txt",  # Eventalign summary filepath
	"--out_dir", methylation_dir_ind,  # Output directory
	# "--genome",  # To run on Genomic coordinates
	# "2>>", os.path.join(pipeline_reports, "methylation_dataprep-report.txt")
	])
	subprocess.run(xpore_dataprep, shell=True)

	# # 3. Output reference sequence around most significantly modified sites
	# tombo_sign = " ".join([
	# "tombo text_output signif_sequence_context",  # Indexing the concat_samples.bam file
	# "--fast5-basedirs", dirs,  # Directory containing fast5 files
	# "--statistics-filename", os.path.join(methylation_dir,"{0}.tombo.stats".format(sample_id)),
	# "--sequences-filename",  os.path.join(methylation_dir,"{0}.tombo_significant_regions.fasta".format(sample_id)),
	# # "2>>", os.path.join(pipeline_reports, "tombo_significants-report.txt")
	# ])
	# subprocess.run(tombo_sign, shell=True)
	return

def methylation_detection_eligos2():
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} METHYLATION DETECTION USING ELIGOS2')


	# methylation_dir_ind = f"{methylation_dir}/{sample_id}"
	if not os.path.exists(eligos_methylation_dir): os.makedirs(eligos_methylation_dir)
	test_samples =  ' '.join(sorted(glob.glob(f'{alignments_dir}/Tumour_*.genome.bam')))
	control_samples = ' '.join(sorted(glob.glob(f'{alignments_dir}/NonTransf_*.genome.bam')))
	ctrl_bams=['/home/stavros/playground/nanodRNA_analysis/batch3_analysis_v3.6.1/alignments/NonTransf_1.genome.bam', '/home/stavros/playground/nanodRNA_analysis/batch3_analysis_v3.6.1/alignments/NonTransf_2.genome.bam', '/home/stavros/playground/nanodRNA_analysis/batch3_analysis_v3.6.1/alignments/NonTransf_3.genome.bam']
	test_bams=['/home/stavros/playground/nanodRNA_analysis/batch3_analysis_v3.6.1/alignments/Tumour_1.genome.bam', '/home/stavros/playground/nanodRNA_analysis/batch3_analysis_v3.6.1/alignments/Tumour_2.genome.bam', '/home/stavros/playground/nanodRNA_analysis/batch3_analysis_v3.6.1/alignments/Tumour_4.genome.bam']

	eligos_pair_diff_mod = " ".join([
	"eligos2 pair_diff_mod",  # Calling eligos2 pair_diff_mod
	"--threads", args.threads,  # Use NUM threads
	"--oddR 5",  # Minimum cut-off for Odd ratio
	"--esb 0.2",  # Minimum cut-off for ratio of error at specific base (ESB)
	"--pval 0.05",  # P-value cut-off
	"--test_bams", '/home/stavros/playground/nanodRNA_analysis/batch3_analysis_v3.6.1/alignments/Tumour_4.genome.bam',  # Test RNA/DNA data in one or more sorted bam files
	"--ctrl_bams", '/home/stavros/playground/nanodRNA_analysis/batch3_analysis_v3.6.1/alignments/NonTransf_3.genome.bam',  # Control RNA/DNA data in one or more sorted bam files
	"--reference", refGenomeGRCh38,  # The genome we are computing a consensus for is in FILE
	"--outdir", eligos_methylation_dir,  # Path of directory name to store output
	"--region", refAnnot_bed12,  # Genes in bed format
	"2>>", os.path.join(pipeline_reports, "methylation_pair_diff_mod-report.T4NT3.txt")])
	subprocess.run(eligos_pair_diff_mod, shell=True)


	eligos_filter = " ".join([
	"eligos2 filter",  # Calling eligos2 filter
	"--adjPval 0.05",  # Adjusted p-value cut-off
	"--homopolymer",  # Filtering out homopolymer
	"--select_base A",  # Select output Base for BedGraph
	"--input", '/home/stavros/playground/nanodRNA_analysis/batch3_analysis_v3.6.1/eligos2_meth_analysis_m6a/Tumour_4.genome_vs_NonTransf_3.genome_on_gencode.v35.primary_assembly.annotation_baseExt0.txt',  # Eligos result in txt file
	"--prefix", '/home/stavros/playground/nanodRNA_analysis/batch3_analysis_v3.6.1/eligos2_meth_analysis_m6a/Tumour_4.genome_vs_NonTransf_3.genome_on_gencode.v35.primary_assembly.annotation_baseExt0',  # Set the output file prefix
	# "2>>", os.path.join(pipeline_reports, "methylation_filter-report.txt")
	])
	subprocess.run(eligos_filter, shell=True)


	output_files = glob.glob(f'{eligos_methylation_dir}/*_combine.txt')
	for file in output_files:
		filename = os.path.basename(file)[:-4]
		eligos_extension = " ".join([
		"table2fa_eligos.mergeextend.sh",  # Calling eligos2 table2fa_eligos.mergeextend.sh
		"--adjPval 0.01",  # Adjusted p-value cut-off
		"--homopolymer",  # Filtering out homopolymer
		"--input", file,  # Eligos result in txt file
		"--prefix", f'{eligos_methylation_dir}/{filename}.fasta',  # Set the output file prefix
		# "2>>", os.path.join(pipeline_reports, "methylation_filter-report.txt")
		])
		subprocess.run(eligos_extension, shell=True)

	return

class expression_analysis:

	def __init__(self):
		self.talon_analysis()
		self.talon_overall_visualisation()
		return

	def talon_analysis(self):
		""" TALON takes transcripts from one or more long read datasets (SAM format) 
		and assigns them transcript and gene identifiers based on a database-bound
		annotation. Novel events are assigned new identifiers """
		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} ANNOTATION AND QUANTIFICATION USING TALON')
		

		talon_database = os.path.join(expression_analysis_dir, 'talon_gencode_v35.db')  ### TALON DB
		if not os.path.exists(differential_expression): os.makedirs(differential_expression)
		temp = os.path.join(expression_analysis_dir, 'transcriptClean_temp')
		if not os.path.exists(temp): os.makedirs(temp)
		filtered_isoforms_final = os.path.join(expression_analysis_dir, "filtered_isoforms_final.csv")
		reference_fasta = os.path.join(expression_analysis_dir, "reference_transcriptome.fasta")

		sample_group = []
		# Create config file that is needed for talon and converting the aligned bam files to sam
		csv_file = os.path.join(expression_analysis_dir,'talon_input.csv')
		if not os.path.exists(csv_file) or os.stat(csv_file).st_size == 0:
			print("{0}  Creating a description file necessary for Talon: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
			with open(csv_file, "w") as talon_out:
				for path, subdir, folder in os.walk(alignments_dir):
					for file in sorted(folder):
						if file.endswith('.genome.bam'):
							sample_group.append(file.split(".")[0].split("_")[0])
							# Output a csv file with 'sample_id, sample_group, technology, input_labeled_sam_file' which is gonna be needed in talon_annotation function
							talon_out.write('{0},{1},ONT,{2}\n'.format(file.split(".")[0], file.split(".")[0].split("_")[0], os.path.join(temp, file).replace(".genome.bam", ".clean_labeled.sam")))
							if not os.path.exists(os.path.join(path, file).replace(".bam", ".sam")):
								print("{0}  Samtools - Converting the genomic bam file ({1}) to sam for Talon: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), file))
								os.system('samtools view -h -@ {0} {1} > {2}'.format(args.threads, os.path.join(path, file), os.path.join(path, file).replace(".bam", ".sam")))


		### Initial step: Correcting mismatches, microindels, and noncanonical splice junctions in long reads that have been mapped to the genome
		sam_files = glob.glob(os.path.join(alignments_dir, "*.genome.sam"))
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} 1/10 TranscriptClean - Correcting mismatches, microindels, and noncanonical splice junctions in long reads: in progress ..')
		for file in sam_files:
			TranscriptClean = " ".join([
			transcriptclean,  # Call TranscriptClean.py
			"--sam", file,  # Input sam file
			"--genome", refGenomeGRCh38_traclean,  # Reference genome fasta file
			"--threads", args.threads,  # Number of threads to be used by the script
			"--canonOnly",  # Output only canonical transcripts and transcripts containing annotated noncanonical junctions to the clean SAM file
			"--spliceJns", sj_ref,  # Splice junction file extracted from the ref. GTF
			"--deleteTmp",  # the temporary directory (TC_tmp) will be removed
			"--outprefix", os.path.join(temp, os.path.basename(file).split(".")[0]),  # Outprefix for the outout file
			"2>>", os.path.join(pipeline_reports, "talon1_transcriptclean-report.txt")])  # Directory where all reports reside
			subprocess.run(TranscriptClean, shell=True)
			os.remove(file)
		
		### Second step: Building the reference database
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} 2/10 TALON - Initiating the database: in progress ..')
		if os.path.exists(talon_database): os.remove(talon_database)
		original = sys.stdout
		sys.stdout = open(f'{pipeline_reports}/talon2_initialize_database-report.txt',"w")
		talon_initialize_database = " ".join([
		"talon_initialize_database",  # Call talon_initialize_database
		"--f", refAnnot,  # GTF annotation containing genes, transcripts, and edges
		"--g", 'hg38',  # Genome build (i.e. hg38) to use
		"--5p 500",  # Maximum allowable distance (bp) at the 5' end during annotation
		"--3p 300",  # Maximum allowable distance (bp) at the 3' end during annotation
		"--a", talon_database[:-3],  # Name of supplied annotation
		"--o", talon_database[:-3],  # Outprefix for the annotation files
		"2>>", os.path.join(pipeline_reports, "talon2_initialize_database-report.txt")])  # Directory where all reports reside
		subprocess.run(talon_initialize_database, shell=True)
		sys.stdout = original

		### Third step: internal priming check
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  3/10 TALON - Run talon_label_reads on each file to compute how likely each read is to be an internal priming product: in progress ..')
		genome_alignments = [sum_file for sum_file in glob.glob(os.path.join(temp, "*clean.sam"))]
		original = sys.stdout
		sys.stdout = open(f'{pipeline_reports}/talon3_talon_priming_check-report.txt',"w")
		for file in genome_alignments:
			talon_priming_check = " ".join([
			"talon_label_reads",  # Call talon talon_label_reads
			"--t", args.threads,  # Number of threads to be used by the script
			"--f", file,  # Input sam file
			"--g", refGenomeGRCh38,  # Reference genome fasta file
			"--tmpDir", os.path.join(expression_analysis_dir, "tmp_label_reads"),  # Path to directory for tmp files
			"--deleteTmp",  # Temporary directory will be removed
			"--o", os.path.join(temp, os.path.basename(file).replace("_clean.sam",".clean")),  # Prefix for output files
			"2>>", os.path.join(pipeline_reports, "talon3_talon_priming_check-report.txt")])  # Directory where all reports reside
			subprocess.run(talon_priming_check, shell=True)
		sys.stdout = original

		### Fourth step: annotating and quantification of the reads
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  4/10 TALON - Annotating and quantification of the reads: in progress ..')
		original = sys.stdout
		sys.stdout = open(f'{pipeline_reports}/talon4_annotationNquantification-report.txt',"w")
		talon_annotation = " ".join([
		"talon",  # Call talon
		"--threads", args.threads,  # Number of threads to be used by the script
		"--db", talon_database,  # TALON database
		"--build", 'hg38',  # Genome build (i.e. hg38) to use
		"--o", os.path.join(expression_analysis_dir, "prefilt"),  # Prefix for output files
		"--f", csv_file,  # Dataset config file: dataset name, sample description, platform, sam file (comma-delimited)
		"2>>", os.path.join(pipeline_reports, "talon4_annotationNquantification-report.txt")])  # Directory where all reports reside
		subprocess.run(talon_annotation, shell=True)
		sys.stdout = original

		### Fifth step: summarising how many of each transcript were found (prior to any filtering)
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  5/10 TALON - Summarising how many of each transcript were found (prior to any filtering): in progress ..')
		talon_summary = " ".join([
		"talon_summarize",  # Call talon_summarize
		"--db", talon_database,  # TALON database
		"--o", os.path.join(expression_analysis_dir, "prefilt"),  # Prefix for output file
		"2>>", os.path.join(pipeline_reports, "talon5_summary-report.txt")])  # Directory where all reports reside
		subprocess.run(talon_summary, shell=True)
		
		### Sixth step: creating an abundance matrix without filtering (for use computing gene expression)
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  6/10 TALON - Creating an abundance matrix without filtering (gene expression): in progress ..')
		talon_abundance = " ".join([
		"talon_abundance",  # Call talon_abundance
		"--db", talon_database,  # TALON database
		"--build", 'hg38',  # Genome build (i.e. hg38) to use
		"--annot", talon_database[:-3],  # Which annotation version to use
		"--o", os.path.join(expression_analysis_dir, "prefilt"),  # Prefix for output file
		"2>>", os.path.join(pipeline_reports, "talon6_abundance-report.txt")])  # Directory where all reports reside
		subprocess.run(talon_abundance, shell=True)
		
		### Seventh step: Applying basic filtering steps and outputting several stats
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  7/10 TALON - Removing internal priming artifacts low abundance transcripts: in progress ..')
		talon_filter = " ".join([
		"talon_filter_transcripts",  # Call talon_filter_transcripts
		"--db", talon_database,  # TALON database
		"--annot", talon_database[:-3],  # Which annotation version to use
		"--minCount 5",  # Number of minimum occurrences required for a novel transcript PER dataset
		"--maxFracA 0.5",  # All of the supporting reads must have 50% or fewer As in the 20 bp interval after alignment
		"--minDatasets", str(3),  # Minimum number of datasets novel transcripts must be found in
		# "--minDatasets", str(min(Counter(sample_group).values())),  # Minimum number of datasets novel transcripts must be found in
		"--o", os.path.join(expression_analysis_dir, "filtered_isoforms.csv"),  # Output
		"2>>", os.path.join(pipeline_reports, "talon7_talon_filter-report.txt")])  # Directory where all reports reside
		subprocess.run(talon_filter, shell=True)
		
		### Eighth step:
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  8/10 TALON - Applying custom filtering of the ISM group based on the polyA estimation: in progress ..')
		prefit_abundance_matrix = os.path.join(expression_analysis_dir, "prefilt_talon_abundance.tsv")
		prefilt_readannot_matrix = os.path.join(expression_analysis_dir, "prefilt_talon_read_annot.tsv")
		filtered_isoforms = os.path.join(expression_analysis_dir, "filtered_isoforms.csv")
		self.polyA_filtering(prefit_abundance_matrix, prefilt_readannot_matrix, filtered_isoforms, filtered_isoforms_final)
		
		### Ninth step: Applying basic filtering steps and outputting several stats
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  9/10 TALON - Removing filtered isoforms based on step 8 and exporting basic statsistics (ISM filters included): in progress ..')
		talon_filter_n_report = " ".join([
		"Rscript",  # Call Rscript
		f"{rscripts}/talon_summarisation.R",  # Calling the talon_summarisation.R script
		os.path.join(expression_analysis_dir, "prefilt_talon_abundance.tsv"),  # Input matrix
		differential_expression,  # Output directory
		os.path.join(expression_analysis_dir, "talon_input.csv"),  # Input annotation matrix
		os.path.join(expression_analysis_dir, "filtered_isoforms_final.csv"),  # Filtered transcripts to maintain, including ISM filtering
		os.path.join(expression_analysis_dir, "filtered_isoforms.csv"),  # Filtered transcripts to maintain, original
		args.ismTheshold,  # # Min. count threshold for the ISM group (either one of the groups)
		"2>>", os.path.join(pipeline_reports, "talon9_summarisation.txt")])  # Directory where all reports reside
		subprocess.run(talon_filter_n_report, shell=True)
		
		### Tenth step: Generating TALON report for each dataset
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  10/10 TALON - Generating TALON report for each dataset (talon filters only): in progress ..')
		for file in glob.glob(os.path.join(alignments_dir, "*.genome.bam")):
			sample_name = os.path.basename(file).split(".")[0]
			output_dir = os.path.join(differential_expression, f"talon_reports/{sample_name}_report")
			if not os.path.exists(output_dir): os.makedirs(output_dir)
			talon_report = " ".join([
			"talon_generate_report",  # Call talon_abundance
			"--db", talon_database,  # TALON database
			"--whitelists", os.path.join(expression_analysis_dir, "filtered_isoforms.csv"),  # Filtered transcripts to be reported
			"--datasets", sample_name,  # Input of the filtered tables produced on the previous step
			"--outdir", output_dir,  # Output dir
			"2>>", os.path.join(pipeline_reports, "talon10_generate_report-report.txt")])  # Directory where all reports reside
			subprocess.run(talon_report, shell=True)

		self.extract_talon_database(talon_database, reference_fasta)  # Extract Talon database in gtf and fasta format
		
		# Examining the antisense isoforms
		self.examine_antisense(talon_database)
		
		### Removing unnecessary directories and files
		subprocess.run(f"rm -r {temp}", shell=True)  
		subprocess.run("rm -r talon_tmp", shell=True)
		return

	def polyA_filtering(self, prefilt_talon_abundance, talon_read_annot, filtered_isoforms, output_file):
		""" Here is where all the magic of the special filtering is taking place. 
		Here we are using the information of the polyA length estimation to filter 
		out novel transcripts in the ISM category """
		
		not_to_keep = {}
		filtered_isolist = {}
		with open(filtered_isoforms) as filtin:
			for line in filtin:
				filtered_isolist[line.strip().split(",")[1]] = line.strip().split(",")[0]

		prefiltered_transcripts = {}
		with open(prefilt_talon_abundance) as filt_tr:
			for line in filt_tr:
				if not line.startswith("gene_ID"):
					if line.strip().split("\t")[1] in filtered_isolist:
						prefiltered_transcripts[line.strip().split("\t")[3]] = line.strip().split("\t")[1]

		read_annot_dict = {}
		# TALON read annotation file and saving it as a dictionary 
		with open(talon_read_annot) as read_annot:
			for line in read_annot:
				if line.strip().split("\t")[12] in prefiltered_transcripts:
					read = line.strip().split("\t")[0]
					sample = line.strip().split("\t")[1]
					transcript_id = line.strip().split("\t")[12]
					read_type = "{0}_{1}_{2}".format(prefiltered_transcripts[transcript_id],line.strip().split("\t")[16],line.strip().split("\t")[17])
					if "ISM" in read_type:  # Removing the ISM Suffix group from the downstream analysis. This will be analysed differently
						read_annot_dict[(sample, read)] = (transcript_id, read_type)


		### Iterating through all "polya_results.tsv" obtained from the Nanopolish polyA function 
		### in order to filter out transcripts in the ISM Prefix group without enough polyA evidence
		filter_dict = {}
		for path, subdir, folder in os.walk(analysis_dir):
			for name in folder:
				if name.endswith("polya_results.tsv"):
					sample = name.split(".")[0]
					polyA_length_est = os.path.join(path, name)
					with open(polyA_length_est) as fin:
						for line in fin:
							if not line.startswith("readname"):
								readname = line.strip().split("\t")[0]
								if (sample ,readname) in read_annot_dict:
									transcript_name = read_annot_dict[(sample ,readname)][0] 
									read_type = read_annot_dict[(sample ,readname)][1]
									# if "ISM_Suffix" in read_type:print(read_type)
									qc = ["OTHER","PASS"][line.strip().split("\t")[-1]=="PASS"]
									transcript_ID = read_type.split("_")[0]
									if "ISM" in read_type:
										if (transcript_name, transcript_ID) in filter_dict:
											if qc ==  "PASS":
												filter_dict[(transcript_name, transcript_ID)][0] += 1
											else:
												filter_dict[(transcript_name, transcript_ID)][1] += 1
										else:
											if qc ==  "PASS":
												filter_dict[(transcript_name, transcript_ID)] = [1,0]
											else:
												filter_dict[(transcript_name, transcript_ID)] = [0,1]
		###621266_ISM_Suffix
		# for a,b in filter_dict.items():
		# 	print(a,b)

		for transcript, qc_tags in filter_dict.items():	# Obtaining the list with the transcripts where the PASS is
			if qc_tags[0] < qc_tags[1]:						# less than the rest
				not_to_keep[transcript[1]] = None

		# If not_to_keep list is empty, raise a warning..
		if len(not_to_keep) == 0: print("WARNING: all ISM transcripts passed the polyA QC filter!")

		with open(output_file, "w") as fout:			# Rewriting the "filtered_isoforms.csv" excluding the
			for trID, gnID in filtered_isolist.items():	# filtered transcripts of the ISM group without enough
				if not trID in not_to_keep:			    # polyA evidence.
					fout.write(f"{gnID},{trID}\n")
		return

	def extract_talon_database(self, talon_database, reference_fasta):

		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Obtaining the transcriptome annotation from the TALON database: in progress ..')
		talon_export_db = " ".join([
		"talon_create_GTF",  # Call talon_create_GTF
		"--db", talon_database,  # TALON database
		"--build hg38",  # Genome build (hg38) to use
		"--annot", talon_database[:-3],  # Which annotation version to use
		"--o", f"{expression_analysis_dir}/database",  # Output
		"--whitelist", f"{expression_analysis_dir}/final_filtered_isoforms_for_db.csv",  # Whitelist file of transcripts to include in the output
		"2>>", os.path.join(pipeline_reports, "talonextract_exportdb-report.txt")])  # Directory where all reports reside
		subprocess.run(talon_export_db, shell=True)

		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Converting TALON database to fasta format: in progress ..')
		convert_db = " ".join([
		"gffread",  # Call gffread
		"--sort-alpha",  # chromosomes (reference sequences) are sorted alphabetically
		"-g", refGenomeGRCh38,  # Fasta file with the genomic sequences for all input mappings
		"-w", reference_fasta,  # Write a fasta file with spliced exons for each transcript
		"-o", f"{expression_analysis_dir}/database_talon.fasta",
		f"{expression_analysis_dir}/database_ALL_talon.gtf",  # Input gtf
		"2>>", os.path.join(pipeline_reports, "talonextract_convertdb-report.txt")])  # Directory where all reports reside
		subprocess.run(convert_db, shell=True)
		return

	def examine_antisense(self, talon_database):
		""" In principal we know that is an antisense isoform is expressed,
		the matching sense isoform is silent. Here we examine if this principal 
		is being supported or not """
		if not os.path.exists(antisense_dir): os.makedirs(antisense_dir)
		output_matrix = f'{antisense_dir}/antisense-sense_abundance.csv'

		# Calling TALON map_antisense_genes_to_sense.py to find the matching sense to the antisense isoforms
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  TALON - Examining the antisense genes: in progress ..')
		run_map_antisense = ' '.join([
		map_antisense,  # Calling TransDecoder.LongOrfs
		"--db", talon_database,  #  TALON database
		"--annot", talon_database[:-3],  # Which annotation version to use
		"--o", f'{antisense_dir}/list_of',  # Path to the intended output directory
		"2>>", os.path.join(pipeline_reports, "antisense-report.txt")])  # Directory where all reports reside
		subprocess.run(run_map_antisense, shell=True)

		antisense_isoforms = {}
		sense_isoforms = {}
		matching_ids = {}
		# Retrieving the isoform expression of the antisense and matching sense isoforms 
		with open(f'{antisense_dir}/list_of_antisense_mapping.csv') as fin:
			for line in fin:
				if not line.startswith("antisense"):
					antisense_isoforms[line.strip().split(",")[0]] = None
					sense_isoforms[line.strip().split(",")[1]] = None
					matching_ids[line.strip().split(",")[0]] = line.strip().split(",")[1]

		header = None
		with open(f'{expression_analysis_dir}/prefilt_talon_abundance.tsv') as refin:
			for line in refin:
				if line.startswith("gene_ID"):
					header = line
				else:
					gene_id = line.strip().split("\t")[0]
					if gene_id in antisense_isoforms:
						if sum([int(i) for i in line.strip().split("\t")[11:]]) > 10:
							antisense_isoforms[gene_id] = line
						else:
							antisense_isoforms.pop(gene_id)
					elif gene_id in sense_isoforms:
						if sum([int(i) for i in line.strip().split("\t")[11:]]) > 10:
							sense_isoforms[gene_id] = line

		# Outputting the final matrix with the antisense (expressed with more than 10 reads)
		# and their matching sense genes
		with open(output_matrix, 'w') as matout:
			matout.write(header)
			for antisense, sense in matching_ids.items():
				if antisense in antisense_isoforms and sense in sense_isoforms:
					if sense_isoforms[sense] != None:
						matout.write(antisense_isoforms[antisense])
						matout.write(sense_isoforms[sense])
						matout.write("\n")
					else:
						matout.write(antisense_isoforms[antisense])
						matout.write("NONE\n")
						matout.write("\n")
		return

	def talon_overall_visualisation(self):
		
		annot = {}
		with open(os.path.join(expression_analysis_dir, "database_talon.gtf")) as ref_in:
			for i, line in enumerate(ref_in):
				if not line.startswith("#"):
					if line.split("\t")[2].strip() == 'transcript':
						transcript_id = line.split("\t")[-1].split(";")[1].split(" ")[-1].strip("\"")
						transcript_type = ["novel", line.split("\t")[-1].split(";")[4].split()[1].strip("\"").strip()][transcript_id.startswith("ENST")]
						if transcript_type.endswith("pseudogene"):
							transcript_type = "pseudogene"
						elif transcript_type in ["Mt_rRNA","miRNA","snoRNA","misc_RNA","snRNA","scaRNA","rRNA","Mt_tRNA"]:
							transcript_type = "ncRNA"
						elif transcript_type == "ribozyme" or transcript_type.startswith("IG_") or transcript_type.startswith("TR_"):
							transcript_type = "protein_coding"
						elif transcript_type == "TEC":
							transcript_type = "To_be_Confirmed"
						annot[transcript_id] = transcript_type

						
		data = {}
		header = []
		with open(f"{expression_analysis_dir}/filt_talon_abundance.csv") as mat_in:
			for i, line in enumerate(mat_in, 1):
				if i == 1:
					header = line.strip().split(",")[9:]
				else:
					transcript = line.strip().split(",")[1]
					values = line.strip().split(",")[9:]
					data[(transcript, annot[transcript])] = values

		# Remove paths from sample names
		header.insert(0,"transcript_id")  # Inserting ttranscript_id in header
		header.insert(1, "transcript_type")  # Inserting transcript_type in header
		
		# Writing output to file 'expression_matrix.csv'
		with open(f"{expression_analysis_dir}/perTranscript_expression_matrix.csv", "w") as fout:
			fout.write("{0}\n".format(','.join(header)))
			for key, values in data.items():
				fout.write("{0},{1}\n".format(','.join(key), ','.join(values)))

		print("{0}  Visualising the RNA categories found in the TALON expression matrix: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		gene_type_sum = " ".join([
		"Rscript",  # Call Rscript
		f"{rscripts}/transcript_type_summary.R",  # Calling the transcript_type_summary.R script
		f"{expression_analysis_dir}/perTranscript_expression_matrix.csv",  # Input matrix
		differential_expression,  # Output dir
		"2>>", os.path.join(pipeline_reports, "R_transcriptType_sum-report.txt")])  # Directory where all reports reside
		subprocess.run(gene_type_sum, shell=True)
		return

class downstream_analysis:

	def __init__(self):
		self.polyA_preprocessing()
		self.differential_polyadelylation_analysis()
		self.differential_expression_analysis()
		self.most_dominant_transcript()
		self.combine_dte_dpa_results()
		self.fusions_events_jaffa()
		return

	def polyA_preprocessing(self):
		""" In this method we are obtaining the output of the Nanopolish polyA function and
		input the transcript annotation nomencalture in the files as well as additional clinical 
		information which will be later be used for the DPA. We are also keeping the transcripts
		that pass the DTE filters. """
		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} POLYA ESTIMATIONS PREPROCESSING')
		


		filtered_transcripts = {}
		filt_talon_abundance = os.path.join(expression_analysis_dir,"filt_talon_abundance.csv")
		with open(filt_talon_abundance) as filt_tr:
			for line in filt_tr:
				if not line.startswith("annot_gene_id"):
					filtered_transcripts[line.strip().split(",")[1]] = None


		read_annot_dict = {}
		# TALON read annotation file and saving it as a dictionary 
		talon_read_annot = os.path.join(expression_analysis_dir,"prefilt_talon_read_annot.tsv")
		with open(talon_read_annot) as read_annot:
			for line in read_annot:
				if line.strip().split("\t")[12] in filtered_transcripts:
					read = line.strip().split("\t")[0]
					sample = line.strip().split("\t")[1]
					transcript_id = line.strip().split("\t")[12]
					gene_id = line.strip().split("\t")[11]
					read_type = "{0}_{1}".format(line.strip().split("\t")[16], line.strip().split("\t")[17])
					read_annot_dict[(sample, read)] = (transcript_id, gene_id, read_type)

		for path, subdir, folder in os.walk(polyA_analysis_dir):
			for name in folder:
				if name.endswith("polya_results.tsv"):
					sample = name.split(".")[0]
					print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Processing the polyadenylated tail lengths of {sample}: in progress ..')
					polyA_length_est = os.path.join(path, name)
					with open(polyA_length_est) as fin, open(polyA_length_est.replace(".tsv",".transcripts.pass.tsv"), "w") as transcript_out:
						for i, line in enumerate(fin, 1):
							if i == 1:
								transcript_out.write("{0}\tsample\tgroup\n".format(line.strip()))
							else:
								readname = line.strip().split("\t")[0]
								if (sample ,readname) in read_annot_dict:  # and line.strip().split("\t")[9] == "PASS":
									transcript_id = read_annot_dict[(sample ,readname)][0] 
									read_type = read_annot_dict[(sample ,readname)][2]
									rest = "\t".join(line.strip().split("\t")[3:])
									group = annotation[sample.split("_")[0]]
									transcript_out.write(f"{readname}\t{transcript_id}\t{read_type}\t{rest}\t{sample}\t{group}\n")	
		return

	def differential_polyadelylation_analysis(self):
		""" Performing differential polyadelylation analysis (DPA) 
		using an edited version of the the polya_diff.py script 
	    from the pipeline-polya-diff repository of ONT. Furthermore
	    we are using custom scripts to visualise the results and 
	    combine with the DTE. """
		
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  PolyA QC - Preprocessing and performing QC in the Nanopolish polyA results: in progress ..')
		polyAqc = " ".join([
		"Rscript",  # Call Rscript
		f"{rscripts}/polyA_prelim_analysis.R",  # Calling the polyA_prelim_analysis.R script
		polyA_analysis_dir,  # Directory of nanopolish output
		args.threads,  # Number of cores to be used
		"2>>", os.path.join(pipeline_reports, "dpa_polyAqc-report.txt")])  # Directory where all reports reside
		subprocess.run(polyAqc, shell=True)


		### Running Nanotail analysis
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Nanopore polyAdiff - Differential polyadenilation analysis: in progress ..')
		polyAdiff = " ".join([
		f"python3 {rscripts}/polya_diff.py",  # Calling the polya_diff.py script
		"-i", f"{polyA_analysis_dir}/dpa_results/all_tails.tsv",  # Input file
		"-g", f"{polyA_analysis_dir}/dpa_results/polya_diff_global.tsv",
		"-t", f"{polyA_analysis_dir}/dpa_results/polya_diff_per_transcript.tsv",
		"-r"  f"{polyA_analysis_dir}/dpa_results/polya_diff_report.pdf",
		"-c", args.minPolyA,  # Min. coverage
		# "-x",  # Plot per-transcript distributions
		"2>>", os.path.join(pipeline_reports, "dpa_polyAdiff-report.txt")])  # Directory where all reports reside
		subprocess.run(polyAdiff, shell=True)
		return

	def differential_expression_analysis(self):
		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} DIFFERENTIAL EXPRESSION ANALYSIS')
		
		
		### First step: Exploratory analysis
		print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")}  1/4 Differential Expression - Exploratory analysis: in progress ..')
		expl_analysis = " ".join([
		"Rscript",  # Call Rscript
		f"{rscripts}/diffExpr_ExplAnalysis.R",  # Calling the diffExpr_ExplAnalysis.R script
		os.path.join(expression_analysis_dir, "prefilt_talon_abundance.tsv"),  # Input filtered matrix
		os.path.join(expression_analysis_dir, "talon_input.csv"),  # Input annotation matrix
		differential_expression,  # Output directory
		args.n_top,  # Top n_top genes for creating the heatmap
		"2>>", os.path.join(pipeline_reports, "diffExpr_exploratory_analysis-report.txt")])  # Directory where all reports reside
		subprocess.run(expl_analysis, shell=True)
		

		### Second step: DGE
		print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")}  2/4 Differential Expression - Differential Gene Expression (DGE) analysis using edgeR: in progress ..')
		dge_analysis = " ".join([
		"Rscript",  # Call Rscript
		f"{rscripts}/diffExpr_DGE.R",  # Calling the diffExpr_DGE.R script
		f"{expression_analysis_dir}/prefilt_talon_abundance.tsv",  # Input filtered matrix
		f"{expression_analysis_dir}/talon_input.csv",  # Input annotation matrix
		differential_expression,  # Output directory
		args.adjPValueThreshold,  # adjPValueThreshold - Adjusted p-value threshold for differential expression
		args.lfcThreshold,  # lfcThreshold - Minimum required log2 fold change for differential expression
		"2>>", os.path.join(pipeline_reports, "diffExpr_dge_analysis-report.txt")])  # Directory where all reports reside
		subprocess.run(dge_analysis, shell=True)
		

		### Third step: DTE
		print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")}  3/4 Differential Expression - Differential Transcript Expression (DTE) analysis using edgeR: in progress ..')
		dte_analysis = " ".join([
		"Rscript",  # Call Rscript
		f"{rscripts}/diffExpr_DTE.R",  # Calling the diffExpr_DTE.R script
		f"{expression_analysis_dir}/filt_talon_abundance.csv",  # Input filtered matrix
		f"{expression_analysis_dir}/talon_input.csv",  # Input annotation matrix
		differential_expression,  # Output directory
		args.adjPValueThreshold,  # adjPValueThreshold - Adjusted p-value threshold for differential expression
		args.lfcThreshold,  # lfcThreshold - Minimum required log2 fold change for differential expression		
		"2>>", os.path.join(pipeline_reports, "diffExpr_dte_analysis-report.txt")])  # Directory where all reports reside
		subprocess.run(dte_analysis, shell=True)
		self.predict_productivity()


		### Fourth step: DTU
		print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")}  4/4 Differential Expression - Differential Transcript Usage (DTU) analysis using IsoformSwitchAnalyzeR: in progress ..')
		dtu_analysis = " ".join([
		"Rscript",  # Call Rscript
		f"{rscripts}/diffExpr_DTU.R",  # Calling the diffExpr_DTU.R script
		f"{expression_analysis_dir}/filt_talon_abundance.csv",  # Input filtered matrix
		f"{expression_analysis_dir}/talon_input.csv",  # Input annotation matrix
		f"{expression_analysis_dir}/reference_transcriptome.fasta",  # Fasta file with spliced exons for each transcript
		f"{expression_analysis_dir}/database_talon.gtf",  # Transcriptome annotation from the TALON database
		differential_expression,  # Output directory
		f"{rscripts}/multi_iupred2a.py",  # Calling the multi_iupred2a.py script
		pfam_dir,  # Pfam database
		args.threads,
		"2>>", os.path.join(pipeline_reports, "diffExpr_dtu_analysis-report.txt")])  # Directory where all reports reside
		subprocess.run(dtu_analysis, shell=True)
		return 

	def most_dominant_transcript(self):
		""" Using the filtered matrix to perform the most dominant transcript analysis (MDT) """
		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} MOST DOMINANT TRANSCRIPT ANALYSIS')


		if not os.path.exists(mdt_dir): os.makedirs(mdt_dir)
		filt_talon_abundance = os.path.join(expression_analysis_dir,"filt_talon_abundance.csv")
		csv_file = os.path.join(expression_analysis_dir,'talon_input.csv')
		

		gr1 = None
		gr2 = None
		groups = {}
		# Obtaining the matching group-sample info
		with open(csv_file) as fin:
		    for line in fin:
		        group = line.strip().split(",")[1]
		        sample = line.strip().split(",")[0]
		        if group in groups:
		            groups[group].append(sample)
		        else:
		            groups[group] = [sample]

		gr1 = sorted(groups.keys())[0]
		gr2 = sorted(groups.keys())[1]

		# Reading the filtered file into a dataframe
		mat = pd.read_csv(filt_talon_abundance, index_col=[2,3], sep=',')
		# print(mat.loc[mat['annot_gene_id'] == 'ENSG00000235098.8'])
		# mat = mat.drop(["annot_gene_name","annot_transcript_name","n_exons","length","gene_novelty","transcript_novelty","ISM_subtype"], axis=1)  # Remove the positional columns 
		mat = mat.drop(["annot_gene_id","annot_transcript_id","n_exons","length","gene_novelty","transcript_novelty","ISM_subtype"], axis=1)  # Remove the positional columns 
		# print(mat)

		# Applying CPM normalisation
		sums = mat.sum(axis=0, skipna=True)
		for col in mat.columns:
		    mat[col] = mat[col] / sums[col] * 1000000

		# Removing lowly expressed transcirpts
		mat = mat.loc[mat.sum(1) > 5]


		# Remove genes with only one transcript
		mat = mat.reset_index()
		mat = mat[mat.groupby('annot_gene_name')['annot_gene_name'].transform('size') > 1].reset_index(drop=True)

		# Calculating the average of each group
		for gr, samples in groups.items():
		    mat[gr] = mat[groups[gr]].mean(axis=1)
		    mat = mat.drop(groups[gr], axis=1)

		# Obtaining the transcript with with max counts per gene per groups
		group1 = mat.loc[mat.reset_index().groupby(['annot_gene_name'])[gr1].idxmax()]
		group1.drop(gr2, inplace=True, axis=1)
		group1 = group1.rename(columns={"annot_gene_name": "GeneID","annot_transcript_name": "TranscriptID_group1"})

		group2 = mat.loc[mat.reset_index().groupby(['annot_gene_name'])[gr2].idxmax()]
		group2.drop(gr1, inplace=True, axis=1)
		group2 = group2.rename(columns={"annot_gene_name": "GeneID_group2","annot_transcript_name": "TranscriptID_group2"})
		group2 = group2[['GeneID_group2', gr2, 'TranscriptID_group2']]  # Rearrange columns
		# print(group1,group2)
		# # print(group1.loc[group1['GeneID'] == 'ENSG00000235098.8'], group2.loc[group2['GeneID_group2'] == 'ENSG00000235098.8'])


		# Adding TranscriptID_group2 column
		group1 = group1.join(group2.set_index('GeneID_group2'), on='GeneID')
		# print(group1)

		# Calling MDT
		group1['MDT'] = np.where(group1['TranscriptID_group1'] == group1['TranscriptID_group2'], 'False', 'True') #create a new column in df1 to check if prices match
		group1.to_csv(f'{mdt_dir}/MDT_results.tsv',sep='\t',index=False)
		return
		
	def predict_productivity(self):
		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} PREDICTING ISOFORM PRODUCTIVITY')
		


		if not os.path.exists(predict_product): os.makedirs(predict_product)
		os.chdir(predict_product)

		# subprocess.run(f'mv {dte}/novel_de_transcripts* {predict_product}', shell=True)
		novel_transcripts_de_seqs = f'{predict_product}/novel_de_transcripts_for_funcAnnotation.fasta'
		novel_transcripts_de = f'{predict_product}/novel_de_transcripts_for_funcAnnotation.tsv'
		novel_transcripts_de_mtg = f'{predict_product}/novel_de_transcripts_for_funcAnnotation.gene_trans_map'
		transdecoder_predictions = f'{predict_product}/novel_de_transcripts_for_funcAnnotation.fasta.transdecoder.pep'


		novel_transcirpts = []
		# Obtaining the list on the novel transcripts that were differentially expressed (DTE)
		with open(novel_transcripts_de) as transcripts_in:
			for line in transcripts_in:
				novel_transcirpts.append(line.strip())
		
		# Isolating their reference sequences
		references = SeqIO.parse(f"{expression_analysis_dir}/reference_transcriptome.fasta", "fasta")
		with open(novel_transcripts_de_seqs, 'w') as ref_out:
			for elm in references:
				if elm.id in novel_transcirpts:
					SeqIO.write(elm, ref_out, "fasta")


		# Calling TransDecoder.LongOrfs to perform ORF prediction
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} 1/ TransDecoder.LongOrfs - TransDecoder identifies candidate coding regions within transcript sequences: in progress ..')
		transdecoder_orfs = ' '.join([
		"TransDecoder.LongOrfs",  # Calling TransDecoder.LongOrfs
		"-m 100",  #  Minimum protein length
		"-t", novel_transcripts_de_seqs,  # Novel isoforms in fasta format to be analysed
		"-O", predict_product,  # Path to intended output directory
		"--gene_trans_map", novel_transcripts_de_mtg,
		"2>>", os.path.join(pipeline_reports, "predproduct1.0_transdecoderORF-report.txt")])  # Directory where all reports reside
		subprocess.run(transdecoder_orfs, shell=True)

		# Use pfam_scan.pl to search the predicted fasta file against a library of Pfam HMMs
		# Search the peptides for protein domains using Pfam
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} 1.1/ Pfam - Searching the predicted fasta file against a library of Pfam HMMs: in progress ..')
		pfamHMMTransd = ' '.join([
		"pfam_scan.pl",  # Calling pfam_scan.pl
		"-cpu", args.threads,  # Number of parallel CPU workers to use for multithreads
		"-fasta", f'{predict_product}/longest_orfs.pep',  # Fasta file
		"-dir", pfam_dir,  # Directory location of Pfam files
		"-outfile", f'{predict_product}/resultsTransdecoder_pfam.txt',  # Output file
		"2>>", os.path.join(pipeline_reports, "predproduct1.1_pfamHMMTransd-report.txt")])  # Directory where all reports reside
		subprocess.run(pfamHMMTransd, shell=True)
		
		# Search a protein database Uniref90 (slow but more comprehensive) using BLASTP
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} 1.2/ BlastP - Search a protein database UniProt: in progress ..')
		blastpTransd = ' '.join([
		"blastp",  # Calling blastp
		"-query", f'{predict_product}/longest_orfs.pep',  # Query fasta file
		"-num_threads", args.threads,  # Number of threads (CPUs) to use in the BLAST search
		"-max_target_seqs 1",  # Maximum number of aligned sequences to keep
		"-evalue 1e-5",  # Expectation value (E) threshold for saving hits
		"-outfmt 6",  # Output formatting option, tabular
		"-db", uniprot,  # BLAST database name (UniRef90 db)
		"-out", f'{predict_product}/resultsTransdecoder_blastp.txt',  # Output file
		"2>>", os.path.join(pipeline_reports, "predproduct1.2_blastpTransd-report.txt")])  # Directory where all reports reside
		subprocess.run(blastpTransd, shell=True)

		# Calling TransDecoder.Predict to perform transcriptome protein prediction
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} 1.3/ TransDecoder.Predict - Performing transcriptome protein prediction: in progress ..')
		transdecoder_predict = ' '.join([
		"TransDecoder.Predict",  # Calling TransDecoder.Predict
		"--single_best_only",  # Retain only the single best orf per transcript (prioritized by homology then orf length)
		"--retain_pfam_hits", f'{predict_product}/resultsTransdecoder_pfam.txt',  # Domain table output file from running hmmscan to search Pfam
		"--retain_blastp_hits", f'{predict_product}/resultsTransdecoder_blastp.txt',  # BlastP output in '-outfmt 6' format
		"-t", novel_transcripts_de_seqs,  # Novel isoforms in fasta format to be analysed
		"-O", predict_product,  # Path to intended output directory
		"2>>", os.path.join(pipeline_reports, "predproduct1.3_transdecoderpredict-report.txt")])  # Directory where all reports reside
		subprocess.run(transdecoder_predict, shell=True)
		os.system(f'rm -r {dte}/*__check*')
		os.system('rm -r {0}/start_* {0}/*.cmds {0}/*.scores {0}/longest_* {0}/resultsTransdecoder*'.format(predict_product))
		

		# Search a protein database UniProt (slow but more comprehensive) using BLASTX
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} 2/ BlastX - Annotating coding regions of the novel nucleotide sequences: in progress ..')
		blastx = ' '.join([
		"blastx",  # Calling blastx
		"-query", novel_transcripts_de_seqs,  # Query fasta file
		"-num_threads", args.threads,  # Number of threads (CPUs) to use in the BLAST search
		"-max_target_seqs 1",  # Maximum number of aligned sequences to keep
		"-evalue 1e-5",  # Expectation value (E) threshold for saving hits
		"-outfmt 6",  # Output formatting option, tabular
		"-db", uniprot,  # BLAST database name (UniRef90 db)
		"-out", f'{predict_product}/results_blastx.txt',  # Output file
		"2>>", os.path.join(pipeline_reports, "predproduct2_blastx-report.txt")])  # Directory where all reports reside
		subprocess.run(blastx, shell=True)

		# Search a protein database UniProt (slow but more comprehensive) using BLASTP
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} 3/ BlastP - Search a protein database UniProt: in progress ..')
		blastp = ' '.join([
		"blastp",  # Calling blastp
		"-query", transdecoder_predictions,  # Query fasta file
		"-num_threads", args.threads,  # Number of threads (CPUs) to use in the BLAST search
		"-max_target_seqs 1",  # Maximum number of aligned sequences to keep
		"-evalue 1e-5",  # Expectation value (E) threshold for saving hits
		"-outfmt 6",  # Output formatting option, tabular
		"-db", uniprot,  # BLAST database name (UniRef90 db)
		"-out", f'{predict_product}/results_blastp.txt',  # Output file
		"2>>", os.path.join(pipeline_reports, "predproduct3_blastp-report.txt")])  # Directory where all reports reside
		subprocess.run(blastp, shell=True)

		# Search sequence(s) against a profile database
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} 4/ HMMER - Searching sequence databases for homologs of protein sequences: in progress ..')
		hmmscan = ' '.join([
		"hmmscan",  # Calling hmmscan
		"--cpu", args.threads,  # Number of parallel CPU workers to use for multithreads
		"--domtblout", f'{predict_product}/results_pfam.txt',  # Output file
		pfam_db,  # Directory location of Pfam database file
		transdecoder_predictions,  # Fasta file
		"2>>", os.path.join(pipeline_reports, "predproduct4_hmmscan-report.txt")])  # Directory where all reports reside
		subprocess.run(hmmscan, shell=True)

		# Predicting signal peptide and cleavage sites in gram+, gram- and eukaryotic amino acid sequences
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} 5/ SignalP - Predicting the presence of signal peptides and the location of their cleavage sites in proteins: in progress ..')
		run_signalp = ' '.join([
		"signalp",  # Calling SignalP
		"-batch 100000",  # Number of sequences that the tool will run simultaneously
		"-format short",  # Output format 'short' for the predictions without plots
		"-fasta", transdecoder_predictions,  # Input file in fasta format
		"-prefix", f'{predict_product}/results_signalp',  # Output file prefix
		"2>>", os.path.join(pipeline_reports, "predproduct5_signalp-report.txt")])  # Directory where all reports reside
		subprocess.run(run_signalp, shell=True)
		subprocess.run(f'mv {predict_product}/results_signalp* {predict_product}/results_signalp.txt', shell=True)

		# Prediction of transmembrane helices in proteins
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} 6/ TMHMM - Predicting of transmembrane helices in proteins: in progress ..')
		tmhmm = ' '.join([
		"tmhmm",  # Calling tmhmm
		"--short", 
		f"< {transdecoder_predictions}",  # Input file in fasta format
		">", f'{predict_product}/results_tmhmm.txt',  # Output file
		"2>>", os.path.join(pipeline_reports, "predproduct6_tmhmm-report.txt")])  # Directory where all reports reside
		subprocess.run(tmhmm, shell=True)
		os.system(f'rm -r {predict_product}/TMHMM*')
		
		# Prediction of ribosomal RNA sub units
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} 7/ RNAmmer - Prediction of ribosomal RNA sub units: in progress ..')
		run_rnammer = ' '.join([
		rnammerTransc,  # Calling RnammerTranscriptome.pl
		"--org_type euk",  # Selected organism
		"--path_to_rnammer", rnammer_exe,  # Path to RNAmmer
		"--transcriptome", novel_transcripts_de_seqs,  # Transcriptome assembly fasta file
		">", f'{predict_product}/results_rnammer.txt',  # Output file
		"2>>", os.path.join(pipeline_reports, "predproduct7_rnammer-report.txt")])  # Directory where all reports reside
		subprocess.run(run_rnammer, shell=True)
		os.system(f'rm {predict_product}/*rnammer.gff {predict_product}/transcriptSuperScaffold*')

		# Functional annotation of the novel transcripts
		self.functional_annotation(novel_transcripts_de_mtg, transdecoder_predictions, novel_transcripts_de_seqs)
		return

	def functional_annotation(self, novel_transcripts_de_mtg, transdecoder_predictions, novel_transcripts_de_seqs):
		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} PREDICTING FUNCTIONAL ANNOTATION OF THE NOVEL ISOFORMS')


		if not os.path.exists(functional_annot): os.makedirs(functional_annot)


		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Trinity - Running Trinotate: in progress ..')
		trinotate = ' '.join([
		f"Trinotate {trinity_sqlite} init",  # Calling Trinotate init
		"--gene_trans_map", novel_transcripts_de_mtg,
		"--transcript_fasta", novel_transcripts_de_seqs,
		"--transdecoder_pep", transdecoder_predictions,  # 
		"2>>", os.path.join(pipeline_reports, "predproduct8_trinotate-report.txt")])  # Directory where all reports reside
		subprocess.run(trinotate, shell=True)
		
		# Loading all results files
		### Transdecoder protein search results:
		subprocess.run(f"Trinotate {trinity_sqlite} LOAD_swissprot_blastp {predict_product}/results_blastp.txt", shell=True)  # Calling Trinotate LOAD_swissprot_blastp
		subprocess.run(f"Trinotate {trinity_sqlite} LOAD_pfam {predict_product}/results_pfam.txt", shell=True)
		subprocess.run(f"Trinotate {trinity_sqlite} LOAD_tmhmm {predict_product}/results_tmhmm.txt", shell=True)
		subprocess.run(f"Trinotate {trinity_sqlite} LOAD_signalp {predict_product}/results_signalp.txt", shell=True)
		### Trinity transcript search results
		subprocess.run(f"Trinotate {trinity_sqlite} LOAD_swissprot_blastx {predict_product}/results_blastx.txt", shell=True)
		subprocess.run(f"Trinotate {trinity_sqlite} LOAD_rnammer {predict_product}/results_rnammer.txt", shell=True)
		
		# Output report and results
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Trinity - Extracting a detailed report and GO Annotation: in progress ..')
		subprocess.run(f"Trinotate {trinity_sqlite} report > {functional_annot}/trinotate_annotation_report.tsv", shell=True)
		subprocess.run(f"extract_GO_assignments_from_Trinotate_xls.pl --trans --Trinotate_xls {functional_annot}/trinotate_annotation_report.tsv > {functional_annot}/trinotate_GO_annot.tsv", shell=True)
		os.system(f'mv {predict_product}/*.cds {predict_product}/*.pep -t {functional_annot}')
		return

	def combine_dte_dpa_results(self):
		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} COMBINING DTE, FUNCTIONAL ANNOTATION AND DPA RESULTS INTO A FINAL MATRIX')
		

		# Output the gene, isoform and ORF status (complete or partial) in a file
		func_annot_dict = {}
		with open(f'{functional_annot}/novel_de_transcripts_for_funcAnnotation.fasta.transdecoder.cds', 'r') as predictin, open(f'{functional_annot}/pred_isoform_status.tsv', 'w') as fout:
			for line in predictin:
				if line.startswith('>'):
					isoform = line.strip().split(" ")[0].split(".")[0].replace(">","")
					gene = line.strip().split(" ")[1].split("~")[0]
					orf_type = line.strip().split("type:")[1].split(" ")[0]
					func_annot_dict[isoform] = orf_type
					fout.write(f'{gene}\t{isoform}\t{orf_type}\n')

		
		polyA_res  = {}
		# Creating a dictionary to save the polyA results
		with open(f'{polyA_analysis_dir}/dpa_results/polya_diff_per_transcript.tsv') as polyin:
			for line in polyin:
				if not line.startswith(("count", "group", "contig", "\t")):
					contig = line.strip().split("\t")[0]
					median_control = line.strip().split("\t")[3]
					median_treatment = line.strip().split("\t")[5]
					median_diff = line.strip().split("\t")[4]
					fdr = line.strip().split("\t")[8]
					polyA_res[contig] = f'{median_control}\t{median_treatment}\t{median_diff}\t{fdr}'

		annot = {}
		with open(os.path.join(expression_analysis_dir, "database_talon.gtf")) as ref_in:
			for i, line in enumerate(ref_in):
				if not line.startswith("#"):
					if line.split("\t")[2].strip() == 'transcript':
						transcript_id = line.split("\t")[-1].split(";")[1].split(" ")[-1].strip("\"")
						transcript_type = ["novel", line.split("\t")[-1].split(";")[4].split()[1].strip("\"").strip()][transcript_id.startswith("ENST")]
						if transcript_type.endswith("pseudogene"):
							transcript_type = "pseudogene"
						elif transcript_type in ["Mt_rRNA","miRNA","snoRNA","misc_RNA","snRNA","scaRNA","rRNA","Mt_tRNA"]:
							transcript_type = "ncRNA"
						elif transcript_type == "ribozyme" or transcript_type.startswith("IG_") or transcript_type.startswith("TR_"):
							transcript_type = "protein_coding"
						elif transcript_type == "TEC":
							transcript_type = "To_be_Confirmed"
						annot[transcript_id] = transcript_type

		
		
		# Creating a new final file of DTE incorporating the polyA findings and functional annotation
		# dte_results = glob.glob(f'{dte}/*_edgeR_topTranscriptsBelow*LFC*.csv')[0]
		dte_results = glob.glob(f'{dte}/*edgeR_allTranscripts*.csv')[0]
		dte_pda_func_mat = dte_results.replace(".csv", ".funcAnnot.polyA.tsv")
		with open(dte_results) as dtein, open(dte_pda_func_mat, 'w')as fout:
			for line in dtein:
				if line.startswith("transcript_id"):
					fout.write(f'{line.strip()}\tpolyA-median_control\tpolyA-median_treatment\tpolyA-median_diff\tpolyA-FDR\tfunct_annot\trna_type\n')
				else:
					transcript = line.strip().split("\t")[0]
					if transcript in polyA_res:
						if transcript in func_annot_dict:
							fout.write(f'{line.strip()}\t{polyA_res[transcript]}\t{func_annot_dict[transcript]}\t{annot[transcript]}\n')
						else:
							if transcript.startswith("EN"):
								fout.write(f'{line.strip()}\t{polyA_res[transcript]}\tknown-functional\t{annot[transcript]}\n')
							else:
								fout.write(f'{line.strip()}\t{polyA_res[transcript]}\tunknown\t{annot[transcript]}\n')
					else:
						if transcript in func_annot_dict:
							fout.write(f'{line.strip()}\t-\t-\t-\t-\t{func_annot_dict[transcript]}\t{annot[transcript]}\n')
						else:
							if transcript.startswith("EN"):
								fout.write(f'{line.strip()}\t-\t-\t-\t-\tknown-functional\t{annot[transcript]}\n')
							else:
								fout.write(f'{line.strip()}\t-\t-\t-\t-\tunknown\t{annot[transcript]}\n')
		return 

	def fusions_events_jaffa(self):
	
		if not os.path.exists(fusions_events_jaffa): os.makedirs(fusions_events_jaffa)

		os.chdir(fusions_events_jaffa)
		jaffa = ' '.join([
		f"{jaffa_pipeline}/tools/bin/bpipe run",  # Calling jaffa
		f"{jaffa_pipeline}/JAFFAL.groovy",
		f'{alignments_dir}/*.filt.fastq',  
		"2>>", os.path.join(pipeline_reports, "fusion_events-jaffa-report.txt")])
		subprocess.run(jaffa, shell=True)
		return

	def annotate_fusions(self, fusion_file, sample_id, reference_annotation):
		""" Using Annotate Gene Fusion (AGFusion) for annotating the detected fusions """
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  AGFusion - Annotating Gene Fusions detected in {sample_id}: in progress ..')
		

		fusion_annot = f'{isoform_fusions}/{sample_id}'
		if not os.path.exists(fusion_annot): os.makedirs(fusion_annot)


		with open(fusion_file) as fusin:
			for line in fusin:
				fusion = line.strip().split("\t")[1].split(" ")[0].replace(":","_")
				agfusion_annotate = " ".join([
				"agfusion annotate",  # Calling AGFusion annotate
				"--gene5prime", reference_annotation[line.strip().split("\t")[1].split(":")[0]],  # 5' gene partner
				"--gene3prime", reference_annotation[line.strip().split("\t")[1].split(" ")[0].split(":")[1]],  #  3' gene partner
				"--junction5prime", line.strip().split("\t")[1].split(" ")[2].split(":")[1],  # Genomic location of predicted fuins for the 5' gene partner
				"--junction3prime", line.strip().split("\t")[1].split(" ")[3].split(":")[1],  # Genomic location of predicted fuins for the 3' gene partner
				"--protein_databases pfam smart superfamily tigrfam pfscan tmhmm seg ncoils prints pirsf signalp",
				"--database", ensembl_db,  # Path to the AGFusion database
				"--out", f"{fusion_annot}/{fusion}",  # Directory to save results
				"--type png",  # Image file type PNG
				"--width 10",  # Image width in inches
				"--height 3",  # Image file height in inches
				"--dpi 900",  # Dots per inch
				"--fontsize 10",  # Fontsize
				"--middlestar",  # Insert a * at the junction position for the cdna, cds, and protein sequences
				"2>>", os.path.join(pipeline_reports, "fusion_annotation-report.txt")])
				subprocess.run(agfusion_annotate, shell=True)
		os.system(f'mv {fusion_file} {fusion_annot}')
		return

def summary():
	
	multiQC = " ".join([
	"multiqc",  # Call MultiQC
	"--quiet",  # Print only log warnings
	"--outdir", postAlignment_reports,  # Create report in the FastQC reports directory
	"--filename", "post-alignment_summarised_report",  # Name of the output report 
	postAlignment_reports,  # Directory where all FastQC and Cutadapt reports reside
	"2>>", os.path.join(pipeline_reports, "summary_postalignment_multiQC-report.txt")])  # Output multiQC report
	subprocess.run(multiQC, shell=True)

	pickle_files = [pk_file for pk_file in glob.glob(os.path.join(postAlignment_reports, "*_qc.pk"))]
	### Wub Compare alignment QC statistics of multiple samples
	bam_multi_qc = ' '.join([
	"bam_multi_qc.py",
	"-r", "{0}/comparison_qc.pdf".format(postAlignment_reports),  # Output pdf file
	' '.join(pickle_files),
	"2>>", os.path.join(pipeline_reports, "summary_multi_qc-report.txt")])
	# subprocess.run(bam_multi_qc, shell=True)

	## Cleaning up reports folder
	figures_pre_dir = os.path.join(initial_qc_reports, "individual_plots")
	if not os.path.exists(figures_pre_dir): os.makedirs(figures_pre_dir)

	os.system('mv {0}/*.png {1}'.format(initial_qc_reports, figures_pre_dir))  # Moving png files to "figures" directory
	
	## REMOVING UNNECESSARY FILES & REPORTS (postanalysis)
	for path, subdir, folder in os.walk(pipeline_reports):
		for name in folder:
			file = os.path.join(path, name)
			if os.stat(file).st_size == 0 or\
			(name.endswith("multiQC-report.txt") and os.stat(file).st_size == 583) or\
			(name.endswith("bamstats-report.txt") and os.stat(file).st_size == 24) or\
	  		(name.endswith("ribution-report.txt") and os.stat(file).st_size == (255 or 256)) or\
	  		(name.endswith("Metrics-report.txt") and os.stat(file).st_size == 3) or\
	  		(name.endswith("bamstats-report.txt") and os.stat(file).st_size == 24) or\
	  		(name.endswith("R_transcriptType_sum-report.txt") and os.stat(file).st_size == 210) or\
	  		(name.endswith("diffExpr_exploratory_analysis-report.txt") and os.stat(file).st_size == 3648) or\
	  		(name.endswith("diffExpr_dge_analysis-report.txt") and os.stat(file).st_size == 347) or\
	  		(name.endswith("strandness-report.txt") and os.stat(file).st_size == 205 or os.stat(file).st_size == 214) or\
	  		 name.endswith(".r") or name.endswith("DupRate.xls") or name.endswith("duplicate_reads-report.txt"):
				os.remove(file)

	## Cleaning up postanalysis folder
	individual_postal_reports = os.path.join(postAlignment_reports, "individual_reports")
	if not os.path.exists(individual_postal_reports): os.makedirs(individual_postal_reports)

	# Moving files to "qc_reports" directory
	# os.system('mv {0}/*_qc.pk {1}'.format(postAlignment_reports, qc_reports))
	os.system('mv {0}/*.txt {1}'.format(postAlignment_reports, individual_postal_reports))
	os.system('mv {0}/*.fragSize {1}'.format(postAlignment_reports, individual_postal_reports))

	### Removing unnecessary files
	# subprocess.run(f'rm {alignments_dir}/*.fastq*', shell=True)  # Removing the extracted fastq
	return



def main():
	


	summary_files = [str(file_path) for file_path in Path(ont_data).glob('**/sequencing_summary.txt') if not "warehouse" in str(file_path)]
	num_of_samples = len(summary_files)

	for sum_file in sorted([s for s in summary_files if os.path.dirname(s).endswith(chosen_samples)]):
		raw_data_dir = os.path.dirname(str(sum_file))
		sample_id = os.path.basename(raw_data_dir)
		fastq_pass = " ".join(glob.glob(os.path.join(raw_data_dir, "pass/*pass.fastq.gz")))
		print(f'\nPROCESSING SAMPLE {sample_id}')
		
		# quality_control(sum_file, sample_id, raw_data_dir)

		alignment_against_ref(fastq_pass, sample_id, raw_data_dir, sum_file)

		# polyA_estimation(sample_id, sum_file, fastq_pass, raw_data_dir)


		# methylation_detection_xpore(sample_id, fastq_pass, raw_data_dir)

	# expression_analysis()

	# downstream_analysis()

	# methylation_detection_eligos2()

	# summary()

	print(f'\t--- The pipeline finisded after {datetime.now() - startTime} ---')

if __name__ == "__main__": main()