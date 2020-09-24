###Stavros Giannoukakos### 
#Version of the program
__version__ = "0.1.7"

import argparse
import subprocess
import numpy as np
from pathlib import Path
from datetime import datetime
import shutil, time, glob, sys, os, re

ont_data =  "/home/stavros/playground/ont_basecalling/guppy_v3_basecalling"

refGenomeGRCh38 = "/home/stavros/references/reference_genome/GRCh38_GencodeV31_primAssembly/GRCh38.primary_assembly.genome.fa"
refGenomeGRCh38_traclean = "/home/stavros/references/reference_genome/GRCh38_GencodeV31_primAssembly/GRCh38.primary_assembly.genome.edited.fa"
refTranscGRCh38_ensembl = "/home/stavros/references/reference_transcriptome/Ensembl/GRCh38.cdna.ncrna.fa"
refAnnot = "/home/stavros/references/reference_annotation/GRCh38_gencode.v31.primary_assembly.annotation.gtf"
refAnnot_ensembl = "/home/stavros/references/reference_annotation/Ensembl/Homo_sapiens.GRCh38.99.gtf"
reference_annotation_bed = "/home/stavros/references/reference_annotation/hg38_gencode.v31.allComprehensive_pseudo.annotation.bed"
rscripts = "{0}/Rscripts".format(os.path.dirname(os.path.realpath(__file__)))

talon_database = "/home/stavros/playground/progs/TALON/talon_db/talon_gencode_v31.db"  ### TALON
transcriptclean = "python3 /home/stavros/playground/progs/TranscriptClean/TranscriptClean.py"  ### TranscriptClean
sj_ref = "/home/stavros/playground/progs/TranscriptClean/GRCh38_SJs.ref"  ### TranscriptClean ref. SJ
chrom_size = "/home/stavros/references/reference_genome/GRCh38_GencodeV31_primAssembly/chrom_lenghts.tsv"   ### FLAIR



usage = "nanodRNA_analysis [options]"
epilog = " -- June 2019 | Stavros Giannoukakos -- "
description = "DESCRIPTION"

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage, description=description, epilog=epilog)
# Number of threads/CPUs to be used
parser.add_argument('-t', '--threads', dest='threads', default=str(40), metavar='', 
                	help="Number of threads to be used in the analysis")
# Genes expressed in minimum this many samples
parser.add_argument('-msge', '--minSampsGeneExpr', dest='minSampsGeneExpr', default=str(3), metavar='', 
                	help="A transcript  must be mapped to a gene in at\nleast this minimum number of samples for the\ngene be included in the analysis")
# Transcripts expressed in minimum this many samples
parser.add_argument('-msfe', '--minSampsFeatureExpr', dest='minSampsFeatureExpr', default=str(2), metavar='', 
                	help="A transcript must be mapped to an isoform at\nleast this minimum number of samples for the\ngene isoform to be included in the analysis")
# Minimum gene counts
parser.add_argument('-mge', '--minGeneExpr', dest='minGeneExpr', default=str(10), metavar='', 
                	help="Min. number of total mapped sequence reads\nfor a gene to be considered expressed")
# Minimum transcript counts
parser.add_argument('-mfe', '--minFeatureExpr', dest='minFeatureExpr', default=str(10), metavar='', 
                	help="Min. number of total mapped sequence reads\nfor a gene isoform to be considered")
# Adjusted p-value threshold for differential expression analysis
parser.add_argument('-p', '--adjPValueThreshold', dest='adjPValueThreshold', default=str(0.05), metavar='', 
                	help="Adjusted p-value threshold for differential\nexpression analysis")
# Minimum required log2 fold change for differential expression analysis
parser.add_argument('-lfc', '--lfcThreshold', dest='lfcThreshold', default=str(1), metavar='', 
                	help="Minimum required log2 fold change for diffe-\nrential expression analysis")
# Top N genes to be used for the heatmap
parser.add_argument('-n', '--n_top', dest='n_top', default=str(50), metavar='', 
                	help="Top N genes to be used for the heatmap")
# Transcripts expressed in minimum this many samples (for DTU)
parser.add_argument('-msfedtu', '--minSampsFeatureExprDTU', dest='minSampsFeatureExprDTU', default=str(1), metavar='', 
                	help="A transcript must be mapped to an isoform at\nleast this minimum number of samples for the\ngene isoform to be included in the analysis\nfor DTU only")
# Minimum transcript counts for DTU
parser.add_argument('-mfedtu', '--minFeatureExprDTU', dest='minFeatureExprDTU', default=str(3), metavar='', 
                	help="Minimum number of total mapped sequence reads\nfor a gene isoform to be considered for DTU")
# Display the version of the pipeline 
parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(__version__))
# Get the options and return them
args = parser.parse_args()

script_dir = os.path.dirname(os.path.realpath(__file__))
startTime = datetime.now()

# Main folder hosting the analysis
analysis_dir = os.path.join(script_dir, "analysis_batch3")
prepr_dir = os.path.join(analysis_dir, "preprocessed_data")
alignments_dir = os.path.join(analysis_dir, "alignments")
reports_dir = os.path.join(analysis_dir, "reports")

initial_qc_reports = os.path.join(analysis_dir, "reports/initial_qc_reports")
postAlignment_reports = os.path.join(analysis_dir, "reports/postAlignment_qc_reports")
pipeline_reports = os.path.join(analysis_dir, "reports/pipeline_reports")

expression_analysis_dir = os.path.join(analysis_dir, "expression_analysis")
polyA_analysis_dir = os.path.join(analysis_dir, "polyA_estimation")

if not os.path.exists(pipeline_reports): os.makedirs(pipeline_reports)

def quality_control(seq_summary_file, sample_id, raw_data_dir):
	# Producing preliminary QC reports 
	if not os.path.exists(initial_qc_reports): os.makedirs(initial_qc_reports)
	print("\t{0} QUALITY CONTROL OF THE INPUT SAMPLES".format(datetime.now().strftime("%d.%m.%Y %H:%M")))

	# Unsing NanoPlot (github.com/wdecoster/NanoPlot) to extract basic information regarding the sequencing
	print("{0}  nanoPlot - Quality Control of the input {1} raw data: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), sample_id))
	nanoPlot = " ".join([
	"NanoPlot",  # Call pycoQC
	"--threads", args.threads,  # Number of threads to be used by the script
	"--summary", seq_summary_file,  # Input of <sequencing_summary> generated by Albacore1.0.0 / Guppy 2.1.3+
	"--color saddlebrown",  # Color for the plots
	"--colormap PuBuGn",  # Colormap for the heatmap
	"--prefix", os.path.join(initial_qc_reports, "{0}_".format(sample_id)),  # Create the report in the reports directory
	"--format png",  # Output format of the plots
	"--dpi 900", # Set the dpi for saving images in high resolution
	"2>>", os.path.join(pipeline_reports, "nanoPlot-prelim_report.txt")])
	subprocess.run(nanoPlot, shell=True)
	return

def alignment_against_ref(fastq_pass, sample_id, raw_data_dir, seq_summary_file):
	if not os.path.exists(alignments_dir): os.makedirs(alignments_dir)
	print("\n\t{0} ALIGNING AGAINSTE THE REF. GENOME AND TRANSCRIPTOME".format(datetime.now().strftime("%d.%m.%Y %H:%M")))

	### ALIGN THE RAW READS AGAINST THE REFERENCE GENOME
	print("{0}  Minimap2 - Mapping {1} against the reference genome: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), sample_id))
	bam_file = os.path.join(alignments_dir, "{0}.genome.bam".format(sample_id))
	minimap2_genome = " ".join([
	"minimap2",  # Call minimap2 (v2.17-r941)
	"-t", args.threads,  # Number of threds to use
	"-ax splice",   # Long-read spliced alignment mode and output in SAM format (-a)
	"-k 14",  # k-mer size
	"-uf",  # Find canonical splicing sites GT-AG - f: transcript strand
	"--secondary=no",  # Do not report any secondary alignments
	"--MD",  # output the MD tag
	# "-o", os.path.join(alignments_dir, "{0}.genome.paf".format(sample_id)),
	refGenomeGRCh38,  # Inputting the reference genome
	fastq_pass,  # Input .fastq.gz file
	"|" "samtools view",
	"--threads", args.threads,  # Number of threads to be used by 'samtools view'
	"-Sb",
	"-q 10",  # Filterring out reads with mapping quality lower than 10
	"|", "samtools sort",  # Calling 'samtools sort' to sort the output alignment file
	"--threads", args.threads,  # Number of threads to be used by 'samtools sort'
	"-",  # Input from standard output
	"-o", bam_file,  # Sorted output  BAM file
	"2>>", os.path.join(pipeline_reports, "minimap2_genome-report.txt")])  # Directory where all reports reside
	subprocess.run(minimap2_genome, shell=True)
	subprocess.run('samtools index -@ {0} {1}'.format(args.threads, bam_file), shell=True)
	
	mapping_qc(sample_id, seq_summary_file, bam_file)
	return

def mapping_qc(sample_id, seq_summary_file, bam_file):
	print("\n\t{0} GENERATING ALIGNMENT STATS".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	if not os.path.exists(postAlignment_reports): os.makedirs(postAlignment_reports)

	# BAM stats
	print("{0}  BAM stats - Generating post-alignment stats of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), sample_id))
	bam_stat = ' '.join([
	"bam_stat.py",
	"-i", bam_file,  # Input BAM file
	"> {0}/{1}.bamstat.txt".format(postAlignment_reports, sample_id),  # Output file
	"2>>", os.path.join(pipeline_reports, "bamstats-report.txt")])
	subprocess.run(bam_stat, shell=True)

	# Picard CollectAlignmentSummaryMetrics
	print("{0}  Picard - Collecting alignment summary metrics of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), sample_id))
	CollectAlignmentSummaryMetrics = ' '.join([
	"picard-tools CollectAlignmentSummaryMetrics",  # Call picard-tools CollectAlignmentSummaryMetrics
	"INPUT= {0}".format(bam_file),  # Input BAM file
	"OUTPUT= {0}/{1}.alignment_metrics.txt".format(postAlignment_reports, sample_id),  # Output
	"REFERENCE_SEQUENCE= {0}".format(refGenomeGRCh38),  # Reference sequence file
	"2>>", os.path.join(pipeline_reports, "collectAlignmentSummaryMetrics-report.txt")])
	subprocess.run(CollectAlignmentSummaryMetrics, shell=True) 

	# # AlignQC 
	# print("{0}  AlignQC - Generating post-alignment stats of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), sample_id))
	# alignqc = ' '.join([
	# "alignqc analyze",  # Calling AlignQC analyze
	# "--threads", str(int(int(args.threads)/2)),
	# "--genome", refGenomeGRCh38,  # Reference in .fasta
	# "--gtf", "{0}.gz".format(refAnnot),  # Input reference annotation in gzipped form
	# "--output", "{0}/{1}.alignqc.xhtml".format(postAlignment_reports, sample_id),  # Output pdf file
	# bam_file,
	# "2>>", os.path.join(pipeline_reports, "alignqc-report.txt")])
	# subprocess.run(alignqc, shell=True)

	# # Wub Alignment based QC plots
	# print("{0}  WUB - Alignment based QC plots of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), sample_id))
	# alignment_qc = ' '.join([
	# "bam_alignment_qc.py",
	# "-f", refGenomeGRCh38,  # Input reference file
	# "-x -Q",  # Do not plot per-reference information/ Keep qiet
	# "-r", "{0}/{1}.bam_alignment_qc.pdf".format(postAlignment_reports, sample_id),  # Output pdf file
	# "-p", "{0}/{1}.bam_alignment_qc.pk".format(postAlignment_reports, sample_id),  # Output pk file
	# bam_file,
	# "2>>", os.path.join(pipeline_reports, "alignment_qc-report.txt")])
	# subprocess.run(alignment_qc, shell=True)

	### EXPORTING ALIGNMENT STATS
	print("{0}  PycoQC - Generating post-alignment stats of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), sample_id))
	pycoQC = " ".join([
	"pycoQC",  # Call pycoQC
	"--summary_file", seq_summary_file, 
	"--bam_file", bam_file,  # Input of bam files from Minimap2
	"--html_outfile", os.path.join(postAlignment_reports, "{0}.pycoQC-report.html".format(sample_id)),  # Create the report in the reports directory
	"--report_title", "\"Post-alignment quality control report\"",  # A title to be used in the html report
	"2>>", os.path.join(pipeline_reports, "pycoQC_postAlignment-report.txt")])
	subprocess.run(pycoQC, shell=True)

	# BAM read distribution
	print("{0}  RSeQC - Generating read distribution stats of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), sample_id))
	read_distribution = ' '.join([
	"read_distribution.py", 
	"-i", bam_file,  # Input BAM file
	"-r", reference_annotation_bed,
	"> {0}/{1}.fragSize".format(postAlignment_reports, sample_id),  # Output file
	"2>>", os.path.join(pipeline_reports, "read_distribution-report.txt")])
	subprocess.run(read_distribution, shell=True)

	# Check the strandness of the reads
	print("{0}  RSeQC - Generating read strandness stats of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), sample_id))
	strandness = ' '.join([
	"infer_experiment.py",  # Call samtools flagstat
	"-i", bam_file,  # Input BAM file
	"-r", reference_annotation_bed,
	"> {0}/{1}.strandness.txt".format(postAlignment_reports, sample_id),  # Output file
	"2>>", os.path.join(pipeline_reports, "strandness-report.txt")])
	subprocess.run(strandness, shell=True)

	# Check duplicate reads
	print("{0}  Picard - Extracting read duplication stats of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), sample_id))
	duplicate_reads = ' '.join([
	"picard-tools MarkDuplicates",  # Call samtools flagstat
	"INPUT= {0}".format(bam_file),  # Input BAM file
	"OUTPUT= {0}/{1}.genome.dedup.bam".format(alignments_dir, sample_id),
	"METRICS_FILE= {0}/{1}.mark_duplicates.txt".format(postAlignment_reports, sample_id),  # Output file
	"2>>", os.path.join(pipeline_reports, "picard_markDuplicate_reads-report.txt")])
	subprocess.run(duplicate_reads, shell=True)
	os.system('rm {0}/{1}.genome.dedup.bam'.format(alignments_dir, sample_id))

	# Number of reads mapped to each chromosome
	print("{0}  RSeQC - Generating mapping stats of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), sample_id))
	mapping_pos = ' '.join([
	"samtools idxstats",  # Call samtools flagstat
	bam_file,  # Input BAM file
	"> {0}/{1}.samtools_idxstats.txt".format(postAlignment_reports, sample_id),  # Output file
	"2>>", os.path.join(pipeline_reports, "samtools_idxstats-report.txt")])
	subprocess.run(mapping_pos, shell=True)
	return

def polyA_estimation(sample_id, sum_file, fastq_pass, raw_data_dir):

	print("\n\t{0} POLYA LENGTH ESTIMATION".format(datetime.now().strftime("%d.%m.%Y %H:%M")))

	polyA_analysis_dir_idv = os.path.join(polyA_analysis_dir, "{0}".format(sample_id))
	if not os.path.exists(polyA_analysis_dir_idv): os.makedirs(polyA_analysis_dir_idv)

	# Extracting the *pass.fastq.gz file 
	extracted_fastq = os.path.join(polyA_analysis_dir_idv, os.path.basename(fastq_pass)[:-3])
	os.system('gzip --decompress --keep --stdout {0} > {1}'.format(fastq_pass, extracted_fastq))
	
	# Nanopolish index
	print("{0}  Nanopolish index - Indexing the output of the guppy basecaller of sample {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), sample_id))
	indexing = ' '.join([
	"nanopolish index",  # Calling nanopolish index
	"--sequencing-summary", sum_file,  # the sequencing summary file from Guppy
	"--directory", "{0}/workspace/fast5_pass".format(raw_data_dir),  # Input BAM file
	extracted_fastq,
	"2>>", os.path.join(pipeline_reports, "nanopolish_index-report.txt")])
	subprocess.run(indexing, shell=True)

	# # Nanopolish polyA
	print("{0}  Nanopolish polyA -  Estimate the polyadenylated tail lengths of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), sample_id))
	polyA_est = ' '.join([
	"nanopolish polya",  # Calling nanopolish polya
	"--threads", args.threads,  # Number of threads to use
	"--genome", refGenomeGRCh38,  # The reference genome assembly that  was used
	"--reads", extracted_fastq,  # The raw 1D ONT direct RNA reads in fastq
	"--bam", os.path.join(alignments_dir, "{0}.genome.bam".format(sample_id)),  # The reads aligned to the genome assembly in BAM format
	"> {0}/{1}.polya_results.tsv".format(polyA_analysis_dir_idv, sample_id),  # Output file
	"2>>", os.path.join(pipeline_reports, "nanopolish_polyA-report.txt")])
	subprocess.run(polyA_est, shell=True)
	os.system("rm {0}".format(extracted_fastq))  # Removing the extracted fastq
	return

class expression_matrix:

	def __init__(self, chosen_samples):
		if not os.path.exists(expression_analysis_dir): os.makedirs(expression_analysis_dir)
		self.genome_perGene_em_featureCounts()
		self.novel_transcripts_detection_talon()
		# self.novel_transcripts_detection_flair(chosen_samples)
		return

	def genome_perGene_em_featureCounts(self):
		print("\n\t{0} GENERATING THE PER-GENE EXPRESSION MATRIX".format(datetime.now().strftime("%d.%m.%Y %H:%M")))

		featureCounts_analysis = os.path.join(expression_analysis_dir, 'featureCounts_analysis')
		if not os.path.exists(featureCounts_analysis): os.makedirs(featureCounts_analysis)
		genome_alignments = [sum_file for sum_file in glob.glob(os.path.join(alignments_dir, "*.genome.bam"))]
		
		print("{0}  FeatureCounts perGene - Counting reads from the genome aligned data: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		featureCounts_gn = " ".join([
		"featureCounts",  # Call featureCounts
		"-T", args.threads,  # Number of threads to be used by the script
		"-a", refAnnot,  # Annotation file in GTF/GFF format
		"-L",  # Count long reads such as Nanopore and PacBio reads
		"-o", os.path.join(featureCounts_analysis, "genome_alignments_perGene_sum.tab"),
		' '.join(genome_alignments),  # Input bam file
		"2>>", os.path.join(pipeline_reports, "featureCounts_genome_summary_perGene-report.txt")])  # Directory where all reports reside
		subprocess.run(featureCounts_gn, shell=True)
		
		print("{0}  Exporting the per gene (genomic) expression matrix: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		subprocess.run("cut -f1,7- {0}/genome_alignments_perGene_sum.tab | sed 1d > {0}/featureCounts_expression_perGene_matrix.txt".format(featureCounts_analysis), shell=True)

		annot = {}
		with open(refAnnot) as ref_in:
			for i, line in enumerate(ref_in):
				if not line.startswith("#"):
					if line.split("\t")[2].strip() == 'transcript' or line.split("\t")[2].strip().endswith('tRNA'):
						gene_id = line.split("\t")[-1].split(";")[0].split(" ")[-1].strip("\"")
						gene_type = [line.split("\t")[-1].split(";")[1].split()[1].strip("\"").strip() ,line.split("\t")[-1].split(";")[2].split()[1].strip("\"").strip()]\
									[line.split("\t")[-1].split(";")[2].split()[0].strip()=="gene_type"]
						annot[gene_id] = gene_type

		data = {}
		header = []
		# nc = ["miRNA", "piRNA", "rRNA", "siRNA", "snRNA", "snoRNA", "tRNA", "vaultRNA"]
		with open("{0}/featureCounts_expression_perGene_matrix.txt".format(featureCounts_analysis)) as mat_in:
			for i, line in enumerate(mat_in, 1):
				if i == 1:
					header = line.split()
				else:
					gene = line.strip().split()[0]
					values = line.strip().split()[1:]
					if not sum(map(int, values)) == 0:
						# Grouping all pseudogenes
						if annot[gene].endswith("pseudogene"):
							data[(gene, "pseudogene")] = values
						# Grouping all immunoglobin genes
						elif annot[gene].startswith("IG_"):
							data[(gene, "Immunoglobulin_gene")] = values
						# Grouping all T-cell receptor genes
						elif annot[gene].startswith("TR_"):
							data[(gene, "Tcell_receptor_gene")] = values
						# Reanming TEC group
						elif annot[gene].startswith("TEC"):
							data[(gene, "To_be_Experimentally_Confirmed")] = values
						# Rest
						else:
							data[(gene, annot[gene])] = values

		# Remove paths from sample names
		header = [elm.split("/")[-1].split(".")[0] for elm in header]
		header.insert(1, "gene_type")  # Inserting gene_type in header
		header[0] = "gene_id"  # Replacing Geneid with gene_id in header

		# Writing output to file 'expression_matrix.csv'
		with open("{0}/perGene_expression_matrix.csv".format(featureCounts_analysis), "a") as fout:
			fout.write("{0}\n".format(','.join(header)))
			for key, values in data.items():
				fout.write("{0},{1}\n".format(','.join(key), ','.join(values)))

		print("{0}  Visualising the RNA categories found in the (genomic) expression matrix: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		gene_type_sum = " ".join([
		"Rscript",  # Call Rscript
		"{0}/gene_type_summary.R".format(rscripts),  # Calling the script
		"{0}/perGene_expression_matrix.csv".format(featureCounts_analysis),  # Input matrix
		featureCounts_analysis,  # Output dir
		"2>>", os.path.join(pipeline_reports, "R_gene_type_sum-report.txt")])  # Directory where all reports reside
		subprocess.run(gene_type_sum, shell=True)
		os.system('rm {0}/featureCounts_expression_perGene_matrix.txt {0}/genome_alignments_perGene_sum.tab'.format(featureCounts_analysis))
		return

	def novel_transcripts_detection_talon(self):
		""" TALON takes transcripts from one or more long read datasets (SAM format) 
		and assigns them transcript and gene identifiers based on a database-bound
		annotation. Novel events are assigned new identifiers """
		print("\n\t{0} ANNOTATION AND QUANTIFICATION USING TALON".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		talon_analysis = os.path.join(expression_analysis_dir, 'talon_analysis')
		temp = os.path.join(talon_analysis, 'transcriptClean_temp')
		if not os.path.exists(temp): os.makedirs(temp)
		R_analysis = os.path.join(talon_analysis, "R_analysis")
		if not os.path.exists(R_analysis): os.makedirs(R_analysis)

		# Create config file that is needed for talon and converting the aligned bam files to sam
		csv_file = os.path.join(talon_analysis,'talon_input.csv')
		if not os.path.exists(csv_file) or os.stat(csv_file).st_size == 0:
			print("{0}  Creating a description file necessary for Talon: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
			with open(csv_file, "w") as talon_out:
				for path, subdir, folder in os.walk(alignments_dir):
					for file in sorted(folder):
						if file.endswith('.genome.bam'):
							# Output a csv file with 'sample_id, sample_group, technology, input_labeled_sam_file' which is gonna be needed in talon_annotation function
							talon_out.write('{0},{1},ONT,{2}\n'.format(file.split(".")[0], file.split(".")[0].split("_")[0], os.path.join(temp, file).replace(".genome.bam", ".clean_labeled.sam")))
							if not os.path.exists(os.path.join(path, file).replace(".bam", ".sam")):
								print("{0}  Samtools - Converting the genomic bam file ({1}) to sam for Talon: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), file))
								os.system('samtools view -h -@ {0} {1} > {2}'.format(args.threads, os.path.join(path, file), os.path.join(path, file).replace(".bam", ".sam")))


		### Initial step: Correcting mismatches, microindels, and noncanonical splice junctions in long reads that have been mapped to the genome
		sam_files = glob.glob(os.path.join(alignments_dir, "*.genome.sam"))
		print("{0} 1/ TranscriptClean - Correcting mismatches, microindels, and noncanonical splice junctions in long reads: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
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
			"2>>", os.path.join(pipeline_reports, "talon2_transcriptclean-report.txt")])  # Directory where all reports reside
			subprocess.run(TranscriptClean, shell=True)
			os.remove(file)
		
		### Second step: Building the reference database
		print("{0} 2/ TALON - Initiating the database: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		if os.path.exists(talon_database): os.remove(talon_database)
		talon_initialize_database = " ".join([
		"talon_initialize_database",  # Call talon_initialize_database
		"--f", refAnnot,  # GTF annotation containing genes, transcripts, and edges
		"--g", 'hg38',  # Genome build (i.e. hg38) to use
		"--5p 500",  # Maximum allowable distance (bp) at the 5' end during annotation
		"--3p 300",  # Maximum allowable distance (bp) at the 3' end during annotation
		"--a", talon_database[:-3],  # Name of supplied annotation
		"--o", talon_database[:-3],  # Outprefix for the annotation files
		"2>>", os.path.join(pipeline_reports, "talon1_initialize_database-report.txt")])  # Directory where all reports reside
		subprocess.run(talon_initialize_database, shell=True)

		### Third step: internal priming check
		print("{0}  3/ TALON - Run talon_label_reads on each file to compute how likely each read is to be an internal priming product: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		genome_alignments = [sum_file for sum_file in glob.glob(os.path.join(temp, "*clean.sam"))]
		for file in genome_alignments:
			talon_priming_check = " ".join([
			"talon_label_reads",  # Call talon talon_label_reads
			"--t", args.threads,  # Number of threads to be used by the script
			"--f", file,  # Input sam file
			"--g", refGenomeGRCh38,  # Reference genome fasta file
			"--tmpDir", os.path.join(talon_analysis, "tmp_label_reads"),  # Path to directory for tmp files
			"--deleteTmp",  # Temporary directory will be removed
			"--o", os.path.join(temp, os.path.basename(file).replace("_clean.sam",".clean")),  # Prefix for output files
			"2>>", os.path.join(pipeline_reports, "talon3_talon_priming_check-report.txt")])  # Directory where all reports reside
			subprocess.run(talon_priming_check, shell=True)

		### Fourth step: annotating and quantification of the reads
		print("{0}  4/ TALON - Annotating and quantification of the reads: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		talon_annotation = " ".join([
		"talon",  # Call talon
		"--threads", args.threads,  # Number of threads to be used by the script
		"--db", talon_database,  # TALON database
		"--build", 'hg38',  # Genome build (i.e. hg38) to use
		"--o", os.path.join(talon_analysis, "prefilt"),  # Prefix for output files
		"--f", csv_file,  # Dataset config file: dataset name, sample description, platform, sam file (comma-delimited)
		"2>>", os.path.join(pipeline_reports, "talon4_annotationNquantification-report.txt")])  # Directory where all reports reside
		subprocess.run(talon_annotation, shell=True)
		
		### Fifth step: summarising how many of each transcript were found (prior to any filtering)
		print("{0}  5/ TALON - Summarising how many of each transcript were found (prior to any filtering): in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		talon_summary = " ".join([
		"talon_summarize",  # Call talon_summarize
		"--db", talon_database,  # TALON database
		"--o", os.path.join(talon_analysis, "prefilt"),  # Prefix for output file
		"2>>", os.path.join(pipeline_reports, "talon5_summary-report.txt")])  # Directory where all reports reside
		subprocess.run(talon_summary, shell=True)
		
		### Sixth step: creating an abundance matrix without filtering (for use computing gene expression)
		print("{0}  6/ TALON - Creating an abundance matrix without filtering (gene expression): in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		talon_abundance = " ".join([
		"talon_abundance",  # Call talon_abundance
		"--db", talon_database,  # TALON database
		"--build", 'hg38',  # Genome build (i.e. hg38) to use
		"--annot", talon_database[:-3],  # Which annotation version to use
		"--o", os.path.join(talon_analysis, "prefilt"),  # Prefix for output file
		"2>>", os.path.join(pipeline_reports, "talon6_abundance-report.txt")])  # Directory where all reports reside
		subprocess.run(talon_abundance, shell=True)

		### Seventh step: Applying basic filtering steps and outputting several stats
		print("{0}  7/ TALON - Removing internal priming artifacts low abundance transcripts: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		talon_filter = " ".join([
		"talon_filter_transcripts",  # Call talon_filter_transcripts
		"--db", talon_database,  # TALON database
		"--annot", talon_database[:-3],  # Which annotation version to use
		"--minCount 1",  # Number of minimum occurrences required for a novel transcript PER dataset
		"--maxFracA 0.5",  # All of the supporting reads must have 50% or fewer As in the 20 bp interval after alignment
		"--minDatasets 1",  # Minimum number of datasets novel transcripts must be found in
		"--o", os.path.join(talon_analysis, "filtered_isoforms.csv"),  # Output
		"2>>", os.path.join(pipeline_reports, "talon7_talon_filter-report.txt")])  # Directory where all reports reside
		subprocess.run(talon_filter, shell=True)

		# ### Eighth step: Applying basic filtering steps and outputting several stats
		print("{0}  8/ TALON - Removing low abundance isoforms based on step 7 and exporting basic statsistics: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		talon_filter_n_report = " ".join([
		"Rscript",  # Call Rscript
		"{0}/talon_summarisation.R".format(rscripts),  # Calling the talon_summarisation.R script
		os.path.join(talon_analysis, "prefilt_talon_abundance.tsv"),  # Input matrix
		R_analysis,  # Output directory
		os.path.join(talon_analysis, "talon_input.csv"),  # Input annotation matrix
		os.path.join(talon_analysis, "filtered_isoforms.csv"),  # Filtered transcripts to maintain
		"2>>", os.path.join(pipeline_reports, "talon8_summarisation.txt")])  # Directory where all reports reside
		subprocess.run(talon_filter_n_report, shell=True)

		### Ninth step: Generating TALON report for each dataset
		print("{0}  9/ TALON - Generating TALON report for each dataset: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		for file in sam_files:
			sample_name = os.path.basename(file).split(".")[0]
			output_dir = os.path.join(R_analysis, f"talon_reports/{sample_name}_report")
			if not os.path.exists(output_dir): os.makedirs(output_dir)
			talon_report = " ".join([
			"talon_generate_report",  # Call talon_abundance
			"--db", talon_database,  # TALON database
			"--whitelists", os.path.join(talon_analysis, "filtered_isoforms.csv"),  # Filtered transcripts to be reported
			"--datasets", sample_name,  # Input of the filtered tables produced on the previous step
			"--outdir", output_dir,  # Output dir
			"2>>", os.path.join(pipeline_reports, "talon9_generate_report-report.txt")])  # Directory where all reports reside
			subprocess.run(talon_report, shell=True)

		### Tenth step: Obtaining the TALON database in GTF format
		print("{0}  10/ TALON - Obtaining the transcriptome annotation from the TALON database: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		talon_export_db = " ".join([
		"talon_create_GTF",  # Call talon_create_GTF
		"--db", talon_database,  # TALON database
		"--build hg38",  # Genome build (hg38) to use
		"--annot", talon_database[:-3],  # Which annotation version to use
		"--whitelist", os.path.join(talon_analysis, "filtered_isoforms.csv"),  # Filtered transcripts to be reported
		"--o", os.path.join(talon_analysis, "database"),  # Output
		"2>>", os.path.join(pipeline_reports, "talon10_export_db-report.txt")])  # Directory where all reports reside
		subprocess.run(talon_export_db, shell=True)
		
		
		############ DIFFERENTIAL ANALYSIS ############

		### Eleventh step: Exploratory analysis
		print("{0}  11/ TALON - Exploratory analysis: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		expl_analysis = " ".join([
		"Rscript",  # Call Rscript
		"{0}/diffExpr_ExplAnalysis.R".format(rscripts),  # Calling the diffExpr_ExplAnalysis.R script
		os.path.join(talon_analysis, "prefilt_talon_abundance.tsv"),  # Input filtered matrix
		os.path.join(talon_analysis, "prefilt_talon_read_annot.tsv"),  # Read annotation matrix
		os.path.join(talon_analysis, "talon_input.csv"),  # Input annotation matrix
		R_analysis,  # Output directory
		args.minGeneExpr,  # minGeneExpr - Minimum number of reads for a gene to be considered expressed
		args.n_top,  # Top n_top genes for creating the heatmap
		"2>>", os.path.join(pipeline_reports, "talon11_exploratory_analysis-report.txt")])  # Directory where all reports reside
		subprocess.run(expl_analysis, shell=True)
		
		# ### Twelfth step: DGE
		# print("{0}  12/ TALON - Differential Gene Expression (DGE) analysis using DRIMSeq/edgeR: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		# dge_analysis = " ".join([
		# "Rscript",  # Call Rscript
		# "{0}/diffExpr_DGE.R".format(rscripts),  # Calling the diffExpr_DGE.R script
		# os.path.join(talon_analysis, "prefilt_talon_abundance.tsv"),  # Input filtered matrix
		# os.path.join(talon_analysis, "talon_input.csv"),  # Input annotation matrix
		# R_analysis,  # Output directory
		# args.minSampsGeneExpr,  # minSampsGeneExpr - Genes expressed in minimum this many samples
		# args.minGeneExpr,  # minGeneExpr - Minimum number of reads for a gene to be considered expressed
		# args.adjPValueThreshold,  # adjPValueThreshold - Adjusted p-value threshold for differential expression
		# args.lfcThreshold,  # lfcThreshold - Minimum required log2 fold change for differential expression
		# "2>>", os.path.join(pipeline_reports, "talon12_dge_analysis-report.txt")])  # Directory where all reports reside
		# subprocess.run(dge_analysis, shell=True)
		
		# ### Thirteenth step: DTE
		# print("{0}  13/ TALON - Differential Transcript Expression (DTE) analysis using DRIMSeq/edgeR: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		# dte_analysis = " ".join([
		# "Rscript",  # Call Rscript
		# "{0}/diffExpr_DTE.R".format(rscripts),  # Calling the diffExpr_DTE.R script
		# os.path.join(talon_analysis, "filt_talon_abundance.csv"),  # Input filtered matrix
		# os.path.join(talon_analysis, "prefilt_talon_read_annot.tsv"),  # Read annotation matrix
		# os.path.join(talon_analysis, "talon_input.csv"),  # Input annotation matrix
		# R_analysis,  # Output directory
		# args.minSampsFeatureExpr,  # minSampsFeatureExpr -  A transcript must be mapped to an isoform at least this minimum number of samples for the gene isoform to be considered
		# args.minFeatureExpr,  # minFeatureExpr - Minimum number of reads for a gene isoform to be considered
		# args.adjPValueThreshold,  # adjPValueThreshold - Adjusted p-value threshold for differential expression
		# args.lfcThreshold,  # lfcThreshold - Minimum required log2 fold change for differential expression		
		# "2>>", os.path.join(pipeline_reports, "talon13_dte_analysis-report.txt")])  # Directory where all reports reside
		# subprocess.run(dte_analysis, shell=True)

		# ### Eleventh step: DTU
		# print("{0}  14/ TALON - Differential Transcript Usage (DTU) analysis using DRIMSeq/DEXSeq/stageR: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		# dtu_analysis = " ".join([
		# "Rscript",  # Call Rscript
		# "{0}/diffExpr_DTU.R".format(rscripts),  # Calling the diffExpr_DTU.R script
		# os.path.join(talon_analysis, "filt_talon_abundance.csv"),  # Input filtered matrix
		# os.path.join(talon_analysis, "talon_input.csv"),  # Input annotation matrix
		# R_analysis,  # Output directory
		# args.minSampsGeneExpr,  # minSampsGeneExpr - Genes expressed in minimum this many samples
		# args.minSampsFeatureExprDTU,  # minSampsFeatureExpr -  A transcript must be mapped to an isoform at least this minimum number of samples for the gene isoform to be considered
		# args.minGeneExpr,  # minGeneExpr - Minimum number of reads for a gene to be considered expressed
		# args.minFeatureExprDTU,  # minFeatureExpr - Minimum number of reads for a gene isoform to be considered
		# args.adjPValueThreshold,  # adjPValueThreshold - Adjusted p-value threshold for differential expression
		# args.lfcThreshold,  # lfcThreshold - Minimum required log2 fold change for differential expression
		# args.threads,  # Num of threads to use
		# "2>>", os.path.join(pipeline_reports, "talon14_dtu_analysis-report.txt")])  # Directory where all reports reside
		# subprocess.run(dtu_analysis, shell=True)


		############ ANNOTATION VISUALISATION ############

		annot = {}
		with open(os.path.join(talon_analysis, "database_talon.gtf")) as ref_in:
			for i, line in enumerate(ref_in):
				if not line.startswith("#"):
					if line.split("\t")[2].strip() == 'transcript':
						transcript_id = line.split("\t")[-1].split(";")[1].split(" ")[-1].strip("\"")
						transcript_type = ["novel", line.split("\t")[-1].split(";")[4].split()[1].strip("\"").strip()][transcript_id.startswith("ENST")]
						annot[transcript_id] = transcript_type
						
		data = {}
		header = []
		nc = ["miRNA", "piRNA", "rRNA", "siRNA", "snRNA", "snoRNA", "tRNA", "vaultRNA"]
		with open("{0}/filt_talon_abundance.csv".format(talon_analysis)) as mat_in:
			for i, line in enumerate(mat_in, 1):
				if i == 1:
					header = line.strip().split(",")[9:]
				else:
					transcript = line.strip().split(",")[1]
					values = line.strip().split(",")[9:]
					if sum(map(int, values)) >= 10:
						# Grouping all pseudogenes
						if annot[transcript].endswith("pseudogene"):
							data[(transcript, "pseudogene")] = values
						# Grouping all immunoglobin genes
						elif annot[transcript].startswith("IG_"):
							data[(transcript, "Immunoglobulin_gene")] = values
						# Grouping all T-cell receptor genes
						elif annot[transcript].startswith("TR_"):
							data[(transcript, "Tcell_receptor_gene")] = values
						# Reanming TEC group
						elif annot[transcript].startswith("TEC"):
							data[(transcript, "To_be_Experimentally_Confirmed")] = values
						# Rest
						else:
							data[(transcript, annot[transcript])] = values

		# Remove paths from sample names
		header.insert(0,"transcript_id")  # Inserting ttranscript_id in header
		header.insert(1, "transcript_type")  # Inserting transcript_type in header
		
		# Writing output to file 'expression_matrix.csv'
		with open("{0}/perTranscript_expression_matrix.csv".format(talon_analysis), "w") as fout:
			fout.write("{0}\n".format(','.join(header)))
			for key, values in data.items():
				fout.write("{0},{1}\n".format(','.join(key), ','.join(values)))

		print("{0}  Visualising the RNA categories found in the TALON expression matrix: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		gene_type_sum = " ".join([
		"Rscript",  # Call Rscript
		"{0}/transcript_type_summary.R".format(rscripts),  # Calling the transcript_type_summary.R script
		"{0}/perTranscript_expression_matrix.csv".format(talon_analysis),  # Input matrix
		R_analysis,  # Output dir
		"2>>", os.path.join(pipeline_reports, "R_transcriptType_sum-report.txt")])  # Directory where all reports reside
		subprocess.run(gene_type_sum, shell=True)

		### Removing unnecessary directories and files
		os.system("rm -r {0}".format(temp))  
		os.system("rm -r {0}".format("talon_tmp"))
		return

	def novel_transcripts_detection_flair(self, chosen_samples):
		print("\n\t{0} ANNOTATION AND QUANTIFICATION USING FLAIR".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		flair_analysis = os.path.join(expression_analysis_dir, 'flair_analysis')
		if not os.path.exists(flair_analysis): os.makedirs(flair_analysis)

		genome_alignments = [sum_file for sum_file in glob.glob(os.path.join(alignments_dir, "*.genome.bam"))]
		print("{0}  1/ FLAIR - Converting BAM files to bed12: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		for file in genome_alignments:
			flair_convert = " ".join([
			"python3 /home/stavros/playground/progs/flair/bin/bam2Bed12.py",  # Calling bam2Bed12.py
			"--input_bam", file,  # Input BAM file
			">", os.path.join(flair_analysis, "{0}.bed12".format(os.path.basename(file).split(".")[0])),  # output file
			"2>>", os.path.join(pipeline_reports, "flair1_convert-report.txt")]) 
			subprocess.run(flair_convert, shell=True)

		genome_alignments_bed12 = [sum_file for sum_file in glob.glob(os.path.join(flair_analysis, "*.bed12"))]
		print("{0}  2/ FLAIR correct - Correcting misaligned splice sites using genome annotations: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		for file in genome_alignments_bed12:
			flair_correct = " ".join([
			"python3 /home/stavros/playground/progs/flair/flair.py correct",  # Calling flair correct
			"--gtf", refAnnot,  # GTF ref. annotation file
			"--threads", args.threads,  # Number of cores to use
			"--genome", refGenomeGRCh38,  # Input ref. genome
			"--output", os.path.join(flair_analysis,'{0}_flair_correct'.format(os.path.basename(file)[:-6])),  # Output name base
			"--query", file,  # Uncorrected bed12 file
			"--chromsizes", chrom_size,  # Chromosome sizes tab-separated file
			"2>>", os.path.join(pipeline_reports, "flair2_correct-report.txt")]) 
			subprocess.run(flair_correct, shell=True)

		print("{0}  3/ FLAIR - combining all samples together by concatenating corrected read psl files together.: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		psl_files = [psl for psl in glob.glob(os.path.join(flair_analysis, "*corrected.psl"))]
		merged_psl = os.path.join(flair_analysis, "merged.psl")
		subprocess.run("cat {0} > {1}".format(' '.join(psl_files), merged_psl), shell=True)

		print("{0}  4/ FLAIR collapse - Defining high-confidence isoforms from corrected reads: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		fastq_files = []
		flair_csv_file = os.path.join(flair_analysis, 'flair_input.csv')  # Create config file that is needed for FLAIR quantify
		summary_files = [str(file_path) for file_path in Path(ont_data).glob('**/batch3*/**/sequencing_summary.txt')]
		with open(flair_csv_file, "w") as flair_out:
			for sum_file in [s for s in sorted(summary_files) if os.path.dirname(s).endswith(chosen_samples)]:
				fastq_file = glob.glob(os.path.join(os.path.dirname(str(sum_file)), "pass/*.fastq.gz"))[0]
				sample_id = os.path.basename(os.path.dirname(str(sum_file)))
				flair_out.write('{0}\t{1}\tMINIONB3\t{2}\n'.format(sample_id.replace("_",""), sample_id.split("_")[0], fastq_file))
				fastq_files.append(fastq_file)
		flair_collapse = " ".join([
		"python3 /home/stavros/playground/progs/flair/flair.py collapse",  # Calling flair collapse
		"--threads", args.threads,  # Number of cores to use
		"--gtf", refAnnot,  # GTF ref. annotation file
		"--genome", refGenomeGRCh38,  # FastA of reference genome
		"--query", merged_psl,  # Corrected psl file
		"--reads", " ".join(fastq_files),
		"--output", os.path.join(flair_analysis, "flair.collapse"),  # output file
		"2>>", os.path.join(pipeline_reports, "flair3_collapse-report.txt")]) 
		subprocess.run(flair_collapse, shell=True)

		print("{0}  5/ FLAIR quantify - Quantification of FLAIR isoform usage across samples: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		collapsed_isoforms = os.path.join(flair_analysis, "flair.collapse.isoforms.fa")
		flair_quantify = " ".join([
		"python3 /home/stavros/playground/progs/flair/flair.py quantify",  # Calling flair quantify
		"--threads", args.threads,  # Number of cores to use
		"--reads_manifest", flair_csv_file,  # Tab delimited file containing: sample id, condition, batch, reads.fq
		"--isoforms", collapsed_isoforms,  # Fasta input from FLAIR collapsed isoforms
		"--tpm",  # TPM column
		"--quality 10",  # Minimum MAPQ of read assignment to an isoform
		"--output", os.path.join(flair_analysis, "flair_expression_matrix.tsv"),  # output file
		"2>>", os.path.join(pipeline_reports, "flair_quantify-report.txt")]) 
		subprocess.run(flair_quantify, shell=True)

		# print("{0}  6/ FLAIR DE - Differential Expression analysis using DRIMSeq/DESeq2: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		# flair_de = " ".join([
		# # "python3 /home/stavros/playground/progs/flair/flair.py diffExp",  # Calling flair diffExp
		# "/home/stavros/playground/progs/anaconda3/envs/flair_env/bin/python /home/stavros/playground/progs/flair/flair.py diffExp",  # Calling flair diffExp
		# "--threads", args.threads,  # Number of cores to use
		# "--counts_matrix", os.path.join(flair_analysis, "flair_expression_matrix.tsv"),  # Tab-delimited isoform count matrix from flair quantify module
		# "--out_dir", os.path.join(flair_analysis, "diffExp"),  # output file
		# "--out_dir_force",
		# "2>>", os.path.join(pipeline_reports, "flair_de-report.txt")]) 
		# subprocess.run(flair_de, shell=True)

		# print("{0}  7/ FLAIR DS - Differential Splicing analysis using DRIMSeq: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		# flair_ds = " ".join([
		# "python3 /home/stavros/playground/progs/flair/flair.py diffSplice",  # Calling flair diffSplice
		# "--threads", args.threads,  # Number of cores to use
		# "--counts_matrix", os.path.join(flair_analysis, "flair_expression_matrix.tsv"),  # Tab-delimited isoform count matrix from flair quantify module
		# "--isoforms", os.path.join(flair_analysis,"flair.collapse.isoforms.psl"),  # Isoforms in bed format
		# "--output", os.path.join(flair_analysis, "flair.diffsplice"),  # output file
		# "2>>", os.path.join(pipeline_reports, "flair_ds-report.txt")]) 
		# subprocess.run(flair_ds, shell=True)

		# # print("{0}  8/ FLAIR predictProductivity - Annotated start codons to identify the longest ORF for each isoform for predicting isoform productivity: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		# # flair_predictProductivity = " ".join([
		# # "python2 /home/stavros/playground/progs/flair/bin/predictProductivity.py",  # Calling flair predictProductivity.py
		# # "--genome_fasta", refGenomeGRCh38,  # Fasta file containing transcript sequences.
		# # "--longestORF",  # Defined ORFs by the longest open reading frame
		# # "--gtf", refAnnot,  # Gencode annotation file
		# # # "--input_isoforms", ,  #
		# # ">", os.path.join(flair_analysis, "productivity.bed"),  # output file
		# # "2>>", os.path.join(pipeline_reports, "flair_predictProductivity-report.txt")]) 
		# # subprocess.run(flair_predictProductivity, shell=True)

		# # print("{0}  7/ FLAIR Visualisation - isoform structures and the percent usage of each isoform in each sample for a given gene: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		# # flair_ds = " ".join([
		# # "python3 /home/stavros/playground/progs/flair/bin/plot_isoform_usage.py",  # Calling plot_isoform_usage.py
		# # ### isoforms.bed file from running predictProductivity.py
		# # "--counts_matrix", os.path.join(flair_analysis, "flair_expression_matrix.tsv"),  # Tab-delimited isoform count matrix from flair quantify module
		# # ### GENENAME,
		# # os.path.join(flair_analysis, "flair.isoform_usage."),  # output file
		# # "2>>", os.path.join(pipeline_reports, "flair_ds-report.txt")]) 
		# # subprocess.run(flair_ds, shell=True)
		return

class special_analysis:

	def __init__(self):
		self.polyA_length_est_analysis()
		# self.methylation_detection()
		return

	def polyA_length_est_analysis(self):
		print("\n\t{0} POLYA LENGTH ESTIMATION ANALYSIS".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		# TALON read annotation file and saving it as a dictionary 
		read_annot_dict = {}
		talon_read_annot = os.path.join(expression_analysis_dir,"talon_analysis/prefilt_talon_read_annot.tsv")
		with open(talon_read_annot) as read_annot:
			for line in read_annot:
				read = line.strip().split("\t")[0]
				sample = line.strip().split("\t")[1]
				transcript_id = line.strip().split("\t")[12]
				gene_id = line.strip().split("\t")[11]
				read_annot_dict[(sample, read)] = (transcript_id,gene_id)

		for path, subdir, folder in os.walk(analysis_dir):
			for name in folder:
				if name.endswith("polya_results.tsv"):
					sample = name.split(".")[0]
					polyA_length_est = os.path.join(path, name)
					with open(polyA_length_est) as fin, \
						 open(polyA_length_est.replace(".tsv",".transcripts.tsv"), "w") as transcript_out,\
						 open(polyA_length_est.replace(".tsv",".genes.tsv"), "w") as gene_out:
						for i, line in enumerate(fin, 1):
							if i == 1:
								transcript_out.write("{0}\n".format(line.strip()))
								gene_out.write("{0}\n".format(line.strip()))
							else:
								readname = line.strip().split("\t")[0]
								contig = "NNNN0000{1}".format(line.strip().split("\t")[1], i)
								transcript_id = read_annot_dict[(sample ,readname)][0] if (sample ,readname) in read_annot_dict else contig
								gene_id = read_annot_dict[(sample ,readname)][1] if (sample ,readname) in read_annot_dict else contig
								rest = "\t".join(line.strip().split("\t")[2:])
								transcript_out.write("{0}\t{1}\t{2}\n".format(readname, transcript_id, rest))
								gene_out.write("{0}\t{1}\t{2}\n".format(readname, gene_id, rest))
					
					# Writing basic info to 'polyA_data_info' for NanoTail analysis in  transcript level
					sample_info_transcripts = "{0}/polyA_transcript_info.csv".format(polyA_analysis_dir)
					with open(sample_info_transcripts, "a") as fout_tr:
						fout_tr.write("{0},{1},{2}/{0}/{0}.polya_results.transcripts.tsv\n".format(sample, sample.split("_")[0], polyA_analysis_dir))
					# nanotail_analysis(sample_info_transcripts, "transcript")  # Transcript level analysis
					
					# # Writing basic info to 'polyA_data_info' for NanoTail analysis in  transcript level
					# sample_info_genes = "{0}/polyA_gene_info.csv".format(polyA_analysis_dir)
					# with open(sample_info_genes, "a") as fout_gene:
					# 	fout_gene.write("{0},{1},{2}/{0}/{0}.polya_results.genes.tsv\n".format(sample, sample.split("_")[0], polyA_analysis_dir))
					# # nanotail_analysis(sample_info_genes, "gene")  # Gene  level analysis
		return

	def nanotail_analysis(self, sample_info, what):
		### Running Nanotail analysis
		print("{0}  NanoTail - Exploratory analysis: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		nanotail_analysis = " ".join([
		"Rscript",  # Call Rscript
		"{0}/polyA_analysis.R".format(rscripts),  # Calling the polyA_analysis.R script
		sample_info,  # Input annotation
		os.path.join(polyA_analysis_dir, "{0}_analysis".format(what)),  # Output directory
		"2>>", os.path.join(pipeline_reports, "nanotail_analysis-report.txt")])  # Directory where all reports reside
		subprocess.run(nanotail_analysis, shell=True)
		return

	def methylation_detection(self):
		if not os.path.exists(methylation_dir): os.makedirs(methylation_dir)
		single_fast5_data = [os.path.dirname(sum_file) for sum_file in glob.glob("/shared/projects/silvia_rna_ont_umc/basecalling/sample*_raw_data/*/")]

		# print(single_fast5_data)
		for dirs in single_fast5_data:
			if dirs.endswith("/0"):

				sample_id = dirs.split("/")[5].split("_")[0]
				# 1. Re-squiggling the raw reads
				tombo_resquiggle = " ".join([
				"tombo resquiggle",
				dirs,
				refGenomeGRCh38,
				"--processes", args.threads,  # Number of threads to be used
				"--num-most-common-errors 5",
				# "2>>", os.path.join(pipeline_reports, "tombo_resquiggle-report.txt")
				])
				subprocess.run(tombo_resquiggle, shell=True)
			
			# 2. Calling tombo to do the methylation analysis
			tombo_methyl = " ".join([
			"tombo detect_modifications de_novo",  # Indexing the concat_samples.bam file
			"--fast5-basedirs", dirs,  # Directory containing fast5 files
			"--statistics-file-basename", os.path.join(methylation_dir,"{0}.tombo.stats".format(sample_id)),
			"--rna",  # Explicitly select canonical RNA mode
			"--processes", args.threads,  # Number of threads to be used
			# "2>>", os.path.join(pipeline_reports, "tombo_methylation-report.txt")
			])
			subprocess.run(tombo_methyl, shell=True)

			# 3. Output reference sequence around most significantly modified sites
			tombo_sign = " ".join([
			"tombo text_output signif_sequence_context",  # Indexing the concat_samples.bam file
			"--fast5-basedirs", dirs,  # Directory containing fast5 files
			"--statistics-filename", os.path.join(methylation_dir,"{0}.tombo.stats".format(sample_id)),
			"--sequences-filename",  os.path.join(methylation_dir,"{0}.tombo_significant_regions.fasta".format(sample_id)),
			# "2>>", os.path.join(pipeline_reports, "tombo_significants-report.txt")
			])
			subprocess.run(tombo_sign, shell=True)

			# # 4. Use line Meme to estimate modified motifs
			# tombo_meme = " ".join([
			# "meme",
	  #  		"-rna",
	  #  		"-mod zoops", 
			# "-oc", os.path.join(methylation_dir,"{0}_de_novo_meme".format(sample_id)),
	  #  		os.path.join(methylation_dir,"{0}.tombo_significant_regions.fasta".format(sample_id)),
	  #  		# "2>>", os.path.join(pipeline_reports, "tombo_meme-report.txt")
			# ])
			# subprocess.run(tombo_meme, shell=True)

			# # 5. This plot will identify the sites in the reference
			# tombo_plot = " ".join([
			# "tombo plot motif_with_stats",
			# " --fast5-basedirs", dirs,
	  #  		"--motif CCWGG",
	  #  		"--genome-fasta", refGenomeGRCh38,
	  #  		"--statistics-filename", os.path.join(methylation_dir,"{0}.tombo.stats".format(sample_id)),
	  #  		"--pdf-filename", os.path.join(methylation_dir,"{0}.tombo.pdf".format(sample_id))
	  #  		# "2>>", os.path.join(pipeline_reports, "tombo_plot-report.txt")
			# ])
			# subprocess.run(tombo_plot, shell=True)
		return

def summary():
	
	multiQC = " ".join([
	"multiqc",  # Call MultiQC
	"--quiet",  # Print only log warnings
	"--outdir", postAlignment_reports,  # Create report in the FastQC reports directory
	"--filename", "post-alignment_summarised_report",  # Name of the output report 
	postAlignment_reports,  # Directory where all FastQC and Cutadapt reports reside
	"2>>", os.path.join(pipeline_reports, "post-alignment_multiQC-report.txt")])  # Output multiQC report
	subprocess.run(multiQC, shell=True)

	pickle_files = [pk_file for pk_file in glob.glob(os.path.join(postAlignment_reports, "*_qc.pk"))]
	### Wub Compare alignment QC statistics of multiple samples
	bam_multi_qc = ' '.join([
	"bam_multi_qc.py",
	"-r", "{0}/comparison_qc.pdf".format(postAlignment_reports),  # Output pdf file
	' '.join(pickle_files),
	"2>>", os.path.join(pipeline_reports, "bam_multi_qc-report.txt")])
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
	return


def main():
	
	chosen_samples = ("NonTransf_1",  "NonTransf_2",  "NonTransf_3", "Tumour_1",  "Tumour_2",  "Tumour_4")
	# chosen_samples = ("Tumour_1")


	summary_files = [str(file_path) for file_path in Path(ont_data).glob('**/sequencing_summary.txt') if not "warehouse" in str(file_path)]
	num_of_samples = len(summary_files)

	for sum_file in [s for s in summary_files if os.path.dirname(s).endswith(chosen_samples)]:
		raw_data_dir = os.path.dirname(str(sum_file))
		sample_id = os.path.basename(raw_data_dir)
		fastq_pass = " ".join(glob.glob(os.path.join(raw_data_dir, "pass/*pass.fastq.gz")))


		quality_control(sum_file, sample_id, raw_data_dir)

		alignment_against_ref(fastq_pass, sample_id, raw_data_dir, sum_file)
			
		# polyA_estimation(sample_id, sum_file, fastq_pass, raw_data_dir)

	expression_matrix(chosen_samples)

	# special_analysis()

	summary()

	print('\t--- The pipeline finisded after {0} ---'.format(datetime.now() - startTime))

if __name__ == "__main__": main()