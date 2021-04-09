### Calling IUPred2A for the downstream analysis of IsoformSwitchAnalyzeR
import os, sys
import subprocess



multifasta_temp_dir = f'{os.path.dirname(sys.argv[1])}/isoformSwitchAnalyzeR_isoform_AA_split_files'

# Splitting the fasta file into multiple fastas containing each sequence
subprocess.run(f'splitfasta {sys.argv[1]}',shell=True)
subprocess.run(f'mv {os.path.dirname(sys.argv[1])}/isoformSwitchAnalyzeR_isoform_AA_*.fasta {multifasta_temp_dir}',shell=True)

for file in os.listdir(multifasta_temp_dir):
	sample_id = os.path.basename(file).split(".")[0]
	run_iupred2a = ' '.join([
	"iupred2a.py",  # Calling iupred2a
	"-a",  # Enable ANCHOR2 prediction
	os.path.join(multifasta_temp_dir, file),  # Input fasta file
	"long",  # IUPred2 long disorder
	">", f'{multifasta_temp_dir}/{sample_id}.txt'])  # Output file
	subprocess.run(run_iupred2a, shell=True)


iupred2a_final = os.path.join(os.path.dirname(sys.argv[1]), 'results_IUPred2A.txt')
with open(iupred2a_final, 'w') as fout:
	for file in os.listdir(multifasta_temp_dir):
		if file.endswith(".txt"):
			seq_id = 0
			fasta_file = os.path.join(multifasta_temp_dir, file.replace('.txt', '.fasta'))
			iupred2a_file = os.path.join(multifasta_temp_dir, file)
			### Obtaining the GENE_ID
			with open(fasta_file) as fin: 
				for line in fin:
					if line.startswith(">"):
						seq_id = line.strip()
			if int(subprocess.check_output(['wc', '-l', iupred2a_file]).split()[0]) > 7:
				fout.write('# IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding\n')
				fout.write('# Balint Meszaros, Gabor Erdos, Zsuzsanna Dosztanyi\n')
				fout.write('# Nucleic Acids Research 2018, Submitted\n')
				fout.write(f'{seq_id}\n')
				fout.write('# POS	AMINO ACID	IUPRED SCORE	ANCHOR SCORE\n')
				with open(iupred2a_file) as res:
					for line in res:
						if not line.startswith("#"):
							fout.write(line)
				fout.write('\n\n################\n')

subprocess.run(f'rm -r {multifasta_temp_dir}', shell=True)