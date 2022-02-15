###############
REQUIREMENTS
##############

1. Python version 3.8+
2. Bowtie2 version 2.4+ (https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/)
3. Primer3 version 0.4.0 (https://sourceforge.net/projects/primer3/)

####################
DOWNLOADS WHEN RUN
####################

When COD_FISH.py program is run, unless they are already present, the following files will be downloaded and organized in a folder titled species.(chosen species):

        1) cDNA.all.fa.gz (source: http://ftp.ensembl.org/pub/current_fasta/+(chosen species)+/cdna/)
        2) ncRNA.fa.gz (source: http://ftp.ensembl.org/pub/current_fasta/+(chosen species)+/ncrna/)
        3) rRNA.SSU.fa,gz (source: https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz)
        4) filtered_rRNA_SSU.fa (same source as #3)
        5) filtered_rRNA_LSU_DNA.fa (same source as #3)

These files are then combined into one large filed named 'transcriptome.fa'

###################
CONFIGURATION FILE
##################

CODFISH comes with a configuration file in which various settings for the probes, such as the melting temperature and the method of probe scoring are speciied.
Here are the variables defined by the configuration file:

1. Species - the name of the species, as specied by the user
2. Numeber of probes (num_probes) - the specified number of probes in a probeset (default is 50)
3. Probe Length (probe_len) - The length of the probes (default is 20)
4. Melting Temperature (Tm) - The minimum melting temperature threshold of the probes (default is 313 Kelvin)
5. Melting Temperature of the hairpin (Tm_hairpin) - The maximum melting temperature threshold of the probe's hairpin (default is 330 Kelvin)
6. Nucleotide Repetition Filter (filter_repeats) - The filtration of probes based on the repetition of nucleotides (default is True)
7. Transcriptome Index (transcriptome_index) - The location and name of the transcriptome index created for Bowtie2 alignment (default is species.(chosen species)/bt2_index)
8. Transcriptome File (transcriptome_file) - The location and name of the transcriptome file created for Bowtie2 alignment (default is species.(chosen species)/transcriptome.fa)
9. Heterodimer Melting Temperature Estimation Method (Heterodimer_melting_temperature_estimation_method) - The method by which the heterodimer melting temeprature estimation can occur. There are two possible methods (default is Primer3):
        a) Primer3 - In this method, the Primer3 module is used to estimate the melting temperature of the probe to the competing RNA, and then from this an specifity score is calculated
        b) Alignment - In this method, the "AS:i:" section of the SAM file created from the bowtie2 alignment is used as the score for the heterodimer of the probe to the competing RNA, and then from this a specifity score is calculated
10. Probe Set Selection Method (Probe_set_selection_method) - The method by which the probe set is selected from a list of probes and their corresponding specifity scores. There are two possible probe set selection methods (default is Greedy Method):
        a) Dynamic Programming - Modeled on Dijkstra's algorithm, this method utilizes a token-passing system and dynamic programming to determine the optimal set of probes
        b) Greedy Method - Chooses the highest scoring probe, and then goes through the list to find the second highest-scoring probe, and so on and so forth, while taking the overlap into account.
11. Temp Directiory (temp_dir) - The location and name of the temp directory (default is ./temp/)
12. Output Directory (output_dir) - The location and name of the output directory (default is ./output_finalProbes/)

######################
RUNNING THE PROGRAM
#####################

In order to create probes from a given mRNA sequence, or from a specific
gene name, use the COD_FISH.py  program.

      $ python COD_FISH.py --method of providing sequence gene symbols/fasta file --species species

Key:
  - "--method of providing sequence" would be replaced either by --mRNA_fasta_file, --gene_symbol_list or --gene_symbol_file
  - "gene symbols/fasta file" would be replaced either by gene symbol (e.g. IFNG), a file containing multiple gene symbols, or a fasta file containing the target  mRNA sequence.
  - "species" would be replaced by whatever species you would like to use (human, mouse, zebrafish, rat, pig, chimpanzee, chicken, cow, fruitfly, celegans)
README.txt

Optional Arguments:
  -a (Advanced flag, would allow user to change the configuration file)
  -bt2 (allows user to change the bowtie2 path)
  -L (seed length, allows the user to choose the length of the seed for the Bowtie2 alignment, default is 15) 
