###############
REQUIREMENTS
##############

1. Python version 3.8+
2. Bowtie2 version 2.4+ (https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/)
3. primer3-py version 0.6.1 (https://pypi.org/project/primer3-py/)

####################
DOWNLOADS WHEN RUN
####################

When COD_FISH.py program is run, unless they are already present, the following files will be downloaded and organized in a folder titled species.(chosen species):
These files are then combined into one large filed named 'transcriptome.fa'

        1) cDNA.fa (from http://ftp.ensembl.org/pub/current_fasta/)
        2) ncRNA.fa (from http://ftp.ensembl.org/pub/current_fasta/)
        3) filtered_rRNA_LSU_DNA.fa (filtered for species from https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz)
        4) filtered_rRNA_SSU_DNA.fa (filtered for species from https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz)
        5) transcriptome.fa

Then bowtie2 builds an index using the transcriptome.fa file which is then used to align candidate probe sequences to.

###################
CONFIGURATION FILE
##################

CODFISH generates a configuration file in which various settings for the probes, such as the melting temperature and the method of probe scoring are speciied. Once generated, this file can be modified to control the probe generation.
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
	
	a) Tm - In this method, the Primer3 package is used to estimate the melting temperature of the probe to the competing RNA, and then from this an specifity score is calculated
	
	b) Alignment - In this method, the "AS:i:" section of the SAM file created from the bowtie2 alignment is used as the score for the heterodimer of the probe to the competing RNA, and then from this a specifity score is calculated
	
10. Probe Set Selection Method (Probe_set_selection_method) - The method by which the probe set is selected from a list of probes and their corresponding specifity scores. There are two possible probe set selection methods (default is Greedy Method):
	
	a) Dynamic - Modeled on Dijkstra's algorithm, this method utilizes a token-passing system and dynamic programming to determine the optimal set of probes
	
	b) Greedy - Chooses the highest scoring probe, and then goes through the list to find the second highest-scoring probe, and so on and so forth, while taking the overlap into account.
	
11. Temp Directiory (temp_dir) - The location and name of the temp directory (default is ./temp/)
12. Output Directory (output_dir) - The location and name of the output directory (default is ./output_finalProbes/)

######################
RUNNING THE PROGRAM
#####################

In order to create probes from a given mRNA sequence, or from a specific gene name, input the sequence as either gene symbols (--gene_symbols/-g), a file with a list of gene symbols (--gene_symbol_file/-gf), a fasta file with raw sequences (--transcripts_file/-fa), or ensembl transcript IDs (ie ENSTXXXXXX) (--ensembl_ids/-ids).

Here is an example when making probes for human PTGS2

        $ python COD_FISH.py -g PTGS2 --species human

Multiple sequences from multiple different formats can be run at the same time. (The species argument (--species/-s) isn't needed once the config.py file and species directory are made)

        $ python COD_FISH.py -g PTGS2 TP53 -ids ENST00000261783 

Additional arguements can be found by running

        $ python COD_FISH.py -h

and are shown below:

  --transcripts_file/-fa        Full path of Fasta file name with target mRNAs, 
                                        write down all target mRNAs in a fasta file
  --ensembl_ids/-ids            List of ensembl gene IDs [List of ensembl gene IDs ...]
                                        list any number of ensembl gene IDs separated by a space (case sensitive!)
  --gene_symbols/-g             List of gene symbols [List of gene symbols ...]
                                        list any number of gene names separated by a space (case sensitive!)
  --gene_symbol_file/-gf        File of gene symbols
                                        file containing a list of genes to process, each gene symbol should be on a separate line
  --species/-s Species name
                                Write down the species you want to use
  --advanced/-a                 Advanced boolean
                                        If specified will let user chose config file criteria
  --bowtie2/-bt2                Bowtie2 path if not in $PATH,  Bowtie2 path if not in $PATH
                                        If Bowtie2 is not a PATH variable then you can specify where the installed bowtie2 is here
  --seed_length/-L Seed Length
                                Length of the seed length used on the probe sequences when aligning to the reanscriptome, lower is more specific but
                                        with diminishing returns, default is 15
