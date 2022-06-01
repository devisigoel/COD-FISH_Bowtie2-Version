# This is the configuration file which defines the various parameters for making probes

time_made = 'Fri Feb 11 19:45:29 2022'
species = 'homo_sapiens'
num_probes = 50 #Number of probes to be created
probe_len = 20 #Length of each probe
Tm = 313 #Minimum melting temperature of probes in Kelvin
Tm_hairpin = 330 #Maximum melting temperature of the hairpin of probes in Kelvin

#Filter repeating nucleotides, AAAAA, GGGG, CCCC, TTTTT
filter_repeats = True
# Complete Human Transcriptome fasta file for determining probe specificity
transcriptome_file = "species.homo_sapiens/transcriptome.fa"
transcriptome_index = "species.homo_sapiens/bt2_index"

#The Heterodimer melting temperature estimation method and the probe set selection methods
Heterodimer_melting_temperature_estimation_method = 'Primer3'
Probe_set_selection_method = 'Dynamic Programming'
#Probe_set_selection_method = 'Greedy Method'

# This is a directory in which the program will create different files, but these files can be deleted 
# once the program has finished running. 
temp_dir = "./temp/"

# This is a directory in which the program will create different files with the output of the program
output_dir = "./output_finalProbes/"
