#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 09:20:13 2021

@author: Ryan 'mentor' Dikdan, Devisi Goel

This is a script which takes gene names and generates smFISH probes for the gene.
The motivation is to find optimal oligos for targeting RNAs given full transcriptome data.
"""

import sys, os, re
import subprocess
import argparse
import time
import multiprocessing as mp

import blast_utils_COD_FISH as utils

# Changes working directory to where the script is so everything lines up well and the directory creation doesn't cause problems
os.chdir(os.path.dirname(os.path.abspath(__file__)))

##############################################################################
#  Input                                                                     #
##############################################################################
parser = argparse.ArgumentParser(description='Generating smFISH probes from given mRNA or gene symbols')
parser.add_argument('--transcripts_file','-fa', metavar='Full path of Fasta file name with target mRNAs', type=str, help='write down all target mRNAs in a fasta file')
parser.add_argument('--ensembl_gene_ids','-gid', metavar='List of ensembl gene IDs',type=str, help='list any number of stable ensembl gene IDs (no decimal) separated by a space (case sensitive!)', nargs='+')
parser.add_argument('--ensembl_transcript_ids','-tid', metavar='List of ensembl transcript IDs',type=str, help='list any number of stable ensembl transcript IDs (no decimal) separated by a space (case sensitive!)', nargs='+')
parser.add_argument('--gene_symbols','-g', metavar='List of gene symbols',type=str, help='list any number of gene names separated by a space (case sensitive!)', nargs='+')
parser.add_argument('--gene_symbol_file','-gf', metavar='File of gene symbols',type=str, help='file containing a list of genes to process, each gene symbol should be on a separate line')

parser.add_argument('--species','-s', metavar='Species name',type=str, help='write down the species you want to use')

parser.add_argument('-a','--advanced', metavar='Advanced bool',type=bool, help='If specified will let user chose config file criteria')
parser.add_argument('-bl', '--blastn', metavar='Bnastn path if not in $PATH',type=str, help='If Blastn is not a PATH variable then you can specify where the installed blastn is here')
parser.add_argument('-w', '--word_size', metavar='Word Size',type=int, help='Length of the word used on the probe sequences when aligning to the transcriptome, lower is more specific but with diminishing returns, default is 7 for blastn')

args = parser.parse_args()

##############################################################################
#  Setup Config file and indexes                                             #
##############################################################################

# Setting the flag for advanced option choices
advanced_flag = args.advanced

# Setting seed length

if args.word_size:
    word_size = args.word_size
else:
    word_size = 7

# Setting the number of threads
threads = mp.cpu_count()

# Setting blastn if path was given
if args.blastn:
    blastn_path = args.blastn
else:
    blastn_path = 'blastn'


# Species selection settings
species_brief = ['human', 'mouse', 'zebrafish', 'rat',
                     'pig', 'chimpanzee', 'chicken', 'cow', 'fruitfly', 'celegans']

ensembl_species = ['acanthochromis_polyacanthus', 'accipiter_nisus', 'ailuropoda_melanoleuca', 'amazona_collaria', 'amphilophus_citrinellus', 'amphiprion_ocellaris', 'amphiprion_percula', 'anabas_testudineus', 'anas_platyrhynchos', 'anas_platyrhynchos_platyrhynchos', 'anas_zonorhyncha', 'ancestral_alleles', 'anolis_carolinensis', 'anser_brachyrhynchus', 'anser_cygnoides', 'aotus_nancymaae', 'apteryx_haastii', 'apteryx_owenii', 'apteryx_rowi', 'aquila_chrysaetos_chrysaetos', 'astatotilapia_calliptera', 'astyanax_mexicanus', 'astyanax_mexicanus_pachon', 'athene_cunicularia', 'balaenoptera_musculus', 'betta_splendens', 'bison_bison_bison', 'bos_grunniens', 'bos_indicus_hybrid', 'bos_mutus', 'bos_taurus', 'bos_taurus_hybrid', 'bubo_bubo', 'buteo_japonicus', 'caenorhabditis_elegans', 'cairina_moschata_domestica', 'calidris_pugnax', 'calidris_pygmaea', 'callithrix_jacchus', 'callorhinchus_milii', 'camarhynchus_parvulus', 'camelus_dromedarius', 'canis_lupus_dingo', 'canis_lupus_familiaris', 'canis_lupus_familiarisbasenji', 'canis_lupus_familiarisgreatdane', 'capra_hircus', 'capra_hircus_blackbengal', 'carassius_auratus', 'carlito_syrichta', 'castor_canadensis', 'catagonus_wagneri', 'catharus_ustulatus', 'cavia_aperea', 'cavia_porcellus', 'cebus_capucinus', 'cercocebus_atys', 'cervus_hanglu_yarkandensis', 'chelonoidis_abingdonii', 'chelydra_serpentina', 'chinchilla_lanigera', 'chlorocebus_sabaeus', 'choloepus_hoffmanni', 'chrysemys_picta_bellii', 'chrysolophus_pictus', 'ciona_intestinalis', 'ciona_savignyi', 'clupea_harengus', 'colobus_angolensis_palliatus', 'corvus_moneduloides', 'cottoperca_gobio', 'coturnix_japonica', 'cricetulus_griseus_chok1gshd', 'cricetulus_griseus_crigri', 'cricetulus_griseus_picr', 'crocodylus_porosus', 'cyanistes_caeruleus', 'cyclopterus_lumpus', 'cynoglossus_semilaevis', 'cyprinodon_variegatus', 'cyprinus_carpio', 'cyprinus_carpio_germanmirror', 'cyprinus_carpio_hebaored', 'cyprinus_carpio_huanghe', 'danio_rerio', 'dasypus_novemcinctus', 'delphinapterus_leucas', 'denticeps_clupeoides', 'dicentrarchus_labrax', 'dipodomys_ordii', 'dromaius_novaehollandiae', 'drosophila_melanogaster', 'echeneis_naucrates', 'echinops_telfairi', 'electrophorus_electricus', 'eptatretus_burgeri', 'equus_asinus_asinus', 'equus_caballus', 'erinaceus_europaeus', 'erpetoichthys_calabaricus', 'erythrura_gouldiae', 'esox_lucius', 'falco_tinnunculus', 'felis_catus', 'ficedula_albicollis', 'fukomys_damarensis', 'fundulus_heteroclitus', 'gadus_morhua', 'gallus_gallus', 'gambusia_affinis', 'gasterosteus_aculeatus', 'geospiza_fortis', 'gopherus_agassizii', 'gopherus_evgoodei', 'gorilla_gorilla', 'gouania_willdenowi', 'haplochromis_burtoni', 'heterocephalus_glaber_female', 'heterocephalus_glaber_male', 'hippocampus_comes', 'homo_sapiens', 'hucho_hucho', 'ictalurus_punctatus', 'ictidomys_tridecemlineatus', 'jaculus_jaculus', 'junco_hyemalis', 'kryptolebias_marmoratus', 'labrus_bergylta', 'larimichthys_crocea', 'lates_calcarifer', 'laticauda_laticaudata', 'latimeria_chalumnae', 'lepidothrix_coronata', 'lepisosteus_oculatus', 'leptobrachium_leishanense', 'lonchura_striata_domestica', 'loxodonta_africana', 'lynx_canadensis', 'macaca_fascicularis', 'macaca_mulatta', 'macaca_nemestrina', 'malurus_cyaneus_samueli', 'manacus_vitellinus', 'mandrillus_leucophaeus', 'marmota_marmota_marmota', 'mastacembelus_armatus', 'maylandia_zebra', 'meleagris_gallopavo', 'melopsittacus_undulatus', 'meriones_unguiculatus', 'mesocricetus_auratus', 'microcebus_murinus', 'microtus_ochrogaster', 'mola_mola', 'monodelphis_domestica', 'monodon_monoceros', 'monopterus_albus', 'moschus_moschiferus', 'mus_caroli', 'mus_musculus', 'mus_musculus_129s1svimj', 'mus_musculus_aj', 'mus_musculus_akrj', 'mus_musculus_balbcj', 'mus_musculus_c3hhej', 'mus_musculus_c57bl6nj', 'mus_musculus_casteij', 'mus_musculus_cbaj', 'mus_musculus_dba2j', 'mus_musculus_fvbnj', 'mus_musculus_lpj', 'mus_musculus_nodshiltj', 'mus_musculus_nzohlltj', 'mus_musculus_pwkphj', 'mus_musculus_wsbeij', 'mus_pahari', 'mus_spicilegus', 'mus_spretus', 'mustela_putorius_furo', 'myotis_lucifugus', 'myripristis_murdjan', 'naja_naja', 'nannospalax_galili', 'neogobius_melanostomus', 'neolamprologus_brichardi', 'neovison_vison', 'nomascus_leucogenys', 'notamacropus_eugenii', 'notechis_scutatus', 'nothobranchius_furzeri', 'nothoprocta_perdicaria', 'numida_meleagris', 'ochotona_princeps', 'octodon_degus', 'oncorhynchus_kisutch', 'oncorhynchus_mykiss', 'oncorhynchus_tshawytscha', 'oreochromis_aureus', 'oreochromis_niloticus', 'ornithorhynchus_anatinus', 'oryctolagus_cuniculus', 'oryzias_javanicus', 'oryzias_latipes', 'oryzias_latipes_hni', 'oryzias_latipes_hsok', 'oryzias_melastigma', 'oryzias_sinensis', 'otolemur_garnettii', 'otus_sunia', 'ovis_aries', 'ovis_aries_rambouillet', 'pan_paniscus', 'pan_troglodytes', 'panthera_leo', 'panthera_pardus', 'panthera_tigris_altaica', 'papio_anubis', 'parambassis_ranga', 'paramormyrops_kingsleyae', 'parus_major', 'pavo_cristatus', 'pelodiscus_sinensis', 'pelusios_castaneus', 'periophthalmus_magnuspinnatus', 'peromyscus_maniculatus_bairdii', 'petromyzon_marinus', 'phascolarctos_cinereus', 'phasianus_colchicus', 'phocoena_sinus', 'physeter_catodon', 'piliocolobus_tephrosceles', 'podarcis_muralis', 'poecilia_formosa', 'poecilia_latipinna', 'poecilia_mexicana', 'poecilia_reticulata', 'pogona_vitticeps', 'pongo_abelii', 'procavia_capensis', 'prolemur_simus', 'propithecus_coquereli', 'pseudonaja_textilis', 'pteropus_vampyrus', 'pundamilia_nyererei', 'pygocentrus_nattereri', 'rattus_norvegicus', 'rhinolophus_ferrumequinum', 'rhinopithecus_bieti', 'rhinopithecus_roxellana', 'saccharomyces_cerevisiae', 'saimiri_boliviensis_boliviensis', 'salarias_fasciatus', 'salmo_salar', 'salmo_trutta', 'salvator_merianae', 'sander_lucioperca', 'sarcophilus_harrisii', 'sciurus_vulgaris', 'scleropages_formosus', 'scophthalmus_maximus', 'serinus_canaria', 'seriola_dumerili', 'seriola_lalandi_dorsalis', 'sinocyclocheilus_anshuiensis', 'sinocyclocheilus_grahami', 'sinocyclocheilus_rhinocerous', 'sorex_araneus', 'sparus_aurata', 'spermophilus_dauricus', 'sphaeramia_orbicularis', 'sphenodon_punctatus', 'stachyris_ruficeps', 'stegastes_partitus', 'strigops_habroptila', 'strix_occidentalis_caurina', 'struthio_camelus_australis', 'suricata_suricatta', 'sus_scrofa', 'sus_scrofa_bamei', 'sus_scrofa_berkshire', 'sus_scrofa_hampshire', 'sus_scrofa_jinhua', 'sus_scrofa_landrace', 'sus_scrofa_largewhite', 'sus_scrofa_meishan', 'sus_scrofa_pietrain', 'sus_scrofa_rongchang', 'sus_scrofa_tibetan', 'sus_scrofa_usmarc', 'sus_scrofa_wuzhishan', 'taeniopygia_guttata', 'takifugu_rubripes', 'terrapene_carolina_triunguis', 'tetraodon_nigroviridis', 'theropithecus_gelada', 'tupaia_belangeri', 'tursiops_truncatus', 'urocitellus_parryii', 'ursus_americanus', 'ursus_maritimus', 'ursus_thibetanus_thibetanus', 'varanus_komodoensis', 'vicugna_pacos', 'vombatus_ursinus', 'vulpes_vulpes', 'xenopus_tropicalis', 'xiphophorus_couchianus', 'xiphophorus_maculatus', 'zalophus_californianus', 'zonotrichia_albicollis', 'zosterops_lateralis_melanops']

brief_to_ensembl_dict = {'human':'homo_sapiens', 'mouse':'mus_musculus', 'zebrafish':'danio_rerio', 'rat':'rattus_norvegicus',
                     'pig':'sus_scrofa', 'chimpanzee':'pan_troglodytes', 'chicken':'gallus_gallus', 'cow':'bos_taurus', 'fruitfly':'drosophila_melanogaster', 'celegans':'caenorhabditis_elegans'}

if os.path.exists("config.py"):
    sys.path.insert(0, os.getcwd())
    import config 
    print('Using config file made at: '+config.time_made+'\nfor the species: '+config.species)
    if args.species in species_brief:
        species = brief_to_ensembl_dict[args.species]
    else: 
        species = args.species
    if species != config.species:
        raise Exception('Error\nInput species is different from loaded config file. Please delete config.py and rerun.')

else:
    if args.species in species_brief:
        species = brief_to_ensembl_dict[args.species]
    else: 
        species = args.species

    if not args.species:
        raise Exception('Error\nPlease make sure to enter a species with the \'--species\' or \'-s\' flags')

    if (args.species not in ensembl_species) and (args.species not in species_brief):
        raise Exception("Error\nThe species that you have entered,", args.species, "cannot be found. \
                        You can choose a species from this list:\n", species_brief, "\nor you can use its full name to access less used species. \
                        A complete list of available species is in the included species.txt file.")
    
    # Making the directory for the transcriptome files
    species_dir = "species." + species + "/"
    if not os.path.exists(species_dir):
        print('Making the species specific directory for',species,'\n')
        os.mkdir(species_dir)
    else:
        print('Species specific directory for',species,'already made\n')
    sys.path.append(os.getcwd()+'/'+species_dir)

    os.chdir(species_dir)

    # Download appropriate files
    raw_cDNA_site_data = subprocess.run(['curl', 'http://ftp.ensembl.org/pub/current_fasta/'+species+'/cdna/'],capture_output=True)
    cDNA_site_data = raw_cDNA_site_data.stdout.decode("utf-8").split("\"")
    for entry in cDNA_site_data: 
        if 'cdna.all.fa.gz' in entry and '</a>' not in entry: 
            cDNA_file_name = entry
    if len(cDNA_file_name)<3:
        raise(Exception('Couldn\'t find the appropriate files on ensemble. Please contact support'))
    if not os.path.exists('cDNA.fa'):
        print("Downloading the cDNA file\n")
        subprocess.run(['curl', 'http://ftp.ensembl.org/pub/current_fasta/'+species+'/cdna/'+cDNA_file_name,
                    '-o', 'cDNA.fa.gz'])
        os.system('gunzip cDNA.fa.gz')
    else:
        print("cDNA file already downloaded")

    raw_ncRNA_site_data = subprocess.run(['curl', 'http://ftp.ensembl.org/pub/current_fasta/'+species+'/ncrna/'],capture_output=True)
    ncRNA_site_data = raw_ncRNA_site_data.stdout.decode("utf-8").split("\"")
    for entry in ncRNA_site_data: 
        if 'ncrna.fa.gz' in entry and '</a>' not in entry: 
            ncRNA_file_name = entry
    if len(ncRNA_file_name)<3:
        raise(Exception('Couldn\'t find the appropriate files on ensemble. Please contact support'))
    if not os.path.exists('ncRNA.fa'):
        print("Downloading the ncRNA file\n")
        subprocess.run(['curl', 'http://ftp.ensembl.org/pub/current_fasta/'+species+'/ncrna/'+ncRNA_file_name,
                    '-o', 'ncRNA.fa.gz'])
        os.system('gunzip ncRNA.fa.gz')
    else:
        print("ncRNA file already downloaded")

    # Combining the cDNA and ncRNA files for simplicity and then extracting it
    if not os.path.exists('transcriptome.fa'):
        print('Combining cDNA and ncRNA files file')
        os.system('cat cDNA.fa ncRNA.fa >> transcriptome.fa')
    else:
        print('\nUsing already made transcriptome file.\n')

    # Using the generated transcriptome file to make the blast database for alignment of probe sequences
    if not os.path.exists('blast_db.nhr'):
        utils.make_blast_db('transcriptome.fa')
    else:
        print('\nUsing already made blast database files.\n')
    
    print("The blast database has successfully been made")

    print("The blast database is stored at blast_db")
    print("The transcriptome file is stored at transcriptome.fa")

    os.chdir('../')

    # Now we will make the config file, which can always be changed manually
    f = open("./config.py", "w")
    f.write("# This is the configuration file which defines the various parameters for making probes\n\n")

    f.write("time_made = \'"+ time.asctime( time.localtime(time.time()) )+'\'\n')
    f.write("species = \'" + species + '\'\n')
    
    if advanced_flag:       # If the advanced flag is chosen
        f.write("num_probes = "+str(int(input('Please enter the number of probes you would like to generate per transcript (default 50)')))+" #Number of probes to be created\n")
    else:
        f.write("num_probes = 50 #Number of probes to be created\n")
    
    if advanced_flag:
        f.write("probe_len = "+str(int(input('Please enter the length of the probes to be made (default 20)')))+" #Length of each probe\n")
    else:
        f.write("probe_len = 20 #Length of each probe\n")
    
    if advanced_flag:
        f.write("probe_dist = "+str(int(input('Please enter the distance between probes (default 2)')))+" #Distance between each probe\n")
    else:
        f.write("probe_dist = 2 #Distance between each probe\n")
    
    if advanced_flag:
        f.write("Tm = "+str(273+int(input('Please enter the minimum melting temperature for each probe in Celsius (default 40 Celsius)')))+" #Minimum melting temperature of probes in Kelvin\n")
    else:
        f.write("Tm = 313 #Minimum melting temperature of probes in Kelvin\n")

    if advanced_flag:
        f.write("Tm_hairpin = "+str(273+int(input('Please enter the maximum melting temperature for probe hairpins in Celsius (default 57 Celsius)')))+" #Maximum melting temperature of the hairpin of probes in Kelvin\n\n")
    else:
        f.write("Tm_hairpin = 330 #Maximum melting temperature of the hairpin of probes in Kelvin\n\n")

    if advanced_flag: 
        f.write("Heterodimer_melting_temperature_estimation_method = \'" +str(input("Please enter Tm or Alignment_score to indicate which method you would like to use to evaluate your probes?")))+'\'\n'
    else:
        f.write("Heterodimer_melting_temperature_estimation_method = \"Tm\"\n")
    
    if advanced_flag:
        f.write("Probe_set_selection_method = \'" + str(input("Please enter Dynamic or Greedy to indicate which probe selection selection you would like to use.")))+'\'\n'
    else:
        f.write("Probe_set_selection_method = \"Greedy\"\n")

    f.write("#Filter repeating nucleotides, AAAAA, GGGG, CCCC, TTTTT\n")
    f.write("filter_repeats = True\n")

    f.write("# Complete Human Transcriptome fasta file for determining probe specificity\n")
    f.write("transcriptome_file = \""+species_dir+"transcriptome.fa\"\n")
    f.write("transcriptome_index = \""+species_dir+"blast_db\"\n\n")

    f.write("# This is a directory in which the program will create different files, but these files can be deleted \n")
    f.write("# once the program has finished running. \n")
    if advanced_flag:
        f.write("temp_dir = \'"+str(input('Please enter the desired temporary directory, which will be deleted after running'))+"\'\n\n")
    else:
        f.write("temp_dir = \"./temp/\"\n\n")

    f.write("# This is a directory in which the program will create different files with the output of the program\n")
    if advanced_flag:
        f.write("output_dir = \'"+str(input('Please enter the desired output directory, which will be deleted after running'))+"\'\n\n")
    else:
        f.write("output_dir = \"./output_finalProbes/\"\n")

    f.close()

    sys.path.insert(0, os.getcwd())
    import config

##############################################################################
#  Configuration Parameters                                                  #
##############################################################################
transcriptome_file = config.transcriptome_file
transcriptome_index = config.transcriptome_index
probe_len = config.probe_len
probe_dist = config.probe_dist
Tm = config.Tm
Tm_hairpin = config.Tm_hairpin
filter_repeats = config.filter_repeats
num_probes = config.num_probes
species = config.species

temp_dir = config.temp_dir
if not os.path.exists(temp_dir):
    os.makedirs(temp_dir)
    
output_dir = config.output_dir
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
print("Designing smFISH Probes...")
print("\t num_probes ",num_probes)
print("\t probe_len ",probe_len)
print("\t Tm ",Tm)
print("\t Writing temporary files in ",temp_dir)
print("\t Writing final probes in ",output_dir)


# Load transcriptome file into a dictionary with keys of the form gene_symbol.transcript_id.gene_id
transcriptome_dict, gene_info_dataframe, rRNA_list = utils.load_transcriptome_dict(transcriptome_file)

# Load rRNA transcript IDs from pybiomart
# Not needed anymore since I can load this data from the transcriptome.fa file to save time
#rRNA_list_mart = utils.load_rRNA_data(species)

# Set up data structures
transcript_targets = {}

ensembl_gene_id_list = []
gene_symbol_list = []
ensembl_transcript_id_list = []


# Add user inputed sequences
if (args.transcripts_file):
    transcript_targets = utils.load_target_transcripts(args.transcripts_file)

# Add gene symbols from the command line
if (args.gene_symbols):
    gene_symbol_list += args.gene_symbols

# Add gene symbols from a file
if (args.gene_symbol_file):
    for line in open(args.gene_symbols_file): 
        gene_symbol = line.strip()
        gene_symbol_list.append(gene_symbol)

# Add gene IDs from the command line
if (args.ensembl_gene_ids):
    for ensembl_gene_id in list(args.ensembl_gene_ids):
        ensembl_gene_id_list.append(ensembl_gene_id)

# If gene symbols or gene IDs were input then use pybiomart to find the primary transcript IDs
if len(gene_symbol_list+ensembl_gene_id_list) > 0:
    primary_ensembl_transcript_ids = utils.load_gene_information(gene_symbol_list, ensembl_gene_id_list, species)
    for ensembl_transcript_id in primary_ensembl_transcript_ids:
        ensembl_transcript_id_list.append(ensembl_transcript_id)

# If transcript IDs were input then they are simply added to the transcript list
if (args.ensembl_transcript_ids):
    for ensembl_transcript_id in list(args.ensembl_transcript_ids):
        ensembl_transcript_id_list.append(ensembl_transcript_id)

# If there are transcript IDs in the transcript list then it loads them from the dictionary
if len(ensembl_transcript_id_list) > 0:
    for transcript_id in ensembl_transcript_id_list:
        transcript_targets[transcript_id] = transcriptome_dict[transcript_id]

# If there were no inputs error out here
elif len(transcript_targets) == 0:
    raise Exception('No genes input, please make sure to use the -g, -fa, -gid, -tid, or -gf \
                    (--gene_symbols, --fasta_file, --ensembl_gene_ids, --ensembl_transcript_ids, or --gene_symbols_file) options')


##############################################################################
# MAIN                                                                       #
##############################################################################

score_dict = {}

for transcript in transcript_targets:
    
    # Pulling gene ID and gene symbol from metadata
    target_gene_id = list(gene_info_dataframe.loc[gene_info_dataframe['transcript_id'] == transcript]['gene_id'])[0]
    target_gene_symbol = list(gene_info_dataframe.loc[gene_info_dataframe['transcript_id'] == transcript]['symbol'])[0]
    
    print("Processing transcript " + target_gene_symbol +"_"+transcript+ " ...\n")
    
    
    # Makes all probe possibilities and filters out ones with low Tm, high hairpin Tm, and repetitive sequences
    filtered_probe_array, filtered_probe_tm_dict = utils.create_and_filter_probe_possibilites(transcript_targets[transcript], probe_len, Tm, Tm_hairpin, num_probes, filter_repeats)
    
    # Write the filtered probe sequences to a fasta file for alignment via blast
    temp_probe_fasta_filename = temp_dir + transcript + ".fa"
    
    probe_file = open(temp_probe_fasta_filename, "w")
    for probe in filtered_probe_array:
        probe_name = probe[0]
        probe_seq = probe[1]
        probe_file.write(">" + probe_name + "\n")
        probe_file.write(probe_seq + "\n")
    probe_file.close()
    
    # Alignment of the filtered probes by blast
    aln_data, blast_run = utils.align_blast(temp_probe_fasta_filename, transcriptome_index, threads, word_size)
    
    print("Alignments: "+str(len(aln_data)))
    print('Blast input args: \n'+' '.join(blast_run.args)+'\n')
    
    # Load the alignment data in, then calculate the off target binding of each probe
    
    if config.Heterodimer_melting_temperature_estimation_method == "Alignment_score": 
        offtarget_scores_list = utils.compute_offtarget_scores(aln_data, transcript, rRNA_list)
    elif config.Heterodimer_melting_temperature_estimation_method == "Tm":
        offtarget_scores_list = utils.compute_offtarget_tm(aln_data, transcript, rRNA_list)
    else: 
        print("Unsupported heterodimer melting temperature estimation method: " , config.Heterodimer_melting_temperature_estimation_method)
        print("Please change it in the config.py file to \"Alignment_score\" or \"Tm\"")
        sys.exit()    
    
    # Makes an output file of all the probes off target scores
    '''
    probe_scores_file_name = output_dir + transcript + "_probe_scores.tsv"
    probe_scores_file = open(probe_scores_file_name, 'w')
    for probe in offtarget_scores_list:
        probe_name = probe[0]
        probe_score = str(probe[1])
        probe_scores_file.write(probe_name + "\t" + probe_score + "\n")
    probe_scores_file.close()
    '''
    
    # Select the final set of probes

    if config.Probe_set_selection_method == "Greedy":
        print('Probe set selection being done by Greedy selection (picking the best probe each time)')
        final_probes, total_score_greedy = utils.select_nonoverlapping_probes(offtarget_scores_list, probe_len, probe_dist)
        #final_probes, total_score_greedy = utils.rand_optimal_set(offtarget_scores_list, probe_len, probe_dist, num_probes, 100000)
        
    elif config.Probe_set_selection_method == "Dynamic":
        print('Probe set selection being done by Dynamic programming (picking the best combination of probes)')
        final_probes_greedy, total_score_greedy = utils.select_nonoverlapping_probes(offtarget_scores_list, probe_len, probe_dist)
        final_probes_ts, total_score_ts = utils.find_optimal_set(offtarget_scores_list, probe_len + probe_dist, num_probes, final_probes_greedy, beam=0.1)
        score_dict[target_gene_symbol] = (sum(tup[1] for tup in final_probes_greedy[0:50]),total_score_ts)
        if total_score_ts < total_score_greedy:
            final_probes = final_probes_ts
        else:
            final_probes = final_probes_greedy
    else: 
        print("Unsupported probe set selection method: " , config.Probe_set_selection_method)
        print("Please change it in the config.py file to \"Greedy\" or \"Dynamic\"")

    # Write the final list of probes to a fasta file and also print
    
    probe_filename = output_dir + target_gene_symbol+"_"+target_gene_id+"_"+transcript + "_probes.fa"
    print("Probe selection is complete. Writing probes to "+ probe_filename)
    
    probe_file = open(probe_filename, "w")
    probes_added = 0
    for probe in final_probes:
        probe_name = probe[0]
        probe_file.write(">" + probe_name + "\n")
        probe_seq = filtered_probe_tm_dict[probe_name][0]
        probe_file.write(probe_seq + "\n")
        probes_added += 1
        if probes_added >= num_probes: 
            break
    probe_file.close()
    
    
print(score_dict)