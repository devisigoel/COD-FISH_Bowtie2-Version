#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 09:20:13 2021

@author: Ryan 'mentor' Dikdan, Devisi Goel

This is a script which takes gene names and generates smFISH probes for the gene.
The motivation is to find optimal oligos for targeting RNAs given full transcriptome data.
"""

##############################################################################
#   Imports                                                                  #
##############################################################################
import subprocess
import primer3
import pandas as pd
import pybiomart
import math
from tqdm import tqdm


##############################################################################
#   Finding primary transcripts and loading fastas                           #
##############################################################################

def load_transcriptome_dict(transcriptome_file):

    print('\nLoading the transcriptome into a python dictionary')
    
    transcriptome_dict = {}
    gene_info_array = []

    for transcriptome_line in tqdm(open(transcriptome_file)):
        transcriptome_line = transcriptome_line.strip()
        if transcriptome_line[0] == ">":
            transcript_full_name = transcriptome_line[1:]
            transcript_id = transcript_full_name.split(' ')[0].split(".")[0]
            
            gene_symbol = transcript_full_name.partition('gene_symbol:')[2].split(' ')[0]
            if gene_symbol == '':
                gene_symbol = 'not_listed'
            transcript_id = transcript_full_name.split(' ')[0].split(".")[0]
            gene_id = transcript_full_name.partition('gene:')[2].split('.')[0]
            gene_biotype = transcript_full_name.partition('gene_biotype:')[2].split(' ')[0]

            gene_info_array.append([gene_symbol, transcript_id, gene_id, gene_biotype])
            
            transcriptome_dict[transcript_id] = ""      # Makes an entry using the name, where the value is blank
        else:
                                                                # Adds each line to the previously identified ID's sequence
            transcriptome_dict[transcript_id] = transcriptome_dict[transcript_id] + transcriptome_line
    gene_info_dataframe = pd.DataFrame(gene_info_array, columns=["symbol","transcript_id","gene_id","gene_biotype"])
    
    # This filters it to make a list of Transcript IDs for rRNAs
    rRNA_list = list(gene_info_dataframe[gene_info_dataframe['gene_biotype'] == 'rRNA']['transcript_id'])
    return(transcriptome_dict, gene_info_dataframe, rRNA_list)

def load_target_transcripts(input_fasta):
    # Parsing, converting fasta to dictionary
    # Fasta should have sequence names of the format:
    # >gene_name.transcript_id.gene_id for example:
    # >PTGS2.
    # This is so that we can exclude binding penalties to the target gene
    transcript_dict = {}
    for line in open(input_fasta):
        line = line.strip()
        if line[0] == ">":
            transcript_full_name = line[1:]
            transcript_dict[transcript_full_name] = ''
        else:
            transcript_dict[transcript_full_name] = transcript_dict[transcript_full_name] + line.upper()
    return(transcript_dict)

def load_gene_information(gene_symbols_list, gene_id_list, species):

    print('\nLoading transcript IDs for genes by symbol and ID through pybiomart')

    # Call the Ensemble server
    server = pybiomart.Server(host='http://www.ensembl.org')
    mart = server['ENSEMBL_MART_ENSEMBL']

    # Convert the species name to the abbreviated name used by the ensembl biomart
    species_abbr = ''
    for name in range(len(species.split('_'))-1):
        species_abbr += species.split('_')[name][0]
    species_abbr += species.split('_')[len(species.split('_'))-1]

    # Pull the dataset from the mart
    dataset = mart[species_abbr+'_gene_ensembl']
    ranking_list_temp = []
    primary_transcript_list = []

    # Pull and filter the appropriate data, using what data is available for the species
    # Multiple try except clauses are used in case the information isn't available for that species
    try:            # Human
        biomart_data = dataset.query(attributes=['ensembl_transcript_id','ensembl_gene_id','external_gene_name','transcript_length', 'transcript_is_canonical'])
        gene_symbol_list_info = biomart_data[biomart_data['Gene name'].isin(gene_symbols_list)]
        gene_symbol_list_info = gene_symbol_list_info[gene_symbol_list_info['Ensembl Canonical'] == 1.0]
        gene_id_list_info = biomart_data[biomart_data['Gene stable ID'].isin(gene_id_list)]
        gene_id_list_info = gene_id_list_info[gene_id_list_info['Ensembl Canonical'] == 1.0]

        for symbol in gene_symbols_list:
            ranking_list_temp = gene_symbol_list_info[gene_symbol_list_info['Gene name'] == symbol]
            if len(ranking_list_temp) > 1:
                ranking_list_temp = ranking_list_temp[ranking_list_temp['Transcript length (including UTRs and CDS)'] == ranking_list_temp['Transcript length (including UTRs and CDS)'].max()]
            primary_transcript_list.append(list(ranking_list_temp.head(1)['Transcript stable ID'])[0])     # Selects 1st entry in case multiple primary transcripts have the same length
        for gene_id in gene_id_list:
            ranking_list_temp = gene_id_list_info[gene_id_list_info['Gene stable ID'] == gene_id]
            if len(ranking_list_temp) > 1:
                ranking_list_temp = ranking_list_temp[ranking_list_temp['Transcript length (including UTRs and CDS)'] == ranking_list_temp['Transcript length (including UTRs and CDS)'].max()]
            primary_transcript_list.append(list(ranking_list_temp.head(1)['Transcript stable ID'])[0])     # Selects 1st entry in case multiple primary transcripts have the same length

    except: 
        try:        # Mouse
            biomart_data = dataset.query(attributes=['ensembl_transcript_id','ensembl_gene_id','external_gene_name','transcript_length', 'transcript_gencode_basic'])
            gene_symbol_list_info = biomart_data[biomart_data['Gene name'].isin(gene_symbols_list)]
            gene_symbol_list_info = gene_symbol_list_info[gene_symbol_list_info['GENCODE basic annotation'] == 'GENCODE basic']
            gene_id_list_info = biomart_data[biomart_data['Gene stable ID'].isin(gene_id_list)]
            gene_id_list_info = gene_id_list_info[gene_id_list_info['GENCODE basic annotation'] == 'GENCODE basic']

            for symbol in gene_symbols_list:
                ranking_list_temp = gene_symbol_list_info[gene_symbol_list_info['Gene name'] == symbol]
                if len(ranking_list_temp) > 1:
                    ranking_list_temp = ranking_list_temp[ranking_list_temp['Transcript length (including UTRs and CDS)'] == ranking_list_temp['Transcript length (including UTRs and CDS)'].max()]
                primary_transcript_list.append(list(ranking_list_temp.head(1)['Transcript stable ID'])[0])      # Selects 1st entry in case multiple primary transcripts have the same length
            for gene_id in gene_id_list:
                ranking_list_temp = gene_id_list_info[gene_id_list_info['Gene stable ID'] == gene_id]
                if len(ranking_list_temp) > 1:
                    ranking_list_temp = ranking_list_temp[ranking_list_temp['Transcript length (including UTRs and CDS)'] == ranking_list_temp['Transcript length (including UTRs and CDS)'].max()]
                primary_transcript_list.append(list(ranking_list_temp.head(1)['Transcript stable ID'])[0])      # Selects 1st entry in case multiple primary transcripts have the same length  

        except:
            try:    # Appris species, cow, chicken, pig, zebrafish, etc. although APPRIS annnotation is for protein coding regions specifically
                    # the longest of the principle transcripts will be used 
                biomart_data = dataset.query(attributes=['ensembl_transcript_id','ensembl_gene_id','external_gene_name','transcript_length', 'transcript_appris'])
                gene_symbol_list_info = biomart_data[biomart_data['Gene name'].isin(gene_symbols_list)]
                gene_symbol_list_info = gene_symbol_list_info[gene_symbol_list_info['APPRIS annotation'] == 'principal1']
                gene_id_list_info = biomart_data[biomart_data['Gene stable ID'].isin(gene_id_list)]
                gene_id_list_info = gene_id_list_info[gene_id_list_info['APPRIS annotation'] == 'principal1']

                for symbol in gene_symbols_list:
                    ranking_list_temp = gene_symbol_list_info[gene_symbol_list_info['Gene name'] == symbol]
                    if len(ranking_list_temp) > 1:
                        ranking_list_temp = ranking_list_temp[ranking_list_temp['Transcript length (including UTRs and CDS)'] == ranking_list_temp['Transcript length (including UTRs and CDS)'].max()]
                    primary_transcript_list.append(list(ranking_list_temp.head(1)['Transcript stable ID'])[0])      # Selects 1st entry in case multiple primary transcripts have the same length
                for gene_id in gene_id_list:
                    ranking_list_temp = gene_id_list_info[gene_id_list_info['Gene stable ID'] == gene_id]
                    if len(ranking_list_temp) > 1:
                        ranking_list_temp = ranking_list_temp[ranking_list_temp['Transcript length (including UTRs and CDS)'] == ranking_list_temp['Transcript length (including UTRs and CDS)'].max()]
                    primary_transcript_list.append(list(ranking_list_temp.head(1)['Transcript stable ID'])[0])      # Selects 1st entry in case multiple primary transcripts have the same length  


            except: # All other species that don't have annotations, the longest transcript will be used which has ~75% accuracy, best we can do given the circumstances
                biomart_data = dataset.query(attributes=['ensembl_transcript_id','external_gene_name','transcript_length'])
                gene_symbol_list_info = biomart_data[biomart_data['Gene name'].isin(gene_symbols_list)]
                for symbol in gene_symbols_list:
                    ranking_list_temp = gene_symbol_list_info[gene_symbol_list_info['Gene name'] == symbol]
                    if len(ranking_list_temp) > 1:
                        ranking_list_temp = ranking_list_temp[ranking_list_temp['Transcript length (including UTRs and CDS)'] == ranking_list_temp['Transcript length (including UTRs and CDS)'].max()]
                    primary_transcript_list.append(list(ranking_list_temp.head(1)['Transcript stable ID'])[0])      # Selects 1st entry in case multiple primary transcripts have the same length

                biomart_data = dataset.query(attributes=['ensembl_transcript_id','ensembl_gene_id','external_gene_name','transcript_length', 'transcript_is_canonical'])
                gene_symbol_list_info = biomart_data[biomart_data['Gene name'].isin(gene_symbols_list)]
                gene_id_list_info = biomart_data[biomart_data['Gene stable ID'].isin(gene_id_list)]

                for symbol in gene_symbols_list:
                    ranking_list_temp = gene_symbol_list_info[gene_symbol_list_info['Gene name'] == symbol]
                    if len(ranking_list_temp) > 1:
                        ranking_list_temp = ranking_list_temp[ranking_list_temp['Transcript length (including UTRs and CDS)'] == ranking_list_temp['Transcript length (including UTRs and CDS)'].max()]
                    primary_transcript_list.append(list(ranking_list_temp.head(1)['Transcript stable ID'])[0])      # Selects 1st entry in case multiple primary transcripts have the same length
                for gene_id in gene_id_list:
                    ranking_list_temp = gene_id_list_info[gene_id_list_info['Gene stable ID'] == gene_id]
                    if len(ranking_list_temp) > 1:
                        ranking_list_temp = ranking_list_temp[ranking_list_temp['Transcript length (including UTRs and CDS)'] == ranking_list_temp['Transcript length (including UTRs and CDS)'].max()]
                    primary_transcript_list.append(list(ranking_list_temp.head(1)['Transcript stable ID'])[0])      # Selects 1st entry in case multiple primary transcripts have the same length  


        if len(primary_transcript_list) < len(gene_symbols_list + gene_id_list):
            print('Some genes listed could not be found, please make sure they were spelt correctly, \n\
            with the correct capitalization, and that the official gene symbols were used. \n\
            A good website to check this on is ensembl.org/index.html\n\
            Search for the gene and then look for the "Name"')

    return(primary_transcript_list)

def load_rRNA_data(species):

    print('\nLoading rRNA transcript IDs through pybiomart')

    # Call the Ensemble server
    server = pybiomart.Server(host='http://www.ensembl.org')
    mart = server['ENSEMBL_MART_ENSEMBL']

    # Convert the species name to the abbreviated name used by the ensembl biomart
    species_abbr = ''
    for name in range(len(species.split('_'))-1):
        species_abbr += species.split('_')[name][0]
    species_abbr += species.split('_')[len(species.split('_'))-1]

    # Pull the dataset from the mart
    dataset = mart[species_abbr+'_gene_ensembl']

    rRNA_query_result = dataset.query(attributes=['ensembl_transcript_id','external_gene_name'], filters={'biotype':'rRNA'})
    rRNA_list = list(rRNA_query_result['Transcript stable ID'])

    return(rRNA_list)
##############################################################################
# Create Probe Possibilities, filter by Tm and Tm_hairpin                    #
##############################################################################
def create_and_filter_probe_possibilites(target_transcript_seq, probe_len, Tm, Tm_hairpin, num_probes, filter_repeats = True):
    
    print("Making and filtering probe possibilities by Tm, repeats, and hairpins\n")

    final_probe_tm_dict = {}
    filtered_by_tm = 0
    filtered_by_repeat = 0
    filtered_by_hairpin = 0

    # Make all probe possibilities
    candidate_probe_list = []
    for location in range(0,len(target_transcript_seq)-probe_len + 1):
        probe_region = target_transcript_seq[location:location+probe_len]
        probe_seq = reverse_complement(probe_region)
        probe_name = "probe_" + str(location+1)
        candidate_probe_list.append((probe_name, probe_seq))

    # Filter out if Tm is too low
    Tm_filtered_probes_list = []
    for candidate_probe in candidate_probe_list:
        probe_name = candidate_probe[0]
        probe_seq = candidate_probe[1]
        probe_tm = 273 + primer3.calcTm(probe_seq)
        if float(probe_tm) > Tm:
            Tm_filtered_probes_list.append((probe_name, probe_seq, probe_tm))
        else:
            filtered_by_tm += 1
    if len(Tm_filtered_probes_list) < num_probes:
        raise Exception("Error: Too many probes filtered out due to low probe binding Tm, \
                        consider lowering the minimum Tm in the config.py file and rerunning.")

    Tm_rep_filtered_probes = []
    if filter_repeats == True:
        for candidate_probe in Tm_filtered_probes_list:
            probe_name = candidate_probe[0]
            probe_seq = candidate_probe[1] 
            probe_tm = candidate_probe[2]
            if ('AAAAA' or 'TTTTT' or 'GGGG' or 'CCCC') not in probe_seq:
                Tm_rep_filtered_probes.append((probe_name,probe_seq,probe_tm))
            else:
                filtered_by_repeat += 1
    if len(Tm_rep_filtered_probes) < num_probes:
        raise Exception("Error: Too many probes filtered out due to repetitive nucleotides (ie. AAAAA or GGGG), \
                        consider changing \"filter_repeats = True\" to \"filter_repeats = False\" in the config.py file and rerunning.")

    # Filter out if high hairpin Tm
    Tm_rep_hairpin_filtered_probes = []
    for candidate_probe in Tm_rep_filtered_probes:
        probe_name = candidate_probe[0]
        probe_seq = candidate_probe[1]
        probe_tm = candidate_probe[2]
        probe_hairpin_tm = primer3.calcHairpinTm(probe_seq) + 273
        if float(probe_hairpin_tm) < Tm_hairpin:
            Tm_rep_hairpin_filtered_probes.append((probe_name,probe_seq,probe_tm,probe_hairpin_tm))
            final_probe_tm_dict[probe_name] = (probe_seq, probe_tm)
        else:
            filtered_by_hairpin += 1
    if len(Tm_rep_hairpin_filtered_probes) < num_probes:
        raise Exception("Error: Too many probes filtered out due to repetitive nucleotides (ie. AAAAA or GGGG), \
                        consider changing \"filter_repeats = True\" to \"filter_repeats = False\" in the config.py file and rerunning.")
    
    print("\nTm filtered: "+str(filtered_by_tm))
    print("Repeat filtered: "+str(filtered_by_repeat))
    print("Hairpin filtered: "+str(filtered_by_hairpin)+'\n')
        
    return(Tm_rep_hairpin_filtered_probes, final_probe_tm_dict)


##############################################################################
#  Scoring functions                                                         #
##############################################################################

def sort_by_score(tuple):
    # This function is used by .sort to sort the tuples based on their score
    # All it does is take a tuple and return the 2nd value which is then used to sort
    return(tuple[1])

def compute_offtarget_scores(alignment_data, target_ensembl_transcript_id, rRNA_list):
    print('Calculating probe candidate off-target scores using blast scores\n')
    offtarget_match_scores = {}

    for aln in alignment_data:
        probe_name = aln[0]
        # Pull ensemble transcript ID from name
        aligned_ensembl_transcript_id = aln[1].split(".")[0]
        #if aligned_ensembl_transcript_id == "*":
            #continue
            # This skips if the probe didn't align to anything
            # I don't believe that this happens with blast
        if aligned_ensembl_transcript_id == target_ensembl_transcript_id:
            continue
            # Doesn't add anything to the probe score if the alignment is an on target alignment

        match_score = int(aln[4])        # This is the score from blast which is the Raw Score

        # Significantly punishes for rRNA alignment. This should make it so these probes aren't selected.
        if aligned_ensembl_transcript_id in rRNA_list:
            match_score = match_score + 10000

        if probe_name not in offtarget_match_scores:
            offtarget_match_scores[probe_name] = 0
            # Adds entry to offtarget scores if probe isn't already in it
        
        # Add the score to the off target score for that probe
        offtarget_match_scores[probe_name] = offtarget_match_scores[probe_name] + match_score

    # Converts dictionary to an array of tuples then sorts it based on the scores
    offtarget_match_scores_tuple_list = list(offtarget_match_scores.items())
    offtarget_match_scores_tuple_list.sort(key=sort_by_score)

    return(offtarget_match_scores_tuple_list)

def compute_offtarget_tm(alignment_data, target_ensembl_id, rRNA_list):
    print('Calculating probe candidate off-target scores using calculated alignment Tm\'s\n')
    offtarget_Tm_scores = {}
    
    for aln in alignment_data:
        probe_name = aln[0]
        aligned_ensembl_transcript_id = aln[1].split(".")[0]
        probe_seq = aln[2]
        aligned_transcript_seq_fragment = reverse_complement(aln[3])

        if aligned_ensembl_transcript_id == target_ensembl_id:
            continue                                                
            # Doesn't add anything to the probe score if the alignment is an on target alignment
        
        aln_tm = (primer3.calcHeterodimer(probe_seq, aligned_transcript_seq_fragment).tm + 273)/(primer3.calcTm(probe_seq) + 273) + 1
        # Punished 1 per alignment as well as 0-1 based on the alignment 
        # (so for a perfect off target alignment penalty of 2 and for a very weak off target alignment penalty of 1)
    
        if aligned_ensembl_transcript_id in rRNA_list:
            aln_tm = aln_tm + (aln_tm + 10000)

        if probe_name not in offtarget_Tm_scores:
            offtarget_Tm_scores[probe_name] = 0
            # Adds entry to offtarget scores if probe isn't already in it
            
        offtarget_Tm_scores[probe_name] = offtarget_Tm_scores[probe_name] + aln_tm
        
    # Converts dictionary to an array of tuples then sorts it based on the cumulative Tm scores
    offtarget_Tm_scores_tuple_list = list(offtarget_Tm_scores.items())
    offtarget_Tm_scores_tuple_list.sort(key=sort_by_score)

    return(offtarget_Tm_scores_tuple_list)
        

##############################################################################
#  Blast functions                                                         #
##############################################################################

def make_blast_db(fasta_file):
    print("\nBlastn will now generate the database with which to align the candidate probe sequences.")
    blast_db_run = subprocess.run(['makeblastdb',
                                '-in', fasta_file,
                                '-dbtype', 'nucl',
                                '-parse_seqids',
                                '-out','blast_db'],
                                capture_output=True)
    if blast_db_run.returncode != 0:
        print('Building the database didn\'t work as planned. Make sure the makeblastdb version you have \
                is up to date and all the necessary files are present and in the right folders.')
    return(blast_db_run)  

def align_blast(probe_fasta_filename, transcriptome_db, threads, word_size = 7):
    print('Aligning probe candidates to transcriptome using blast\n')
    # To understand how the scoring and alignment works look at this site: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastn_application_options/
    # In summary, for the raw score 'score' in task blastn-short (which we're doing) the reward is 1 for each match, and -3 for each mismatch. 
    # If a gap is present then it penalizes -5 for gap opening and -2 for extending the gap.
    blast_run = subprocess.run(['blastn',
                                '-task', 'blastn-short',
                                '-db', transcriptome_db,
                                '-query', probe_fasta_filename, 
                                '-num_threads', str(threads), 
                                '-word_size', str(word_size),
                                '-strand', 'minus',
                                '-outfmt','6 qacc sacc qseq sseq score'],
                                capture_output=True)

    raw_data = blast_run.stdout.decode("utf-8").split('\n')    # Converts to string then to a list where each element is a line
    raw_data.remove('') 
    alignment_data = []
    for line in raw_data:
        alignment_data.append(line.split('\t'))               # Splits each alignment lines
    return(alignment_data, blast_run)

##############################################################################
# "Greedy" Set Selection Method                                              #
##############################################################################
def select_nonoverlapping_probes(offtarget_scores_tuple_list, probe_len, probe_dist):
    # Takes sorted probe tuple list
    nucleotide_idx_covered_already = []
    final_probe_set = []

    total_score = 0.0
    for probe_tuple in offtarget_scores_tuple_list:
        probe_name = probe_tuple[0]
        probe_idx = int(probe_name.split('_')[1])
        probe_score = probe_tuple[1]
        
        if probe_idx in nucleotide_idx_covered_already or probe_idx+probe_len in nucleotide_idx_covered_already:
            # Won't add probe if already covered by a probe
            continue
        
        final_probe_set.append(probe_tuple)
        total_score += probe_score
        nucleotide_idx_covered_already.extend(range(probe_idx-probe_dist,probe_idx+probe_len+probe_dist))
        #if len(final_probe_set) == num_probes:
            #break
    final_probe_set.sort(key=sort_by_score)
    #print(final_probe_set[0:50])
    #print(sum(tup[1] for tup in final_probe_set[0:50]))
    return(final_probe_set, sum(tup[1] for tup in final_probe_set[0:50]))


##############################################################################
#  DNA helper functions                                                      #
##############################################################################
def complement(seq):
    complement_dict = {'A':'T','T':'A','C':'G','G':'C'}
    complement_list = []
    for nt in seq.upper():
        complement_list.append(complement_dict[nt])
    complement_seq = ''.join(complement_list)
    return(complement_seq)

def reverse_complement(seq):
    return(complement(seq)[::-1])


##############################################################################
#  Token passing dynamic programming selection function                      #
##############################################################################

def find_optimal_set(probe_score_list, overlap, output_set_size, final_probes_greedy, beam=5.0):

    #STEP 1: Create an ordered list of probe numbers and score, sorted by probe numbers. 
    #   We assume that the probes are numbered based on their nucleotide position in the mRNA 
    probe_dat = []
    for ii in probe_score_list:
        p_num = int(ii[0].split("_")[1])
        p_score = ii[1]
        assert(p_score >= 0.), 'all scores need to be >=0.'
        probe_dat.append((p_num, p_score))

    # sort by probe position ... required for the algorithm to work
    probe_dat.sort(key=lambda p: p[0], reverse=False)


    #STEP 2: Make a score profile based on greedy selection approach.  This is used in pruning
    greedy_probeidx_score = []
    for ii in final_probes_greedy:
        probe_idx = int(ii[0].split("_")[1])
        probe_scr = ii[1]
        greedy_probeidx_score.append( (probe_idx, probe_scr) )
    greedy_probeidx_score.sort(key=lambda p: p[0], reverse=False)
    
    search_score_profile = []
    search_score_profile.append(greedy_probeidx_score[0][1])
    for ii in range(1,len(greedy_probeidx_score)):
        search_score_profile.append(search_score_profile[-1] + greedy_probeidx_score[ii][1])
    if (len(search_score_profile) < output_set_size):
        lastscr = search_score_profile[-1]
        for ii in range(len(search_score_profile), output_set_size):
            search_score_profile.append(lastscr)


    #STEP 3: Define a token class 
    class probe_set:
        def __init__(self, Tm=0):
            self.back_ptr = -1
            self.num_probes = 0
            self.Tm = Tm
            

    #STEP 4: initialize token list
    token_list = []
    for aa in probe_dat:
        #add a list with set_size tokens
        t_list = []
        for jj in range(output_set_size):
            t_list.append(probe_set())
        if (aa[1] <= search_score_profile[0] + beam):
            t = probe_set()
            t.num_probes = 1
            t.Tm = aa[1]
            t_list[0] = t
        token_list.append(t_list)


    #STEP 5: token passing algorithm
    tot_work = len(token_list)*(len(token_list) + 1)/2.0
    done_work = 0
    for ii in range(len(token_list)):
        probe_ii = probe_dat[ii][0]
        t_list_ii = token_list[ii]

        for jj in range(ii+1,len(token_list)):
            probe_jj = probe_dat[jj][0]
            jj_tm = probe_dat[jj][1]
 
            #check if probe_jj is allowed to be in the same set as probe_ii
            if probe_ii + overlap > probe_jj: continue

            #now propagate all tokens in t_list to token list at jj
            t_list_jj = token_list[jj]
            for tt in t_list_ii:
                if tt.num_probes == 0:
                    continue
                
                #let's try to extend tt by adding probe_jj to it
                new_score = tt.Tm + jj_tm
                newPos = tt.num_probes
           
                #are we going beyond set_size?
                if (newPos >= output_set_size): continue

                #beam pruning
                if (new_score > search_score_profile[newPos] + beam): continue
                if (new_score < search_score_profile[newPos]): search_score_profile[newPos] = new_score
            
                if (t_list_jj[newPos].num_probes == 0) or (t_list_jj[newPos].Tm > new_score): 
                    #Make a new token which extends tt by adding probe jj:
                    x = probe_set()
                    x.back_ptr = ii
                    x.num_probes = tt.num_probes + 1
                    x.Tm = tt.Tm + jj_tm
                    t_list_jj[newPos] = x

            done_work += 1
        if ii%100 == 0: 
            print('optimal probe set search.  done % : ', round(done_work*100.0/tot_work,2), end='\r')
    
    #STEP 6: all token propagation is done! now let's find best probe set of set_size
    #first find out what's the max number of probes we can generate
    num_probes_available = 0
    for tL in token_list:
        for tt in tL:
            if (tt.num_probes > num_probes_available):
                num_probes_available = tt.num_probes
            
    bestTm = math.inf
    bestii = -1
    for ii in range(len(token_list)):
        if (token_list[ii][num_probes_available - 1].num_probes == 0): continue
        if token_list[ii][num_probes_available - 1].Tm < bestTm:
            bestTm = token_list[ii][num_probes_available - 1].Tm
            bestii = ii
                
    finaltoken = token_list[bestii][num_probes_available - 1] 

    final_probes = [probe_dat[bestii]]
    while( finaltoken.back_ptr != -1 ):
        final_probes.append( probe_dat[finaltoken.back_ptr] )
        finaltoken = token_list[finaltoken.back_ptr][finaltoken.num_probes - 2]

    final_probes.sort(key=lambda p: p[1], reverse=False)

    retL = []
    total_score = 0.
    for ii in final_probes:
        nm = "probe_" + str(ii[0])
        sc = ii[1]
        total_score += sc
        retL.append( (nm, sc) )
    print('optimal probe set search.  done % : ', 100.0)
    #print(retL)
    #print(total_score)
    
    return( retL, total_score )