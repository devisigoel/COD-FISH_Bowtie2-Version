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
import math, copy
import re
import functools
import sys, json


##############################################################################
#   Finding primary transcripts and loading fastas                           #
##############################################################################

def load_transcriptome_dict(transcriptome_file):

    print('\nLoading the transcriptome into a python dictionary')
    
    transcriptome_dict = {}
    for transcriptome_line in open(transcriptome_file):
        transcriptome_line = transcriptome_line.strip()
        if transcriptome_line[0] == ">":
            transcript_full_name = transcriptome_line[1:]
            transcript_ensembl_id = transcript_full_name.split(".")[0]
            transcriptome_dict[transcript_ensembl_id] = ""      # Makes an entry using the name, where the value is blank
        else:
                                                                # Adds each line to the previously identified ID's sequence
            transcriptome_dict[transcript_ensembl_id] = transcriptome_dict[transcript_ensembl_id] + transcriptome_line
    return(transcriptome_dict)

def load_target_transcripts(input_fasta):
    # Parsing, converting fasta to dictionary
    transcript_dict = {}
    for line in open(input_fasta):
        line = line.strip()
        if line[0] == ">":
            target_name = line[1:]
            target_ensembl_id = target_name.split(".")[0] 
            transcript_dict[target_ensembl_id] = ''
        else:
            transcript_dict[target_ensembl_id] = transcript_dict[target_ensembl_id] + line.upper()
    return(transcript_dict)

def load_gene_information(list_of_genes, species):

    print('\nLoading transcript targets and genes by symbol and ID through pybiomart')

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
    primary_transcript_list = pd.DataFrame()

    # Pull and filter the appropriate data, using what data is available for the species
    # Multiple try except clauses are used in case the information isn't available for that species
    try:            # Human
        biomart_data = dataset.query(attributes=['ensembl_transcript_id','external_gene_name','transcript_length', 'transcript_is_canonical'])
        gene_list_info = biomart_data[biomart_data['Gene name'].isin(list_of_genes)]
        print(gene_list_info)
        gene_list_info = gene_list_info[gene_list_info['Ensembl Canonical'] == 1.0]
        for gene in list_of_genes:
            ranking_list_temp = gene_list_info[gene_list_info['Gene name'] == gene]
            if len(ranking_list_temp) > 1:
                ranking_list_temp = ranking_list_temp[ranking_list_temp['Transcript length (including UTRs and CDS)'] == ranking_list_temp['Transcript length (including UTRs and CDS)'].max()]
            primary_transcript_list = primary_transcript_list.append(ranking_list_temp.head(1))     # Selects 1st entry in case multiple primary transcripts have the same length
    except: 
        try:        # Mouse
            biomart_data = dataset.query(attributes=['ensembl_transcript_id','external_gene_name','transcript_length', 'transcript_gencode_basic'])
            gene_list_info = biomart_data[biomart_data['Gene name'].isin(list_of_genes)]
            gene_list_info = gene_list_info[gene_list_info['GENCODE basic annotation'] == 'GENCODE basic']
            for gene in list_of_genes:
                ranking_list_temp = gene_list_info[gene_list_info['Gene name'] == gene]
                if len(ranking_list_temp) > 1:
                    ranking_list_temp = ranking_list_temp[ranking_list_temp['Transcript length (including UTRs and CDS)'] == ranking_list_temp['Transcript length (including UTRs and CDS)'].max()]
                primary_transcript_list = primary_transcript_list.append(ranking_list_temp.head(1))     # Selects 1st entry in case multiple primary transcripts have the same length 
        except:
            try:    # Appris species, cow, chicken, pig, zebrafish, etc. although APPRIS annnotation is for protein coding regions specifically
                    # the longest of the principle transcripts will be used 
                biomart_data = dataset.query(attributes=['ensembl_transcript_id','external_gene_name','transcript_length','transcript_appris'])
                gene_list_info = biomart_data[biomart_data['Gene name'].isin(list_of_genes)]
                gene_list_info = gene_list_info[gene_list_info['APPRIS annotation'] == 'principal1']
                for gene in list_of_genes:
                    ranking_list_temp = gene_list_info[gene_list_info['Gene name'] == gene]
                    if len(ranking_list_temp) > 1:
                        ranking_list_temp = ranking_list_temp[ranking_list_temp['Transcript length (including UTRs and CDS)'] == ranking_list_temp['Transcript length (including UTRs and CDS)'].max()]
                    primary_transcript_list = primary_transcript_list.append(ranking_list_temp.head(1))     # Selects 1st entry in case multiple primary transcripts have the same length 
            except: # All other species that don't have annotations, the longest transcript will be used which has ~75% accuracy, best we can do given the circumstances
                biomart_data = dataset.query(attributes=['ensembl_transcript_id','external_gene_name','transcript_length'])
                print(biomart_data[100:200])
                gene_list_info = biomart_data[biomart_data['Gene name'].isin(list_of_genes)]
                for gene in list_of_genes:
                    ranking_list_temp = gene_list_info[gene_list_info['Gene name'] == gene]
                    if len(ranking_list_temp) > 1:
                        ranking_list_temp = ranking_list_temp[ranking_list_temp['Transcript length (including UTRs and CDS)'] == ranking_list_temp['Transcript length (including UTRs and CDS)'].max()]
                    primary_transcript_list = primary_transcript_list.append(ranking_list_temp.head(1))     # Selects 1st entry in case multiple primary transcripts have the same length 
        if len(primary_transcript_list) < len(list_of_genes):
            print('Some genes listed could not be found, please make sure they were spelt correctly, \n\
            with the correct capitalization, and that the official gene symbols were used. \n\
            A good website to check this on is ensembl.org/index.html\n\
            Search for the gene and then look for the "Name"')
    
    return(primary_transcript_list)
##############################################################################
# Create Probe Possibilities, filter by Tm and Tm_hairpin                    #
##############################################################################
def create_and_filter_probe_possibilites(target_transcript_seq, probe_len, Tm, Tm_hairpin, num_probes, filter_repeats = True):
    
    print("Making and filtering probe possibilities by Tm, repeats, and hairpins\n")


    probe_Tm_dict = {}
    filtered_by_tm = 0
    filtered_by_repeat = 0
    filtered_by_hairpin = 0

    # Make all probe possibilities
    candidate_probe_list = []
    for location in range(0,len(target_transcript_seq)-probe_len + 1):
        probe_region = target_transcript_seq[location:location+probe_len]
        probe_seq = reverse_Complement(probe_region)
        probe_name = "probe_" + str(location+1)
        candidate_probe_list.append((probe_name, probe_seq))

    # Filter out if Tm is too low
    Tm_filtered_probes = []
    for candidate_probe in candidate_probe_list:
        probe_name = candidate_probe[0]
        probe_seq = candidate_probe[1]
        probe_tm = 273 + primer3.calcTm(probe_seq)
        if float(probe_tm) > Tm:
            Tm_filtered_probes.append((probe_name, probe_seq, probe_tm))
        else:
            filtered_by_tm += 1
    if len(Tm_filtered_probes) < num_probes:
        raise Exception("Error: Too many probes filtered out due to low probe binding Tm, consider lowering the minimum Tm in the config.py file and rerunning.")

    Tm_rep_filtered_probes = []
    if filter_repeats == True:
        for candidate_probe in Tm_filtered_probes:
            probe_name = candidate_probe[0]
            probe_seq = candidate_probe[1] 
            probe_tm = candidate_probe[2]
            if ('AAAAA' or 'TTTTT' or 'GGGG' or 'CCCC') not in probe_seq:
                Tm_rep_filtered_probes.append((probe_name,probe_seq,probe_tm))
            else:
                filtered_by_repeat += 1
    if len(Tm_rep_filtered_probes) < num_probes:
        raise Exception("Error: Too many probes filtered out due to repetitive nucleotides (ie. AAAAA or GGGG), consider changing \"filter_repeats = True\" to \"filter_repeats = False\" in the config.py file and rerunning.")

    # Filter out if high hairpin Tm
    Tm_rep_hairpin_filtered_probes = []
    for candidate_probe in Tm_rep_filtered_probes:
        probe_name = candidate_probe[0]
        probe_seq = candidate_probe[1]
        probe_tm = candidate_probe[2]
        probe_hairpin_tm = primer3.calcHairpinTm(probe_seq) + 273
        if float(probe_hairpin_tm) < Tm_hairpin:
            Tm_rep_hairpin_filtered_probes.append((probe_name,probe_seq,probe_tm,probe_hairpin_tm))
            probe_Tm_dict[probe_name] = (probe_seq, probe_tm)
        else:
            filtered_by_hairpin += 1
    if len(Tm_rep_hairpin_filtered_probes) < num_probes:
        raise Exception("Error: Too many probes filtered out due to repetitive nucleotides (ie. AAAAA or GGGG), consider changing \"filter_repeats = True\" to \"filter_repeats = False\" in the config.py file and rerunning.")
    
    print("\nTm filtered: "+str(filtered_by_tm))
    print("Repeat filtered: "+str(filtered_by_repeat))
    print("Hairpin filtered: "+str(filtered_by_hairpin))
        
    return(Tm_rep_hairpin_filtered_probes, probe_Tm_dict)



##############################################################################
#  Scoring functions                                                         #
##############################################################################

def compare(i1, i2):
    if i1[1] < i2[1]:
        return(-1)
    if i1[1] > i2[1]:
        return(1)
    if i1[1] == i2[1] and i1[0] > i2[0]:
        return(1)
    if i1[1] == i2[1] and i1[0] < i2[0]:
        return(-1)

def compute_offtarget_scores_matches(sam_data, target_ensembl_id, probe_len):
    print('Calculating probe candidate off-target scores\n')
    offtarget_match_scores = {}

    for aln in sam_data:
        probe_name = aln[0]
        # Pull ensemble transcript ID from name
        aligned_ensembl_transcript_id = aln[2].split(".")[0]
        if aligned_ensembl_transcript_id == "*":
            continue
            # This skips if the probe didn't align to anything
        if aligned_ensembl_transcript_id == target_ensembl_id:
            continue                                                
            # Doesn't add anything to the probe score if the alignment is an on target alignment

        match_score = int(aln[11].split(':')[2])/(2*probe_len)+1     # divided by the maximum alignment score. Then add 1 to punish for aligning at all
        # Note that this normalization only works for the --local mode and must be adjusted for --end-to-end mode which has a different maximum alignment score.

        # Punishes 20x for rRNA alignment. Should make it so these probes aren't selected
        if aligned_ensembl_transcript_id.startswith('rRNA'):
            match_score = match_score + int(aln[11].split(':')[2])/(2*probe_len)+1000000

        if probe_name not in offtarget_match_scores:
            offtarget_match_scores[probe_name] = 0
            # Adds entry to offtarget scores if probe isn't already in it
        
        # Add the score of mismatches and indels to the off target score for that probe
        offtarget_match_scores[probe_name] = match_score + offtarget_match_scores[probe_name]

    offTarget_match_scoresList = list(offtarget_match_scores.items())
     
    keyfun =  functools.cmp_to_key(compare)
    offTarget_match_scoresList.sort(key=keyfun)

    return(offTarget_match_scoresList)

def compute_offtarget_scores_tm(sam_data, target_ensembl_id, transcriptome_dict, probe_len):
    print('Calculating probe candidate off-target scores\n')
    offtarget_TM_scores = {}
    for aln in sam_data:
        probe_name = aln[0]
        if probe_name not in offtarget_TM_scores:
            offtarget_TM_scores[probe_name] = 0
        probe_seq = reverse_Complement(aln[9])
        alignment_pos = int(aln[3]) - 1 # since it's one-based counting
        # Pull ensemble transcript ID from name      
        aligned_ensembl_transcript_id = aln[2].split(".")[0]
        # Adds entry to offtarget scores if probe isn't already in it
        if aligned_ensembl_transcript_id == "*":
            continue
            # This skips if the probe didn't align to anything
        if aligned_ensembl_transcript_id == target_ensembl_id:
            continue                                                
            # Doesn't add anything to the probe score if the alignment is an on target alignment
        aligned_transcript_seq = transcriptome_dict[aligned_ensembl_transcript_id]

        CIGAR_score = aln[5]
        split_CIGAR = re.split('(\d+)',CIGAR_score)[1:]
        dict_CIGAR = {}
        dict_CIGAR['S1'] = 0
        dict_CIGAR['S2'] = 0
        dict_CIGAR['I'] = 0
        dict_CIGAR['D'] = 0
        dict_CIGAR['M'] = 0
        for place in range(len(split_CIGAR)):
            if split_CIGAR[place] == 'I':
                dict_CIGAR['I'] = dict_CIGAR['I'] + int(split_CIGAR[place-1])
            if split_CIGAR[place] == 'M':
                dict_CIGAR['M'] = dict_CIGAR['M'] + int(split_CIGAR[place-1])
            if split_CIGAR[place] == 'D':
                dict_CIGAR['D'] = dict_CIGAR['D'] + int(split_CIGAR[place-1])
            if split_CIGAR[place] == 'S' and place == 1:
                dict_CIGAR['S1'] = int(split_CIGAR[place-1])
            if split_CIGAR[place] == 'S' and place == len(split_CIGAR)-1:
                dict_CIGAR['S2'] = int(split_CIGAR[place-1])

        # Use position to pull out the fragment of the transcript that aligns to the probe
        # accounts for insertions and deletions as well as soft clipping, then 2 nucleotides 
        # are added on each end for a more precise estimate of binding
        aligned_transcript_seq_fragment = aligned_transcript_seq[alignment_pos-dict_CIGAR['S1']-2:alignment_pos+probe_len+dict_CIGAR['S2']+dict_CIGAR['I']-dict_CIGAR['D']+2]
        if aligned_transcript_seq_fragment == '':
            print('\n')
            print(aln)
            print(alignment_pos-dict_CIGAR['S1']-2)
            print(alignment_pos+probe_len+dict_CIGAR['S2']+dict_CIGAR['I']-dict_CIGAR['D']+2)
            print(aligned_transcript_seq)
              
        aln_tm = (primer3.calcHeterodimer(probe_seq, aligned_transcript_seq_fragment).tm + 273)/(primer3.calcTm(probe_seq) + 273) + 1
        if aligned_ensembl_transcript_id.startswith('rRNA'):
            aln_tm = 1000000*aln_tm
        offtarget_TM_scores[probe_name] = aln_tm + offtarget_TM_scores[probe_name]

    offTarget_TM_scoresList = list(offtarget_TM_scores.items())
    #offTarget_TM_scoresList.sort(key=lambda p: p[1], reverse=False)
    
    keyfun =  functools.cmp_to_key(compare)
    offTarget_TM_scoresList.sort(key=keyfun)

    return(offTarget_TM_scoresList)
        


##############################################################################
#  Bowtie2 functions                                                         #
##############################################################################
def make_bowtie2_index(fasta_file, bowtie2_path, threads):
    print("\nBowtie2 will now generate the index with which to align the candidate probe sequences\nUsing",threads,"threads.")
    bt2_run = subprocess.run([bowtie2_path +'-build',
                                '--threads',str(threads),
                                '-f', fasta_file, 
                                'bt2_index'],
                                capture_output=True)
    return(bt2_run)   

# Align the probes with the gencode transcriptome
def align_bowtie2(testprobeFile, transcriptomeIndex, bowtie2_path, threads, seed_length):
    print('Aligning probe candidates to transcriptome using bowtie2\n')
    bt2_run = subprocess.run([bowtie2_path, '-f', '-a',
                              '-x', transcriptomeIndex,
                              '-U', testprobeFile, 
                              '-p', str(threads), 
                              '-L', str(seed_length),
                              '--nofw',
                              '--local',
                              '--score-min', 'L,20,0'],
                             capture_output=True)

    sam_raw = bt2_run.stdout.decode("utf-8").split('\n')    # Converts to string then to a list where each element is a line
    sam_raw.remove('') 
    sam_data = []
    for line in sam_raw:
        if line[0] != '@':                                  # Removes header information which isn't informative
            sam_data.append(line.split('\t'))               # Splits each alignment lines
    return(sam_data, bt2_run)


##############################################################################
# "Greedy" Set Selection Method                                              #
##############################################################################
def select_nonoverlapping_probes(offTarget_TmscoresList, probe_len, num_probes):

    print('Selecting nonoverlapping probes based on scoring rank\n')

    overlapping_probes = []
    final_probe_set = []

    total_score = 0.0
    for probeandScore in offTarget_TmscoresList:
        if probeandScore[0] in overlapping_probes:
            continue
        final_probe_set.append(probeandScore)
        total_score += probeandScore[1]
        probeName_list = probeandScore[0].split("_")
        probe_idx = int(probeName_list[1])
        addingOverlap = probe_idx + probe_len + 2
        subtractingOverlap = probe_idx - (probe_len + 2)
        if (subtractingOverlap < 0): subtractingOverlap = 0
        for bb in range(probe_idx + 1, addingOverlap):
            overlapProbeName = "probe_" + str(bb)
            overlapping_probes.append(overlapProbeName)
        for cc in range(subtractingOverlap, probe_idx):
            overlapProbeName2 = "probe_" + str(cc)
            overlapping_probes.append(overlapProbeName2) 

    final_probe_set.sort(key=lambda p: p[1], reverse=False)

    return(final_probe_set, total_score)


##############################################################################
#  helper functions                                                          #
##############################################################################
def complement(seq):
    compl = {'A':'T','T':'A','C':'G','G':'C'}
    compList = []
    for nt in seq.upper():
        compList.append(compl[nt])
    compSeq = ''.join(compList)
    return(compSeq)

def reverse_Complement(seq):
    compl = {'A':'T','T':'A','C':'G','G':'C'}   
    compList = []
    for nt in seq.upper():
        compList.append(compl[nt])
    compSeq = ''.join(compList)
    return(compSeq[::-1])


##############################################################################
#  optimal set selection with token passing                                  #
##############################################################################
def find_optimal_set(probe_score_list, overlap, output_set_size, final_probes_greedy, beam=5.0):
    
    print('Selecting optimal nonoverlapping probe set based on scoring rank\n')

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


    #STEP 3: Define a token class 
    class token:
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
            t_list.append(token())
        if (aa[1] <= search_score_profile[0] + beam):
            t = token()
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
                newTm = tt.Tm + jj_tm
                newPos = tt.num_probes
           
                #are we going beyond set_size?
                if (newPos >= output_set_size): continue

                #beam pruning
                if (newTm > search_score_profile[newPos] + beam): continue
                if (newTm < search_score_profile[newPos]): search_score_profile[newPos] = newTm
            
                if (t_list_jj[newPos].num_probes == 0) or (t_list_jj[newPos].Tm > newTm): 
                    #Make a new token which extends tt by adding probe jj:
                    x = token()
                    x.back_ptr = ii
                    x.num_probes = tt.num_probes + 1
                    x.Tm = tt.Tm + jj_tm
                    t_list_jj[newPos] = x

            done_work += 1
        if ii%100 == 0: print('optimal probe set search.  done % : ', done_work*100.0/tot_work)
    
    #STEP 6: all token propagation is done! now let's find best probe set of set_size
    bestTm = math.inf
    bestii = -1
    for ii in range(len(token_list)):
        if (token_list[ii][output_set_size - 1].num_probes == 0): continue
        if token_list[ii][output_set_size - 1].Tm < bestTm:
            bestTm = token_list[ii][output_set_size - 1].Tm
            bestii = ii
                
    finaltoken = token_list[bestii][output_set_size - 1] 

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

    return( retL, total_score )

