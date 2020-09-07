"""
add up all PTMs found from peptide.tsv file and return all PTM for a given protein
"""
import re
import numpy as np
from collections import defaultdict
from tsv_reader import peptide_phospho_reader


def ptm_reader(protein_seq_dict,uniprot_id,peptide_tsv,mod=79.9663):
    """
    main function to combine PTM from peptide.tsv onto protein seq
    :param protein_seq_dict:
    :param uniprot_id:
    :param peptide_tsv:
    :param mod:
    :return: a dictionary with protein id as key, a list of modification locations as value
    """

    protein_ptm_loc_dict = defaultdict(list)

    pep_w_mod_dict = peptide_phospho_reader(peptide_tsv,mod)
    protein_seq = protein_seq_dict[uniprot_id]
    prot_ptm_loc_list = []

    # check each peptide with mod
    for each in pep_w_mod_dict:

        # if peptide with mod is found in the protein
        if each in protein_seq:
            pep_loc = protein_seq.find(each)
            # positive lookahead for each ptm found from peptide.tsv
            for each_mod in pep_w_mod_dict[each]:
                ptm_loc_list = [int(each_loc)+pep_loc for each_loc in re.findall(r'\d+(?=[A-Z])',each_mod)]
                prot_ptm_loc_list+=ptm_loc_list
            protein_ptm_loc_dict[uniprot_id] = prot_ptm_loc_list

        else:
            continue
    return protein_ptm_loc_dict


def uniprot_id_matched_pep_getter(protein_dict,pep_list):
    """
    get a uniprot id with list of matched peptides dictionary
    :param protein_dict:
    :param pep_list:
    :return:
    """
    import seq_operation, aho_corasick
    from multiprocessing_naive_algorithym import creat_ID_pep_dict

    uni_id_list, seq_list = seq_operation.extract_UNID_and_seq(protein_dict)
    seq_line = seq_operation.creat_total_seq_line(seq_list)
    aho_automaton = aho_corasick.automaton_trie(pep_list)
    aho_result = aho_corasick.automaton_matching(aho_automaton, seq_line)
    pos_id_dict = seq_operation.read_position_ID_into_dict(uni_id_list, seq_list, seq_line)
    id_pep_dict = creat_ID_pep_dict(aho_result, pos_id_dict)

    return id_pep_dict


def freq_array_generator(peptide_list, protein_seq_string):
    """

    :param peptide_list: matched psm list from aho
    :param protein_seq_string: one protein sequence
    :return:
    """
    freq_array = np.zeros(len(protein_seq_string))

    # calculation
    for pep in peptide_list:

        start_pos = protein_seq_string.find(pep)
        end_pos = start_pos + len(pep) -1
        freq_array[start_pos:end_pos + 1] += 1

    return freq_array


def ptm_sites_counting(matched_pep_list,peptide_tsv,psm_tsv,mod,regex=r'\w\[\d+\]'):
    """
    count the ptm sites in one protein given matched peptides
    :param matched_pep_list: peptide list that mapped to one protein seq from aho
    :param peptide_tsv:
    :param psm_tsv:
    :param mod:
    :return:
    """
    from tsv_reader import pep_mod_pep_dict_gen

    PTM_sites_counting = defaultdict(int)
    pep_mod_pep_dict = pep_mod_pep_dict_gen(psm_tsv,fragpipe_ver=13.0)
    pep_w_mod_dict = peptide_phospho_reader(peptide_tsv, mod)
    # changing matched peptide list to all peptide list and a condition would be faster
    for each_pep in matched_pep_list:
        # if peptide has modification

        if each_pep in pep_w_mod_dict:
            aa_w_mods_list = np.array([re.findall(regex,each_mod_pep)
                                        for each_mod_pep in pep_mod_pep_dict[each_pep]]).flatten()
            print(aa_w_mods_list)
            for each_aa in aa_w_mods_list:
                PTM_sites_counting[each_aa]+=1
        else:
            continue
    return PTM_sites_counting


#def html_vis(matched_pep_list,protein_seq,uniprot_id,gene_name,output_path,PTM_dict={})

if __name__=='__main__':
    from protein_coverage import fasta_reader,read_fasta_info_dict2
    from tsv_reader import pep_mod_pep_dict_gen,peptide_counting
    fasta_file = 'D:/data/proteome_fasta/UniProt_Mouse_04-17-2019_reversed.fasta'
    pep_tsv = 'D:/data/phospho_wang/2020-09-06/result/B_phos/peptide.tsv'
    psm_tsv = 'D:/data/phospho_wang/2020-09-06/result/B_phos/psm.tsv'
    peptide_list = peptide_counting(pep_tsv)
    protein_dict = read_fasta_info_dict2(fasta_file)
    print (ptm_reader(protein_dict,'P28652',pep_tsv))
    print (pep_mod_pep_dict_gen(psm_tsv)['NSSAITSPK'])

    matched_pep_list = uniprot_id_matched_pep_getter(protein_dict,peptide_list)['P28652']
    print (ptm_sites_counting(matched_pep_list,pep_tsv,psm_tsv,mod=79.9663))