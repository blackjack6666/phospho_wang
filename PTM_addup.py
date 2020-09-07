"""
add up all PTMs found from peptide.tsv file and return all PTM for a given protein
"""
import re
from collections import defaultdict


def ptm_reader(protein_seq_dict,uniprot_id,peptide_tsv,mod=79.9663):
    """
    main function
    :param protein_seq_dict:
    :param uniprot_id:
    :param peptide_tsv:
    :param mod:
    :return: a dictionary with protein id as key, a list of modification locations as value
    """
    from tsv_reader import peptide_phospho_reader

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


if __name__=='__main__':
    from protein_coverage import fasta_reader,read_fasta_info_dict2
    fasta_file = 'D:/data/proteome_fasta/UniProt_Mouse_04-17-2019_reversed.fasta'
    pep_tsv = 'D:/data/phospho_wang/2020-09-06/result/B_phos/peptide.tsv'
    protein_dict = read_fasta_info_dict2(fasta_file)
    print (ptm_reader(protein_dict,'P28652',pep_tsv))