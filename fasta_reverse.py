from protein_coverage import fasta_reader, read_description
import numpy as np
from MS_tools_in_silico_digestion import peptide_generator
from MS_tools_parameters import max_len, min_len, missed
import re


def my_replace(match_obj):
    match_obj = match_obj.group()
    return match_obj[::-1]  # gives back the matched object backwards

def fasta_reverse_generator(uniprot_id_list, fasta_file_in, fasta_file_out, reverse_algorithm='null'):
    print ('reverse_algorithm = null, or other protease')
    # read protein sequence into dic
    protein_dict = fasta_reader(fasta_file_in)

    # read description into dic
    ID_descrip_dict = read_description(fasta_file_in)  # dictionary format: {'ID': ('sp'or'tr', description)}

    # write id and reverse sequence into fasta_file_out
    with open(fasta_file_out, 'w', newline='\n') as file_open:
        for id in uniprot_id_list:
            forward_seq = protein_dict[id]
            if reverse_algorithm == 'null':  # reverse all aa in protein sequence
                rev_seq = protein_dict[id][::-1]

            else:  # only reverse digested peptides and remain the original sequence as same
                protein_seq_str = protein_dict[id]

                peptide_list = list(peptide_generator((id, protein_seq_str), reverse_algorithm)[id])
                new_peptide_list = [pep for pep in peptide_list if min_len<= len(pep) <= max_len]  # filter the peptides based on min_len and max_len

                for pep in new_peptide_list:
                    protein_seq_str = re.sub(pep, my_replace, protein_seq_str)  # replace all digested peptides backward
                rev_seq = protein_seq_str

            block = range(0,len(forward_seq)+60,60)

            # write forward
            file_open.write('>'+ ID_descrip_dict[id][0]+'|'+id+'|'+ID_descrip_dict[id][1]+'\n')
            for i in range(len(block)-1):
                file_open.write(forward_seq[block[i]:block[i+1]]+'\n')

            # write reverse
            file_open.write('>Rev_'+ ID_descrip_dict[id][0]+'|'+id+'|'+ID_descrip_dict[id][1]+'\n')
            for i in range(len(block)-1):
                file_open.write(rev_seq[block[i]:block[i+1]]+'\n')
    return fasta_file_out

fasta_file_input = 'D:/data/ext_evo_pj/matched_protein_extended_dbKRW_6_23.fasta'
fasta_file_out = 'D:/data/ext_evo_pj/matched_protein_extended_dbKRW_6_23_rev.fasta'
ID_list = [i for i in fasta_reader(fasta_file_input).keys()]
fasta_reverse_generator(ID_list, fasta_file_input, fasta_file_out, 'null')