import re
from collections import defaultdict
import pandas as pd

def protein_tsv_reader(protein_tsv_file):
    protein_list = []
    with open(protein_tsv_file, 'r') as file_open:
        next(file_open)
        for line in file_open:
            line_split = line.split("\t")
            protein_ID = line_split[3]
            protein_list.append(protein_ID)
    return protein_list

def protein_phospho_counting(protein_csv_file):
    phospho_protein_dict = {}
    with open(protein_csv_file, 'r') as file_open:
        next(file_open)
        pattern = re.compile('\d+\w{1}\(79\.9663\)')  # look for phosphorylation
        for line in file_open:
            protein_ID = line.split('\t')[3]
            regex = re.findall(pattern, line)
            for ele in regex:
                if ele != '':
                    phospho_protein_dict[protein_ID] = regex
    return phospho_protein_dict


def info_getter(protein_tsv_file):
    info_dict = {}
    with open(protein_tsv_file, 'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            protein_ID = line_split[3]
            length = int(line_split[6])
            coverage = float(line_split[7])
            total_spec_count = int(line_split[-9])
            total_intensity = line_split[-5]
            info_dict[protein_ID] = (length,coverage,total_spec_count,total_intensity)
    return info_dict


def info_getter1(protein_tsv_file):
    df = pd.read_csv(protein_tsv_file, delimiter='\t', header=0)
    return df

def peptide_counting(peptide_tsv_file):
    peptide_list = []
    with open(peptide_tsv_file, 'r') as file_open:
        next(file_open)
        for line in file_open:
            peptide_seq = line.split("\t")[0]
            peptide_list.append(peptide_seq)
    return peptide_list


def peptide_phospho_reader(peptide_tsv_file, mod=79.9663): # 79.9663 is the delta mass of phosphorylation on STY
    pep_phos_dict = defaultdict()
    with open(peptide_tsv_file) as file_open:
        for i in range(1):
            next(file_open)
        for line in file_open:
            pep_seq = line.split('\t')[0]

            pattern = re.compile('\d+\w{1}\('+str(mod)+'\)')
            regex = re.findall(pattern, line)
            for ele in regex:
                if ele != '':
                    pep_phos_dict[pep_seq]=regex
    return pep_phos_dict

def venn_diagram_gen(dictionary, title=''): # parameter could be a dictionary of proteins or peptides from different samples {'sample1': [], 'sample2': []}
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2, venn3

    value_list_of_sets = [set(l) for l in dictionary.values()]
    sample_name_list = [n for n in dictionary.keys()]
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    if len(dictionary) == 2:  # two samples for venn diagram
        venn2(value_list_of_sets, set_labels=sample_name_list)

    elif len(dictionary) == 3:  # 3 samples for venn diagram
        venn3(value_list_of_sets, set_labels=sample_name_list)

    else:
        print ('Error: only 2 or 3 comparison for venn diagram are accepted in this script.')

    plt.title(title)
    plt.show()

def psm_reader(psm_path):
    pep_spec_count_dict = defaultdict(int)
    ret_pep_dict = {}
    with open(psm_path, 'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            pep_seq = line_split[1]
            retention_time = float(line_split[4])/60  # in minute
            pep_spec_count_dict[pep_seq]+=1
            ret_pep_dict[retention_time] = pep_seq
    return pep_spec_count_dict, ret_pep_dict

def spectra_num_counting(peptide_tsv, psm_tsv, fasta_path, reverse=0):
    import aho_corasick
    import multiprocessing_naive_algorithym
    from protein_coverage import fasta_reader, read_fasta_info_dict2

    protein_matched_pep_dict = defaultdict(list)
    protein_total_spec_count_dict = defaultdict(int)
    if reverse == 0:
        protein_dict = fasta_reader(fasta_path)
    else:
        protein_dict = read_fasta_info_dict2(fasta_path)

    # exclude the contaminant proteins


    ID_list,seq_list = multiprocessing_naive_algorithym.extract_UNID_and_seq(protein_dict)
    seq_line = multiprocessing_naive_algorithym.creat_total_seq_line(seq_list)

    peptide_list = peptide_counting(peptide_tsv)

    # ahocorasick trie building and iterate seq_line with trie structure
    A = aho_corasick.automaton_trie(peptide_list)
    aho_result = aho_corasick.automaton_matching(A, seq_line)

    # get a dictionary that has each seq line position as key, and corresponding protein ID as value
    pos_ID_dict = multiprocessing_naive_algorithym.read_position_ID_into_dict(ID_list,seq_list, seq_line)

    # build up a dictionary with protein as key, matched peptides list as value
    for each in aho_result:
        protein_matched_pep_dict[pos_ID_dict[each[0]]].append(each[2])

    # get the peptide spectrum count dict
    psm_count_dict = psm_reader(psm_tsv)[0]

    # calculate total spectrum count for each protein
    for protein in protein_matched_pep_dict:
        for pep in protein_matched_pep_dict[protein]:
            protein_total_spec_count_dict[protein] += psm_count_dict[pep]

    return protein_total_spec_count_dict

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn3, venn3_circles
    import pandas as pd
    from protein_coverage import read_fasta_info_dict2, fasta_reader
    import seq_operation, aho_corasick
    from multiprocessing_naive_algorithym import creat_pep_ID_dict,creat_ID_pep_dict
    import pickle as ppp
    '''
    SC_1_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_1/protein.tsv'
    SC_1_protein_set = (protein_tsv_reader(SC_1_tsv_path))  

    SC_2_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_2_3/SC_2/protein.tsv'
    SC_2_protein_set = (protein_tsv_reader(SC_2_tsv_path))
    SC_3_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_2_3/SC_3/protein.tsv'
    SC_3_protein_set = (protein_tsv_reader(SC_3_tsv_path))
    SC_combined_set = set(SC_1_protein_set+SC_2_protein_set+SC_3_protein_set)
    print (len(SC_combined_set))
    '''

    fasta_path = 'D:/data/proteome_fasta/uniprot-proteome_UP000005640.fasta'
    protein_dict=fasta_reader(fasta_path)

    pep_tsv = 'D:/data/deep_proteome/20200716/T_5min_search/peptide.tsv'
    prot_tsv = 'D:/data/deep_proteome/20200716/T_5min_search/protein.tsv'
    psm_tsv = 'D:/data/deep_proteome/20200716/T_5min_search/psm.tsv'
    print (len(protein_tsv_reader(prot_tsv)))
    print (len(spectra_num_counting(pep_tsv,psm_tsv,fasta_path,reverse=0)))
    peptide_list = peptide_counting(pep_tsv)
    phos_peptide_list = [key for key in peptide_phospho_reader(pep_tsv)]
    print (phos_peptide_list)

    venn_dict = {'PXD001364': [key for key in ppp.load(open('C:/Users/gao lab computer/PycharmProjects/ms2_spectra/PXD001364_ext_pep_cosScore_dict.p', 'rb'))],
                 'PXD001723':[key for key in ppp.load(open('C:/Users/gao lab computer/PycharmProjects/ms2_spectra/PXD001723_ext_pep_cosScore_dict.p', 'rb'))]}
    venn_diagram_gen(venn_dict,title='C_elegans Ext peptides overlap between two PRIDE projects')
    """
    uni_id_list, seq_list = seq_operation.extract_UNID_and_seq(protein_dict)
    seq_line = seq_operation.creat_total_seq_line(seq_list)
    aho_automaton = aho_corasick.automaton_trie(phos_peptide_list)
    aho_result = aho_corasick.automaton_matching(aho_automaton, seq_line)
    pos_id_dict = seq_operation.read_position_ID_into_dict(uni_id_list, seq_list, seq_line)
    id_pep_dict = creat_ID_pep_dict(aho_result,pos_id_dict)
    print(len(id_pep_dict))
    psm_count_dict = psm_reader(psm_tsv)[0]
    protein_spec_dict = spectra_num_counting(pep_tsv,psm_tsv,fasta_path,reverse=1) 
    phosphoprotein_spec_dict = defaultdict(int)
    for phospho_protein in id_pep_dict:
        for pep in id_pep_dict[phospho_protein]:
            phosphoprotein_spec_dict[phospho_protein]+=psm_count_dict[pep]

    print (phosphoprotein_spec_dict['P28652'],phosphoprotein_spec_dict['P11798'])

    df = pd.DataFrame([])
    """

    '''
    fig, ax = plt.subplots(1,1,figsize=(10,10))
    venn3([SC_1_protein_set, SC_2_protein_set, SC_3_protein_set],set_labels=('SC_1', 'SC_2', 'SC_3'))
    c = venn3_circles([SC_1_protein_set, SC_2_protein_set, SC_3_protein_set],linewidth=1)
    c[0].set_color('magenta')
    c[0].set_edgecolor('none')
    c[0].set_alpha(0.4)
    plt.title('Venn diagram for proteins identified in spinal cord without phosphopeptide enrichment')
    plt.show()
    '''
    '''
    SC_phos_1_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_phos_2nd_search/SC_phos_1/protein.tsv'
    SC_phos_2_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_phos_2nd_search/SC_phos_2/protein.tsv'
    SC_phos_3_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_phos_2nd_search/SC_phos_3/protein.tsv'

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    venn3([set(protein_tsv_reader(SC_phos_1_tsv_path)), set(protein_tsv_reader(SC_phos_2_tsv_path)),
           set(protein_tsv_reader(SC_phos_3_tsv_path))], set_labels=('SC_1', "SC_2", "SC_3"))
    plt.title('Venn diagram for proteins identified in spinal cord after phospho-peptide enrichment')
    plt.show()
    '''

    '''
    SC_1_peptide_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_1/peptide.tsv'
    pep_phos_dict = peptide_phospho_reader(SC_1_peptide_tsv_path)
    print (len(pep_phos_dict))
    '''