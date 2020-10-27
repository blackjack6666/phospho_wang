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
                ptm_loc_list = [int(each_loc)+pep_loc-1 for each_loc in re.findall(r'\d+(?=[A-Z])',each_mod)]
                prot_ptm_loc_list += ptm_loc_list
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
            for each_ele in aa_w_mods_list:
                try:
                    PTM_sites_counting[each_ele]+=1
                except TypeError: # some elements in aa_w_mod_list are list
                    for each_aa in each_ele:
                        PTM_sites_counting[each_aa]+=1
        else:
            continue
    return PTM_sites_counting


def gen_cov_graph(peptide_tsv,
                  psm_tsv,
                  fasta_path,
                  Uniprot_name,
                  Gene_name,
                  output_filename,
                  PTM_dict=None,
                  fasta_rev=0,
                  mod=79.9663):
    """

    :param peptide_tsv:
    :param psm_tsv:
    :param fasta_path:
    :param Uniprot_name:
    :param Gene_name:
    :param output_filename:
    :param PTM_dict: example: {"S[167]":"Serine phosphorylation",
                            "T[181]":"Threonine phosphorylation",
                            'Y[243]':'Tyrosine phosphorylation'}
    :param fasta_rev:
    :return:
    """
    from tsv_reader import psm_reader, peptide_counting
    from protein_coverage import fasta_reader, read_fasta_info_dict2

    if fasta_rev == 0:
        uniprot_dict = fasta_reader(fasta_path)
    else:
        uniprot_dict = read_fasta_info_dict2(fasta_path)

    pep_list = peptide_counting(peptide_tsv)
    psm_dict = psm_reader(psm_tsv)[0]

    string_of_full_protein_sequence = uniprot_dict[Uniprot_name]

    matched_list_of_peptides = uniprot_id_matched_pep_getter(uniprot_dict, pep_list)[Uniprot_name]
    matched_list_of_psms = [each_pep for each_pep in matched_list_of_peptides for i in range(psm_dict[each_pep])]

    aa_frequency, PTM_loc_list, PTM_sites_counting = freq_array_generator(matched_list_of_psms,
                                                                          string_of_full_protein_sequence), \
                                                     ptm_reader(uniprot_dict, Uniprot_name, peptide_tsv,mod=mod)[Uniprot_name], \
                                                     ptm_sites_counting(matched_list_of_peptides, peptide_tsv, psm_tsv,
                                                                        mod=mod)
    print (PTM_sites_counting)
    # some required variables from other function


    # coverage calculation
    non_zero = np.count_nonzero(aa_frequency)
    coverage = float(non_zero)/len(aa_frequency)*100

    # split sequence
    #split_seq = range(0, len(string_of_full_protein_sequence)+50, 50) # 50 amino acids per line on html file
    split_seq = np.arange(0,len(string_of_full_protein_sequence),50)
    split_seq = np.append(split_seq,len(string_of_full_protein_sequence))

    # required stat for color labeling
    max_freq = np.max(aa_frequency)

    # HTML contents, color markers
    content = '<!DOCTYPE HTML>\n<html>\n<body>\n<style>\n'
    div = '#outPopUp { \n\tmargin-right: auto;\n\tmargin-left: auto;\n\twidth: 500px;\n\t}\n'
    mark = 'mark {\n\tbackground-color: #2C3E50;\n\tcolor: white;\n}\n'
    mark += 'mark1 {\n\tbackground-color: #2471A3;\n\tcolor: white;\n}\n'
    mark += 'mark2 {\n\tbackground-color: #AED6F1;\n\tcolor: white;\n}\n'
    mark += 'mark3 {\n\tbackground-color: #95A5A6;\n\tcolor: white;\n}\n'
    mark += 'mark4 {\n\tbackground-color: #D7DBDD;\n\tcolor: white;\n}\n'
    mark += 'mark5 {\n\tbackground-color: white;\n\tcolor: red;\n}\n'

    # HTML content, protein info and color legend
    protein_info = '<p>The protein\t' + Uniprot_name + '(Gene name:' + Gene_name + ')' + '\thas sequence coverage of: '
    note = '<p><font size="2">Note*<br>frequency: <mark>color</mark>><mark1>color</mark1>><mark2>color</mark2>><mark3>' \
           'color</mark3>><mark4>color</mark4></font></p>\n'

    # start writing html script
    with open(output_filename, 'w') as f:
        f.write(content+div+mark+'</style>\n\t<div id="outPopUp">\n'+protein_info+'<strong>'+str(coverage)+'%</strong></p>\n'+note+'\n\n<pre>\n\n\n')

        for i in range(len(split_seq)-1):
            for j in range(split_seq[i],split_seq[i+1]):
                if j not in PTM_loc_list:

                    if aa_frequency[j] == 0:
                        f.write(string_of_full_protein_sequence[j])
                    elif 1 <= aa_frequency[j] < 0.2*max_freq:
                        f.write('<mark4>'+string_of_full_protein_sequence[j]+'</mark4>')
                    elif 0.2*max_freq <= aa_frequency[j] < 0.4*max_freq:
                        f.write('<mark3>'+string_of_full_protein_sequence[j]+'</mark3>')
                    elif 0.4*max_freq <= aa_frequency[j] < 0.6*max_freq: # color legend changeable
                        f.write('<mark2>'+string_of_full_protein_sequence[j]+'</mark2>')
                    elif 0.6*max_freq <= aa_frequency[j] < 0.8*max_freq: # color legend changeable
                        f.write('<mark1>'+string_of_full_protein_sequence[j]+'</mark1>')
                    else:
                        f.write('<mark>'+string_of_full_protein_sequence[j]+'</mark>')
                else:
                    f.write('<mark5>'+string_of_full_protein_sequence[j]+'</mark5>')

            f.write('\n')
        f.write('\n\n</pre>\n')
        if PTM_dict:
            f.write('<pre>PTM sites stats:\n')
            for PTM in PTM_dict:
                number_of_each_site = '%s' % PTM_sites_counting[PTM]
                f.write(PTM+'\t'+PTM_dict[PTM]+':<strong>'+ number_of_each_site+'</strong> times identified\n')
        f.write('</pre>\n</div>\n</body>\n</html>\n')
    f.close()
    return output_filename


def gen_cov_graph2(pep_list,
                  psm_dict,
                  fasta_path,
                  Uniprot_name,
                  Gene_name,
                  output_filename,
                  PTM_dict=None,
                  fasta_rev=0,
                  mod=79.9663):
    """
    same as gen_cov_graph, change file path input to pep list and psm reader, no ptm visul function
    :param peptide_tsv:
    :param psm_tsv:
    :param fasta_path:
    :param Uniprot_name:
    :param Gene_name:
    :param output_filename:
    :param PTM_dict:
    :param fasta_rev:
    :param mod:
    :return:
    """
    from tsv_reader import psm_reader, peptide_counting
    from protein_coverage import fasta_reader, read_fasta_info_dict2

    if fasta_rev == 0:
        uniprot_dict = fasta_reader(fasta_path)
    else:
        uniprot_dict = read_fasta_info_dict2(fasta_path)


    string_of_full_protein_sequence = uniprot_dict[Uniprot_name]

    matched_list_of_peptides = uniprot_id_matched_pep_getter(uniprot_dict, pep_list)[Uniprot_name]
    matched_list_of_psms = [each_pep for each_pep in matched_list_of_peptides for i in range(psm_dict[each_pep])]

    aa_frequency = freq_array_generator(matched_list_of_psms,string_of_full_protein_sequence)
    # aa_frequency, PTM_loc_list, PTM_sites_counting = freq_array_generator(matched_list_of_psms,
    #                                                                       string_of_full_protein_sequence), \
    #                                                  ptm_reader(uniprot_dict, Uniprot_name, peptide_tsv, mod=mod)[
    #                                                      Uniprot_name], \
    #                                                  ptm_sites_counting(matched_list_of_peptides, peptide_tsv, psm_tsv,
    #                                                                     mod=mod)

    # some required variables from other function

    # coverage calculation
    non_zero = np.count_nonzero(aa_frequency)
    coverage = float(non_zero) / len(aa_frequency) * 100

    # split sequence
    # split_seq = range(0, len(string_of_full_protein_sequence)+50, 50) # 50 amino acids per line on html file
    split_seq = np.arange(0, len(string_of_full_protein_sequence), 50)
    split_seq = np.append(split_seq, len(string_of_full_protein_sequence))

    # required stat for color labeling
    max_freq = np.max(aa_frequency)

    # HTML contents, color markers
    content = '<!DOCTYPE HTML>\n<html>\n<body>\n<style>\n'
    div = '#outPopUp { \n\tmargin-right: auto;\n\tmargin-left: auto;\n\twidth: 500px;\n\t}\n'
    mark = 'mark {\n\tbackground-color: #2C3E50;\n\tcolor: white;\n}\n'
    mark += 'mark1 {\n\tbackground-color: #2471A3;\n\tcolor: white;\n}\n'
    mark += 'mark2 {\n\tbackground-color: #AED6F1;\n\tcolor: white;\n}\n'
    mark += 'mark3 {\n\tbackground-color: #95A5A6;\n\tcolor: white;\n}\n'
    mark += 'mark4 {\n\tbackground-color: #D7DBDD;\n\tcolor: white;\n}\n'
    mark += 'mark5 {\n\tbackground-color: white;\n\tcolor: red;\n}\n'

    # HTML content, protein info and color legend
    protein_info = '<p>The protein\t' + Uniprot_name + '(Gene name:' + Gene_name + ')' + '\thas sequence coverage of: '
    note = '<p><font size="2">Note*<br>frequency: <mark>color</mark>><mark1>color</mark1>><mark2>color</mark2>><mark3>' \
           'color</mark3>><mark4>color</mark4></font></p>\n'

    # start writing html script
    with open(output_filename, 'w') as f:
        f.write(content + div + mark + '</style>\n\t<div id="outPopUp">\n' + protein_info + '<strong>' + str(
            coverage) + '%</strong></p>\n' + note + '\n\n<pre>\n\n\n')

        for i in range(len(split_seq) - 1):
            for j in range(split_seq[i], split_seq[i + 1]):
                # if j not in PTM_loc_list:

                if aa_frequency[j] == 0:
                    f.write(string_of_full_protein_sequence[j])
                elif 1 <= aa_frequency[j] < 0.2 * max_freq:
                    f.write('<mark4>' + string_of_full_protein_sequence[j] + '</mark4>')
                elif 0.2 * max_freq <= aa_frequency[j] < 0.4 * max_freq:
                    f.write('<mark3>' + string_of_full_protein_sequence[j] + '</mark3>')
                elif 0.4 * max_freq <= aa_frequency[j] < 0.6 * max_freq:  # color legend changeable
                    f.write('<mark2>' + string_of_full_protein_sequence[j] + '</mark2>')
                elif 0.6 * max_freq <= aa_frequency[j] < 0.8 * max_freq:  # color legend changeable
                    f.write('<mark1>' + string_of_full_protein_sequence[j] + '</mark1>')
                else:
                    f.write('<mark>' + string_of_full_protein_sequence[j] + '</mark>')
            # else:
                    #f.write('<mark5>' + string_of_full_protein_sequence[j] + '</mark5>')

            f.write('\n')
        f.write('\n\n</pre>\n')
        # if PTM_dict:
        #     f.write('<pre>PTM sites stats:\n')
        #     for PTM in PTM_dict:
        #         number_of_each_site = '%s' % PTM_sites_counting[PTM]
        #         f.write(PTM + '\t' + PTM_dict[PTM] + ':<strong>' + number_of_each_site + '</strong> times identified\n')
        f.write('</pre>\n</div>\n</body>\n</html>\n')
    f.close()
    return output_filename

if __name__=='__main__':
    from protein_coverage import fasta_reader,read_fasta_info_dict2
    from tsv_reader import pep_mod_pep_dict_gen,peptide_counting, psm_reader
    fasta_file = 'D:/data/proteome_fasta/uniprot-proteome_UP000005640.fasta'
    tryp_pep_tsv = 'D:/data/deep_proteome/20200915_tryp_37C_240min/peptide.tsv'
    tryp_psm_tsv = 'D:/data/deep_proteome/20200915_tryp_37C_240min/psm.tsv'

    chymo_pep_tsv = 'D:/data/deep_proteome/20200915_ct_37C_240min/peptide.tsv'
    chymo_psm_tsv = 'D:/data/deep_proteome/20200915_ct_37C_240min/psm.tsv'
    tryp_pep_list,ct_pep_list = peptide_counting(tryp_pep_tsv),peptide_counting(chymo_pep_tsv)
    combined_pep_list = tryp_pep_list+ct_pep_list
    psm_dict = psm_reader(tryp_psm_tsv)[0]
    psm_dict.update(psm_reader(chymo_psm_tsv)[0])

    gen_cov_graph2(combined_pep_list,psm_dict,fasta_file,'P05386','RPLP1','P05386_comb_240min.html')

    #print (ptm_reader(protein_dict,'P28652',pep_tsv))
    #print (pep_mod_pep_dict_gen(psm_tsv)['NSSAITSPK'])

    # matched_pep_list = uniprot_id_matched_pep_getter(protein_dict,peptide_list)['P28652']
    #print (ptm_sites_counting(matched_pep_list,pep_tsv,psm_tsv,mod=79.9663))

    # pep_tsv = 'D:/data/phospho_wang/2020-09-06/result/B_phos/peptide.tsv'
    # psm_tsv = 'D:/data/phospho_wang/2020-09-06/result/B_phos/psm.tsv'
    #
    # gen_cov_graph(pep_tsv,
    #               psm_tsv,
    #               fasta_file,
    #               'P11798',
    #               'Camk2a',
    #               'P11798.html',
    #               PTM_dict={"S[167]":"Serine phosphorylation","T[181]":"Threonine phosphorylation",'Y[243]':'Tyrosine phosphorylation'},
    #               fasta_rev=0)