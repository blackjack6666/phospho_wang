def fasta_reader(fasta_file_path):
    protein_seq_dict = {}
    with open(fasta_file_path, 'r') as file_open:
        file_split = file_open.read().split('\n>')
        for each in file_split:
            split_line = each.split('\n')
            protein_seq_dict[split_line[0].split('|')[1]] = ''.join(split_line[1:])
    return protein_seq_dict

def read_fasta_info_dict2(filename):
    """
    -----
    read fasta file including reverse seq, exclude rev seq in returned dictionary
    -----
    :param filename:
    :return:
    """
    protein_seq_dict = {}
    with open(filename, 'r') as file_open:
        file_split = file_open.read().split('\n>')
        # read first entry
        protein_seq_dict[file_split[0].split('\n')[0].split('|')[1]] = ''.join(file_split[0].split('\n')[1:])
        # read the rest entry
        for each in file_split[1:]:
            if each.startswith('Rev'):
                continue
            else:
                split_line = each.split('\n')
                protein_seq_dict[split_line[0].split('|')[1]] = ''.join(split_line[1:])

    return protein_seq_dict

def read_description_into_dict(fasta_file):
    ID_description_dict = {}
    with open(fasta_file, 'r') as file_open:
        for line in file_open:
            if line.startswith('>'):
                ID_description_dict[line.split('|')[1]] = line.split('|')[2].rstrip('\n')
    return ID_description_dict

def read_description(fasta_file):
    """
       read description and prefix in fasta file
       :param fasta_file:
       :return:
       """

    # dictionary format: {'ID': ('sp', 'description')}
    ID_description_dict = {}
    with open(fasta_file, 'r') as file_open:
        file_split = file_open.read().split('\n>')
        # first entry special case
        ID_description_dict[file_split[0].split('\n')[0].split('|')[1]] = (
        file_split[0].split('\n')[0].split('|')[0][1:], file_split[0].split('\n')[0].split('|')[2])
        for each in file_split[1:]:
            split_line = each.split('\n')
            ID_description_dict[split_line[0].split('|')[1]] = (
            split_line[0].split('|')[0], split_line[0].split('|')[2])
    return ID_description_dict


def read_fasta_into_dict(filename):

    # only read forward sequence, do not include reverse seq
    dict_human_ref = dict()
    #key = ""
    value = ""
    aa_frequency_all_zero_dict = {}
    with open(filename, 'r') as file_open:
        Reverse = 0
        for line in file_open:
            if line.startswith('>sp') and value == '':
                key = line.split('|')[1]
                Reverse = 0
                #dict_human_ref[key] = value
                #value = ""
                #key = line.split('|')[1]
            elif line.startswith('>sp') and value != '':
                dict_human_ref[key] = value
                key = line.split('|')[1]
                value = ''
                Reverse = 0
            elif line.startswith('>Rev'):
                Reverse = 1
            elif Reverse == 0:
                value += line.rstrip('\n')

    if key not in dict_human_ref.keys():
        dict_human_ref[key] = value
    for ID in [ID for ID in dict_human_ref.keys()][1:]:
        sequence = dict_human_ref[ID]
        aa_frequency_all_zero_dict[ID] = [0]*len(sequence)


    return dict_human_ref, aa_frequency_all_zero_dict
if __name__=='__main__':
    fasta_path = 'D:/data/ext_evo_pj/c_elegans_gbff_reference.fasta'
    protein_dict = fasta_reader(fasta_path)
    print (len(protein_dict))
