import numpy as np
import re


def extract_UNID_and_seq(protein_dict):
    UNID_list = [key for key in protein_dict.keys()]
    seq_list = [value for value in protein_dict.values()]
    return UNID_list, seq_list


def read_seq_length_into_nparray(seq_list):
    seq_length_array = np.array([])
    for each_seq in seq_list:
        seq_length_array = np.append(seq_length_array, len(each_seq))
    return seq_length_array


def creat_total_seq_line(seq_list):
    seq_line = '|'.join(seq_list)
    return seq_line


def zero_line_for_seq(seq_line):
    zero_line = np.zeros(len(seq_line))
    return zero_line


# the following function create a dictionary that read each position in long sequence line as key, corresponding UNIPORT ID as value.
def read_position_ID_into_dict(UNID_list, seq_list, zero_line):
    m = 0
    j = 0
    seq_line_ID_dict = dict()
    for i in range(len(zero_line)):
        if j < len(seq_list[m]):
            seq_line_ID_dict[i] = UNID_list[m]
            j += 1
        else:
            j = 0

            m += 1
    return seq_line_ID_dict


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def separator_pos(seq_line):
    sep_pos_array = np.array([m.start() for m in re.finditer('\|', seq_line)])
    sep_pos_array = np.insert(sep_pos_array, 0, 0)
    sep_pos_array = np.append(sep_pos_array, len(seq_line))
    return sep_pos_array

