import numpy as np
from collections import defaultdict
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


def find_all(a_str, sub):  # a generator
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)  # use start += 1 to find overlapping matches


def find_pep_start(param):
    seq_line, peptide_list = param
    res_dict = {}
    for peptide in peptide_list:
        res_dict[peptide] = [m for m in find_all(seq_line, peptide)]
    return res_dict


# the following function returns a dictionary with each peptide as key and corresponding start position list as value.
def start_end_pos_dict(res_dicts):
    start_pos_dict = {}
    end_pos_dict = {}
    for res_dict in res_dicts:
        for peptide in res_dict:
            start_pos_dict[peptide] = res_dict[peptide]
            end_pos_dict[peptide] = [i + len(peptide) - 1 for i in res_dict[peptide]]
    return start_pos_dict, end_pos_dict


# the following function returns a number_line after matching peptides to long_seq_line.
def adding_numbers_to_zero_line(zero_line, start_dict, end_dict):
    for peptide in start_dict:
        start_list_for_each = start_dict[peptide]
        end_list_for_each = end_dict[peptide]
        for i, j in zip(start_list_for_each, end_list_for_each):
            zero_line[i:j] += 1
    return zero_line


def separator_pos(seq_line):
    sep_pos_array = np.array([m.start() for m in re.finditer('\|', seq_line)])
    sep_pos_array = np.insert(sep_pos_array, 0, 0)
    sep_pos_array = np.append(sep_pos_array, len(seq_line))
    return sep_pos_array


# the following function use aho-corasick result and pos_ID_dict as parameter to generate a dict with pep sequence as key and uniprotID as value.
def creat_pep_ID_dict(aho_result, pos_ID_dict):
    pep_ID_dict = defaultdict(set)
    for i in aho_result:
        pep_ID_dict[i[2]].add(pos_ID_dict[i[0]])
    return pep_ID_dict


# a dictionary that has Uniport ID as key, identified peptides for that ID as value
def creat_ID_pep_dict(aho_result, pos_ID_dict):
    ID_pep_dict = defaultdict(set)
    for i in aho_result:
        ID_pep_dict[pos_ID_dict[i[0]]].add(i[2])
    return ID_pep_dict
