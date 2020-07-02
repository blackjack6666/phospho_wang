import ahocorasick



def automaton_trie(peptide_list):
    A = ahocorasick.Automaton()
    for idx, peptide in enumerate(peptide_list):
        A.add_word(peptide, (idx, peptide))
    A.make_automaton()
    return A

def automaton_matching(A, seq_line):
    result = []
    for end_idx, (insert_order, original_value) in A.iter(seq_line):
        start_idx = end_idx - len(original_value) + 1
        result.append((start_idx, end_idx, original_value))
        assert seq_line[start_idx:start_idx+len(original_value)] == original_value
    return result
