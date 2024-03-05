import random
import numpy as np
import pandas as pd
from ete3 import PhyloTree  # in-order for this package to run you need to install PyQt5 package
from Bio import SeqIO, AlignIO
from os.path import isdir, isfile


def scoring_matrix():
    # here we use a substitution matrix for DNA sequences,
    # the global alignment to be conducted in class Alignment
    score = [['', 'A', 'C', 'T', 'G'], ['A', '1', '-1', '-1', '-1'], ['C', '-1', '1', '-1', '-1'],
             ['T', '-1', '-1', '1', '-1'], ['G', '-1', '-1', '-1', '1']]

    # convert to dataframe
    score = pd.DataFrame(score)
    new_header = score.iloc[0]
    score = score[1:]
    score.columns = new_header
    score = score.set_index([''])

    return score  # return simpli score matrix for DNA alignment


def extract(file1, file2):
    """ This function extracts sequences from files for alignment """
    x = ''
    y = ''

    for i in SeqIO.parse(file1, 'fasta'):
        x = x + str(i.seq)

    for j in SeqIO.parse(file2, 'fasta'):
        y = y + str(j.seq)
        break

    score = scoring_matrix()
    test = Alignment(score, -1)
    test.get_alignment(x, y)
    test.probability_likelihood(x, y)
    return f"the first sequence is {x} \nwith length {len(x)}" \
        f"\nthe second is sequence is {y} \nwith length {len(y)}"


class Alignment:

    """
    sequence alignment tool
    requires a scoring matrix
    """

    def __init__(self, scoring_matrix, gap_penality):
        self.scoring_matrix = scoring_matrix
        self.gap_penality = gap_penality

    def get_alignment(self, seq1, seq2):
        # returns alignment string and score
        score, path_matrix = self._alignment_helper(seq1, seq2)
        align_top = ""
        matches = ""
        align_bottom = ""
        i = len(seq1) - 1
        j = len(seq2) - 1
        while (i, j) != (-1, -1):
            path = path_matrix[i+1][j+1]
            species1 = seq1[i]
            species2 = seq2[j]
            match = " "
            if path == 0:
                i -= 1
                j -= 1
                if species1 == species2:
                    match = "|"
            elif path == 1:
                species1 = "-"
                j -= 1
            else:
                species2 = "-"
                i -= 1
            align_top = species1 + align_top
            align_bottom = species2 + align_bottom
            matches = match + matches
        return print(f"{align_top}\n{matches}\n{align_bottom}\nscore: {score}")

    def _alignment_helper(self, seq1, seq2):
        # sets scoring matrix
        matrix, path_matrix = self.run_needle_wunsch(seq1, seq2)
        alignment_score = matrix[-1][-1]
        return alignment_score, path_matrix

    def run_needle_wunsch(self, seq1, seq2):
        # returns scoring matrix based on input sequence
        """
        [1,0]
        [2,x]
        x: our target
        """
        n = len(seq1)
        m = len(seq2)
        matrix = [[0 for i in range(m + 1)] for j in range(n + 1)]
        path_matrix = [[0 for i in range(m + 1)] for j in range(n + 1)]

        for i in range(1, n + 1):
            matrix[i][0] = self.gap_penality * i
            path_matrix[i][0] = 2

        for j in range(1, m + 1):
            matrix[0][j] = self.gap_penality * j
            path_matrix[0][j] = 1

        for i in range(1, n+1):
            for j in range(1, m+1):

                species1 = seq1[i-1]
                species2 = seq2[j-1]

                diag_score = matrix[i-1][j-1] + int(self.scoring_matrix[species1][species2])
                right_score = matrix[i-1][j] + self.gap_penality
                down_score = matrix[i][j-1] + self.gap_penality
                score = max(diag_score, max(right_score, down_score))

                matrix[i][j] = score
                if score == diag_score:
                    path_matrix[i][j] = 0
                elif score == down_score:
                    path_matrix[i][j] = 1
                else:
                    path_matrix[i][j] = 2

        return matrix, path_matrix

    def get_p_value(self, seq1, seq2, n=1e4):
        # n is the num. of possibilities that both sequences are actually related
        # following a bootstraping pseudo-method

        """ a function that gives a value corresponding to how much of a fluke the alignment could be """
        score, _ = self._alignment_helper(seq1, seq2)
        k = 0
        for i in range(int(n)):
            seq1_permutation = ''.join(random.sample(seq1, len(seq1)))
            score_permutation, _ = self._alignment_helper(seq1_permutation, seq2)
            if score_permutation >= score:
                k += 1
        return (k + 1) / (n + 1)

    def probability_likelihood(self, seq1, seq2):

        score, _ = self._alignment_helper(seq1, seq2)
        n = 2
        p = self.get_p_value(seq1, seq2)

        mean = x / n
        sigma = sum(x - math.square(mean)) / n - 1

        part_1 = (n * math.log10(float(1 / sigma * math.sqrt(2 * math.pi))))
        part_2 = float(p) * (sum(score ** 2 - mean) / 2 * sigma)
        mle = part_1 - part_2
        return mle


def phylo(sequences):
    # in this function, the gene_tree variable has to be filled by the user

    gene_tree = "((((P10, NCTC_8325), M48), (NRS1, AR465)), R50);"

    t = PhyloTree(gene_tree)
    t.link_to_alignment(alignment=fasta_txt, alg_format="fasta")

    return t


if __name__ == "__main__":
    #  first enter the names of the files as variables
    # or copy the absolute path of the files, and it will do the job
    file1 = "P10.fasta"
    file2 = "R50.fasta"
    file3 = "AR465.fasta"
    file4 = "M48.fasta"
    file5 = "NCTC_8325.fasta"
    file6 = "NRS1.fasta"
    #  then call the extract function which calls the alignment class
    #  this is repeated multiple times for however many sequences you enter
    seq_extract = extract(file1, file2)
    seq_extract1 = extract(file3, file4)
    seq_extract2 = extract(file5, file6)

    print(seq_extract)
    print(seq_extract1)
    print(seq_extract2)
