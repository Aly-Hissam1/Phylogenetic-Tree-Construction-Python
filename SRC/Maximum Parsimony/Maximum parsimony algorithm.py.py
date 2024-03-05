import numpy as np
from Bio.Seq import Seq
#input function to parse the files (fasta,seq,dna)
def input_function(fasta_file1):      
        fasta_open=open(fasta_file1)
        fasta1=fasta_open.readlines()
        fasta_seq1=""
        specie=""
        dict={}
        for line in fasta1:
            line=line.upper()
            if line[0]!=">":
                fasta_seq1+=line.rstrip()
            if line.startswith(">"):
                specie=line[1:]
            dict.update({specie:fasta_seq1})
        fasta_open.close()
        return dict

def traceback_and_Traceback_matrices(seq1,seq2):
    
    gap_open_penalty=8
    gap_extend_penalty=8
    match_reward=1
    mismatch_penalty=8
    
    #Initilaizing score matrix
    num_rows = len(seq2) + 1
    num_cols = len(seq1) + 1
    score_matrix = np.zeros(shape=(num_rows, num_cols), dtype=np.int)
    
    
    #Initilaizing Traceback matrix
    Traceback = np.full(shape=(num_rows, num_cols), fill_value=" ", dtype=np.str)
    
#    fillig the 1st row and 1st column in score_matrix
    score_matrix[0][0] = 0
    for i in range(1, num_rows):
        score_matrix[i][0] = score_matrix[i-1][0]-gap_open_penalty
    
    for j in range(1, num_cols):
        score_matrix[0][j] = score_matrix[0][j-1]-gap_open_penalty
    #print(F)
    
#    fillig the 1st row and 1st column in Traceback_matrix    
    Traceback[0][0] = "•"
    for i in range(1, num_rows):
        Traceback[i][0] = "↑"
    
    for j in range(1, num_cols):
        Traceback[0][j] = "←"
    
    
    match_mismatch="↖"
    vgap="↑"
    hgap="←"
  
        
    for seq2_pos, seq2_char in enumerate(seq2):
        seq2_pos +=1
# print(seq2_pos, seq2_char)
        for seq1_pos, seq1_char in enumerate(seq1):
            seq1_pos +=1
            
# compute the score for a match/mismatch
            if seq1_char==seq2_char:
                diag_score=(score_matrix[seq2_pos-1][seq1_pos-1]+match_reward,match_mismatch)
            else:
                diag_score=(score_matrix[seq2_pos-1][seq1_pos-1]-mismatch_penalty,match_mismatch)
            
# compute the score for adding a gap in seq2 (vertical)
            if seq1_pos==len(seq1):
                up_score=(score_matrix[seq2_pos-1][seq1_pos],vgap)
            elif Traceback[seq2_pos-1][seq1_pos]==vgap:
                up_score=(score_matrix[seq2_pos-1][seq1_pos]-gap_extend_penalty,vgap)
            else:
                up_score=(score_matrix[seq2_pos-1][seq1_pos]-gap_open_penalty,vgap)
    
                
# compute the score for adding a gap in seq1 (horizontal)
            if seq2_pos==len(seq2):
                left_score=(score_matrix[seq2_pos][seq1_pos-1],hgap)
            elif Traceback[seq2_pos-1][seq1_pos]==hgap:
                left_score=(score_matrix[seq2_pos][seq1_pos-1]-gap_extend_penalty,hgap)
            else:
                left_score=(score_matrix[seq2_pos][seq1_pos-1]-gap_open_penalty,hgap)
         
            
# identify the largest score, and use that information to populate
# the score and traceback matrices
            best_score=max(up_score,left_score,diag_score)
            score_matrix[seq2_pos][seq1_pos]=best_score[0]
            Traceback[seq2_pos][seq1_pos]=best_score[1]
    return (score_matrix,Traceback,best_score)
#    return score_matrix
#    return best_score


##Initializing Traceback function
def Traceback_function(score_matrix,Traceback,seq1,seq2):

    match_mismatch="↖"
    vgap="↑"
    hgap="←"
    aligned_seq1=""
    aligned_seq2=""
    
    x,y=Traceback.shape
    current_row = x-1
    current_col = y-1
    
    current_value = None
    gap_character="_"
    
    while current_value !=Traceback[0][0]:
        current_value = Traceback[current_row, current_col]
        if current_value == match_mismatch:
            aligned_seq1 =(seq1[current_col-1])+aligned_seq1
            aligned_seq2 =(seq2[current_row-1])+aligned_seq2
            current_row -= 1
            current_col -= 1
                
        elif current_value == vgap:
            aligned_seq1=(gap_character)+aligned_seq1
            aligned_seq2=(seq2[current_row-1])+aligned_seq2
            current_row -= 1
            
        elif current_value == hgap:
            aligned_seq1=(seq1[current_col-1])+aligned_seq1
            aligned_seq2=(gap_character)+aligned_seq2
            current_col -= 1
            
        elif current_value == Traceback[0][0]:
                continue
    return(aligned_seq1,aligned_seq2)

# Calling the alignment functions
def call_and_apply():
    seq_list=[]
    sorted_organism=[]
    dic = {}
    c=0
    aligned_seq_list=[]
    Entry=input(" Please choose the expected number for the input fasta files\n")
    for i in range(int(Entry)):
        fasta_input=input("Please Enter the path file " + str (i+1) + ":   ")
        dic=input_function(fasta_input)
        for key,value in dic.items():
            #organism.append(key)
            seq_list.append(value)
        #seq_list.append(input_function(fasta_input))
        del fasta_input
    seq_list.sort(key=len)
    #make a list with the name of sorted organisms based on the sorted list of sequences
    for i in seq_list:
        if i == dic.values():
            x = dic.keys()
            sorted_organism.append(x)

    for y in range(len(seq_list)-1):
        if c==0:
            x=traceback_and_Traceback_matrices(seq_list[-1],seq_list[y])
            seq1= seq_list[-1]
            seq2= seq_list[y]
            T=Traceback_function(x[0],x[1],seq1,seq2)
            del seq1
            del seq2
            aligned_seq_list.append(T[0])
            aligned_seq_list.append(T[1])
            c+=1
        else:
            x=traceback_and_Traceback_matrices(seq_list[-1],seq_list[y])
            seq1= seq_list[-1]
            seq2= seq_list[y]
            T=Traceback_function(x[0],x[1],seq1,seq2)
            del seq1
            del seq2
            aligned_seq_list.append(T[1])
    return (aligned_seq_list)


def informative (column):
    c=[]
    for i in range(len(column)-1):
        c_num=0
        for y in range(len(column)-(i+1)):
            if column[i] == column[-(y+1)]:
                c_num += 1
            else:
                 continue
        c.append(str(c_num))
        del c_num
    result = 'True'
    for num in c:
        if int(num) >1:
            del result
            result = "False"
            break
        else:
            continue
    return (result)

def informative_filter(aligned_seq_list):
    list_of_informatives=[]
    R = [list(i) for i in zip(*aligned_seq_list)] # flip the matrix to get the columns
    for lists in R:
        if informative(lists) == "True":
            list_of_informatives.append(lists)
        else:
            continue
    return (list_of_informatives)

#x=["TTA","CTA","CTG","CTG"]
#print(informative_filter(x))

alignment= call_and_apply()
print (alignment)
def _nni(self, list_of_informatives, alignment):
    """Search for the best parsimony tree using the NNI algorithm (PRIVATE)."""
    best_tree = list_of_informatives
    while True:
        best_score = self.scorer.get_score(best_tree, alignment)
        temp = best_score
        for t in self._get_neighbors(best_tree):
            score = self.scorer.get_score(t, alignment)
            if score < best_score:
                best_score = score
                best_tree = t
        # stop if no smaller score exist
        if best_score >= temp:
            break
    return best_tree


x= _nni(informative_filter(alignment),alignment)
print(x)


