# -*- coding: utf-8 -*-
"""
Created on Sat May 15 03:57:45 2021

Nussinov RNA folding algorithm + recursive backtrack. 
Implemented by: Muhammad Elsadany
"""
import numpy as np
import time




def init_matrix(rna):
    """
    initialize matrix with zeros where can't have pairings.
    NxN matrix that stores the scores of the optimal pairings.
    """
    M = len(rna)
    # init matrix
    nm = np.empty([M, M])
    nm[:] = np.NAN
    # init diaganols to 0
    # other ways to do this: np.fill_diaganol(), np.diag(), nested loop, ...
    nm[range(M), range(M)] = 0
    nm[range(1, len(rna)), range(len(rna) - 1)] = 0
    return nm


def couple(pair):
    """
    Return True if RNA nucleotides are Watson-Crick base pairs.
    """
    pairs = {"A": "U", "U": "A", "G": "C", "C": "G"} # or a list of tuples
    # check if pair is couplable
    if pair in pairs.items():
        return True
    
    return False

def fill(nm, rna):
    """
    Fill the matrix as per the Nussinov algorithm.
    """
    minimal_loop_length = 0
    for k in range(1, len(rna)):
        for i in range(len(rna) - k):
            j = i + k
            if j - i >= minimal_loop_length:
                down = nm[i + 1][j] # 1st rule
                left = nm[i][j - 1] # 2nd rule
                diag = nm[i + 1][j - 1] + couple((rna[i], rna[j])) # 3rd rule
                rc = max([nm[i][t] + nm[t + 1][j] for t in range(i, j)]) # 4th rule
                nm[i][j] = max(down, left, diag, rc) # max of all
                # returns the score of the optimal pairing between indices i and j
                            
            else:
                nm[i][j] = 0
    return nm


def traceback(nm, rna, fold, i, L):
    """
    Traceback through complete Nussinov matrix to find optimial RNA 
    secondary structure solution through max base-pairs.
    """
    j = L

    if i < j:
        if nm[i][j] == nm[i + 1][j]: # 1st rule
            traceback(nm, rna, fold, i + 1, j)
        elif nm[i][j] == nm[i][j - 1]: # 2nd rule
            traceback(nm, rna, fold, i, j - 1)
        elif nm[i][j] == nm[i + 1][j - 1] + couple((rna[i], rna[j])): # 3rd rule
            fold.append((i, j))
            traceback(nm, rna, fold, i + 1, j - 1)
        else:
            for k in range(i + 1, j - 1):
                if nm[i][j] == nm[i, k] + nm[k + 1][j]: # 4th rule
                    traceback(nm, rna, fold, i, k)
                    traceback(nm, rna, fold, k + 1, j)
                    break

    return fold



def dot_write(rna, fold):
    """ 
    returns a written form of pairing between indices from generated fold structure.
    """
    dot = ["." for i in range(len(rna))]

    for s in fold:
        dot[min(s)] = "("
        dot[max(s)] = ")"

    return "".join(dot)


def Nussinov(sequence):
    """ 
    returns the secondary structure of entered RNA sequence.
    """
    nm = init_matrix(sequence)
    nm = fill(nm, sequence)
    # print(nm)
    fold = traceback(nm, sequence, [], 0, len(sequence)-1)
    # print(fold)
    # print(dot_write(sequence, fold))
    return(fold, dot_write(sequence, fold))






##############################################################################
# main


            
# enter the file path for your RNA sequence            
File = open("G:/BMS 422/Project/sequence1kbp.txt", "r")
seq_ = File.readlines()
Seq = seq_[0]
Seq_Capitalized = Seq.upper()
# print(Seq_Capitalized)


Start_time = time.time()
print(Nussinov(Seq_Capitalized))
Finish_time = time.time()-Start_time
print("time taken: ", Finish_time)