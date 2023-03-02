import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
from sklearn.decomposition import TruncatedSVD
import random
# import seaborn as sns
import os.path as path
import os
# import matplotlib
# import matplotlib.font_manager
# import matplotlib.pyplot as plt # graphs plotting
# import Bio
from Bio import SeqIO # some BioPython that will come in handy
#matplotlib inline
import timeit
from itertools import cycle

from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
from sklearn.metrics import roc_auc_score


# from matplotlib import rc
# # for Arial typefont
# matplotlib.rcParams['font.family'] = 'Arial'

from sklearn.model_selection import train_test_split
from sklearn.linear_model import Lasso, LogisticRegression
from sklearn.feature_selection import SelectFromModel
from sklearn.preprocessing import StandardScaler

from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn import metrics
from sklearn import svm

from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn import metrics

import pandas as pd
import sklearn
from sklearn import preprocessing
from sklearn.model_selection import train_test_split 
from sklearn.preprocessing import StandardScaler  
from sklearn.neural_network import MLPClassifier 
from sklearn.metrics import classification_report, confusion_matrix 

from sklearn.neighbors import KNeighborsClassifier

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix

from pandas import DataFrame

from sklearn.model_selection import KFold 
from sklearn.model_selection import RepeatedKFold

from sklearn.metrics import confusion_matrix

from numpy import mean

import seaborn as sns

import itertools
from itertools import product

import math

filename = 'filesubmited'
input_file_name = os.path.join('uploads', filename)
sequences_dictionary = {sequence.id : sequence.seq for sequence in SeqIO.parse(input_file_name,'fasta')}
attribute_data = [s.id for s in SeqIO.parse(input_file_name,'fasta')]
seq_data = [s.seq for s in SeqIO.parse(input_file_name,'fasta')]


def pwm(sequence): 
    print('task: '+str(sequence))
    seq = seq_data[sequence]
    ################ Generate k-mers (Start) #########################
    L = len(seq)
    Kmer = 9
    k_mers_final = []
    for i in np.arange(0, L-Kmer+1):
        sub_f=seq[i:i+Kmer]
        k_mers_final.append(sub_f)
    #     print(sub_f)

    ################ Generate k-mers (end) #########################


    ################ Generate PWM (Start) #########################
    # To create a new txt file for writing "EI_nine_pwm.txt"
    # Initialize the PWM with four rows and nine columns [i.e., 4 lists of zeros]
    # a = [0]*9
    # c = [0]*9
    # g = [0]*9
    # t = [0]*9

    a_val = [0]*Kmer
    b_val = [0]*Kmer
    c_val = [0]*Kmer
    d_val = [0]*Kmer
    e_val = [0]*Kmer
    f_val = [0]*Kmer
    g_val = [0]*Kmer
    h_val = [0]*Kmer
    i_val = [0]*Kmer
    j_val = [0]*Kmer
    k_val = [0]*Kmer
    l_val = [0]*Kmer
    m_val = [0]*Kmer
    n_val = [0]*Kmer
    p_val = [0]*Kmer
    q_val = [0]*Kmer
    r_val = [0]*Kmer
    s_val = [0]*Kmer
    t_val = [0]*Kmer
    v_val = [0]*Kmer
    w_val = [0]*Kmer
    x_val = [0]*Kmer
    y_val = [0]*Kmer
    z_val = [0]*Kmer

    # input_file = open("E:/RA/Position Weight Matrix/Code/EI_nine.txt","r")   
    count_lines = 0 # Initialize the total number of sequences to 0
    # Read line by line, stripping the end of line character and
    # updating the PWM with the frequencies of each base at the 9 positions
    for ii in range(0,len(k_mers_final)):
        line = k_mers_final[ii]
        count_lines += 1 # Keep counting the sequences
        for i in range(0,len(line)):
            if line[i] == 'A':
                a_val[i] = a_val[i]+1
            elif line[i] == 'B':
                b_val[i] = b_val[i]+1
            elif line[i] == 'C':
                c_val[i] = c_val[i]+1
            elif line[i] == 'D':
                d_val[i] = d_val[i]+1
            elif line[i] == 'E':
                e_val[i] = e_val[i]+1
            elif line[i] == 'F':
                f_val[i] = f_val[i]+1
            elif line[i] == 'G':
                g_val[i] = g_val[i]+1       
            elif line[i] == 'H':
                h_val[i] = h_val[i]+1
            elif line[i] == 'I':
                i_val[i] = i_val[i]+1
            elif line[i] == 'J':
                j_val[i] = j_val[i]+1
            elif line[i] == 'K':
                k_val[i] = k_val[i]+1
            elif line[i] == 'L':
                l_val[i] = l_val[i]+1
            elif line[i] == 'M':
                m_val[i] = m_val[i]+1
            elif line[i] == 'N':
                n_val[i] = n_val[i]+1
            elif line[i] == 'P':
                p_val[i] = p_val[i]+1
            elif line[i] == 'Q':
                q_val[i] = q_val[i]+1
            elif line[i] == 'R':
                r_val[i] = r_val[i]+1
            elif line[i] == 'S':
                s_val[i] = s_val[i]+1
            elif line[i] == 'T':
                t_val[i] = t_val[i]+1
            elif line[i] == 'V':
                v_val[i] = v_val[i]+1            
            elif line[i] == 'W':
                w_val[i] = w_val[i]+1
            elif line[i] == 'X':
                x_val[i] = x_val[i]+1
            elif line[i] == 'Y':
                y_val[i] = y_val[i]+1
            elif line[i] == 'Z':
                z_val[i] = z_val[i]+1

    LaPlace_pseudocount = 0.1
    equal_prob_nucleotide = 0.04

    for i in np.arange(0,len(k_mers_final[0])):
        a_val[i] = round(math.log(((a_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        b_val[i] = round(math.log(((b_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        c_val[i] = round(math.log(((c_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        d_val[i] = round(math.log(((d_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        e_val[i] = round(math.log(((e_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        f_val[i] = round(math.log(((f_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        g_val[i] = round(math.log(((g_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        h_val[i] = round(math.log(((h_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        i_val[i] = round(math.log(((i_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        j_val[i] = round(math.log(((j_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        k_val[i] = round(math.log(((k_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        l_val[i] = round(math.log(((l_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        m_val[i] = round(math.log(((m_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        n_val[i] = round(math.log(((n_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        p_val[i] = round(math.log(((p_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        q_val[i] = round(math.log(((q_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        r_val[i] = round(math.log(((r_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        s_val[i] = round(math.log(((s_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        t_val[i] = round(math.log(((t_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        v_val[i] = round(math.log(((v_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        w_val[i] = round(math.log(((w_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        x_val[i] = round(math.log(((x_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        y_val[i] = round(math.log(((y_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
        z_val[i] = round(math.log(((z_val[i] + LaPlace_pseudocount)/(len(k_mers_final) + 0.4))/equal_prob_nucleotide,2),3)
#         a_val[i] = round(math.log((a_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         b_val[i] = round(math.log((b_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         c_val[i] = round(math.log((c_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         d_val[i] = round(math.log((d_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         e_val[i] = round(math.log((e_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         f_val[i] = round(math.log((f_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         g_val[i] = round(math.log((g_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         h_val[i] = round(math.log((h_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         i_val[i] = round(math.log((i_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         j_val[i] = round(math.log((j_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         k_val[i] = round(math.log((k_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         l_val[i] = round(math.log((l_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         m_val[i] = round(math.log((m_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         n_val[i] = round(math.log((n_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         p_val[i] = round(math.log((p_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         q_val[i] = round(math.log((q_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         r_val[i] = round(math.log((r_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         s_val[i] = round(math.log((s_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         t_val[i] = round(math.log((t_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         v_val[i] = round(math.log((v_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         w_val[i] = round(math.log((w_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         x_val[i] = round(math.log((x_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         y_val[i] = round(math.log((y_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)
#         z_val[i] = round(math.log((z_val[i] + LaPlace_pseudocount)/(count_lines + 0.4)/equal_prob_nucleotide,2),3)

    ################ Generate PWM (End) #########################

    ################ Assign Individual k-mers Score (Start) #########################
    each_k_mer_score = []
    for ii in range(len(k_mers_final)):
        line = k_mers_final[ii]
        score = 0
        for i in range(0,len(line)):
            if line[i] == 'A':
                score += a_val[i]
            elif line[i] == 'B':
                score +=  b_val[i]
            elif line[i] == 'C':
                score += c_val[i]
            elif line[i] == 'D':
                score += d_val[i]
            elif line[i] == 'E':
                score += e_val[i]
            elif line[i] == 'F':
                score += f_val[i]
            elif line[i] == 'G':
                score += g_val[i]       
            elif line[i] == 'H':
                score += h_val[i]
            elif line[i] == 'I':
                score += i_val[i]
            elif line[i] == 'J':
                score += j_val[i]
            elif line[i] == 'K':
                score += k_val[i]
            elif line[i] == 'L':
                score += l_val[i]
            elif line[i] == 'M':
                score += m_val[i]
            elif line[i] == 'N':
                score += n_val[i]
            elif line[i] == 'P':
                score += p_val[i]
            elif line[i] == 'Q':
                score += q_val[i]
            elif line[i] == 'R':
                score += r_val[i]
            elif line[i] == 'S':
                score += s_val[i]
            elif line[i] == 'T':
                score += t_val[i]
            elif line[i] == 'V':
                score += v_val[i]            
            elif line[i] == 'W':
                score += w_val[i]
            elif line[i] == 'X':
                score += x_val[i]
            elif line[i] == 'Y':
                score += y_val[i]
            elif line[i] == 'Z':
                score += z_val[i]
        each_k_mer_score.append(round(score, 3))
    return sequence,each_k_mer_score