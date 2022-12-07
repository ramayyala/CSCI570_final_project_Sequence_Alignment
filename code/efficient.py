import sys 
import numpy as np
from resource import * 
import time
import psutil
import os
import tracemalloc


#define gloabl variables
mismatch_matrix = [[0,110,48,94],
                        [110,0,118,48],
                        [48,118,0,110],
                        [94,48,110,0]]
nuc_dict={
        'A' : 0,
        'C' : 1,
        'G' : 2,
        'T' : 3
    }

gap_pen=30

#Parsing Input File 
def input_read(file):
    base_strings = []
    index_1 = []
    index_2 = []
    j=0
    k=0
    string_counter=0
    with open(file,"r") as f:
        for line in f.readlines():
            #remove all whitespaces in text file from reading in line by line 
            line=line.strip()
            # select only the base sequence to append to the sequence list 
            if not line.isnumeric():
                base_strings.append(line)
                #counter to make sure we are only appending strings not not the indices 
                string_counter+=1
            ##if the line is not a base sequence, then it is an index
            else:
                #if the string counter ==2, then that means are appending the second set of indices
                if string_counter == 2:
                    index_2.append(int(line))
                    k+=1
                #if the string !=2, then that means we are appending the first set of indices 
                else:
                    index_1.append(int(line))
                    j+=1
    return base_strings,index_1,index_2,k,j

#Generating Sequences from input_read output 
def generate_sequences(base_strings,index_1,index_2,k,j):
    sequences=[]
    seq1=base_strings[0]
    seq2=base_strings[1]
    ##check if its string 1 or string 2 we are using 
    for i in index_1:
        seq1 = seq1[0:i+1] + seq1 + seq1[i+1:]
    if len(seq1) == (2**j)*len(base_strings[0]):
        sequences.append(seq1)
    for i in index_2:
        seq2 = seq2[0:i+1] + seq2 + seq2[i+1:]
    if len(seq2) == (2**k)*len(base_strings[1]):   
        sequences.append(seq2)
    return sequences

#mismatch penalty function to determine mismatch penalty from matrix 
def mismatch_penalty(base_1,base_2,nuc_lib, mis_matrix):
    return mis_matrix[nuc_lib[base_1]][nuc_lib[base_2]]

#alignment function that fills dp table and outputs the alignment for both sequences and the alignment score 
def alignment(seq1,seq2):
    # build dp table
    m=len(seq1)
    n=len(seq2)
    dp = np.zeros([m+1,n+1], dtype=int)
    #initialize dp table 
    dp[0:(m+1),0] = [i*gap_pen for i in range(m+1)]
    dp[0,0:(n+1)] = [j*gap_pen for j in range(n+1)]
    #fill dp table with penalty scores 
    for i in range(1,m+1):
        for j in range(1,n+1):
            # if the nucs are equal to each other 
            if seq1[i-1] == seq2[j-1]:
                dp[i][j] = dp[i-1][j-1]
            # if the nucs don't match each other, then find min of the potential penalties 
            else:
                dp[i][j] = min(dp[i-1][j-1] + mismatch_penalty(seq1[i-1],seq2[j-1], nuc_dict, mismatch_matrix),
                                dp[i-1][j] + gap_pen,
                                dp[i][j-1] + gap_pen)
    
    # find optimal path from dp table 
    i= m
    j= n
    #tracing path 
    seq1_list=[]
    seq2_list=[]
    while i!=0 or j!=0:
        #case 1: seq2 has gap 
        if dp[i][j]==dp[i-1][j]+gap_pen:
            seq1_list.append(seq1[i-1])
            seq2_list.append('_')
            i-=1
        #case 2: seq1 has gap 
        elif dp[i][j]==dp[i][j-1]+gap_pen:
            seq1_list.append('_')
            seq2_list.append(seq2[j-1])
            j-=1

        #case 3: nucs match, so no penalty or nucs mismatch, but penalty is already accounted for so same path 
        else:
            seq1_list.append(seq1[i-1])
            seq2_list.append(seq2[j-1])
            i-=1
            j-=1


    s1_final=''
    s2_final=''
    for k in range(len(seq1_list)):
        s1_final += seq1_list.pop()
        s2_final += seq2_list.pop()
    seq = [s1_final, s2_final]
    return seq

# dc alignment algo 
def dc_alignment(seq1,seq2):
 
    m=len(seq1)
    n=len(seq2)
    #if the length of the sequences is less than 2, then just use the dp algorithm as a bound 
    if n<=2 or m<=2:
        return alignment(seq1,seq2)
    f = np.zeros((m + 1, 2))
    g = np.zeros((m + 1, 2))

    if n%2==0:
        front_len=n//2
        back_len=n//2
    else:
        front_len=n//2+1
        back_len=n//2


    for i in range(m + 1):
        f[i, 0] = i * gap_pen

    for j in range(1, front_len+1):
        f[0, 1] = j * gap_pen
        for i in range(1, m + 1):
            f[i, 1] = min(f[i - 1, 0] + mismatch_penalty(seq1[i-1],seq2[j-1], nuc_dict, mismatch_matrix),
                            f[i - 1, 1] + gap_pen,
                            f[i, 0] + gap_pen)
        f=np.flip(f,1)
    
    #reverse the strings
    seq1_inv=seq1[::-1]
    seq2_inv=seq2[::-1]

    for i in range(m + 1):
        g[i, 0] = i * gap_pen

    for j in range(1, back_len+1):
        g[0, 1] = j * gap_pen
        for i in range(1, m + 1):
            g[i, 1] = min(g[i - 1, 0] + mismatch_penalty(seq1_inv[i-1],seq2_inv[j-1],nuc_dict, mismatch_matrix),
                            g[i - 1, 1] + gap_pen,
                            g[i, 0] + gap_pen)
        g=np.flip(g,1)
    q = 0
    for i in range(m + 1):
        if f[i, 0] + g[m - i, 0] < f[q, 0] + g[m - q, 0]:
            q = i
    # print(f[q, 0] + g[m - q, 0])
    #print("F table:")
    #print(f)
    #print("G table:")
    #print(g)
    alignment_x1, alignment_y1 = dc_alignment(seq1[:q], seq2[:front_len])
    alignment_x2, alignment_y2 = dc_alignment(seq1[q:], seq2[front_len:])
    return alignment_x1 + alignment_x2, alignment_y1 + alignment_y2

#score function to calculate alginment score 
def score(output):
    mismatch=0
    gap=0
    for i in range(len(output[0])):
        if output[0][i]!=output[1][i]:
            if output[0][i]=='_' or output[1][i]=='_':
                gap+=30
            else:
                mismatch+=mismatch_penalty(output[0][i],output[1][i],nuc_dict, mismatch_matrix)
    score=mismatch+gap
    return score 


#memory function to track memory 
def process_memory():
    pid=os.getpid()
    process = psutil.Process(pid) 
    memory_info = process.memory_full_info() 
    memory_consumed = int(memory_info.rss/1024)
    return memory_consumed

#main function that orders the workflow of the functions to finding alignment 
def main():
    file = sys.argv[1]
    output_name=sys.argv[2]
    start_time = time.time() 
    tracemalloc.start(25)
    inputs=input_read(file)
    sequences=generate_sequences(inputs[0],inputs[1],inputs[2],inputs[3],inputs[4])
    alignment_output=dc_alignment(sequences[0],sequences[1])
    alignment_score=score(alignment_output)
    end_time = time.time()
    size, peak = tracemalloc.get_traced_memory()
    time_taken = (end_time - start_time)*1000

    #print(alignment_output[0],"\n",alignment_output[1][0],"\n",alignment_output[1][1],"\n",time_taken,"\n",memory_used)
    output= (alignment_score,alignment_output[0],alignment_output[1],time_taken,(peak/1024))
    output=np.array(output)
    np.savetxt(output_name, output, fmt='%s', newline='\n')

if __name__ == "__main__":
    main()
    
    