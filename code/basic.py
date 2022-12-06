import sys 
import numpy as np
from resource import * 
import time
import psutil
import os 

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
def mismatch_penalty(base_1,base_2):
    nucs="ACGT"
    index_1=nucs.index(base_1)
    index_2=nucs.index(base_2)
    mismatch_matrix = [[0,110,48,94],
                        [110,0,118,48],
                        [48,118,0,110],
                        [94,48,110,0]]
    return mismatch_matrix[index_1][index_2]

#alignment function that fills dp table and outputs the alignment for both sequences and the alignment score 
def alignment(seq1,seq2):
    # set value of gaps 
    gap_pen=30
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
                dp[i][j] = min(dp[i-1][j-1] + mismatch_penalty(seq1[i-1],seq2[j-1]),
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
    return dp[m][n],seq

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
    start_memory=process_memory()
    inputs=input_read(file)
    sequences=generate_sequences(inputs[0],inputs[1],inputs[2],inputs[3],inputs[4])
    alignment_output=alignment(sequences[0],sequences[1])
    end_time = time.time()
    end_memory=process_memory()
    time_taken = (end_time - start_time)*1000
    memory_used=(end_memory - start_memory)
    print(alignment_output[0],"\n",alignment_output[1][0],"\n",alignment_output[1][1],"\n",time_taken,"\n",memory_used)
    output= (alignment_output[0],alignment_output[1][0],alignment_output[1][1],time_taken,memory_used)
    output=np.array(output)
    np.savetxt(output_name, output, fmt='%s', newline='\n')

if __name__ == "__main__":
    main()