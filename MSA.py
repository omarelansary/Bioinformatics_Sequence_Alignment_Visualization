from Bio.Align.Applications import MuscleCommandline
from io import StringIO 
from Bio.Align.Applications import ClustalwCommandline
from plotly.offline import *
import subprocess
from Bio import SeqIO
import numpy as np
import numpy
import plotly.graph_objects as go
from math import log
dic={'A': 10, 'B': 20, 'C': 30, 'D': 40, 'E': 50, 'F': 60, 'G': 70, 'H': 80, 'I': 90, 'J':100, 'K': 110, 'L': 120, 'M': 130, 'N': 140, 'O': 150, 'P': 160, 'Q': 170, 'R': 180, 'S': 190, 'T': 200, 'U': 210, 'V': 220, 'W': 230, 'X': 235, 'Y': 240, 'Z': 245,'-':100}
# # from StringIO import StringIO
# from Bio import AlignIO
# #####
# muscle_exe = r"C:\Muscle.exe" #  path
# # MSA using Clustal2 and Generating output FASTA file of result
# # cmd = MuscleCommandline(muscle_exe, infile=filepath,
# #                             type="DNA", output='FASTA', outfile=MSA_outputfile)
# # std_out, std_err = cmd()
# ######
# cline = MuscleCommandline(input="msa.fasta", out="msa.txt",clw=True)
# print(cline)

#=================
# muscle_exe= r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\muscle5.1.win64.exe" # path
# in_file = r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\msa.fasta"
# out_file = r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\msa_output.txt"
# muscle_cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
# print(muscle_cline)
#===================
# cmd = ClustalwCommandline("clustalw2",
# infile=r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\msa.fasta")
# print(cmd)
# stdout, stderr = cmd()
#========================
#Align sequences with MUSCLE (using parameters to make the alignment
#process as fast as possible)
#==============================================================
# muscle_cline = MuscleCommandline(input="E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\msa.fasta", 
#                                  out="SARS-CoV-2_aligned.fasta", 
#                                  diags = True, 
#                                  maxiters = 1, 
#                                  log="../../data/raw/align_log.txt")
# muscle_cline()
# muscle_exe= r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\muscle5.1.win64.exe" # path
# in_file = r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\msa.fasta"
# out_file = r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\msa_output.fasta"
# muscle_result = subprocess.check_output([muscle_exe, "-in", in_file, "-out", out_file])

import subprocess

#define compute function
def Compute_Scores(a,b,match,mismatch,gap):
    scores_list=[]
    colors_seq1=[]
    colors_seq2=[]
    for x, y in zip(a, b):
        col1=int(dic[x])
        col2=int(dic[y])
        colors_seq1.append(col1)
        colors_seq2.append(col2)
        if x == y:
            scores_list.append(match)
        elif x!=y:
            if x=='-' or y=='-':
                scores_list.append(gap)
            else:
                scores_list.append(mismatch)
    tot_score=np.sum(scores_list)
    scores_list.append(tot_score)
    Matrix_colors=[
    colors_seq2,
    colors_seq1,
    ]
    # Matrix_colors=np.row_stack((colors_seq2,colors_seq1))
    
    # print("after matrix cols:")
    # print(Matrix_colors)
    return scores_list,Matrix_colors

class Frequency_dictionary(dict): 
    def __init__(self): 
        self = dict()   
    def add(self, key, value): 
        self[key] = value 

def compare(a,b):
        identical_count = 0
        all_pairs_count=0
        sum_of_pairs=0
        for x, y in zip(a, b):
            all_pairs_count+=1
            if ((x == y) and (x!='-') and(y!='-')):
                identical_count += 1
                sum_of_pairs+=1
            elif(((x=='-')and (y!='-')) or ((y=='-')and(x!='-'))):
                sum_of_pairs+=-1
            elif((x=='-')and (y=='-')):
                sum_of_pairs+=0
            else:  #mismatch       
                sum_of_pairs+=-1           
            

        return identical_count,all_pairs_count,sum_of_pairs
        
def Mutual_Information(records):

    freq_dic = Frequency_dictionary() 

    for i in range(len(records)):
        Curr_seq=records[i].seq
        for residue in Curr_seq:
            if residue in freq_dic:
               freq_dic[residue] +=1
            else:
               freq_dic.add(residue,1)

    for i in range(len(records)):
        for j in range(i+1,len(records)):
            seq1=records[i].seq
            seq2=records[j].seq
            for residue1,residue2 in zip(seq1, seq2):
                if residue1+residue2 in freq_dic:
                   freq_dic[residue1+residue2] +=1
                else:
                   freq_dic.add(residue1+residue2,1)


    total_residues=sum(freq_dic.values())
    Mutual_Information=0
    for i in range(len(records)):
        for j in range(i+1,len(records)):
            seq1=records[i].seq
            seq2=records[j].seq
            for residue1,residue2 in zip(seq1, seq2):
                prob_residue1_residue2=(freq_dic[residue1+residue2])/total_residues
                prob_residue1=freq_dic[residue1]/total_residues
                prob_residue2=freq_dic[residue2]/total_residues
                Mutual_Information+=log(prob_residue1_residue2/(prob_residue1*prob_residue2))
    print("Mutual Information")
    print(Mutual_Information)
    Normalized_mutual_information=Mutual_Information/log(total_residues)
    print("Normalized Mutual Information ")
    print(Normalized_mutual_information)
    return Mutual_Information,Normalized_mutual_information

def MSA(input):
    # output = subprocess.check_output(
    #     [ r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\muscle5.1.win64.exe",
    #     "-align", r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\Group_5.fasta",
    #     "-output", r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\Group_5_msa_output.fasta"],
    #     text=True)
    output = subprocess.check_output(
        [ r"muscle5.1.win64.exe",
        "-align", input,
        "-output", r"Group_5_msa_output.fasta"],
        text=True)#ha
    #input="E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\"+input
    # output = subprocess.check_output(
    #     [ r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\muscle5.1.win64.exe",
    #     "-align", r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\Group_5.fasta",
    #     "-output", r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\Group_5_msa_output.fasta"],
    #     text=True)
    records = list(SeqIO.parse("Group_5_msa_output.fasta", "fasta"))
    # subprocess.call('muscle -align %s -output %s'%(in_file,out_file))
    #Open fasta file
    # recs = SeqIO.parse("E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\Group_5_msa_output.fasta", 'fasta')
    # print(recs)
    # rec = next(recs)
    # print(rec)  
    seq1=[]
    seq2=[]
    total_pairs=0
    count=0
    countall=0
    print(len(records))
    matrix_seqs=[]
        #yesssss
    for i in range(len(records)):
        l1=records[i].seq
        print("MSA records")
        print(l1)
        seq_2l=[i for i in l1]
        matrix_seqs.append(seq_2l)
        for j in range(i+1,len(records)):
            seq1=records[i].seq
            seq2=records[j].seq
            #For percent identity analysis
            identicalcount,countall,sum_of_pairs=compare(seq1,seq2)
            total_pairs+=countall
            # sop=sop+identicalcount
    # print("sop= ",sop)
    try:
        percent_Identity=identicalcount/total_pairs*100
        print("Percent Identity:")
        print(percent_Identity)
    except:
        print("No pairs")
    #Mutual Info
    MutualInfo,NormalizedMutualInfo=Mutual_Information(records)

    
    #================================================================
    # ScoresList,Matrix_colors=Compute_Scores(,match,mismatch,gap)

    Matrix_colors=np.zeros((len(records),len(records[0].seq)))
    for i in range(len(records)):
        for j in range(len(records[0].seq)):
            Matrix_colors[i][j]+=dic[matrix_seqs[i][j]]
    Matrix_colors.tolist()
    print("Matrix after call:")
    print(Matrix_colors)
    fig = go.Figure(data=go.Heatmap(
                    z=Matrix_colors
                                    ,
                    text=matrix_seqs,
                    texttemplate="%{text}",
                    textfont={"size":5}))

    fig.write_image("MultipleSeq_Align.png")
    #fig.show()
    path=plot(fig)
    return path,percent_Identity,sum_of_pairs,MutualInfo,NormalizedMutualInfo
    

        
# MSA("E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Final_Integration\\Group_5.fasta")  
  
# print(records[0].seq)  # first record
# print(records[-1].seq)  # last record
#sop=  2936
#sop=  8764