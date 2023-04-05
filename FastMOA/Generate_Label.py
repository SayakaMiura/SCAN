import sys
import os
from Bio import SeqIO


path = ''#os.path.join(os.path.expanduser('~'), os.getcwd())+os.sep


def find(S, P):
    count = ['1' for x in range(len(S)) if S[x][P] =='T']
    return P, count

filename = str(sys.argv[1])
outDir = str(sys.argv[2])
Remove_Doublet = str(sys.argv[3])

if Remove_Doublet == 'true':
    fasta_sequences = SeqIO.parse(open(path + outDir + os.sep + filename),'fasta')
if Remove_Doublet == 'false':
    fasta_sequences = SeqIO.parse(open(path +filename),'fasta')

SEQ = []
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    SEQ.append(sequence)

Normal = SEQ[0]
label_names = []
days = []

for n in range(len(Normal)):
    value = find(SEQ, n)
    label_names.append(value[0])
    days.append(value[1])


label_file = open(path+outDir + os.sep +'Labels.txt', 'w')
for i in range(len(label_names)):
    
    label_file.write(str(label_names[i]+1) + '\t1\n')	
label_file.close()


