#!/usr/bin/python
import sys
import argparse
import glob
import os

#create needed arguments
parser = argparse.ArgumentParser(description='Subset NCBI nt library to alignments found in sublibraries')
parser.add_argument("--ml", type=str, help="Location of full NCBI nt library")
parser.add_argument("--sa", type=str, help="Location of MALT aligned NCBI library subsets")
parser.add_argument("--o", type=str, help="Full output file path")

args=parser.parse_args()
ml = args.ml
sa = args.sa
o = args.o

print("Your ")

#define needed functions
def fasta2dict(fil):
    """
    Read fasta-format file fil, return dict of form scaffold:sequence.
    Note: Uses only the unique identifier of each sequence, rather than the
    entire header, for dict keys.
    """
    dic = {}
    cur_scaf = ''
    cur_seq = []
    for line in open(fil):
        if line.startswith(">") and cur_scaf == '':
            cur_scaf = line.split(' ')[0]
        elif line.startswith(">") and cur_scaf != '':
            dic[cur_scaf] = ''.join(cur_seq)
            cur_scaf = line.split(' ')[0]
            cur_seq = []
        else:
            cur_seq.append(line.rstrip())
    dic[cur_scaf] = ''.join(cur_seq)
    return dic

def get_aligned_keys(blast_txt):
    """
    Read in a BLAST text format file & export a list of species ID's matched by the alignment
    """
    unset_keys = []
    for line in open(blast_txt):
        check_line = line.split('|')[0]
        if len(check_line) > 20:
            check_line = line.split(' ')[0]
        unset_keys.append(check_line)
    set_keys = set(unset_keys)
    return set_keys

fasta_dict = fasta2dict(ml)
file_list = glob.glob(os.path.join(sa, "*.found"))

all_alignments = []
for i in len(file_list):
    d = get_aligned_keys(file_list[i])
    all_alignments.extend(d) #now all_alignments is all of the potential matches

sys.stdout = open(o, "w")
for key in all_aligned:
    print(key,'\n',fasta_dict[key])
sys.stdout.close()
